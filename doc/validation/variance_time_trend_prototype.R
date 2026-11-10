# Source the three unchanged preceding diagnostic modules first.
tv_base_kernel <- iv_kernel
tv_base_inputs <- iv_inputs
tv_base_parameters <- iv_parameters

tv_kernel <- function(source, destination) {
  stage <- tempfile(fileext = ".stan"); on.exit(unlink(stage), add = TRUE)
  tv_base_kernel(source, stage)
  text <- paste(readLines(stage, warn = FALSE), collapse = "\n")
  once <- function(old, new) {
    stopifnot(length(gregexpr(old, text, fixed = TRUE)[[1L]]) == 1L,
      grepl(old, text, fixed = TRUE))
    text <<- sub(old, new, text, fixed = TRUE)
  }
  once("  int<lower=0, upper=1> asset_variance_active;", paste(
    "  int<lower=0, upper=1> variance_time_trend_active;",
    "  vector[N] variance_time01;",
    "  real<lower=0> variance_time_prior_sd;",
    "  int<lower=0, upper=1> asset_variance_active;", sep = "\n"))
  once("  real log_sigma_share;", paste(
    "  vector[variance_time_trend_active] variance_time_gamma;",
    "  real log_sigma_share;", sep = "\n"))
  once("* exp(asset_variance_active == 1 ? asset_gamma * institution_log_assets_z[i] : 0.0);",
    paste("* exp((asset_variance_active == 1 ? asset_gamma * institution_log_assets_z[i] : 0.0)",
      "      + (variance_time_trend_active == 1 ? 0.5 * variance_time_gamma[1] * variance_time01[i] : 0.0));", sep = "\n"))
  once("  target += normal_lpdf(b_group_centered | 0, sigma_b_strat_economic);", paste(
    "  target += normal_lpdf(variance_time_gamma | 0, variance_time_prior_sd);",
    "  target += normal_lpdf(b_group_centered | 0, sigma_b_strat_economic);", sep = "\n"))
  writeLines(c("// History-only linear log-VARIANCE trend; base scale is defined at 1960.", text), destination)
  invisible(destination)
}

tv_inputs <- function(paper, template, config) {
  prepared <- tv_base_inputs(paper, template, config)
  for (name in names(prepared$inputs)) {
    years <- as.integer(prepared$inputs[[name]]$rows$year)
    time <- (years - config$trend_base_year) / config$trend_scale_years
    stopifnot(all(is.finite(time)), all(time >= 0), all(time <= 1))
    prepared$inputs[[name]]$data$variance_time_trend_active <- 0L
    prepared$inputs[[name]]$data$variance_time01 <- time
    prepared$inputs[[name]]$data$variance_time_prior_sd <- config$variance_time_prior_sd
  }
  prepared
}

tv_parameters <- function(fit, input, candidate, output) {
  tv_base_parameters(fit, input, candidate, output)
  if (candidate) {
    d <- ds_private(".bayes_backend_extract")(fit, variables = "variance_time_gamma")
    gamma <- as.numeric(d$variance_time_gamma)
    summarize <- function(x) c(Mean = mean(x), Lower = unname(quantile(x, .025)), Upper = unname(quantile(x, .975)))
    p <- read.csv(file.path(output, "parameter_summary.csv"))
    for (name in c("variance_time_gamma", "variance_ratio_1965_1960", "sd_ratio_1965_1960")) {
      x <- switch(name, variance_time_gamma = gamma,
        variance_ratio_1965_1960 = exp(gamma), sd_ratio_1965_1960 = exp(.5 * gamma))
      p <- rbind(p, data.frame(Parameter = name, t(summarize(x))))
    }
    write.csv(p, file.path(output, "parameter_summary.csv"), row.names = FALSE)
    scales <- do.call(rbind, lapply(0:5, function(t) data.frame(Year = 1960 + t,
      t(summarize(exp(gamma * t / 5))))))
    write.csv(scales, file.path(output, "variance_time_multipliers.csv"), row.names = FALSE)
  }
}

tv_strip_data <- function(data) {
  data[c("variance_time_trend_active", "variance_time01", "variance_time_prior_sd")] <- NULL
  data
}

tv_reference <- function(previous, source_cell, input, production, config) {
  stopifnot(source_cell %in% c("history_reference", "history_candidate"))
  hash <- function(path) digest::digest(file = path, algo = "sha256")
  seals <- read.csv(file.path(previous, "artifact_seals.csv"), stringsAsFactors = FALSE)
  verify <- function(path, records = seals) {
    record <- records[records$Path == normalizePath(path, mustWork = TRUE), , drop = FALSE]
    stopifnot(nrow(record) == 1L, identical(hash(path), record$SHA256[[1L]]))
  }
  verify(file.path(previous, "manifest.rds"))
  manifest <- readRDS(file.path(previous, "manifest.rds"))
  stopifnot(identical(manifest$production_kernel_sha256, hash(production)),
    all(vapply(names(manifest$source_hash), function(path)
      identical(hash(path), manifest$source_hash[[path]]), logical(1))))
  if (source_cell == "history_reference") {
    fit_path <- manifest$reused_reference_paths[[source_cell]]
    stopifnot(identical(hash(fit_path), manifest$reused_references[[source_cell]]$fit_sha256))
    original_seals <- read.csv(file.path(dirname(dirname(fit_path)), "completed_artifact_seals.csv"))
    verify(fit_path, original_seals)
  } else {
    fit_path <- file.path(previous, source_cell, "fit.rds")
    verify(fit_path)
    original_seals <- seals
  }
  saved <- readRDS(fit_path)
  stopifnot(identical(iv_base_data(saved$data), iv_base_data(tv_strip_data(input$data))),
    identical(saved$config$chains, config$chains), identical(saved$config$warmup, config$warmup),
    identical(saved$config$retained_per_chain, config$retained_per_chain))
  for (path in saved$native_fit$output_files()) verify(path, original_seals)
  if (source_cell == "history_reference") {
    stopifnot(identical(saved$kernel_sha256, hash(production)))
  } else {
    expected <- tv_strip_data(input$data); expected$asset_variance_active <- 1L
    stopifnot(identical(saved$data, expected), isTRUE(saved$candidate),
      identical(saved$config$fixed_variance_exponent, config$fixed_variance_exponent))
    kernel <- file.path(previous, "build/inverse_asset_variance_logit_supply.stan")
    verify(kernel)
    generated <- tempfile(fileext = ".stan"); on.exit(unlink(generated), add = TRUE)
    tv_base_kernel(production, generated)
    stopifnot(identical(hash(kernel), saved$kernel_sha256), identical(hash(generated), saved$kernel_sha256))
  }
  files <- file.path(previous, source_cell, c("policy_draws.rds", "policy_validity.csv",
    "comparison.csv", "sampler_diagnostics.csv", "parameter_summary.csv"))
  for (path in files) verify(path)
  policy <- readRDS(files[[1L]])
  stopifnot(length(policy$indices) == config$policy_draws)
  list(fit_path = normalizePath(fit_path), fit_sha256 = hash(fit_path), files = files,
    copied_policy_seed = config$policy_seed, policy_indices = policy$indices)
}
