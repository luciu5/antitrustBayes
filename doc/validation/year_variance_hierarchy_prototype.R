# Source unchanged demand-shock, asset, inverse-asset and trend helpers first.
yh_base_kernel <- iv_kernel
yh_trend_kernel <- tv_kernel
yh_base_inputs <- iv_inputs
yh_base_parameters <- iv_parameters

yh_kernel <- function(source, destination) {
  stage <- tempfile(fileext = ".stan"); on.exit(unlink(stage), add = TRUE)
  yh_base_kernel(source, stage)
  text <- paste(readLines(stage, warn = FALSE), collapse = "\n")
  once <- function(old, new) {
    stopifnot(length(gregexpr(old, text, fixed = TRUE)[[1L]]) == 1L,
      grepl(old, text, fixed = TRUE))
    text <<- sub(old, new, text, fixed = TRUE)
  }
  once("  int<lower=0, upper=1> asset_variance_active;", paste(
    "  int<lower=0, upper=1> year_variance_active;",
    "  int<lower=2> N_variance_year;",
    "  array[N] int<lower=1, upper=N_variance_year> variance_year_idx;",
    "  matrix[N_variance_year, N_variance_year-1] variance_year_basis;",
    "  real<lower=0> year_logvar_prior_scale;",
    "  int<lower=0, upper=1> asset_variance_active;", sep = "\n"))
  once("  real log_sigma_share;", paste(
    "  vector[year_variance_active * (N_variance_year-1)] year_variance_raw;",
    "  vector<lower=0>[year_variance_active] year_logvar_sd;",
    "  real log_sigma_share;", sep = "\n"))
  once("transformed parameters {", paste("transformed parameters {",
    "  vector[N_variance_year] year_log_variance_effect = rep_vector(0.0, N_variance_year);",
    "  vector[N_variance_year] year_log_variance_offset = rep_vector(0.0, N_variance_year);", sep = "\n"))
  once("  for (i in 1:N) sigma_share_observation[i] = sigma_share_market[market[i]]", paste(
    "  if (year_variance_active == 1) {",
    "    year_log_variance_effect = year_logvar_sd[1] * variance_year_basis * year_variance_raw;",
    "    year_log_variance_offset = year_log_variance_effect - rep_vector(year_log_variance_effect[1], N_variance_year);",
    "  }",
    "  for (i in 1:N) sigma_share_observation[i] = sigma_share_market[market[i]]", sep = "\n"))
  once("* exp(asset_variance_active == 1 ? asset_gamma * institution_log_assets_z[i] : 0.0);",
    paste("* exp((asset_variance_active == 1 ? asset_gamma * institution_log_assets_z[i] : 0.0)",
      "      + 0.5 * year_log_variance_offset[variance_year_idx[i]]);", sep = "\n"))
  once("  target += normal_lpdf(b_group_centered | 0, sigma_b_strat_economic);", paste(
    "  target += std_normal_lpdf(year_variance_raw);",
    "  target += normal_lpdf(year_logvar_sd | 0, year_logvar_prior_scale)",
    "    - num_elements(year_logvar_sd) * normal_lccdf(0 | 0, year_logvar_prior_scale);",
    "  target += normal_lpdf(b_group_centered | 0, sigma_b_strat_economic);", sep = "\n"))
  writeLines(c("// Noncentered hierarchical year log-variance effects; base scale at 1960.", text), destination)
  invisible(destination)
}

yh_inputs <- function(paper, template, config) {
  prepared <- yh_base_inputs(paper, template, config)
  prepared$inputs <- prepared$inputs["history"]
  years <- as.character(config$history_years); n <- length(years)
  stopifnot(n == 6L, years[1L] == config$target_year)
  h <- contr.helmert(n)
  h <- sweep(h, 2L, sqrt(colSums(h^2)), "/")
  stopifnot(max(abs(crossprod(h) - diag(n - 1L))) < 1e-12,
    max(abs(colSums(h))) < 1e-12)
  input <- prepared$inputs$history
  index <- match(input$rows$year, years)
  stopifnot(!anyNA(index), setequal(index, seq_len(n)))
  input$data$year_variance_active <- 0L
  input$data$N_variance_year <- n
  input$data$variance_year_idx <- index
  input$data$variance_year_basis <- h
  input$data$year_logvar_prior_scale <- config$year_logvar_prior_scale
  prepared$inputs$history <- input
  prepared$variance_year_labels <- years
  prepared
}

yh_parameters <- function(fit, input, candidate, output) {
  yh_base_parameters(fit, input, candidate, output)
  if (candidate) {
    d <- ds_private(".bayes_backend_extract")(fit,
      variables = c("year_logvar_sd", "year_log_variance_offset"))
    summarize <- function(x) c(Mean = mean(x), Lower = unname(quantile(x, .025)), Upper = unname(quantile(x, .975)))
    p <- read.csv(file.path(output, "parameter_summary.csv"))
    p <- rbind(p, data.frame(Parameter = "year_logvar_sd", t(summarize(as.numeric(d$year_logvar_sd)))))
    write.csv(p, file.path(output, "parameter_summary.csv"), row.names = FALSE)
    labels <- sort(unique(as.integer(input$rows$year)))
    scales <- do.call(rbind, lapply(seq_along(labels), function(j) data.frame(Year = labels[j],
      t(summarize(exp(d$year_log_variance_offset[, j]))))))
    write.csv(scales, file.path(output, "variance_year_multipliers.csv"), row.names = FALSE)
  }
}

yh_strip_data <- function(data) {
  data[c("year_variance_active", "N_variance_year", "variance_year_idx",
    "variance_year_basis", "year_logvar_prior_scale")] <- NULL
  data
}

yh_reference <- function(previous, cell, input, production, config) {
  stopifnot(cell %in% c("history_reference", "history_inverse", "history_trend"))
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
  if (cell %in% names(manifest$reused_reference_paths)) {
    fit_path <- manifest$reused_reference_paths[[cell]]
    stopifnot(identical(hash(fit_path), manifest$reused_references[[cell]]$fit_sha256))
    ledger <- if (cell == "history_reference") "completed_artifact_seals.csv" else "artifact_seals.csv"
    original_seals <- read.csv(file.path(dirname(dirname(fit_path)), ledger))
    verify(fit_path, original_seals)
  } else {
    fit_path <- file.path(previous, cell, "fit.rds"); verify(fit_path)
    original_seals <- seals
  }
  saved <- readRDS(fit_path)
  expected <- yh_strip_data(input$data)
  stopifnot(identical(iv_base_data(tv_strip_data(saved$data)), iv_base_data(expected)),
    identical(saved$config$chains, config$chains), identical(saved$config$warmup, config$warmup),
    identical(saved$config$retained_per_chain, config$retained_per_chain))
  for (path in saved$native_fit$output_files()) verify(path, original_seals)
  if (cell == "history_reference") {
    stopifnot(identical(saved$kernel_sha256, hash(production)))
  } else {
    expected$asset_variance_active <- 1L
    if (cell == "history_trend") {
      expected$variance_time_trend_active <- 1L
      expected$variance_time01 <- (as.integer(input$rows$year) - 1960) / 5
      expected$variance_time_prior_sd <- 1
    }
    stopifnot(identical(saved$data, expected), isTRUE(saved$candidate))
    generated <- tempfile(fileext = ".stan"); on.exit(unlink(generated), add = TRUE)
    if (cell == "history_trend") yh_trend_kernel(production, generated) else yh_base_kernel(production, generated)
    stopifnot(identical(hash(generated), saved$kernel_sha256))
  }
  files <- file.path(previous, cell, c("policy_draws.rds", "policy_validity.csv",
    "comparison.csv", "sampler_diagnostics.csv", "parameter_summary.csv"))
  for (path in files) verify(path)
  policy <- readRDS(files[[1L]])
  stopifnot(length(policy$indices) == config$policy_draws)
  list(fit_path = normalizePath(fit_path), fit_sha256 = hash(fit_path), files = files,
    copied_policy_seed = config$policy_seed, policy_indices = policy$indices)
}
