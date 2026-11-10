# Source the unchanged demand_shock_prototype.R and asset_variance_prototype.R
# first. All modifications below are separate research-only functions.
iv_base_kernel <- av_kernel
iv_base_inputs <- av_inputs
iv_base_parameters <- av_parameters

iv_kernel <- function(source, destination) {
  stage <- tempfile(fileext = ".stan")
  on.exit(unlink(stage), add = TRUE)
  iv_base_kernel(source, stage)
  text <- paste(readLines(stage, warn = FALSE), collapse = "\n")
  once <- function(old, new) {
    stopifnot(length(gregexpr(old, text, fixed = TRUE)[[1L]]) == 1L,
      grepl(old, text, fixed = TRUE))
    text <<- sub(old, new, text, fixed = TRUE)
  }
  once("  real<lower=0> asset_gamma_prior_sd;", "  real<lower=0> institution_log_assets_log_sd;")
  once("  vector[asset_variance_active] asset_gamma;", "")
  once("transformed parameters {", paste("transformed parameters {",
    "  real asset_gamma = asset_variance_active == 1 ? -0.5 * institution_log_assets_log_sd : 0.0;", sep = "\n"))
  once("asset_gamma[1] * institution_log_assets_z[i]", "asset_gamma * institution_log_assets_z[i]")
  once("  target += normal_lpdf(asset_gamma | 0, asset_gamma_prior_sd);", "")
  once("// Separate diagnostic: log sigma_i = log sigma + gamma z(log institution assets).",
    "// Fixed inverse asset VARIANCE: sigma_i^2 = sigma_base^2 * A_reference/A_i; no additional sampler parameter.")
  writeLines(text, destination)
  invisible(destination)
}

iv_inputs <- function(paper, template, config) {
  stopifnot(identical(config$fixed_variance_exponent, 1))
  prepared <- iv_base_inputs(paper, template, config)
  for (name in names(prepared$inputs)) {
    prepared$inputs[[name]]$data$asset_gamma_prior_sd <- NULL
    prepared$inputs[[name]]$data$institution_log_assets_log_sd <- prepared$assets$scale
    z <- prepared$inputs[[name]]$rows$institution_log_assets_z
    assets <- prepared$inputs[[name]]$rows$total_assets
    stopifnot(max(abs(exp(-prepared$assets$scale * z) - exp(prepared$assets$center) / assets)) < 1e-12)
  }
  prepared
}

iv_parameters <- function(fit, input, candidate, output) {
  # Keep scalar base-prior parameter summaries separate from the fixed exponent.
  iv_base_parameters(fit, input, FALSE, output)
  if (candidate) {
    draw <- ds_private(".bayes_backend_extract")(fit, variables = "sigma_share_observation")
    target <- which(input$rows$year == "1960")
    size <- input$rows[target, c("group", "bank", "total_assets", "institution_log_assets_z")]
    size$Mean_Residual_SD <- colMeans(draw$sigma_share_observation[, target, drop = FALSE])
    write.csv(size, file.path(output, "asset_residual_scales.csv"), row.names = FALSE)
    write.csv(data.frame(Fixed_Variance_Exponent = 1,
      Fixed_Standardized_Log_SD_Slope = -.5 * input$data$institution_log_assets_log_sd),
      file.path(output, "fixed_asset_rule.csv"), row.names = FALSE)
  }
}

iv_base_data <- function(data) {
  data[c("asset_variance_active", "institution_log_assets_z", "asset_gamma_prior_sd",
    "institution_log_assets_log_sd")] <- NULL
  data
}

iv_reference <- function(previous, cell, input, production, config) {
  stopifnot(cell %in% c("snapshot_reference", "history_reference"))
  hash <- function(path) digest::digest(file = path, algo = "sha256")
  seals <- read.csv(file.path(previous, "artifact_seals.csv"), stringsAsFactors = FALSE)
  verify <- function(path) {
    record <- seals[seals$Path == normalizePath(path, mustWork = TRUE), , drop = FALSE]
    stopifnot(nrow(record) == 1L, identical(hash(path), record$SHA256[[1L]]))
    invisible(path)
  }
  verify(file.path(previous, "manifest.rds"))
  previous_manifest <- readRDS(file.path(previous, "manifest.rds"))
  stopifnot(identical(previous_manifest$production_kernel_sha256, hash(production)),
    all(vapply(names(previous_manifest$source_hash), function(path)
      identical(hash(path), previous_manifest$source_hash[[path]]), logical(1))))
  if (cell == "snapshot_reference") {
    fit_path <- file.path(previous, cell, "fit.rds")
  } else {
    fit_path <- previous_manifest$reused_history_path
    # History lives in the preceding experiment and retains its original seal.
    original_seals <- read.csv(file.path(dirname(dirname(fit_path)), "completed_artifact_seals.csv"))
    record <- original_seals[original_seals$Path == fit_path, , drop = FALSE]
    stopifnot(nrow(record) == 1L, identical(hash(fit_path), record$SHA256[[1L]]))
  }
  if (cell == "snapshot_reference") verify(fit_path)
  saved <- readRDS(fit_path)
  stopifnot(identical(iv_base_data(saved$data), iv_base_data(input$data)),
    identical(saved$config$chains, config$chains),
    identical(saved$config$warmup, config$warmup),
    identical(saved$config$retained_per_chain, config$retained_per_chain))
  if (cell == "snapshot_reference") {
    stopifnot(identical(saved$candidate, FALSE), saved$data$asset_variance_active == 0L)
    kernel <- file.path(previous, "build/asset_variance_logit_supply.stan")
    verify(kernel)
    stopifnot(identical(hash(kernel), saved$kernel_sha256))
    generated <- tempfile(fileext = ".stan"); on.exit(unlink(generated), add = TRUE)
    iv_base_kernel(production, generated)
    stopifnot(identical(hash(generated), saved$kernel_sha256))
    for (path in saved$native_fit$output_files()) verify(path)
  } else {
    stopifnot(identical(saved$kernel_sha256, hash(production)))
    for (path in saved$native_fit$output_files()) {
      record <- original_seals[original_seals$Path == path, , drop = FALSE]
      stopifnot(nrow(record) == 1L, identical(hash(path), record$SHA256[[1L]]))
    }
  }
  files <- file.path(previous, cell, c("policy_draws.rds", "policy_validity.csv",
    "comparison.csv", "sampler_diagnostics.csv", "parameter_summary.csv"))
  for (path in files) verify(path)
  policy <- readRDS(files[[1L]])
  stopifnot(length(policy$indices) == config$policy_draws)
  list(fit_path = normalizePath(fit_path), fit_sha256 = hash(fit_path),
    files = files, native_fit = saved$native_fit,
    copied_policy_seed = 10902L, policy_indices = policy$indices)
}
