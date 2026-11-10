# Two new fits maximum. Reuse sealed reference/fixed-.10 comparisons unchanged.
args <- commandArgs(TRUE)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
root <- normalizePath(file.path(dirname(script), "../.."), mustWork = TRUE)
option <- function(key, default) {
  x <- args[startsWith(args, paste0("--", key, "="))]
  if (length(x) > 1L) stop("Repeated option: ", key)
  if (!length(x)) default else substring(x, nchar(key) + 4L)
}
if (any(!args %in% c("--plan", "--execute") & !grepl("^--(output|paper-root|worker)=", args)))
  stop("Use --plan or --execute [--output=PATH] [--paper-root=PATH].")
if (all(c("--plan", "--execute") %in% args)) stop("Choose one mode.")
config <- source(file.path(root, "doc/validation/pnb_measurement_prior_config.R"))$value
paper <- normalizePath(option("paper-root", file.path(root, "../mergerBayes")), mustWork = TRUE)
output <- option("output", file.path(root, ".scratch/measurement_prior_diagnostic_20261109"))
worker <- option("worker", "")
plan <- expand.grid(Arm = c("reference", "fixed", "observed", "estimated"),
  Information = c("snapshot", "history"), stringsAsFactors = FALSE)
plan$Cell <- paste(plan$Information, plan$Arm, sep = "_")
plan$Status <- "not_run"; plan$Exit <- NA_integer_; plan$Error <- NA_character_
if (!"--execute" %in% args && !nzchar(worker)) {
  print(plan); print(config)
  cat("Two new jointly estimated measurement-SD fits; four sealed comparisons reused; two observed-anchor policy replays. No writes, compilation, sampling, bridges or BMA in plan mode.\n")
  quit(save = "no")
}
if (!identical(Sys.getenv("IOVRSE_MEASUREMENT_PRIOR_WRAPPED"), "1"))
  stop("Execute through run_pnb_measurement_prior_diagnostic.py for the absolute cap.")
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1", STAN_NUM_THREADS = "1", MAKEFLAGS = "-j1")
.libPaths(c(file.path(paper, ".scratch/paper-library"), .libPaths()))
suppressPackageStartupMessages(library(dplyr))
pkgload::load_all(root, quiet = TRUE)
for (module in c("demand_shock", "asset_variance", "inverse_asset_variance", "measurement_prior"))
  source(file.path(root, "doc/validation", paste0(module, "_prototype.R")))
evaluate <- mp_make_evaluator()
observed_evaluate <- mp_make_evaluator(observed = TRUE)
hash <- function(path) digest::digest(file = path, algo = "sha256")
clock_path <- file.path(output, "clock.json")
remaining <- function() {
  clock <- jsonlite::read_json(clock_path)
  as.numeric(as.POSIXct(sub("\\+00:00$", "Z", clock$deadline_utc),
    format = "%Y-%m-%dT%H:%M:%OSZ", tz = "UTC")) - as.numeric(Sys.time())
}
if (nzchar(worker)) {
  stopifnot(worker %in% c("snapshot_estimated", "history_estimated"), remaining() > 15)
  manifest <- readRDS(file.path(output, "manifest.rds"))
  stopifnot(identical(config, manifest$config), all(vapply(names(manifest$source_hash),
    function(path) identical(hash(path), manifest$source_hash[[path]]), logical(1))))
  input <- readRDS(file.path(output, "inputs.rds"))[[sub("_.*", "", worker)]]
  ds_validate_data(input$data)
  destination <- file.path(output, worker); dir.create(destination, showWarnings = FALSE)
  if (file.exists(file.path(destination, "fit.rds"))) stop("Completed fit is preserved.")
  kernel <- file.path(output, "build/estimated_measurement_logit_supply.stan")
  model <- cmdstanr::cmdstan_model(kernel, quiet = TRUE, cpp_options = list(stan_threads = FALSE))
  stopifnot(remaining() > 15)
  fit <- model$sample(data = input$data, chains = config$chains,
    parallel_chains = config$parallel_chains, threads_per_chain = 1L,
    iter_warmup = config$warmup, iter_sampling = config$retained_per_chain,
    seed = config$seed + match(worker, plan$Cell), adapt_delta = config$adapt_delta,
    max_treedepth = config$max_treedepth, metric = "diag_e", init = 0,
    sig_figs = 18L, refresh = 100L, output_dir = destination)
  saveRDS(list(native_fit = fit, data = input$data, profile = config$profile,
    candidate = TRUE, kernel_sha256 = hash(kernel), data_sha256 = digest::digest(input$data, algo = "sha256"),
    config = config), file.path(destination, "fit.rds"))
  # Inspect every varying native parameter/transformed quantity, not just m.
  diagnostic <- av_diagnose(fit, FALSE)
  energy <- fit$diagnostic_summary()$ebfmi
  diagnostic$Min_EBFMI <- if (length(energy)) min(energy) else NA_real_
  write.csv(diagnostic, file.path(destination, "sampler_diagnostics.csv"), row.names = FALSE)
  write.csv(fit$summary(variables = c("log_sigma_measurement", "diagnostic_measurement_sd",
    "sigma_demand_shock", "rho_demand_shock", "sigma_share")),
    file.path(destination, "variance_split_sampler_diagnostics.csv"), row.names = FALSE)
  mp_parameters(fit, input, destination, config)
  evaluate(fit, input, TRUE, destination, n_draws = config$policy_draws, seed = config$policy_seed)
  quit(save = "no")
}
dir.create(output, recursive = TRUE, showWarnings = FALSE)
output <- normalizePath(output, mustWork = TRUE)
if (file.exists(file.path(output, "manifest.rds"))) stop("Existing execution preserved; use a distinct output.")
stopifnot(file.exists(clock_path), remaining() > 15)
production <- file.path(root, "inst/stan/logit_supply.stan")
seal <- file.path(paper, "results/banking/share_residual_pilot_20261108",
  "banking_operating_v3_robustness_share_residual_reference_pnb_1960_20261108",
  "pnb_1960/seed_19091/pilot_artifacts_b.csv")
saved <- ds_saved_fit(seal)
input <- mp_inputs(paper, saved$fit, config)
saveRDS(input, file.path(output, "inputs.rds"))
dir.create(file.path(output, "build"), showWarnings = FALSE)
mp_kernel(production, file.path(output, "build/estimated_measurement_logit_supply.stan"))
writeLines(deparse(evaluate), file.path(output, "executed_policy_evaluator.R"))
writeLines(deparse(observed_evaluate), file.path(output, "executed_observed_evaluator.R"))
sources <- c(production, file.path(root, "doc/validation", c("demand_shock_prototype.R",
  "asset_variance_prototype.R", "inverse_asset_variance_prototype.R", "measurement_prior_prototype.R",
  "measurement_prior_prototype_checks.R", "pnb_measurement_prior_config.R",
  "run_pnb_measurement_prior_diagnostic.R", "run_pnb_measurement_prior_diagnostic.py")),
  file.path(root, "R", c("params.R", "latent_quality_simulation.R", "relative_shares.R")),
  system.file("R/antitrust.rdb", package = "antitrust"),
  file.path(paper, "data", c("analysis/supreme_data_operating_v3.RData", "raw/event_overlaps.RData")),
  file.path(paper, "code/banking", c("banking_ladder_common.R", "supreme_antitrustbayes.R")))
manifest <- list(profile = config$profile, config = config,
  source_hash = setNames(vapply(sources, hash, character(1)), sources),
  package_source_revision = system2("git", c("-C", shQuote(root), "rev-parse", "HEAD"), stdout = TRUE),
  dirty_source = TRUE, source_is_development_checkout = TRUE,
  packages = setNames(lapply(c("antitrustBayes", "antitrust", "coordination", "cmdstanr"),
    function(p) list(version = as.character(packageVersion(p)), path = find.package(p))),
    c("antitrustBayes", "antitrust", "coordination", "cmdstanr")),
  cmdstan_version = as.character(cmdstanr::cmdstan_version()), cmdstan_path = cmdstanr::cmdstan_path(),
  production_kernel_sha256 = hash(production), reference_template_path = saved$record$Path,
  reference_template_sha256 = saved$record$SHA256, new_hmc_attempts = 0L,
  inference_provenance = "Separate jointly fitted measurement-prior diagnostic. One added primitive log measurement SD; independent proper lognormal prior; persistent quality and demand-shock priors and total_fraction Jacobian unchanged from fixed-.10 prototype. No production-fit/evidence/BMA eligibility. Full Bayesian posterior, not a saved-fit split.",
  identification = "Likelihood identifies total residual variance and shock-margin covariance, not the separate transitory variances with free rho_shock. Prior-to-posterior changes remain conditional on both variance priors and correlation support.")
saveRDS(manifest, file.path(output, "manifest.rds"))
write.csv(plan, file.path(output, "execution_status.csv"), row.names = FALSE)
collect <- function() {
  pieces <- lapply(seq_len(nrow(plan)), function(j) {
    folder <- file.path(output, plan$Cell[j]); path <- file.path(folder, "comparison.csv")
    if (!file.exists(path)) return(NULL)
    x <- read.csv(path); d <- read.csv(file.path(folder, "sampler_diagnostics.csv"))
    x$Information <- plan$Information[j]; x$Arm <- plan$Arm[j]; x$Execution_Status <- plan$Status[j]
    for (n in c("Max_Rhat", "Min_Bulk_ESS", "Min_Tail_ESS", "Divergences", "Treedepth_Hits")) x[[n]] <- d[[n]]
    x$Convergence <- if (d$Max_Rhat <= 1.01 && d$Min_Bulk_ESS >= 400 &&
      d$Divergences == 0 && d$Treedepth_Hits == 0) "adequate_for_bounded_pilot" else "inconclusive_convergence"
    x
  })
  if (any(lengths(pieces) > 0)) write.csv(do.call(rbind, pieces), file.path(output, "comparison.csv"), row.names = FALSE)
  write.csv(plan, file.path(output, "execution_status.csv"), row.names = FALSE)
}
# Reuse stronger 4x500 original references, including byte-identical policy states.
references <- list()
previous <- file.path(root, ".scratch/asset_variance_diagnostic_20261109")
for (information in c("snapshot", "history")) {
  cell <- paste0(information, "_reference")
  reference <- iv_reference(previous, cell,
    list(data = mp_base_data(input[[information]]$data)), production, config)
  references[[cell]] <- reference[c("fit_path", "fit_sha256", "copied_policy_seed", "policy_indices")]
  destination <- file.path(output, cell); dir.create(destination, showWarnings = FALSE)
  stopifnot(all(file.copy(reference$files, destination, overwrite = FALSE)))
}
# Fixed-.10 prototype fits and all 100 original policy states retain their seals.
fixed_root <- file.path(root, ".scratch/demand_shock_diagnostic_20261109")
fixed_seals <- read.csv(file.path(fixed_root, "completed_artifact_seals.csv"))
verify_fixed <- function(path) {
  record <- fixed_seals[fixed_seals$Path == normalizePath(path, mustWork = TRUE), , drop = FALSE]
  stopifnot(nrow(record) == 1L, identical(hash(path), record$SHA256[[1L]]))
}
fixed_kernel <- file.path(fixed_root, "build/period_demand_shock.stan")
verify_fixed(fixed_kernel)
generated <- tempfile(fileext = ".stan"); ds_kernel(production, generated)
stopifnot(identical(hash(generated), hash(fixed_kernel))); unlink(generated)
for (information in c("snapshot", "history")) {
  source_cell <- paste0(information, "_candidate"); cell <- paste0(information, "_fixed")
  path <- file.path(fixed_root, source_cell, "fit.rds"); verify_fixed(path)
  fit <- readRDS(path)
  stopifnot(identical(fit$data, mp_base_data(input[[information]]$data)),
    isTRUE(fit$prototype), identical(fit$kernel_sha256, hash(fixed_kernel)),
    fit$config$chains == config$chains, fit$config$retained_per_chain == config$retained_per_chain,
    fit$config$warmup == config$warmup, identical(fit$config$measurement_sd, .1))
  for (csv in fit$native_fit$output_files()) verify_fixed(csv)
  files <- file.path(fixed_root, source_cell, c("policy_draws.rds", "policy_validity.csv", "comparison.csv", "sampler_diagnostics.csv"))
  for (p in files) verify_fixed(p)
  destination <- file.path(output, cell); dir.create(destination, showWarnings = FALSE)
  stopifnot(all(file.copy(files, destination, overwrite = FALSE)))
  references[[cell]] <- list(fit_path = normalizePath(path), fit_sha256 = hash(path),
    copied_policy_seed = 10702L, policy_indices = readRDS(files[[1L]])$indices)
}
manifest$reused_references <- references
manifest$reused_reference_paths <- vapply(references, `[[`, character(1), "fit_path")
saveRDS(manifest, file.path(output, "manifest.rds"))
for (j in seq_len(nrow(plan))) {
  cell <- plan$Cell[j]; information <- plan$Information[j]
  if (cell %in% names(references)) {
    plan$Status[j] <- "reused_sealed_fit_and_identical_policy_draws"; plan$Exit[j] <- 0L
    collect(); next
  }
  if (remaining() < 120) { plan$Status[j] <- "not_run_budget"; collect(); next }
  if (information == "history" && plan$Status[plan$Cell == "snapshot_estimated"] != "completed") {
    plan$Status[j] <- "not_run_snapshot_incomplete"; collect(); next
  }
  if (information == "history") {
    d <- read.csv(file.path(output, "snapshot_estimated/sampler_diagnostics.csv"))
    if (d$Divergences > 0 || d$Treedepth_Hits > 0 || d$Max_Rhat > 1.01 || d$Min_Bulk_ESS < 400) {
      plan$Status[j] <- "not_run_snapshot_sampling_inconclusive"; collect(); next
    }
  }
  plan$Status[j] <- "running"; collect()
  if (plan$Arm[j] == "observed") {
    path <- references[[paste0(information, "_reference")]]$fit_path
    destination <- file.path(output, cell); dir.create(destination, showWarnings = FALSE)
    observed_evaluate(readRDS(path)$native_fit, input[[information]], FALSE,
      destination, n_draws = config$policy_draws, seed = config$policy_seed)
    stopifnot(file.copy(file.path(output, paste0(information, "_reference"), "sampler_diagnostics.csv"), destination))
    plan$Status[j] <- "completed_observed_anchor_replay"; plan$Exit[j] <- 0L
    collect(); next
  }
  manifest$new_hmc_attempts <- manifest$new_hmc_attempts + 1L
  stopifnot(manifest$new_hmc_attempts <= 2L)
  saveRDS(manifest, file.path(output, "manifest.rds"))
  result <- processx::run(Sys.which("timeout"), c("--signal=TERM", "--kill-after=5s",
    paste0(floor(remaining()) - 10L, "s"), file.path(R.home("bin"), "Rscript"), "--vanilla", normalizePath(script),
    paste0("--worker=", cell), paste0("--output=", output), paste0("--paper-root=", paper)),
    env = "current", stdout = file.path(output, paste0(cell, ".log")), stderr_to_stdout = TRUE,
    cleanup_tree = TRUE, error_on_status = FALSE)
  plan$Exit[j] <- result$status
  plan$Status[j] <- if (result$status == 0) "completed" else if (result$status == 124) "incomplete_time_cap" else "failed"
  if (result$status != 0) plan$Error[j] <- paste("Preserved worker log; exit", result$status)
  collect()
}
stopifnot(identical(hash(production), manifest$production_kernel_sha256),
  identical(hash(saved$record$Path), manifest$reference_template_sha256))
if (any(plan$Status == "failed")) quit(save = "no", status = 1L)
if (any(grepl("not_run|incomplete", plan$Status))) quit(save = "no", status = 124L)
cat("Measurement-prior diagnostic complete. No bridges, BMA, production default changes or pushes.\n")
