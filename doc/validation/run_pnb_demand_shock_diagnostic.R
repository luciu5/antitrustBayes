# Plan by default; no production controller, inference cache or selection is resumed.
args <- commandArgs(TRUE)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
root <- normalizePath(file.path(dirname(script), "../.."), mustWork = TRUE)
option <- function(key, default) {
  value <- args[startsWith(args, paste0("--", key, "="))]
  if (length(value) > 1L) stop("Repeated option: ", key)
  if (!length(value)) default else substring(value, nchar(key) + 4L)
}
if (any(!args %in% c("--plan", "--execute") & !grepl("^--(output|paper-root|worker)=", args)))
  stop("Use --plan or --execute [--output=PATH] [--paper-root=PATH].")
if (all(c("--plan", "--execute") %in% args)) stop("Choose either --plan or --execute.")
config <- source(file.path(root, "doc/validation/pnb_demand_shock_config.R"))$value
paper <- normalizePath(option("paper-root", file.path(root, "../mergerBayes")), mustWork = TRUE)
output <- option("output", file.path(root, ".scratch/demand_shock_diagnostic_20261109"))
cell <- option("worker", "")
plan <- data.frame(Cell = c("snapshot_reference", "snapshot_candidate", "history_reference", "history_candidate"),
  Information = rep(c("snapshot", "history"), each = 2L),
  Arm = rep(c("reference", "candidate"), 2L), Status = "not_run",
  Exit = NA_integer_, Error = NA_character_, stringsAsFactors = FALSE)
if (!"--execute" %in% args && !nzchar(cell)) {
  print(plan); print(config)
  cat("At most four fits; reuse compatible saved references. No compilation, sampling, bridges, BMA or writes.\n")
  quit(save = "no")
}
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1", STAN_NUM_THREADS = "1", MAKEFLAGS = "-j1")
.libPaths(c(file.path(paper, ".scratch/paper-library"), .libPaths()))
suppressPackageStartupMessages(library(dplyr))
pkgload::load_all(root, quiet = TRUE)
source(file.path(root, "doc/validation/demand_shock_prototype.R"))
clock_path <- file.path(output, "clock.json")
remaining <- function() {
  clock <- jsonlite::read_json(clock_path)
  stamp <- sub("\\+00:00$", "Z", clock$deadline_utc)
  as.numeric(as.POSIXct(stamp, format = "%Y-%m-%dT%H:%M:%OSZ", tz = "UTC")) - as.numeric(Sys.time())
}
diagnose <- function(fit, candidate) {
  variables <- c("a_market", "s0_prob", "sigma_b_strat_economic", "sigma_share", "sigma_margin", "rho")
  if (candidate) variables <- c(variables, "sigma_demand_shock", "rho_demand_shock")
  s <- fit$summary(variables = variables)
  diagnostic <- fit$diagnostic_summary()
  shape <- dim(fit$draws(variables = "a_market", format = "draws_array"))
  data.frame(Chains = shape[[2L]], Retained_Draws = shape[[1L]] * shape[[2L]],
    Max_Rhat = max(s$rhat, na.rm = TRUE), Min_Bulk_ESS = min(s$ess_bulk, na.rm = TRUE),
    Min_Tail_ESS = min(s$ess_tail, na.rm = TRUE),
    Divergences = sum(diagnostic$num_divergent), Treedepth_Hits = sum(diagnostic$num_max_treedepth))
}
if (nzchar(cell)) {
  stopifnot(cell %in% plan$Cell, remaining() > 10)
  manifest <- readRDS(file.path(output, "manifest.rds"))
  stopifnot(identical(config, manifest$config), all(vapply(names(manifest$source_hash),
    function(path) digest::digest(file = path, algo = "sha256") == manifest$source_hash[[path]], logical(1))))
  input <- readRDS(file.path(output, "inputs.rds"))[[sub("_.*", "", cell)]]
  candidate <- grepl("candidate$", cell)
  destination <- file.path(output, cell); dir.create(destination, showWarnings = FALSE)
  if (file.exists(file.path(destination, "fit.rds"))) stop("Existing fit is preserved; this worker never refits a completed cell.")
  kernel <- file.path(output, "build", if (candidate) "period_demand_shock.stan" else "reference_logit_supply.stan")
  model <- cmdstanr::cmdstan_model(kernel, quiet = TRUE, cpp_options = list(stan_threads = FALSE))
  stopifnot(remaining() > 10)
  fit <- model$sample(data = input$data, chains = 4L, parallel_chains = 4L, threads_per_chain = 1L,
    iter_warmup = 500L, iter_sampling = 500L, seed = 10710L + match(cell, plan$Cell),
    adapt_delta = .999, max_treedepth = 13L, metric = "diag_e", init = 0,
    sig_figs = 18L, refresh = 100L, output_dir = destination)
  saveRDS(list(native_fit = fit, data = input$data, prototype = candidate,
    profile = config$profile, kernel_sha256 = digest::digest(file = kernel, algo = "sha256"),
    data_sha256 = digest::digest(input$data, algo = "sha256"), config = config), file.path(destination, "fit.rds"))
  write.csv(diagnose(fit, candidate), file.path(destination, "sampler_diagnostics.csv"), row.names = FALSE)
  ds_evaluate(fit, input, candidate, destination, seed = 10702L)
  quit(save = "no")
}
dir.create(output, recursive = TRUE, showWarnings = FALSE)
output <- normalizePath(output, mustWork = TRUE)
lock <- paste0(output, ".lock")
if (!dir.create(lock, showWarnings = FALSE)) stop("Existing diagnostic controller lock is preserved: ", lock)
run_controller <- function() {
  on.exit(file.remove(lock), add = TRUE)
  if (!file.exists(clock_path)) jsonlite::write_json(list(
    started_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    deadline_utc = format(Sys.time() + config$budget_seconds, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    budget_seconds = config$budget_seconds), clock_path, auto_unbox = TRUE, pretty = TRUE)
  if (file.exists(file.path(output, "manifest.rds"))) stop("Existing execution is preserved; choose a distinct output.")
  stopifnot(remaining() > 15)
  seal <- file.path(paper, "results/banking/share_residual_pilot_20261108",
    "banking_operating_v3_robustness_share_residual_reference_pnb_1960_20261108",
    "pnb_1960/seed_19091/pilot_artifacts_b.csv")
  saved <- ds_saved_fit(seal); ds_alias_check(seal, output)
  input <- ds_inputs(paper, saved$fit); saveRDS(input, file.path(output, "inputs.rds"))
  dir.create(file.path(output, "build"), showWarnings = FALSE)
  production <- file.path(root, "inst/stan/logit_supply.stan")
  stopifnot(file.copy(production, file.path(output, "build/reference_logit_supply.stan"), overwrite = FALSE))
  ds_kernel(production, file.path(output, "build/period_demand_shock.stan"))
  sources <- c(production, file.path(root, "doc/validation",
    c("demand_shock_prototype.R", "run_pnb_demand_shock_diagnostic.R", "pnb_demand_shock_config.R",
      "run_pnb_demand_shock_diagnostic.py")),
    file.path(root, "R", c("params.R", "latent_quality_simulation.R", "relative_shares.R")),
    system.file("R/antitrust.rdb", package = "antitrust"),
    file.path(paper, "data", c("analysis/supreme_data_operating_v3.RData", "raw/event_overlaps.RData")))
  manifest <- list(profile = config$profile, config = config,
    source_hash = setNames(vapply(sources, function(p) digest::digest(file = p, algo = "sha256"), character(1)), sources),
    packages = setNames(lapply(c("antitrustBayes", "antitrust", "coordination", "cmdstanr"),
      function(p) list(version = as.character(packageVersion(p)), path = find.package(p))),
      c("antitrustBayes", "antitrust", "coordination", "cmdstanr")),
    production_kernel_unchanged = digest::digest(file = production, algo = "sha256"),
    saved_reference_path = saved$record$Path, saved_reference_sha256 = saved$record$SHA256,
    new_hmc_attempts = 0L,
    measurement_assumption = "Independent constant .10 log-scale measurement; shock-margin rho retains existing prior.",
    inference_provenance = "Separate diagnostic kernel; candidate is not an unmodified specification-14 BayesFit, and has no bridge/BMA eligibility.")
  saveRDS(manifest, file.path(output, "manifest.rds"))
  comparison <- list()
  collect <- function(j) {
    destination <- file.path(output, plan$Cell[j]); path <- file.path(destination, "comparison.csv")
    if (file.exists(path)) {
      x <- read.csv(path, stringsAsFactors = FALSE)
      x$Information <- plan$Information[j]; x$Arm <- plan$Arm[j]
      x$Execution_Status <- plan$Status[j]
      diagnostic <- read.csv(file.path(destination, "sampler_diagnostics.csv"))
      x$Max_Rhat <- diagnostic$Max_Rhat; x$Min_Bulk_ESS <- diagnostic$Min_Bulk_ESS
      x$Divergences <- diagnostic$Divergences; x$Treedepth_Hits <- diagnostic$Treedepth_Hits
      x$Convergence <- if (x$Max_Rhat[[1L]] <= 1.01 && x$Min_Bulk_ESS[[1L]] >= 400 &&
        diagnostic$Divergences == 0 && diagnostic$Treedepth_Hits == 0) "adequate_for_bounded_pilot" else "inconclusive_convergence"
      comparison[[plan$Cell[j]]] <<- x
    } else if (length(comparison)) {
      x <- comparison[[1L]]
      numerical <- vapply(x, is.numeric, logical(1))
      numerical[names(numerical) == "Observed_Share_PP"] <- FALSE
      x[numerical] <- NA_real_
      x$Information <- plan$Information[j]; x$Arm <- plan$Arm[j]
      x$Execution_Status <- plan$Status[j]; x$Convergence <- "not_available"
      comparison[[plan$Cell[j]]] <<- x
    }
    if (length(comparison)) write.csv(do.call(rbind, comparison), file.path(output, "comparison.csv"), row.names = FALSE)
    write.csv(plan, file.path(output, "execution_status.csv"), row.names = FALSE)
  }
  destination <- file.path(output, "snapshot_reference"); dir.create(destination, showWarnings = FALSE)
  ds_evaluate(saved$fit@stanfit, input$snapshot, FALSE, destination)
  write.csv(diagnose(saved$fit@stanfit, FALSE), file.path(destination, "sampler_diagnostics.csv"), row.names = FALSE)
  plan$Status[1L] <- "reused_compatible_saved_fit_2_chains_250_retained"
  plan$Exit[1L] <- 0L; collect(1L)
  for (j in 2:4) {
    if (j > 2L && plan$Status[2L] != "completed") {
      plan$Status[j] <- "not_run_snapshot_incomplete"; collect(j); next
    }
    if (remaining() < 120) { plan$Status[j] <- "not_run_budget"; collect(j); next }
    plan$Status[j] <- "running"; collect(j)
    manifest$new_hmc_attempts <- manifest$new_hmc_attempts + 1L
    saveRDS(manifest, file.path(output, "manifest.rds"))
    result <- processx::run(Sys.which("timeout"), c("--signal=TERM", "--kill-after=5s",
      paste0(floor(remaining()) - 10L, "s"), file.path(R.home("bin"), "Rscript"), "--vanilla",
      normalizePath(script), paste0("--worker=", plan$Cell[j]), paste0("--output=", output), paste0("--paper-root=", paper)),
      env = "current", stdout = file.path(output, paste0(plan$Cell[j], ".log")),
      stderr_to_stdout = TRUE, cleanup_tree = TRUE, error_on_status = FALSE)
    plan$Exit[j] <- result$status
    plan$Status[j] <- if (result$status == 0) "completed" else if (result$status == 124) "incomplete_time_cap" else "failed"
    if (result$status != 0) plan$Error[j] <- paste("See preserved worker log; exit", result$status)
    collect(j)
  }
  stopifnot(identical(digest::digest(file = production, algo = "sha256"), manifest$production_kernel_unchanged),
    manifest$new_hmc_attempts <= 3L)
  cat("Finished. All planned cells/invalid draws retained; no bridges or BMA.\n")
  if (any(plan$Status == "failed")) quit(save = "no", status = 1L)
  if (any(plan$Status %in% c("not_run_budget", "incomplete_time_cap")))
    quit(save = "no", status = 124L)
}
run_controller()
