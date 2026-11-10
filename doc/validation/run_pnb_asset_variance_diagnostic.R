# Bounded institution-assets variance test; plan-only unless --execute.
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
if (all(c("--plan", "--execute") %in% args)) stop("Choose either --plan or --execute.")
config <- source(file.path(root, "doc/validation/pnb_asset_variance_config.R"))$value
paper <- normalizePath(option("paper-root", file.path(root, "../mergerBayes")), mustWork = TRUE)
output <- option("output", file.path(root, ".scratch/asset_variance_diagnostic_20261109"))
worker <- option("worker", "")
plan <- data.frame(Cell = c("snapshot_reference", "snapshot_candidate", "history_reference", "history_candidate"),
  Information = rep(c("snapshot", "history"), each = 2L), Arm = rep(c("reference", "candidate"), 2L),
  Status = "not_run", Exit = NA_integer_, Error = NA_character_, stringsAsFactors = FALSE)
if (!"--execute" %in% args && !nzchar(worker)) {
  print(plan); print(config)
  cat("Institution totals and global bank-year asset scaling. At most four fits; reuse compatible history. No writes or sampling.\n")
  quit(save = "no")
}
if (!identical(Sys.getenv("IOVRSE_ASSET_PILOT_WRAPPED"), "1"))
  stop("Execute through run_pnb_asset_variance_diagnostic.py, which enforces the total deadline.")
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1", STAN_NUM_THREADS = "1", MAKEFLAGS = "-j1")
.libPaths(c(file.path(paper, ".scratch/paper-library"), .libPaths()))
suppressPackageStartupMessages(library(dplyr))
pkgload::load_all(root, quiet = TRUE)
source(file.path(root, "doc/validation/demand_shock_prototype.R"))
source(file.path(root, "doc/validation/asset_variance_prototype.R"))
evaluate <- av_make_evaluator()
clock_path <- file.path(output, "clock.json")
remaining <- function() {
  clock <- jsonlite::read_json(clock_path)
  as.numeric(as.POSIXct(sub("\\+00:00$", "Z", clock$deadline_utc),
    format = "%Y-%m-%dT%H:%M:%OSZ", tz = "UTC")) - as.numeric(Sys.time())
}
hash <- function(path) digest::digest(file = path, algo = "sha256")
if (nzchar(worker)) {
  stopifnot(worker %in% plan$Cell, remaining() > 15)
  manifest <- readRDS(file.path(output, "manifest.rds"))
  stopifnot(identical(config, manifest$config), all(vapply(names(manifest$source_hash),
    function(path) hash(path) == manifest$source_hash[[path]], logical(1))))
  input <- readRDS(file.path(output, "inputs.rds"))[[sub("_.*", "", worker)]]
  candidate <- grepl("candidate$", worker)
  input$data$asset_variance_active <- as.integer(candidate)
  destination <- file.path(output, worker); dir.create(destination, showWarnings = FALSE)
  if (file.exists(file.path(destination, "fit.rds"))) stop("Completed fit is preserved.")
  kernel <- file.path(output, "build/asset_variance_logit_supply.stan")
  model <- cmdstanr::cmdstan_model(kernel, quiet = TRUE, cpp_options = list(stan_threads = FALSE))
  stopifnot(remaining() > 15)
  fit <- model$sample(data = input$data, chains = config$chains,
    parallel_chains = config$parallel_chains, threads_per_chain = 1L,
    iter_warmup = config$warmup, iter_sampling = config$retained_per_chain,
    seed = config$seed + match(worker, plan$Cell), adapt_delta = config$adapt_delta,
    max_treedepth = config$max_treedepth, metric = "diag_e", init = 0,
    sig_figs = 18L, refresh = 100L, output_dir = destination)
  saveRDS(list(native_fit = fit, data = input$data, profile = config$profile,
    candidate = candidate, kernel_sha256 = hash(kernel), data_sha256 = digest::digest(input$data, algo = "sha256"),
    config = config), file.path(destination, "fit.rds"))
  write.csv(av_diagnose(fit, candidate), file.path(destination, "sampler_diagnostics.csv"), row.names = FALSE)
  av_parameters(fit, input, candidate, destination)
  evaluate(fit, input, candidate, destination, n_draws = config$policy_draws, seed = 10902L)
  quit(save = "no")
}
dir.create(output, recursive = TRUE, showWarnings = FALSE)
output <- normalizePath(output, mustWork = TRUE)
if (file.exists(file.path(output, "manifest.rds"))) stop("Existing execution is preserved; choose a distinct output.")
if (!file.exists(clock_path)) jsonlite::write_json(list(
  started_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  deadline_utc = format(Sys.time() + config$budget_seconds, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  budget_seconds = config$budget_seconds), clock_path, auto_unbox = TRUE, pretty = TRUE)
stopifnot(remaining() > 15)
seal <- file.path(paper, "results/banking/share_residual_pilot_20261108",
  "banking_operating_v3_robustness_share_residual_reference_pnb_1960_20261108",
  "pnb_1960/seed_19091/pilot_artifacts_b.csv")
saved <- ds_saved_fit(seal)
prepared <- av_inputs(paper, saved$fit, config)
input <- prepared$inputs
saveRDS(input, file.path(output, "inputs.rds"))
saveRDS(prepared$assets, file.path(output, "asset_scaling.rds"))
write.csv(prepared$assets$rows, file.path(output, "institution_asset_index.csv"), row.names = FALSE)
dir.create(file.path(output, "build"), showWarnings = FALSE)
production <- file.path(root, "inst/stan/logit_supply.stan")
av_kernel(production, file.path(output, "build/asset_variance_logit_supply.stan"))
writeLines(deparse(evaluate), file.path(output, "executed_policy_evaluator.R"))
sources <- c(production, file.path(root, "doc/validation", c("asset_variance_prototype.R",
  "demand_shock_prototype.R", "pnb_asset_variance_config.R", "run_pnb_asset_variance_diagnostic.R",
  "run_pnb_asset_variance_diagnostic.py")),
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
  asset_scaling = prepared$assets[c("center", "scale", "definition")],
  inference_provenance = "Separate research kernel; no production-fit, bridge or BMA eligibility. Existing quality, total residual and rho priors retained; new gamma ~ N(0,.35).")
saveRDS(manifest, file.path(output, "manifest.rds"))
write.csv(plan, file.path(output, "execution_status.csv"), row.names = FALSE)
collect <- function() {
  pieces <- lapply(seq_len(nrow(plan)), function(j) {
    folder <- file.path(output, plan$Cell[j]); path <- file.path(folder, "comparison.csv")
    if (!file.exists(path)) return(NULL)
    x <- read.csv(path); d <- read.csv(file.path(folder, "sampler_diagnostics.csv"))
    x$Information <- plan$Information[j]; x$Arm <- plan$Arm[j]; x$Execution_Status <- plan$Status[j]
    for (n in names(d)) x[[n]] <- d[[n]]
    x$Convergence <- if (d$Max_Rhat <= 1.01 && d$Min_Bulk_ESS >= 400 &&
      d$Divergences == 0 && d$Treedepth_Hits == 0) "adequate_for_bounded_pilot" else "inconclusive_convergence"
    x
  })
  if (any(lengths(pieces) > 0L)) write.csv(do.call(rbind, pieces), file.path(output, "comparison.csv"), row.names = FALSE)
  write.csv(plan, file.path(output, "execution_status.csv"), row.names = FALSE)
}
# The preceding experiment's history reference is reusable only after exact
# input/kernel equality and the original artifact seals pass; no seals change.
old_root <- file.path(root, ".scratch/demand_shock_diagnostic_20261109")
history_path <- file.path(old_root, "history_reference/fit.rds")
history_reference <- NULL
if (file.exists(history_path)) {
  sealed <- read.csv(file.path(old_root, "completed_artifact_seals.csv"), stringsAsFactors = FALSE)
  path_column <- intersect(c("Path", "path"), names(sealed)); sha_column <- intersect(c("SHA256", "sha256"), names(sealed))
  if (length(path_column) != 1L || length(sha_column) != 1L) stop("Unrecognized existing seal schema.")
  record <- sealed[sealed[[path_column]] == history_path, , drop = FALSE]
  stopifnot(nrow(record) == 1L, identical(hash(history_path), record[[sha_column]][[1L]]))
  old <- readRDS(history_path)
  base <- input$history$data
  base[c("asset_variance_active", "institution_log_assets_z", "asset_gamma_prior_sd")] <- NULL
  stopifnot(identical(old$data, base), identical(old$kernel_sha256, hash(production)))
  history_reference <- old$native_fit
  for (path in history_reference$output_files()) {
    record <- sealed[sealed[[path_column]] == path, , drop = FALSE]
    stopifnot(nrow(record) == 1L, identical(hash(path), record[[sha_column]][[1L]]))
  }
  manifest$reused_history_path <- history_path; manifest$reused_history_sha256 <- hash(history_path)
  saveRDS(manifest, file.path(output, "manifest.rds"))
}
for (j in seq_len(nrow(plan))) {
  if (j > 2L && plan$Status[2L] != "completed") {
    plan$Status[j] <- "not_run_snapshot_incomplete"; collect(); next
  }
  if (remaining() < 120) { plan$Status[j] <- "not_run_budget"; collect(); next }
  if (j == 3L && !is.null(history_reference)) {
    destination <- file.path(output, plan$Cell[j]); dir.create(destination, showWarnings = FALSE)
    evaluate(history_reference, input$history, FALSE, destination, n_draws = config$policy_draws, seed = 10902L)
    av_parameters(history_reference, input$history, FALSE, destination)
    write.csv(av_diagnose(history_reference, FALSE), file.path(destination, "sampler_diagnostics.csv"), row.names = FALSE)
    plan$Status[j] <- "reused_compatible_saved_fit_4_chains_500_retained"; plan$Exit[j] <- 0L
    collect(); next
  }
  plan$Status[j] <- "running"; collect()
  manifest$new_hmc_attempts <- manifest$new_hmc_attempts + 1L
  stopifnot(manifest$new_hmc_attempts <= 4L)
  saveRDS(manifest, file.path(output, "manifest.rds"))
  result <- processx::run(Sys.which("timeout"), c("--signal=TERM", "--kill-after=5s",
    paste0(floor(remaining()) - 10L, "s"), file.path(R.home("bin"), "Rscript"), "--vanilla", normalizePath(script),
    paste0("--worker=", plan$Cell[j]), paste0("--output=", output), paste0("--paper-root=", paper)),
    env = "current", stdout = file.path(output, paste0(plan$Cell[j], ".log")), stderr_to_stdout = TRUE,
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
cat("Four cells collected; production source and saved references unchanged.\n")
