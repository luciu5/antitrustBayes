#!/usr/bin/env Rscript

args0 <- commandArgs(FALSE)
script <- sub("^--file=", "", grep("^--file=", args0, value = TRUE)[[1L]])
source(file.path(dirname(normalizePath(script)), "..", "..", "inst", "benchmarks",
                 "market_intercept_enterprise", "common.R"))

run_dir <- path.expand(mic_env("MIC_ENTERPRISE_RUN_DIR",
  file.path(tempdir(), "market-intercept-enterprise")))
manifest_path <- path.expand(mic_env("MIC_ENTERPRISE_MANIFEST",
                                     file.path(run_dir, "manifest.csv")))
task_value <- mic_env("MIC_ENTERPRISE_TASK_ID", mic_env("SLURM_ARRAY_TASK_ID"))
task_id <- suppressWarnings(as.integer(task_value))
if (length(task_id) != 1L || is.na(task_id) || task_id < 1L)
  stop("Set MIC_ENTERPRISE_TASK_ID or SLURM_ARRAY_TASK_ID to a positive integer.",
       call. = FALSE)
if (!file.exists(manifest_path)) stop("Missing manifest: ", manifest_path, call. = FALSE)

manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE,
                            na.strings = c("", "NA"))
mic_validate_manifest(manifest)
task <- manifest[manifest$task_id == task_id, , drop = FALSE]
if (nrow(task) != 1L) stop("Task id is absent or duplicated: ", task_id, call. = FALSE)
manifest_hash <- unique(manifest$manifest_hash)
if (length(manifest_hash) != 1L || !nzchar(manifest_hash))
  stop("Manifest hash is absent or inconsistent.", call. = FALSE)
if (!identical(as.character(mic_manifest_hash(manifest)), as.character(manifest_hash)))
  stop("Manifest content does not match its frozen hash.", call. = FALSE)

mode <- tolower(mic_env("MIC_ENTERPRISE_MODE", "smoke"))
if (!mode %in% c("smoke", "run")) stop("Mode must be smoke or run.", call. = FALSE)
task_dir <- file.path(run_dir, "tasks", sprintf("task-%04d", task_id))
dir.create(task_dir, recursive = TRUE, showWarnings = FALSE)
attempts <- list.files(task_dir, pattern = "^attempt-[0-9]{3}\\.csv$", full.names = TRUE)
if (length(attempts)) {
  prior <- do.call(rbind, lapply(attempts, utils::read.csv, stringsAsFactors = FALSE))
  desired <- if (mode == "run") "completed" else c("completed", "smoke_passed")
  good <- prior$status %in% desired &
    prior$run_key == task$run_key & prior$manifest_hash == manifest_hash
  if (any(good)) {
    cat("Task", task_id, "already completed with matching provenance.\n")
    quit(save = "no", status = 0L)
  }
  substantive_failure <- prior$status == "error" &
    prior$run_key == task$run_key & prior$manifest_hash == manifest_hash
  if (any(substantive_failure) &&
      !mic_truthy(mic_env("MIC_ENTERPRISE_RETRY", "false"))) {
    stop("A failed attempt is preserved for this task. Set MIC_ENTERPRISE_RETRY=true ",
         "for a separately recorded retry.", call. = FALSE)
  }
}
attempt <- length(attempts) + 1L
result_path <- file.path(task_dir, sprintf("attempt-%03d.csv", attempt))
summary_path <- file.path(task_dir, sprintf("attempt-%03d-posterior.csv", attempt))

base_row <- function(status, error = "") data.frame(
  schema_version = 1L, manifest_hash = manifest_hash, task_id = task_id,
  run_key = task$run_key, attempt = attempt, status = status, error = error,
  tier = task$tier, source = task$source, fixture = task$fixture,
  demand = task$demand, model = task$model, intercept = task$intercept,
  scenario = task$scenario, dgp_seed = task$dgp_seed,
  benchmark_kind = task$benchmark_kind,
  true_k = task$true_k, true_role_scope = task$true_role_scope,
  sampler_seed = task$sampler_seed, chains = task$chains, iter = task$iter,
  warmup = task$warmup, cores = task$cores, adapt_delta = task$adapt_delta,
  max_treedepth = task$max_treedepth, slurm_job_id = mic_env("SLURM_JOB_ID"),
  slurm_array_task_id = mic_env("SLURM_ARRAY_TASK_ID"), host = Sys.info()[["nodename"]],
  started_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  stringsAsFactors = FALSE)

if (mode == "smoke") {
  row <- tryCatch({
    checked <- mic_smoke_task(task)
    cbind(base_row("smoke_passed"), checked[, setdiff(names(checked),
      c("task_id", "run_key", "source", "fixture", "demand", "model", "intercept",
        "true_k", "true_role_scope", "status")), drop = FALSE])
  }, error = function(e) base_row("error", conditionMessage(e)))
  mic_atomic_csv(row, result_path)
  cat("Task", task_id, row$status, "\n")
  quit(save = "no", status = if (row$status == "smoke_passed") 0L else 1L)
}

lib <- mic_env("MIC_ENTERPRISE_LIB")
if (nzchar(lib)) .libPaths(c(path.expand(lib), .libPaths()))
early_failure <- function(message) {
  row <- base_row("error", message)
  row$wall_seconds <- 0
  mic_atomic_csv(row, result_path)
  cat("Task", task_id, "preflight error:", message, "\n")
  quit(save = "no", status = 1L)
}
needed <- c("antitrustBayes", "rstan", "posterior", "digest")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) early_failure(paste("Missing run dependencies:",
                                         paste(missing, collapse = ", ")))
spec_version <- tryCatch(get(".antitrust_bayes_model_spec_version",
                            envir = asNamespace("antitrustBayes")),
                         error = function(e) NA_integer_)
if (!identical(as.integer(spec_version), 11L)) {
  early_failure(paste0("Enterprise benchmark requires antitrustBayes model ",
    "specification 11; installed ", spec_version, "."))
}
required_api <- c("share_likelihood", "margin_likelihood", "quality_parameterization",
  "quality_scope", "quality_scale_parameterization", "market_intercept_parameterization",
  "s0_pooling", "market_s0_parameterization", "partition_prior", "k_support", "k_prior")
missing_api <- setdiff(required_api,
  names(formals(getExportedValue("antitrustBayes", "bayes_est"))))
if (length(missing_api)) early_failure(paste("Installed bayes_est() lacks:",
                                             paste(missing_api, collapse = ", ")))
expected_kernel_hashes <- mic_expected_kernel_hashes()
installed_models <- get("stanmodels", envir = asNamespace("antitrustBayes"))
if (!all(names(expected_kernel_hashes) %in% names(installed_models)))
  early_failure("Installed antitrustBayes lacks one or more required Stan kernels.")
installed_kernel_hashes <- vapply(names(expected_kernel_hashes), function(name)
  digest::digest(installed_models[[name]]@model_code, algo = "sha256",
                 serialize = FALSE), character(1))
if (!identical(unname(installed_kernel_hashes), unname(expected_kernel_hashes))) {
  bad <- names(expected_kernel_hashes)[installed_kernel_hashes != expected_kernel_hashes]
  early_failure(paste("Installed Stan kernel hash mismatch:", paste(bad, collapse = ", ")))
}

## One array element owns one four-chain fit.  Prevent accidental nested BLAS
## oversubscription; chain parallelism remains controlled by bayes_est().
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
repo <- mic_enterprise_root()
repo_info <- mic_git_info(repo)
metadata_path <- file.path(run_dir, "manifest_metadata.csv")
if (!file.exists(metadata_path)) early_failure("Missing manifest_metadata.csv.")
manifest_metadata <- tryCatch(utils::read.csv(metadata_path, stringsAsFactors = FALSE),
                              error = function(e) NULL)
if (is.null(manifest_metadata) || nrow(manifest_metadata) != 1L ||
    !identical(as.character(manifest_metadata$manifest_hash),
               as.character(manifest_hash))) {
  early_failure("Manifest metadata is missing or does not match the manifest hash.")
}
if (!identical(as.character(repo_info$sha),
               as.character(manifest_metadata$antitrustBayes_git))) {
  early_failure("The checkout SHA differs from the SHA that froze the manifest.")
}
if (isTRUE(repo_info$dirty) &&
    !mic_truthy(mic_env("MIC_ENTERPRISE_ALLOW_DIRTY", "false"))) {
  early_failure("The antitrustBayes checkout is dirty; commit or use a new frozen run.")
}

## Independent preflight that the local DGP schedules have not drifted from
## the package equations used by estimation.
for (d in c("logit", "ces")) for (m in c("B", "C", "MC", "A")) {
  q <- c(.03, .11, .24); s0 <- .27; supply <- mic_model_contract(m)$supply
  local <- mic_conduct_schedule(q, m, d, s0, gamma = 3.5)
  public <- antitrustBayes::conduct_schedule(q, supply, s0 = s0,
    demand_model = d, gamma = 3.5)
  if (!isTRUE(all.equal(local, public, tolerance = 1e-13))) {
    early_failure(paste("DGP/public conduct schedule mismatch for", d, m))
  }
}
fixture <- NULL; fit <- NULL; warnings <- character(); started <- proc.time()[["elapsed"]]

summarise_fit <- function(fit, task, fixture, wall_seconds) {
  stanfit <- fit@stanfit
  sampler <- rstan::get_sampler_params(stanfit, inc_warmup = FALSE)
  divergences <- sum(vapply(sampler, function(x) sum(x[, "divergent__"]), numeric(1)))
  depth_hits <- sum(vapply(sampler, function(x)
    sum(x[, "treedepth__"] >= task$max_treedepth), numeric(1)))
  bfmi <- vapply(sampler, function(x) {
    e <- x[, "energy__"]
    if (length(e) < 2L || !is.finite(stats::var(e)) || stats::var(e) <= 0) return(NA_real_)
    mean(diff(e)^2) / stats::var(e)
  }, numeric(1))
  pars <- intersect(c("alpha_structural_market", "mu_log_a", "log_kappa", "s0_prob",
    "mu_b_market_vec", "sigma_b_market_vec", "mu_s0_pool_vec", "sigma_s0_pool_vec",
    "sigma_b_strat_economic", "sigma_share", "sigma_margin", "rho",
    "share_variance_fraction"), unique(sub("\\[.*$", "", stanfit@sim$fnames_oi)))
  draws <- posterior::as_draws_array(rstan::extract(stanfit, pars = pars,
    permuted = FALSE, inc_warmup = FALSE))
  sm <- posterior::summarise_draws(draws, "mean", "sd", "mcse_mean", "rhat",
                                   "ess_bulk", "ess_tail")
  mic_atomic_csv(cbind(data.frame(schema_version = 1L, manifest_hash = manifest_hash,
    task_id = task$task_id, run_key = task$run_key, attempt = attempt,
    stringsAsFactors = FALSE), sm), summary_path)
  finite <- function(x) x[is.finite(x)]
  safe_max <- function(x) { x <- finite(x); if (length(x)) max(x) else NA_real_ }
  safe_min <- function(x) { x <- finite(x); if (length(x)) min(x) else NA_real_ }
  max_rhat <- safe_max(sm$rhat)
  min_bulk <- safe_min(sm$ess_bulk)
  min_tail <- safe_min(sm$ess_tail)
  native <- tryCatch(rstan::get_elapsed_time(stanfit), error = function(e) NULL)
  k_error <- NA_real_
  if (mic_model_contract(task$model)$structure != "simultaneous") {
    pk <- tryCatch(antitrustBayes::posterior_k(fit), error = function(e) NULL)
    if (is.data.frame(pk) && all(c("market", "posterior") %in% names(pk))) {
      sums <- tapply(pk$posterior, pk$market, sum)
      k_error <- max(abs(sums - 1))
    }
  }
  dat <- fixture$data
  row <- base_row("completed")
  row$N <- nrow(dat); row$N_market <- length(unique(dat$market))
  row$N_period <- length(unique(interaction(dat$market, dat$year, drop = TRUE)))
  row$observed_margins <- sum(dat$margin_observed); row$fixture_hash <- fixture$source_hash
  row$fixed_rank_hash <- mic_rank_hash(fixture)
  row$wall_seconds <- wall_seconds
  row$native_warmup_seconds_mean <- if (is.null(native)) NA_real_ else mean(native[, "warmup"])
  row$native_sampling_seconds_mean <- if (is.null(native)) NA_real_ else mean(native[, "sampling"])
  row$peak_rss_kb <- {
    z <- readLines("/proc/self/status", warn = FALSE)
    hit <- z[grepl("^VmHWM:", z)]
    if (length(hit)) as.numeric(sub("^VmHWM:\\s+([0-9]+).*$", "\\1", hit[[1L]])) else NA_real_
  }
  row$divergences <- divergences; row$treedepth_hits <- depth_hits
  row$min_bfmi <- if (all(is.na(bfmi))) NA_real_ else min(bfmi, na.rm = TRUE)
  row$max_rank_rhat <- max_rhat; row$min_bulk_ess <- min_bulk; row$min_tail_ess <- min_tail
  row$min_bulk_ess_per_second <- min_bulk / wall_seconds
  row$min_tail_ess_per_second <- min_tail / wall_seconds
  row$posterior_k_sum_error <- k_error
  row$usable <- if (task$benchmark_kind == "timing") NA else
    divergences == 0L && depth_hits == 0L && all(is.finite(bfmi)) &&
      min(bfmi) >= .30 && is.finite(max_rhat) && max_rhat <= 1.01 &&
      is.finite(min_bulk) && min_bulk >= 100 && is.finite(min_tail) &&
      min_tail >= 100 && (is.na(k_error) || k_error <= 1e-10)
  row$warnings <- paste(unique(warnings), collapse = " | ")
  row$antitrustBayes_version <- as.character(utils::packageVersion("antitrustBayes"))
  row$model_specification_version <- spec_version
  row$antitrustBayes_git <- repo_info$sha; row$antitrustBayes_dirty <- repo_info$dirty
  row$installed_kernel_hashes <- paste(installed_kernel_hashes, collapse = ";")
  row$R_version <- R.version.string
  for (i in seq_along(bfmi)) row[[paste0("bfmi_chain", i)]] <- bfmi[[i]]
  row
}

row <- tryCatch({
  fixture <- mic_load_fixture(task)
  invisible(mic_smoke_task(task))
  fit_args <- mic_fit_args(task, fixture)
  fit <- withCallingHandlers(do.call(antitrustBayes::bayes_est, fit_args),
    warning = function(w) { warnings <<- c(warnings, conditionMessage(w));
      invokeRestart("muffleWarning") })
  elapsed <- proc.time()[["elapsed"]] - started
  summarise_fit(fit, task, fixture, elapsed)
}, error = function(e) {
  out <- base_row("error", conditionMessage(e))
  out$wall_seconds <- proc.time()[["elapsed"]] - started
  out$warnings <- paste(unique(warnings), collapse = " | ")
  out$antitrustBayes_git <- repo_info$sha; out$antitrustBayes_dirty <- repo_info$dirty
  out
})
mic_atomic_csv(row, result_path)
if (!is.null(fit)) rm(fit)
invisible(gc())
cat("Task", task_id, "attempt", attempt, row$status,
    if ("usable" %in% names(row)) paste("usable", row$usable) else "", "\n")
quit(save = "no", status = if (row$status == "completed") 0L else 1L)
