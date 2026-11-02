# Benchmark-only overlay: no package installation or production source edits.
# Freeze both the source and plan before execution; checkpoint complete calls.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: runner.R MODE SEED")
mode <- args[[1L]]
seed <- as.integer(args[[2L]])
allowed <- c("legacy_default", "candidate_default", "candidate_seeded_serial",
             "candidate_seeded_parallel")
stopifnot(mode %in% allowed, seed %in% c(707101L, 707102L))
library(antitrustBayes)
recovery <- "doc/statistical_toolkit_runs/20260916/performance-audit/recovery"
compact <- "doc/performance_audit/recovery/bridge"
source("doc/performance_audit/recovery/checkpoint.R")
ns <- asNamespace("antitrustBayes")
overlay <- new.env(parent = globalenv())
for (n in ls(ns, all.names = TRUE)) assign(n, get(n, ns), overlay)
candidate <- !identical(mode, "legacy_default")
bma <- if (candidate) file.path(recovery, "worktrees/bridge/R/bma.R") else "R/bma.R"
sys.source(bma, overlay)
if (candidate) sys.source(file.path(recovery, "worktrees/bridge/R/model_universe.R"), overlay)
stopifnot(identical(as.numeric(formals(overlay$bayes_compare)$bridge_maxiter), 5000),
          identical(as.numeric(formals(overlay$bayes_compare)$bridge_repetitions), 3))
paths <- setNames(file.path("doc/statistical_toolkit_runs/20260916/hybrid/fits",
                   paste0("logit_", c("BF", "BS"), "_noisy_271828__D.rds")), c("BF", "BS"))
fits <- lapply(paths, readRDS)
invisible(lapply(fits, function(fit) overlay$validate_comparison_contract(list(fit))))
stopifnot(all(vapply(fits, overlay$.bayes_is_collapsed_quality, logical(1))))
fit_hashes <- vapply(fits, digest::digest, character(1), algo = "sha256")
cpp_path <- file.path(find.package("antitrustBayes"), "libs/antitrustBayes.so")
source_files <- c(bma, if (candidate) file.path(recovery, "worktrees/bridge/R/model_universe.R"),
                  "doc/performance_audit/recovery/checkpoint.R",
                  "doc/performance_audit/recovery/bridge/run_archived_bridge_benchmark.R")
spec <- list(plan_version = 2L, mode = mode, seed = seed, cores = if (grepl("parallel$", mode)) 2L else 1L,
  bridge_seed = if (grepl("seeded", mode)) seed else NULL,
  model_seeds = setNames(seed + c(0L, 10000L), names(fits)),
  bridge_maxiter = 5000L, bridge_repetitions = 3L,
  fit_file_sha256 = vapply(paths, function(p) digest::digest(file = p, algo = "sha256"), character(1)),
  source_sha256 = vapply(source_files, function(p) digest::digest(file = p, algo = "sha256"), character(1)),
  native_sha256 = digest::digest(file = cpp_path, algo = "sha256"),
  R_version = R.version.string, package_versions = vapply(c("rstan", "bridgesampling"),
    function(p) as.character(packageVersion(p)), character(1)),
  scope = "two separate complete public comparisons; different DGPs, no new HMC")
run <- function() {
  set.seed(seed)
  rng_before <- .Random.seed
  comparison <- lapply(seq_along(fits), function(i) {
    inputs <- c(fits[i], list(bridge = TRUE, bridge_maxiter = 5000L, bridge_repetitions = 3L))
    if (candidate) inputs <- c(inputs, list(bridge_cores = spec$cores,
      bridge_seed = if (is.null(spec$bridge_seed)) NULL else unname(spec$model_seeds[[i]])))
    do.call(overlay$bayes_compare, inputs)
  })
  names(comparison) <- names(fits)
  immutable <- identical(fit_hashes, vapply(fits, digest::digest, character(1), algo = "sha256"))
  if (!immutable) stop("Bridge comparison mutated an input fit.")
  if (!is.null(spec$bridge_seed) && !identical(rng_before, .Random.seed))
    stop("Explicit bridge seed failed to restore caller RNG.")
  list(comparison = comparison, rng_after = .Random.seed, fit_immutable = immutable)
}
state <- audit_checkpoint_run(paste(mode, seed, sep = "-"), spec, run,
  artifact_root = file.path(recovery, "bridge-archived"))
row <- data.frame(mode = mode, seed = seed, status = state$status, cached = state$cached,
  elapsed_seconds = unname(state$compute_time[["elapsed"]]),
  warning_count = length(state$warnings), error = if (is.null(state$error)) "" else state$error,
  native_sha256 = spec$native_sha256, artifact_fingerprint = state$fingerprint)
if (identical(state$status, "COMPLETE")) {
  rows <- lapply(state$result$comparison, function(cmp) {
    bridges <- attr(cmp, "bridge_list")[[1L]]
    do.call(rbind, lapply(seq_along(bridges), function(j) {
      b <- bridges[[j]]
      data.frame(mode = mode, seed = seed, model = cmp$model[[1L]], repetition = j,
        logml = if (is.null(b)) NA_real_ else b$logml,
        niter = if (is.null(b) || is.null(b$niter)) NA_real_ else sum(b$niter),
        method = if (is.null(b)) "FAILED" else as.character(b$method),
        fit_immutable = state$result$fit_immutable)
    }))
  })
  write.csv(do.call(rbind, rows), file.path(compact, paste0(mode, "-", seed, "-evidence.csv")), row.names = FALSE)
  row$failure_count <- sum(vapply(state$result$comparison, function(cmp) {
    if (candidate) {
      failures <- attr(cmp, "bridge_failure_records")
      if (is.null(failures)) 0L else nrow(failures)
    } else sum(vapply(attr(cmp, "bridge_list"),
      function(reps) sum(vapply(reps, is.null, logical(1))), integer(1)))
  }, integer(1)))
} else row$failure_count <- NA_integer_
write.csv(row, file.path(compact, paste0(mode, "-", seed, "-timing.csv")), row.names = FALSE)
cat(mode, seed, state$status, "elapsed", row$elapsed_seconds, "warnings", row$warning_count, "\n")
if (!identical(state$status, "COMPLETE")) quit(status = 1L)
