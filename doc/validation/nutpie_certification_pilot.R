#!/usr/bin/env Rscript
# A new, bounded sampler-control experiment; the original validation stays sealed.
# Use run_nutpie_certification.sh --prepare, --pilot, or --collect.
# Sourcing defines functions only. This runner never enables public evidence/BMA.

nutpie_pilot_design <- function() {
  tasks <- expand.grid(fixture = c("B", "CF_pooled", "CES_CS"),
    variant = c("cmdstanr", "nutpie_diag", "nutpie_low_rank"), stringsAsFactors = FALSE)
  tasks <- tasks[order(match(tasks$fixture, c("B", "CF_pooled", "CES_CS")),
    match(tasks$variant, c("cmdstanr", "nutpie_diag", "nutpie_low_rank"))), ]
  rownames(tasks) <- NULL
  tasks$id <- paste(tasks$fixture, tasks$variant, sep = "__")
  list(tasks = tasks,
    controls = list(chains = 4L, warmup = 2000L, iter = 3000L, cores = 2L,
      seed = 419977L, init = 0, refresh = 500L,
      control = list(adapt_delta = .99, max_treedepth = 13L)),
    diagnostic_thresholds = c(max_divergences = 0, max_depth_hits = 0, min_bfmi = .3,
      max_rhat = 1.01, min_bulk_ess = 400, min_tail_ess = 400),
    confirmation = list(seeds = c(420101L, 420211L), chains = 4L, warmup = 2000L,
      retained = 2000L, target_accept = .99, maximum_depth = 13L,
      evidence_models = c("B", "MC", "C"), policy_model = "CES_CS",
      bridge_repetitions = 3L, max_logml_tolerance = .25, max_probability_difference = .05,
      policy_four_mcse_caps = c(mean_or_quantile_price_pp = .25,
        sd_relative_to_reference_sd = .10, threshold_probability = .03)),
    decision = paste("Pilot only: report every attempt; no retries, replacement seeds,",
      "automatic control escalation, public gate changes, or automatic confirmation.",
      "Only a backend/adaptation passing all three cases is eligible for this",
      "cross-case confirmation design. Confirmation uses independent fixed seeds.",
      "A failed pilot does not establish an upstream bug or invalidate the target."))
}

nutpie_pilot_hash <- function(path) digest::digest(file = path, algo = "sha256")

nutpie_pilot_runtime <- function(repository) {
  list(runtime = nutpie_validation_runtime(repository),
    runner = nutpie_pilot_hash(file.path(repository, "doc/validation/nutpie_certification_pilot.R")),
    wrapper = nutpie_pilot_hash(file.path(repository, "doc/validation/run_nutpie_certification.sh")))
}

nutpie_pilot_assert_manifest <- function(manifest, repository) {
  if (!identical(manifest$design, nutpie_pilot_design()) ||
      !identical(manifest$runtime, nutpie_pilot_runtime(repository)))
    stop("Pilot design, runtime or source changed after preparation; retain this run and choose a new output root.")
  if (!identical(nutpie_pilot_hash(manifest$original_manifest), manifest$original_manifest_hash))
    stop("Original validation manifest changed; no original artifacts may be overwritten.")
  original <- readRDS(manifest$original_manifest)
  if (!identical(original$fixture_hash, digest::digest(original$fixtures, algo = "sha256")) ||
      !identical(manifest$fixtures, original$fixtures[unique(manifest$design$tasks$fixture)]))
    stop("Pilot fixtures differ from the sealed original observations/specifications.")
  for (reference in c(manifest$reference_fits, manifest$constant_reference_fits)) {
    if (!identical(nutpie_pilot_hash(reference$path), reference$sha256))
      stop("Original reference fit changed: ", reference$path)
  }
  for (fixture in names(manifest$reference_fits)) {
    reference <- readRDS(manifest$reference_fits[[fixture]]$path)
    native <- readRDS(manifest$constant_reference_fits[[fixture]]$path)
    if (!identical(nutpie_validation_identity(reference), nutpie_validation_identity(native)))
      stop("Original native constant reference has a different target: ", fixture)
  }
  sources <- vapply(names(manifest$stan_source_hash), function(kernel)
    nutpie_pilot_hash(antitrustBayes:::.cmdstan_kernel_file(kernel)), character(1))
  if (!identical(sources, manifest$stan_source_hash)) stop("Stan target source changed.")
  invisible(TRUE)
}

nutpie_pilot_constant_reference <- function(draws) {
  values <- lapply(seq_len(dim(draws)[3L]), function(i) {
    x <- as.numeric(draws[, , i])
    if (length(x) && all(is.finite(x)) && all(x == x[1L])) x[1L] else NULL
  })
  names(values) <- dimnames(draws)[[3L]]
  values[!vapply(values, is.null, logical(1))]
}

nutpie_pilot_posterior_availability <- function(draws, summary, reference_constants) {
  variables <- dimnames(draws)[[3L]]
  if (anyDuplicated(summary$variable) || !setequal(variables, summary$variable))
    stop("Posterior diagnostic columns do not match retained draw columns.")
  summary <- summary[match(variables, summary$variable), , drop = FALSE]
  do.call(rbind, lapply(seq_along(variables), function(i) {
    values <- as.numeric(draws[, , i])
    finite <- all(is.finite(values))
    constant <- finite && length(values) > 0L && all(values == values[1L])
    known_constant <- constant && variables[i] %in% names(reference_constants) &&
      identical(values[1L], as.numeric(reference_constants[[variables[i]]]))
    available <- finite && if (constant) known_constant else
      all(is.finite(as.numeric(summary[i, c("rhat", "ess_bulk", "ess_tail")])) )
    data.frame(variable = variables[i], constant = constant,
      constant_matches_reference = known_constant, diagnostics_available = available)
  }))
}

nutpie_pilot_diagnostic_pass <- function(d, thresholds) {
  needed <- c("divergences", "depth_hits", "bfmi", "max_rhat", "min_bulk_ess", "min_tail_ess",
    "unavailable_posterior_diagnostics")
  if (!is.data.frame(d) || nrow(d) != 4L || !all(needed %in% names(d)) ||
      !all(is.finite(as.matrix(d[needed])))) return(FALSE)
  all(d$unavailable_posterior_diagnostics == 0) &&
    all(d$divergences <= thresholds[["max_divergences"]]) &&
    all(d$depth_hits <= thresholds[["max_depth_hits"]]) &&
    all(d$bfmi >= thresholds[["min_bfmi"]]) && all(d$max_rhat <= thresholds[["max_rhat"]]) &&
    all(d$min_bulk_ess >= thresholds[["min_bulk_ess"]]) &&
    all(d$min_tail_ess >= thresholds[["min_tail_ess"]])
}

nutpie_pilot_collect <- function(manifest, root, report_dir) {
  results <- lapply(manifest$design$tasks$id, function(id) {
    path <- file.path(root, id, "result.rds")
    if (file.exists(path)) readRDS(path) else list(status = "not_run", reason = "No execution artifact.")
  })
  status <- do.call(rbind, lapply(seq_along(results), function(i) {
    task <- manifest$design$tasks[i, ]; value <- results[[i]]
    data.frame(task, status = value$status,
      diagnostics_passed = isTRUE(value$diagnostics_passed),
      identity_equal = isTRUE(value$identity_equal),
      reason = if (is.null(value$reason)) "" else value$reason)
  }))
  utils::write.csv(status, file.path(report_dir, "nutpie_certification_pilot_status.csv"), row.names = FALSE)
  rows <- lapply(seq_along(results), function(i) {
    result <- results[[i]]
    if (is.null(result[["diagnostics"]])) return(NULL)
    cbind(manifest$design$tasks[i, , drop = FALSE], result[["diagnostics"]]$summary,
      compile_lookup_seconds = result$compile_lookup_seconds,
      estimation_seconds = result$estimation_seconds,
      chains = manifest$design$controls$chains, warmup = manifest$design$controls$warmup,
      retained = manifest$design$controls$iter - manifest$design$controls$warmup,
      seed = manifest$design$controls$seed, target_accept = manifest$design$controls$control$adapt_delta)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows)) utils::write.csv(do.call(rbind, rows),
    file.path(report_dir, "nutpie_certification_pilot_diagnostics.csv"), row.names = FALSE)
  eligible <- vapply(c("cmdstanr", "nutpie_diag", "nutpie_low_rank"), function(variant) {
    selected <- status[status$variant == variant, ]
    nrow(selected) == 3L && all(selected$status == "completed") &&
      all(selected$identity_equal) && all(selected$diagnostics_passed)
  }, logical(1))
  outcome <- list(status = if (all(status$status == "completed")) "completed" else "incomplete",
    eligible = eligible, public_evidence_gate_changed = FALSE, tasks = status)
  saveRDS(outcome, file.path(root, "outcome.rds"))
  invisible(outcome)
}

nutpie_pilot_worker <- function(task_path) {
  if (!identical(Sys.getenv("AB_NUTPIE_PILOT_EXECUTE"), "true"))
    stop("Pilot inference requires AB_NUTPIE_PILOT_EXECUTE=true.")
  work <- readRDS(task_path)
  .libPaths(work$library_paths)
  pkgload::load_all(work$repository, quiet = TRUE)
  source(file.path(work$repository, "doc/validation/nutpie_backend_validation.R"))
  manifest <- readRDS(file.path(work$root, "manifest.rds"))
  nutpie_pilot_assert_manifest(manifest, work$repository)
  options(antitrustBayes.cmdstan_cache = file.path(manifest$original_root, "cmdstan-cache"),
    antitrustBayes.nutpie_cache = file.path(manifest$original_root, "nutpie-cache"))
  if (nzchar(Sys.getenv("AB_NUTPIE_CMDSTAN"))) cmdstanr::set_cmdstan_path(Sys.getenv("AB_NUTPIE_CMDSTAN"))
  task <- work$task
  result <- tryCatch({
    args <- manifest$fixtures[[task$fixture]]$arguments
    args$stan_options <- manifest$design$controls
    backend <- if (task$variant == "cmdstanr") "cmdstanr" else "nutpie"
    args$backend <- backend
    if (backend == "cmdstanr") {
      args$stan_options$output_dir <- file.path(work$directory, "cmdstan")
      dir.create(args$stan_options$output_dir, recursive = TRUE, showWarnings = FALSE)
      if (isTRUE(manifest$fixtures[[task$fixture]]$pooled)) args$stan_options$control$metric <- "dense_e"
    } else args$nutpie_options <- list(adaptation = if (task$variant == "nutpie_diag") "diag" else "low_rank")
    saveRDS(args, file.path(work$directory, "arguments.rds"))
    reference <- readRDS(manifest$reference_fits[[task$fixture]]$path)
    kernel <- reference@metadata$stan_kernel
    start <- proc.time()[["elapsed"]]
    if (backend == "cmdstanr") antitrustBayes:::.cmdstan_model_for_kernel(kernel) else
      antitrustBayes:::.nutpie_model_for_kernel(kernel)
    compile_lookup_seconds <- proc.time()[["elapsed"]] - start
    start <- proc.time()[["elapsed"]]
    fit <- do.call(antitrustBayes::bayes_est, args)
    estimation_seconds <- proc.time()[["elapsed"]] - start
    fit_path <- file.path(work$directory, "fit.rds")
    saveRDS(fit, fit_path)
    reloaded <- readRDS(fit_path)
    expected <- c(manifest$design$controls$iter - manifest$design$controls$warmup,
      manifest$design$controls$chains)
    draws <- antitrustBayes:::.bayes_backend_draws(reloaded@stanfit)
    stopifnot(identical(as.integer(dim(draws)[1:2]), as.integer(expected)),
      identical(draws, antitrustBayes:::.bayes_backend_draws(fit@stanfit)))
    identity_equal <- identical(nutpie_validation_identity(reloaded), nutpie_validation_identity(reference))
    if (!identity_equal) stop("Pilot target/data/prior/model identity differs from the sealed reference.")
    diagnostics <- nutpie_validation_diagnostics(reloaded)
    constant_reference <- if (backend == "cmdstanr") reference else
      readRDS(manifest$constant_reference_fits[[task$fixture]]$path)
    reference_constants <- nutpie_pilot_constant_reference(antitrustBayes:::.bayes_backend_draws(constant_reference@stanfit))
    availability <- nutpie_pilot_posterior_availability(draws, diagnostics$common$posterior, reference_constants)
    diagnostics$summary$unavailable_posterior_diagnostics <- sum(!availability$diagnostics_available)
    list(status = "completed", reason = "", identity_equal = identity_equal,
      diagnostics = diagnostics, parameters = nutpie_validation_parameters(reloaded),
      posterior_diagnostic_availability = availability,
      diagnostics_passed = nutpie_pilot_diagnostic_pass(diagnostics$summary, manifest$design$diagnostic_thresholds),
      compile_lookup_seconds = compile_lookup_seconds, estimation_seconds = estimation_seconds,
      fit_sha256 = nutpie_pilot_hash(fit_path), provenance = reloaded@metadata,
      posterior_sha256 = digest::digest(draws, algo = "sha256"))
  }, error = function(e) list(status = "failed", reason = conditionMessage(e)))
  saveRDS(result, file.path(work$directory, "result.rds"))
}

nutpie_pilot_main <- function(action) {
  if (!action %in% c("--prepare", "--pilot", "--collect")) stop("Use --prepare, --pilot, or --collect.")
  repository <- normalizePath(".", mustWork = TRUE)
  pkgload::load_all(repository, quiet = TRUE)
  source(file.path(repository, "doc/validation/nutpie_backend_validation.R"))
  root <- Sys.getenv("AB_NUTPIE_PILOT_OUTPUT", file.path(repository, ".scratch/nutpie-certification/pilot-accept99"))
  original <- normalizePath(Sys.getenv("AB_NUTPIE_ORIGINAL_OUTPUT", file.path(repository, ".scratch/nutpie-validation")), mustWork = TRUE)
  if (identical(normalizePath(root, mustWork = FALSE), original)) stop("Pilot requires a separate output root.")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  root <- normalizePath(root)
  path <- file.path(root, "manifest.rds")
  if (!file.exists(path)) {
    old_path <- file.path(original, "manifest.rds"); old <- readRDS(old_path)
    if (!identical(old$fixture_hash, digest::digest(old$fixtures, algo = "sha256"))) stop("Original fixture seal failed.")
    design <- nutpie_pilot_design()
    fixtures <- old$fixtures[unique(design$tasks$fixture)]
    if (any(vapply(fixtures, is.null, logical(1)))) stop("Missing original pilot fixture.")
    references <- lapply(names(fixtures), function(id) {
      fit <- file.path(original, paste0(id, "__cmdstanr"), "fit.rds")
      list(path = fit, sha256 = nutpie_pilot_hash(fit))
    }); names(references) <- names(fixtures)
    constant_references <- lapply(names(fixtures), function(id) {
      fit <- file.path(original, paste0(id, "__nutpie_diag"), "fit.rds")
      list(path = fit, sha256 = nutpie_pilot_hash(fit))
    }); names(constant_references) <- names(fixtures)
    manifest <- list(design = design, fixtures = fixtures, runtime = nutpie_pilot_runtime(repository),
      original_root = original, original_manifest = old_path, original_manifest_hash = nutpie_pilot_hash(old_path),
      stan_source_hash = old$source_hash, reference_fits = references,
      constant_reference_fits = constant_references,
      prepared_at = Sys.time(), starting_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE))
    saveRDS(manifest, path)
  } else manifest <- readRDS(path)
  nutpie_pilot_assert_manifest(manifest, repository)
  # Keep a newly prepared run from overwriting an earlier run's tracked evidence.
  report_dir <- file.path(root, "reports")
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  if (action != "--pilot") {
    value <- nutpie_pilot_collect(manifest, root, report_dir)
    message("Pilot ", value$status, "; no inference. Manifest: ", path)
    return(invisible(value))
  }
  if (!identical(Sys.getenv("AB_NUTPIE_PILOT_EXECUTE"), "true")) stop("Explicit pilot execution opt-in required.")
  if (!nzchar(Sys.which("timeout"))) stop("GNU timeout is required.")
  for (i in seq_len(nrow(manifest$design$tasks))) {
    task <- manifest$design$tasks[i, ]; directory <- file.path(root, task$id)
    if (file.exists(file.path(directory, "task.rds"))) next # Preserve all prior attempts, including failures.
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    task_path <- file.path(directory, "task.rds")
    saveRDS(list(task = task, directory = directory, repository = repository, root = root,
      library_paths = .libPaths()), task_path)
    code <- system2("timeout", c("--signal=TERM", "--kill-after=5s", "1800",
      shQuote(file.path(R.home("bin"), "Rscript")), "--vanilla",
      shQuote(file.path(repository, "doc/validation/nutpie_certification_pilot.R")), "--worker", shQuote(task_path)),
      stdout = file.path(directory, "worker.log"), stderr = file.path(directory, "worker.log"))
    if (code != 0L || !file.exists(file.path(directory, "result.rds")))
      saveRDS(list(status = "worker_failed", reason = paste("Worker exit", code, "; artifacts preserved; no automatic retry.")),
        file.path(directory, "result.rds"))
    nutpie_pilot_collect(manifest, root, report_dir)
  }
  nutpie_pilot_collect(manifest, root, report_dir)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 2L && identical(args[[1L]], "--worker")) nutpie_pilot_worker(args[[2L]]) else {
    if (length(args) != 1L) stop("Use --prepare, --pilot, or --collect.")
    nutpie_pilot_main(args[[1L]])
  }
}
