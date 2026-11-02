#!/usr/bin/env Rscript
# Sourcing defines functions only. Execution has a durable per-invocation
# budget; completed tasks are reused. Never changes the public evidence gate.
# Prepare/collect: Rscript --vanilla doc/validation/nutpie_evidence_validation.R
# Saved-fit checks: AB_NUTPIE_EVIDENCE_EXECUTE=true Rscript --vanilla ...
# Additional fits ALSO require AB_NUTPIE_EVIDENCE_SAMPLE=true.

nutpie_evidence_pairs <- function() data.frame(
  kernel = c("logit_supply", "logit_stackelberg", "ces_supply", "ces_stackelberg"),
  fixture = c("B", "CS", "CES_B", "CES_CS"), stringsAsFactors = FALSE)

nutpie_evidence_tasks <- function() {
  pairs <- nutpie_evidence_pairs()
  data.frame(id = c(paste0("target_", pairs$kernel), "normalized_reference",
    paste0("sample_", rep(c("MC", "C"), each = 2L), "_", rep(c("cmdstanr", "nutpie"), 2L)),
    paste0("bridge_", rep(c("B", "MC", "C"), each = 2L), "_", rep(c("cmdstanr", "nutpie"), 3L))),
    stage = c(rep("target", 4L), "reference", rep("sample", 4L), rep("bridge", 6L)),
    model = c(pairs$fixture, NA, rep(c("MC", "C"), each = 2L), rep(c("B", "MC", "C"), each = 2L)),
    backend = c(rep(NA, 5L), rep(c("cmdstanr", "nutpie"), 5L)),
    kernel = c(pairs$kernel, rep(NA, 11L)), stringsAsFactors = FALSE)
}

nutpie_evidence_budget <- function(invocation, now = Sys.time()) {
  if (is.null(invocation$started_at) || is.null(invocation$deadline))
    stop("A durable evidence invocation clock is required before execution.")
  remaining <- as.numeric(difftime(invocation$deadline, now, units = "secs"))
  if (!is.finite(remaining) || remaining <= 5)
    stop("The current evidence invocation budget is exhausted; completed artifacts are retained.")
  floor(remaining - 5)
}

nutpie_evidence_invocation <- function(directory, seconds = 14400, now = Sys.time()) {
  if (!is.numeric(seconds) || length(seconds) != 1L || !is.finite(seconds) ||
      seconds < 10 || seconds > 14400 || seconds != as.integer(seconds))
    stop("AB_NUTPIE_EVIDENCE_SECONDS must be an integer between 10 and 14400.")
  path <- file.path(directory, "invocations")
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  filename <- tempfile(paste0(format(now, "%Y%m%dT%H%M%S"), "-"), tmpdir = path, fileext = ".rds")
  invocation <- list(started_at = now, deadline = now + seconds,
    seconds = as.integer(seconds), pid = Sys.getpid(),
    scope = "This invocation only; saved scientific manifest and previous attempts are unchanged.")
  saveRDS(invocation, filename)
  filename
}

nutpie_evidence_manifest_contract <- function(manifest) {
  if (!identical(manifest$fixture_hash, digest::digest(manifest$fixtures, algo = "sha256")))
    stop("The frozen backend fixture manifest has changed.")
  kernels <- nutpie_evidence_pairs()$kernel
  if (!setequal(names(manifest$source_hash), kernels)) stop("The full four-kernel source contract is required.")
  current <- vapply(names(manifest$source_hash), function(kernel)
    digest::digest(file = antitrustBayes:::.cmdstan_kernel_file(kernel), algo = "sha256"), character(1))
  if (!identical(current, manifest$source_hash)) stop("Frozen and current Stan source hashes differ.")
  invisible(TRUE)
}

nutpie_evidence_fit_path <- function(root, model, backend, shared = FALSE) {
  if (shared && model != "B")
    return(file.path(root, "evidence", paste0("sample_", model, "_", backend), "fit.rds"))
  variant <- if (backend == "cmdstanr") "cmdstanr" else "nutpie_diag"
  file.path(root, paste0(model, "__", variant), "fit.rds")
}

nutpie_evidence_fit_ready <- function(path) {
  if (!file.exists(path)) return(FALSE)
  directory <- dirname(path)
  state_path <- if (startsWith(basename(directory), "sample_"))
    file.path(directory, "result.rds") else file.path(directory, "status.rds")
  if (!file.exists(state_path)) return(FALSE)
  state <- readRDS(state_path)$status
  state %in% c("fit_saved", "completed", "postfit_failed", "passed")
}

nutpie_evidence_load_fit <- function(path, seal_path) {
  if (!file.exists(path)) stop("Required saved fit is unavailable: ", path)
  hash <- digest::digest(file = path, algo = "sha256")
  if (file.exists(seal_path)) {
    seal <- readRDS(seal_path)
    if (!identical(seal$fit_hash, hash) || !identical(seal$path, normalizePath(path)))
      stop("A previously sealed validation fit changed: ", path)
  }
  fit <- readRDS(path)
  antitrustBayes:::.bayes_evidence_contract(fit)
  antitrustBayes:::.bayes_validate_bridge(list(fit))
  identity <- nutpie_validation_identity(fit)
  seal <- list(path = normalizePath(path), fit_hash = hash,
    identity_hash = digest::digest(identity, algo = "sha256"),
    backend = fit@metadata$backend,
    source_hash = fit@metadata$stan_source_hash, data_hash = fit@metadata$stan_data_hash)
  if (file.exists(seal_path)) {
    if (!identical(seal, readRDS(seal_path))) stop("Saved fit identity changed after sealing.")
  } else saveRDS(seal, seal_path)
  fit
}

nutpie_evidence_point_checks <- function(normalized, proportional, cmdstan,
                                          no_jacobian, proportional_no_jacobian,
                                          cmdstan_no_jacobian) {
  inputs <- list(normalized, proportional, cmdstan, no_jacobian,
                 proportional_no_jacobian, cmdstan_no_jacobian)
  if (length(unique(lengths(inputs))) != 1L || !length(normalized) ||
      any(!is.finite(unlist(inputs)))) stop("Target checks need finite paired evaluations.")
  tolerance <- 1e-7 + 1e-9 * pmax(1, abs(normalized), abs(cmdstan))
  data.frame(point = seq_along(normalized), normalized, proportional, cmdstan,
    no_jacobian, proportional_no_jacobian, cmdstan_no_jacobian,
    normalized_error = normalized - cmdstan,
    proportional_error = normalized - proportional,
    no_jacobian_error = no_jacobian - cmdstan_no_jacobian,
    jacobian = normalized - no_jacobian,
    jacobian_error = (normalized - no_jacobian) - (cmdstan - cmdstan_no_jacobian),
    proportional_no_jacobian_error = no_jacobian - proportional_no_jacobian,
    tolerance,
    passed = abs(normalized - cmdstan) <= tolerance &
      abs(normalized - proportional) <= tolerance &
      abs(no_jacobian - cmdstan_no_jacobian) <= tolerance &
      abs(no_jacobian - proportional_no_jacobian) <= tolerance &
      abs((normalized - no_jacobian) - (cmdstan - cmdstan_no_jacobian)) <= tolerance)
}

nutpie_evidence_target <- function(task, root, directory) {
  paths <- lapply(c("cmdstanr", "nutpie"), function(backend)
    nutpie_evidence_fit_path(root, task$model, backend))
  fits <- lapply(seq_along(paths), function(i) nutpie_evidence_load_fit(paths[[i]],
    file.path(directory, paste0("fit-seal-", i, ".rds"))))
  if (!identical(nutpie_validation_identity(fits[[1L]]), nutpie_validation_identity(fits[[2L]])))
    stop("Paired target fits have different source/data/prior/support identities.")
  if (!all(vapply(fits, function(fit) identical(fit@metadata$stan_kernel, task$kernel), logical(1))))
    stop("Target task and saved fit kernels differ.")
  cmdstan <- antitrustBayes:::.bayes_bridge_restore_context(fits[[1L]])
  on.exit(cmdstan$cleanup(), add = TRUE)
  nutpie <- antitrustBayes:::.bayes_nutpie_bridge_context(fits[[2L]])
  on.exit(nutpie$cleanup(), add = TRUE)
  model <- environment(nutpie$log_posterior)$model
  canonical_names <- function(x) vapply(strsplit(x, ".", fixed = TRUE), function(parts) {
    if (length(parts) == 1L) parts else paste0(parts[[1L]], "[", paste(parts[-1L], collapse = ","), "]")
  }, character(1))
  cmdstan_coordinates <- colnames(cmdstan$samples$unconstrain_draws(format = "matrix"))
  if (!identical(canonical_names(cmdstan_coordinates), canonical_names(colnames(nutpie$samples))))
    stop("CmdStan and BridgeStan unconstrained parameter names/order differ.")
  selected <- unique(as.integer(round(seq(1, nrow(nutpie$samples), length.out = 20L))))
  points <- nutpie$samples[selected, , drop = FALSE]
  ## Nonposterior interior points detect a constant offset or order error that
  ## coincidentally disappears at a single saved state.
  center <- colMeans(nutpie$samples)
  extra <- t(vapply(c(.05, -.05, .1), function(scale)
    center + scale * sin(seq_along(center)), numeric(length(center))))
  points <- rbind(points, extra)
  evaluate <- function(fun) vapply(seq_len(nrow(points)), function(i)
    as.numeric(fun(points[i, ])), numeric(1))
  checks <- nutpie_evidence_point_checks(
    evaluate(function(q) model$log_density(q, propto = FALSE, jacobian = TRUE)),
    evaluate(function(q) model$log_density(q, propto = TRUE, jacobian = TRUE)),
    evaluate(function(q) cmdstan$samples$log_prob(q, jacobian = TRUE)),
    evaluate(function(q) model$log_density(q, propto = FALSE, jacobian = FALSE)),
    evaluate(function(q) model$log_density(q, propto = TRUE, jacobian = FALSE)),
    evaluate(function(q) cmdstan$samples$log_prob(q, jacobian = FALSE)))
  checks$kernel <- task$kernel
  list(status = if (all(checks$passed) && any(abs(checks$jacobian) > 1e-8)) "passed" else "failed",
    checks = checks, points = points, coordinates = colnames(nutpie$samples),
    cmdstan_coordinates = cmdstan_coordinates,
    constrained_coordinates = nutpie$constrained_parameter_names,
    identity = nutpie_validation_identity(fits[[1L]]),
    model_info = model$model_info(), model_version = model$model_version(),
    bridgestan_client = as.character(utils::packageVersion("bridgestan")),
    conventions = list(propto = FALSE, jacobian = TRUE, bounds = "unbounded",
      order = "Stan unconstrained order; iterations within chains"))
}

nutpie_evidence_summary <- function(repetitions) {
  if (length(repetitions) != 3L || any(!vapply(repetitions, function(x)
    is.numeric(x$logml) && length(x$logml) == 1L && is.finite(x$logml), logical(1))))
    stop("Every one of the three planned bridge repetitions must succeed.")
  values <- vapply(repetitions, `[[`, numeric(1), "logml")
  mcse <- vapply(repetitions, function(x) {
    value <- x$mcse_logml
    if (is.numeric(value) && length(value) == 1L && is.finite(value) && value >= 0)
      as.numeric(value) else NA_real_
  }, numeric(1))
  spread <- stats::sd(values)
  list(logml = stats::median(values), logml_sd = spread,
    uncertainty = max(c(spread, mcse[is.finite(mcse)])),
    mcse_available = all(is.finite(mcse)), repetitions = values, mcse = mcse)
}

nutpie_evidence_reference_density <- function(q, data = NULL) {
  ## A normalized Exp(1) x Uniform(0,1) reference on unconstrained coordinates.
  ## This is an integration check, never an alternate economic-model target.
  -exp(q[[1L]]) + q[[1L]] + stats::plogis(q[[2L]], log.p = TRUE) +
    stats::plogis(q[[2L]], lower.tail = FALSE, log.p = TRUE)
}

nutpie_evidence_reference <- function(directory) {
  set.seed(419983L)
  samples <- cbind(scale = log(stats::rexp(4000L)), share = stats::qlogis(stats::runif(4000L)))
  bounds <- stats::setNames(c(Inf, Inf), colnames(samples))
  repetitions <- lapply(1:3, function(i) {
    set.seed(419984L + i)
    value <- bridgesampling::bridge_sampler(samples, log_posterior = nutpie_evidence_reference_density,
      data = NULL, lb = -bounds, ub = bounds, repetitions = 1L, method = "normal",
      cores = 1L, silent = TRUE, maxiter = 5000L)
    saveRDS(value, file.path(directory, paste0("bridge-", i, ".rds")))
    value
  })
  summary <- nutpie_evidence_summary(repetitions)
  tolerance <- max(.02, 4 * summary$uncertainty)
  list(status = if (abs(summary$logml) <= tolerance) "passed" else "failed",
    expected_logml = 0, tolerance = tolerance, summary = summary,
    scope = "Analytic normalized reference, iid exact draws; no Stan or HMC.")
}

nutpie_evidence_shared_arguments <- function(backend_manifest, model, backend, directory) {
  fixture <- backend_manifest$fixtures$B
  if (is.null(fixture) || !model %in% c("MC", "C")) stop("Unknown shared-data evidence model.")
  args <- fixture$arguments
  args$supply_model <- if (model == "MC") "monopolistic" else "cournot"
  args$structure <- "simultaneous"
  args$backend <- backend
  args$stan_options <- list(chains = 4L, iter = 2000L, warmup = 1000L, cores = 2L,
    seed = 419771L, refresh = 100L, init = 0,
    control = list(adapt_delta = .95, max_treedepth = 13L))
  if (backend == "cmdstanr") args$stan_options$output_dir <- file.path(directory, "cmdstan")
  if (backend == "nutpie") args$nutpie_options <- list(adaptation = "diag")
  args
}

nutpie_evidence_sample <- function(task, root, directory, backend_manifest) {
  if (!identical(Sys.getenv("AB_NUTPIE_EVIDENCE_SAMPLE", "false"), "true"))
    stop("Additional inference requires explicit AB_NUTPIE_EVIDENCE_SAMPLE=true opt-in.")
  args <- nutpie_evidence_shared_arguments(backend_manifest, task$model, task$backend, directory)
  dir.create(file.path(directory, "cmdstan"), recursive = TRUE, showWarnings = FALSE)
  saveRDS(args, file.path(directory, "arguments.rds"))
  fit <- do.call(antitrustBayes::bayes_est, args)
  saveRDS(fit, file.path(directory, "fit.rds"))
  reloaded <- nutpie_evidence_load_fit(file.path(directory, "fit.rds"), file.path(directory, "fit-seal.rds"))
  stopifnot(identical(antitrustBayes:::.bayes_backend_draws(fit@stanfit),
                      antitrustBayes:::.bayes_backend_draws(reloaded@stanfit)))
  list(status = "passed", identity = nutpie_validation_identity(reloaded),
    diagnostics = nutpie_validation_diagnostics(reloaded),
    fit_hash = digest::digest(file = file.path(directory, "fit.rds"), algo = "sha256"))
}

nutpie_evidence_require_targets <- function(root) {
  paths <- file.path(root, "evidence", paste0("target_", nutpie_evidence_pairs()$kernel), "result.rds")
  if (!all(file.exists(paths)) || !all(vapply(paths, function(path)
    identical(readRDS(path)$status, "passed"), logical(1))))
    stop("All four native target checks must pass before any economic-model bridge calculation.")
  invisible(TRUE)
}

nutpie_evidence_bridge <- function(task, root, directory) {
  nutpie_evidence_require_targets(root)
  models <- c("B", "MC", "C")
  fits <- stats::setNames(lapply(models, function(model) nutpie_evidence_load_fit(
    nutpie_evidence_fit_path(root, model, task$backend, shared = TRUE),
    file.path(directory, paste0("fit-seal-", model, ".rds")))), models)
  antitrustBayes:::validate_comparison_contract(fits)
  antitrustBayes:::.bayes_validate_bridge(fits)
  fit <- fits[[task$model]]
  context <- if (task$backend == "nutpie") antitrustBayes:::.bayes_nutpie_bridge_context(fit) else
    antitrustBayes:::.bayes_bridge_restore_context(fit)
  on.exit(context$cleanup(), add = TRUE)
  repetitions <- lapply(1:3, function(i) {
    set.seed(419997L + i)
    value <- do.call(bridgesampling::bridge_sampler, c(list(samples = context$samples,
      repetitions = 1L, method = "normal", maxiter = 5000L, silent = TRUE), context$bridge_args))
    saveRDS(value, file.path(directory, paste0("bridge-", i, ".rds")))
    value
  })
  list(status = "passed", model = task$model, backend = task$backend,
    summary = nutpie_evidence_summary(repetitions),
    identity = nutpie_validation_identity(fit),
    diagnostics = nutpie_validation_diagnostics(fit),
    common_observation_hash = digest::digest(fit@data, algo = "sha256"))
}

nutpie_evidence_compare <- function(reference, alternative) {
  if (!identical(reference$identity, alternative$identity))
    stop("Bridge results describe different source/data/prior/support identities.")
  difference <- alternative$summary$logml - reference$summary$logml
  uncertainty <- sqrt(reference$summary$uncertainty^2 + alternative$summary$uncertainty^2)
  tolerance <- max(1e-6, 4 * uncertainty)
  mcse_available <- isTRUE(reference$summary$mcse_available) &&
    isTRUE(alternative$summary$mcse_available)
  data.frame(model = reference$model, cmdstan_logml = reference$summary$logml,
    nutpie_logml = alternative$summary$logml, difference, combined_uncertainty = uncertainty,
    tolerance, mcse_available, precision_passed = is.finite(tolerance) && tolerance <= .25,
    passed = is.finite(difference) && is.finite(tolerance) &&
      mcse_available && tolerance <= .25 && abs(difference) <= tolerance)
}

nutpie_evidence_diagnostic_pass <- function(result) {
  d <- result$diagnostics$summary
  needed <- c("divergences", "depth_hits", "bfmi", "max_rhat", "min_bulk_ess", "min_tail_ess")
  all(needed %in% names(d)) && nrow(d) == 4L &&
    all(is.finite(as.matrix(d[needed]))) && all(d$divergences == 0) &&
    all(d$depth_hits == 0) && all(d$bfmi >= .3) && all(d$max_rhat <= 1.01) &&
    all(d$min_bulk_ess >= 100) && all(d$min_tail_ess >= 100)
}

nutpie_evidence_collect <- function(manifest, root, report_dir) {
  results <- list()
  status <- lapply(seq_len(nrow(manifest$tasks)), function(i) {
    task <- manifest$tasks[i, ]
    path <- file.path(root, "evidence", task$id, "result.rds")
    value <- if (file.exists(path)) readRDS(path) else list(status = "not_run", reason = "No execution artifact.")
    results[[task$id]] <<- value
    data.frame(task = task$id, stage = task$stage, status = value$status,
      reason = if (is.null(value$reason)) "" else value$reason, stringsAsFactors = FALSE)
  })
  status <- do.call(rbind, status)
  utils::write.csv(status, file.path(report_dir, "nutpie_evidence_status.csv"), row.names = FALSE)
  target_results <- results[paste0("target_", nutpie_evidence_pairs()$kernel)]
  checks <- lapply(target_results, `[[`, "checks")
  checks <- checks[vapply(checks, is.data.frame, logical(1))]
  if (length(checks)) utils::write.csv(do.call(rbind, checks),
    file.path(report_dir, "nutpie_evidence_targets.csv"), row.names = FALSE)
  complete <- all(status$status == "passed")
  protocol <- list(status = "incomplete", public_gate_changed = FALSE,
    shared_prior = c(B = 1/3, MC = 1/3, C = 1/3),
    diagnostic_thresholds = c(max_rhat = 1.01, min_bulk_ess = 100, min_tail_ess = 100,
      min_bfmi = .3, max_divergences = 0, max_depth_hits = 0),
    scope = "Four native kernels; same-B-data atomic B/MC/C universe; mixed backends remain unsupported.")
  if (complete) {
    models <- c("B", "MC", "C")
    paired <- lapply(models, function(model) nutpie_evidence_compare(
      results[[paste0("bridge_", model, "_cmdstanr")]], results[[paste0("bridge_", model, "_nutpie")]]))
    comparison <- do.call(rbind, paired)
    diagnostic_results <- results[manifest$tasks$id[manifest$tasks$stage == "bridge"]]
    diagnostic_pass <- all(vapply(diagnostic_results, nutpie_evidence_diagnostic_pass, logical(1)))
    universe <- lapply(c("cmdstanr", "nutpie"), function(backend) {
      fits <- stats::setNames(lapply(models, function(model) readRDS(
        nutpie_evidence_fit_path(root, model, backend, shared = TRUE))), models)
      antitrustBayes:::validate_comparison_contract(fits)
      source_evidence <- data.frame(model = models,
        logml = vapply(models, function(model) results[[paste0("bridge_", model, "_", backend)]]$summary$logml, numeric(1)),
        logml_sd = vapply(models, function(model) results[[paste0("bridge_", model, "_", backend)]]$summary$logml_sd, numeric(1)))
      descriptor <- antitrustBayes::bayes_model_universe(as.list(stats::setNames(models, models)), protocol$shared_prior)
      rows <- antitrustBayes:::.mu_fit_descriptor_rows(fits, descriptor)
      value <- antitrustBayes:::.mu_result_from_evidence(rows, source_evidence, descriptor$prior_weights)
      stopifnot(identical(value$status, "complete"), nrow(value$table) == 3L,
                  all(is.finite(value$table$posterior_probability)))
      value$table$backend <- backend
      value$table
    })
    universe <- do.call(rbind, universe)
    left <- universe[universe$backend == "cmdstanr", ]
    right <- universe[universe$backend == "nutpie", ]
    right <- right[match(left$hypothesis, right$hypothesis), ]
    ## The softmax derivative has l1 norm <= 1/2, so a uniform log-evidence
    ## error bound d gives a conservative probability-error bound d/2.
    probability_tolerance <- min(.05, max(comparison$tolerance) / 2 + 1e-6)
    probability_pass <- all(abs(left$posterior_probability - right$posterior_probability) <= probability_tolerance)
    protocol$status <- if (all(comparison$passed) && diagnostic_pass && probability_pass) "passed" else "inconclusive"
    protocol$bridge_comparison <- comparison
    protocol$diagnostics_passed <- diagnostic_pass
    protocol$probability_tolerance <- probability_tolerance
    protocol$probabilities_passed <- probability_pass
    protocol$universe <- universe
    utils::write.csv(comparison, file.path(report_dir, "nutpie_evidence_comparison.csv"), row.names = FALSE)
    utils::write.csv(universe, file.path(report_dir, "nutpie_evidence_universe.csv"), row.names = FALSE)
  }
  saveRDS(protocol, file.path(root, "evidence", "protocol.rds"))
  invisible(protocol)
}

nutpie_evidence_worker <- function(task_path) {
  if (!identical(Sys.getenv("AB_NUTPIE_EVIDENCE_EXECUTE", "false"), "true"))
    stop("Evidence workers require AB_NUTPIE_EVIDENCE_EXECUTE=true.")
  work <- readRDS(task_path)
  .libPaths(work$library_paths)
  pkgload::load_all(work$repository, quiet = TRUE)
  source(file.path(work$repository, "doc", "validation", "nutpie_backend_validation.R"), local = .GlobalEnv)
  options(antitrustBayes.cmdstan_cache = file.path(work$root, "cmdstan-cache"),
          antitrustBayes.nutpie_cache = file.path(work$root, "nutpie-cache"))
  cmdstan <- Sys.getenv("AB_NUTPIE_CMDSTAN", "")
  if (nzchar(cmdstan)) cmdstanr::set_cmdstan_path(cmdstan)
  backend_manifest <- readRDS(file.path(work$root, "manifest.rds"))
  nutpie_evidence_manifest_contract(backend_manifest)
  evidence_manifest <- readRDS(file.path(work$root, "evidence", "manifest.rds"))
  current_r_source <- vapply(names(evidence_manifest$r_source_hash), function(path)
    digest::digest(file = file.path(work$repository, path), algo = "sha256"), character(1))
  if (!identical(current_r_source, evidence_manifest$r_source_hash))
    stop("Evidence R source changed after preparing its immutable manifest.")
  nutpie_evidence_budget(readRDS(work$invocation_path))
  result <- tryCatch(switch(work$task$stage,
    target = nutpie_evidence_target(work$task, work$root, work$directory),
    reference = nutpie_evidence_reference(work$directory),
    sample = nutpie_evidence_sample(work$task, work$root, work$directory, backend_manifest),
    bridge = nutpie_evidence_bridge(work$task, work$root, work$directory)),
    error = function(e) list(status = "failed", reason = conditionMessage(e)))
  result$versions <- vapply(c("antitrustBayes", "cmdstanr", "nutpieR", "bridgestan", "bridgesampling"),
    function(p) if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else NA_character_, character(1))
  saveRDS(result, file.path(work$directory, "result.rds"))
}

nutpie_evidence_main <- function() {
  repository <- normalizePath(".", mustWork = TRUE)
  pkgload::load_all(repository, quiet = TRUE)
  source(file.path(repository, "doc", "validation", "nutpie_backend_validation.R"), local = .GlobalEnv)
  root <- Sys.getenv("AB_NUTPIE_OUTPUT", file.path(repository, ".scratch", "nutpie-validation"))
  if (!file.exists(file.path(root, "manifest.rds"))) stop("Prepare the backend benchmark manifest first.")
  root <- normalizePath(root)
  backend_manifest <- readRDS(file.path(root, "manifest.rds"))
  nutpie_evidence_manifest_contract(backend_manifest)
  directory <- file.path(root, "evidence")
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  manifest_path <- file.path(directory, "manifest.rds")
  identity <- list(script_hash = digest::digest(file = file.path(repository, "doc", "validation", "nutpie_evidence_validation.R"), algo = "sha256"),
    fixture_hash = backend_manifest$fixture_hash, source_hash = backend_manifest$source_hash,
    r_source_hash = vapply(c("R/bma.R", "R/estimate.R", "R/nutpie_backend.R", "R/sampler_backend.R", "R/model_universe.R",
      "doc/validation/nutpie_backend_validation.R"),
      function(path) digest::digest(file = file.path(repository, path), algo = "sha256"), character(1)),
    tasks = nutpie_evidence_tasks(), repetitions = 3L, max_cores = 2L,
    inference = "Four explicitly opted-in additional same-B-data fits only; four chains, 1000 warmup and retained draws each; all other fits reused.")
  if (file.exists(manifest_path)) {
    manifest <- readRDS(manifest_path)
    if (!identical(manifest, identity)) stop("Evidence script/design changed after preparing its immutable manifest.")
  } else { manifest <- identity; saveRDS(manifest, manifest_path) }
  report_dir <- file.path(repository, "doc", "validation")
  execute <- identical(Sys.getenv("AB_NUTPIE_EVIDENCE_EXECUTE", "false"), "true")
  if (!execute) {
    result <- nutpie_evidence_collect(manifest, root, report_dir)
    message("Prepared/collected evidence validation; no compilation or inference. Status: ", result$status)
    return(invisible(manifest))
  }
  if (!nzchar(Sys.which("timeout"))) stop("GNU timeout is required for the invocation deadline.")
  invocation_path <- nutpie_evidence_invocation(directory,
    seconds = suppressWarnings(as.numeric(Sys.getenv("AB_NUTPIE_EVIDENCE_SECONDS", "14400"))))
  retry <- identical(Sys.getenv("AB_NUTPIE_EVIDENCE_RETRY", "false"), "true")
  for (i in seq_len(nrow(manifest$tasks))) {
    task <- manifest$tasks[i, ]
    task_directory <- file.path(directory, task$id)
    replace_attempt <- FALSE
    if (file.exists(file.path(task_directory, "task.rds"))) {
      result_path <- file.path(task_directory, "result.rds")
      if (file.exists(result_path) && identical(readRDS(result_path)$status, "passed")) next
      if (!retry) next
      replace_attempt <- TRUE
    }
    if (task$stage == "sample" && !identical(Sys.getenv("AB_NUTPIE_EVIDENCE_SAMPLE", "false"), "true")) next
    if (task$stage == "sample" && !isTRUE(tryCatch({
      nutpie_evidence_require_targets(root); TRUE
    }, error = function(e) FALSE))) next
    if (task$stage == "target") {
      paths <- vapply(c("cmdstanr", "nutpie"), function(backend)
        nutpie_evidence_fit_path(root, task$model, backend), character(1))
      if (!all(vapply(paths, nutpie_evidence_fit_ready, logical(1)))) next
    }
    if (task$stage == "bridge") {
      paths <- vapply(c("B", "MC", "C"), function(model)
        nutpie_evidence_fit_path(root, model, task$backend, shared = TRUE), character(1))
      if (!all(vapply(paths, nutpie_evidence_fit_ready, logical(1)))) next
    }
    remaining <- nutpie_evidence_budget(readRDS(invocation_path))
    if (replace_attempt) {
      archive <- file.path(directory, "attempts")
      dir.create(archive, recursive = TRUE, showWarnings = FALSE)
      previous <- tempfile(paste0(task$id, "-"), tmpdir = archive)
      if (!file.rename(task_directory, previous)) stop("Cannot preserve previous task attempt: ", task$id)
    }
    dir.create(task_directory, recursive = TRUE, showWarnings = FALSE)
    task_path <- file.path(task_directory, "task.rds")
    saveRDS(list(task = task, repository = repository, root = root, directory = task_directory,
      library_paths = .libPaths(), invocation_path = invocation_path), task_path)
    code <- system2("timeout", c("--signal=TERM", "--kill-after=5s", as.character(remaining),
      shQuote(file.path(R.home("bin"), "Rscript")), "--vanilla",
      shQuote(file.path(repository, "doc", "validation", "nutpie_evidence_validation.R")), "--worker", shQuote(task_path)),
      stdout = file.path(task_directory, "worker.log"), stderr = file.path(task_directory, "worker.log"),
      env = c("OMP_NUM_THREADS=1", "OPENBLAS_NUM_THREADS=1", "MKL_NUM_THREADS=1",
        "CARGO_BUILD_JOBS=2", "MAKEFLAGS=-j2", "CMAKE_BUILD_PARALLEL_LEVEL=2"))
    if (code != 0L || !file.exists(file.path(task_directory, "result.rds")))
      saveRDS(list(status = if (code %in% c(124L, 137L)) "timed_out" else "worker_failed",
        reason = paste("Worker exit", code, "; partial artifacts retained; no automatic retry.")),
        file.path(task_directory, "result.rds"))
    nutpie_evidence_collect(manifest, root, report_dir)
  }
  nutpie_evidence_collect(manifest, root, report_dir)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 2L && identical(args[[1L]], "--worker")) nutpie_evidence_worker(args[[2L]]) else
    nutpie_evidence_main()
}
