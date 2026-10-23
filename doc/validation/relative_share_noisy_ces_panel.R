# Separate finite-observation-error CES M4T4 backend check.
# Sourcing defines functions only: no sampling, compilation, or file writes.
# After the other sampling workload has stopped, with compatible development
# source and the already compiled kernel caches selected:
# source('doc/validation/relative_share_noisy_ces_panel.R')
# relative_noisy_ces_panel_run(ORIGINAL_SAMPLING_ROOT, NEW_OUTPUT_ROOT)
# The original noiseless stress run is neither overwritten nor relabeled.
# max_seconds is checked between stages; the caller owns the external process
# timeout and the shared sampling budget. No compilation is requested here.

relative_noisy_ces_panel_arguments <- function(original_arguments) {
  args <- original_arguments
  required <- c('prices', 'quantities', 'margins', 'market', 'year', 'group')
  if (!all(required %in% names(args))) stop('Original CES panel arguments are incomplete.')
  n <- length(args$prices)
  stopifnot(n > 0L, all(vapply(args[required], length, integer(1)) == n),
    identical(args$demand_model, 'ces'), identical(args$share_basis, 'relative'),
    all(is.finite(args$quantities) & args$quantities > 0),
    all(is.finite(args$margins) & args$margins > 0 & args$margins < 1))
  market_index <- match(as.character(args$market), unique(as.character(args$market)))
  period_index <- match(as.character(args$year), unique(as.character(args$year)))
  cell <- interaction(market_index, period_index, drop = TRUE)
  stopifnot(length(unique(market_index)) == 4L, length(unique(period_index)) == 4L,
    length(levels(cell)) == 16L)
  product_index <- integer(n)
  for (index in split(seq_len(n), cell)) {
    stopifnot(!anyDuplicated(as.character(args$group[index])))
    product_index[index] <- seq_along(index)
  }
  share_error <- .015 * sin(product_index + 2 * period_index + market_index)
  margin_error <- .02 * cos(2 * product_index + period_index + market_index)
  observed_quantity <- args$quantities * exp(share_error)
  for (index in split(seq_len(n), cell))
    observed_quantity[index] <- observed_quantity[index] *
      sum(args$quantities[index]) / sum(observed_quantity[index])
  args$quantities <- observed_quantity
  args$margins <- original_arguments$margins * exp(margin_error)
  stopifnot(all(is.finite(args$quantities) & args$quantities > 0),
    all(is.finite(args$margins) & args$margins > 0 & args$margins < 1),
    identical(args$prices, original_arguments$prices),
    identical(args[setdiff(names(args), c('quantities', 'margins'))],
      original_arguments[setdiff(names(original_arguments), c('quantities', 'margins'))]))
  cells <- split(seq_len(n), cell)
  total_error <- max(vapply(cells, function(i)
    abs(sum(args$quantities[i]) - sum(original_arguments$quantities[i])), numeric(1)))
  stopifnot(total_error <= 1e-12 * max(original_arguments$quantities))
  list(arguments = args, observation_table = data.frame(row = seq_len(n),
    market = args$market, year = args$year, group = args$group,
    market_index, period_index, product_index,
    price = args$prices, original_quantity = original_arguments$quantities,
    observed_quantity = args$quantities, original_margin = original_arguments$margins,
    observed_margin = args$margins, prescribed_log_quantity_error = share_error,
    actual_log_quantity_error = log(args$quantities / original_arguments$quantities),
    prescribed_log_margin_error = margin_error), maximum_cell_total_error = total_error,
    specification = list(share_amplitude = .015, margin_amplitude = .02,
      product_index = 'One-based original observation order within each market-period cell.',
      market_period_index = 'One-based first appearance in the original argument vectors.',
      quantity_transform = 'q*exp(.015*sin(j+2*t+m)); renormalize to the original cell total.',
      margin_transform = 'margin*exp(.02*cos(2*j+t+m)).',
      interpretation = 'Prescribed deterministic observation errors; structural truth, prices, priors and pooling unchanged.'))
}

relative_noisy_ces_panel_source_status <- function(original, repository = '.') {
  expected <- original$runtime$source_hashes
  paths <- names(expected)
  if (!is.character(expected) || !length(expected) || is.null(paths) ||
      anyNA(expected) || anyNA(paths) || any(!nzchar(paths)) || anyDuplicated(paths) ||
      any(grepl('(^|/)\\.\\.(/|$)|^/|^~|^[A-Za-z]:', paths)))
    stop('Original runtime must contain named source hashes within the repository.')
  current <- vapply(file.path(repository, paths), function(path)
    if (file.exists(path)) digest::digest(file = path, algo = 'sha256') else NA_character_, character(1))
  data.frame(path = paths, original_sha256 = unname(expected), current_sha256 = unname(current),
    matches = !is.na(current) & unname(current) == unname(expected), stringsAsFactors = FALSE)
}

relative_noisy_ces_panel_original_state <- function(original_sampling_root) {
  backend <- c('cmdstanr', 'nutpie')
  paths <- file.path(original_sampling_root, paste0('ces_C_M4T4_', backend), 'fit.rds')
  present <- file.exists(paths)
  fits <- data.frame(backend, path = paths, fit_file_present = present,
    sha256 = vapply(paths, function(path) if (file.exists(path))
      digest::digest(file = path, algo = 'sha256') else NA_character_, character(1)))
  files <- list.files(file.path(original_sampling_root, 'ces_C_M4T4_cmdstanr'),
    pattern = '[.]csv$', recursive = TRUE, full.names = TRUE)
  inspect <- function(path) {
    before <- digest::digest(file = path, algo = 'sha256')
    lines <- readLines(path, warn = FALSE)
    if (!identical(before, digest::digest(file = path, algo = 'sha256')))
      stop('Original CmdStan CSV changed while being read; retry after that workload stops.')
    header <- which(grepl('^lp__(,|$)', trimws(lines)))
    numeric_rows <- 0L; incomplete_rows <- 0L
    if (length(header)) {
      width <- length(strsplit(lines[header[1L]], ',', fixed = TRUE)[[1L]])
      candidate <- lines[seq_along(lines) > header[1L] & nzchar(trimws(lines)) &
        !grepl('^#', trimws(lines))]
      candidate <- candidate[!grepl('^lp__(,|$)', trimws(candidate))]
      numeric_pattern <- '^[+-]?(([0-9]+[.]?[0-9]*|[.][0-9]+)([eE][+-]?[0-9]+)?|nan|inf(inity)?)$'
      valid <- vapply(strsplit(candidate, ',', fixed = TRUE), function(row)
        length(row) == width && all(grepl(numeric_pattern, trimws(row), ignore.case = TRUE)), logical(1))
      numeric_rows <- sum(valid); incomplete_rows <- sum(!valid)
    }
    data.frame(path, sha256 = before, sampler_csv = length(header) > 0L,
      numeric_rows = if (length(header)) numeric_rows else NA_integer_,
      incomplete_or_nonnumeric_rows = if (length(header)) incomplete_rows else NA_integer_)
  }
  csv <- if (length(files)) do.call(rbind, lapply(files, inspect)) else data.frame(
    path = character(), sha256 = character(), sampler_csv = logical(), numeric_rows = integer(),
    incomplete_or_nonnumeric_rows = integer())
  sampler <- csv[csv$sampler_csv, , drop = FALSE]
  rows <- if (nrow(sampler)) sum(sampler$numeric_rows) else NA_integer_
  note <- paste0(if (!any(present)) 'No complete retained fit is saved (no fit.rds files); ' else
      paste0(sum(present), ' saved-fit file(s) present, without certifying their completeness; '),
    if (is.na(rows)) 'no identified CmdStan sampler CSVs available for row counting' else
      paste0(rows, ' numeric row(s) in ', nrow(sampler), ' identified CmdStan sampler CSV file(s)'),
    '. Numeric rows describe saved CSV contents, including warmup if requested; they do not certify fit completion.')
  list(inspected_at = Sys.time(), saved_fit_files = fits, cmdstan_csv_files = csv,
    total_cmdstan_numeric_rows = rows, status_note = note)
}

relative_noisy_ces_panel_run <- function(original_sampling_root, output_root,
    max_seconds = 2400L, repository = '.') {
  stopifnot(length(max_seconds) == 1L, is.finite(max_seconds),
    max_seconds > 0, max_seconds <= 14400)
  repository <- normalizePath(repository, mustWork = TRUE)
  original_sampling_root <- normalizePath(original_sampling_root, mustWork = TRUE)
  original_manifest_path <- file.path(original_sampling_root, 'sampling-manifest.rds')
  original_manifest_sha256 <- digest::digest(file = original_manifest_path, algo = 'sha256')
  original <- readRDS(original_manifest_path)
  if (!identical(original_manifest_sha256, digest::digest(file = original_manifest_path, algo = 'sha256')))
    stop('Original sampling manifest changed while being read; retry after that workload stops.')
  original_source_status <- relative_noisy_ces_panel_source_status(original, repository)
  if (!all(original_source_status$matches)) stop('Original runtime source mismatch: ',
    paste(original_source_status$path[!original_source_status$matches], collapse = ', '))
  if (!is.list(original$fixtures$ces_C_M4T4)) stop('Original manifest lacks ces_C_M4T4 arguments.')
  prepared <- relative_noisy_ces_panel_arguments(original$fixtures$ces_C_M4T4)
  # Controls are held fixed, with separately preregistered seeds for this new
  # observation-error experiment and independent streams for the two backends.
  controls <- original$controls
  stopifnot(identical(as.integer(controls$chains), 4L), identical(as.integer(controls$iter), 1000L),
    identical(as.integer(controls$warmup), 500L), identical(as.integer(controls$cores), 4L),
    identical(controls$control$adapt_delta, .99), identical(as.integer(controls$control$max_treedepth), 13L))
  sampling_seeds <- c(cmdstanr = 453593L, nutpie = 453594L)
  bridge_seeds <- c(cmdstanr = 453693L, nutpie = 453694L)
  helper <- new.env(parent = environment())
  sys.source(file.path(repository, 'doc/validation/nutpie_backend_validation.R'), helper)
  sys.source(file.path(repository, 'doc/validation/relative_share_validation_collect.R'), helper)
  source_paths <- c(file.path(repository, 'inst/stan/ces_supply.stan'),
    file.path(repository, 'doc/validation/relative_share_noisy_ces_panel.R'),
    file.path(repository, 'doc/validation/relative_share_validation_collect.R'))
  source_hashes <- vapply(source_paths, digest::digest, character(1), file = TRUE, algo = 'sha256')
  names(source_hashes) <- substring(source_paths, nchar(repository) + 2L)
  kernel_hash <- unname(source_hashes['inst/stan/ces_supply.stan'])
  if (!identical(kernel_hash, unname(original$kernel_hashes['inst/stan/ces_supply.stan'])))
    stop('Original fixture manifest and current CES kernel differ; do not relabel targets.')
  if (file.exists(output_root)) stop('Select a new, nonexistent output root; existing attempts are preserved.')
  original_state <- relative_noisy_ces_panel_original_state(original_sampling_root)
  stopifnot(dir.create(output_root, recursive = TRUE, showWarnings = FALSE))
  output_root <- normalizePath(output_root, mustWork = TRUE)
  case <- 'ces_C_M4T4'
  calls <- lapply(names(sampling_seeds), function(backend) {
    args <- prepared$arguments; args$backend <- backend; args$stan_options <- controls
    args$stan_options$seed <- unname(sampling_seeds[backend])
    if (backend == 'cmdstanr') {
      args$stan_options$output_dir <- file.path(output_root, paste(case, backend, sep = '_'))
      args$stan_options$sig_figs <- 18L
    }
    args
  })
  names(calls) <- names(sampling_seeds)
  manifest <- list(experiment_id = 'ces_C_M4T4_finite_observation_error_v1',
    started = Sys.time(), deadline = Sys.time() + max_seconds, status = 'prepared',
    original_noiseless = list(root = original_sampling_root,
      manifest_path = original_manifest_path, manifest_sha256 = original_manifest_sha256,
      status_note = paste(original_state$status_note,
        'This is an additional observation-error experiment; the original result is not replaced.'),
      artifact_state = original_state,
      attempts = original$attempts[grep('^ces_C_M4T4_', names(original$attempts))]),
    original_runtime_source_verification = original_source_status,
    original_arguments = original$fixtures$ces_C_M4T4,
    original_arguments_sha256 = digest::digest(original$fixtures$ces_C_M4T4, algo = 'sha256'),
    fixtures = stats::setNames(list(prepared$arguments), case),
    observed_arguments_sha256 = digest::digest(prepared$arguments, algo = 'sha256'),
    observation_error = prepared$specification, maximum_cell_total_error = prepared$maximum_cell_total_error,
    controls = controls, sampling_seeds = sampling_seeds, bridge_seeds = bridge_seeds,
    bridge_repetitions = 3L, scheduled_calls = calls,
    kernel_hashes = stats::setNames(kernel_hash, 'inst/stan/ces_supply.stan'),
    validation_source_hashes = source_hashes,
    runtime = helper$nutpie_validation_runtime(repository), attempts = list(), warnings = list())
  manifest_path <- file.path(output_root, 'sampling-manifest.rds')
  save <- function() saveRDS(manifest, manifest_path)
  saveRDS(original, file.path(output_root, 'original-noiseless-manifest-snapshot.rds'))
  write.csv(original_source_status, file.path(output_root, 'original-source-verification.csv'), row.names = FALSE)
  write.csv(original_state$saved_fit_files, file.path(output_root, 'original-saved-fit-files.csv'), row.names = FALSE)
  write.csv(original_state$cmdstan_csv_files, file.path(output_root, 'original-cmdstan-csv-rows.csv'), row.names = FALSE)
  write.csv(prepared$observation_table, file.path(output_root, 'observation-errors.csv'), row.names = FALSE)
  manifest$observation_table_sha256 <- digest::digest(file = file.path(output_root, 'observation-errors.csv'), algo = 'sha256')
  save()
  active <- NULL; ended <- FALSE
  on.exit({
    if (!ended) {
      manifest$status <- 'interrupted'; manifest$ended <- Sys.time()
      if (!is.null(active)) manifest$attempts[[active]]$interrupted <- TRUE
      save()
    }
  }, add = TRUE)
  keep <- function(expr, backend, stage) tryCatch(withCallingHandlers(expr, warning = function(w) {
    manifest$warnings[[length(manifest$warnings) + 1L]] <<- list(backend = backend, stage = stage,
      at = Sys.time(), message = conditionMessage(w)); save(); invokeRestart('muffleWarning')
  }), error = identity)
  fits <- parameters <- diagnostics <- evidence <- evidence_summary <- evidence_attempts <- list()
  for (backend in names(calls)) {
    if (Sys.time() > manifest$deadline) {manifest$status <- 'budget_exhausted'; break}
    active <- paste(case, backend, sep = '_'); directory <- file.path(output_root, active)
    stopifnot(dir.create(directory, showWarnings = FALSE))
    manifest$attempts[[active]] <- list(status = 'fitting', started = Sys.time(), arguments = calls[[backend]])
    manifest$status <- 'running'; save()
    cat('FIT finite-observation-error', active, format(Sys.time()), '\n'); flush.console()
    started <- proc.time()[['elapsed']]
    fit <- keep(do.call(antitrustBayes::bayes_est, calls[[backend]]), backend, 'sampling')
    manifest$attempts[[active]]$fit_elapsed_seconds <- proc.time()[['elapsed']] - started
    if (inherits(fit, 'error')) {
      saveRDS(fit, file.path(directory, 'fit-error.rds'))
      manifest$attempts[[active]]$status <- 'failed'; manifest$attempts[[active]]$reason <- conditionMessage(fit)
      save(); next
    }
    path <- file.path(directory, 'fit.rds'); saveRDS(fit, path)
    manifest$attempts[[active]]$fit_sha256 <- digest::digest(file = path, algo = 'sha256')
    identity <- helper$nutpie_validation_identity(fit)
    data_fields <- fit@metadata$stan_data_fields
    observed_verified <- all(vapply(c('prices', 'quantities', 'margins'), function(field)
      identical(as.numeric(fit@data[[field]]), as.numeric(prepared$arguments[[field]])), logical(1))) &&
      all(vapply(c('market', 'year', 'group'), function(field)
        identical(as.character(fit@data[[field]]), as.character(prepared$arguments[[field]])), logical(1)))
    verified <- identical(fit@metadata$backend, backend) && identical(fit@metadata$stan_kernel, 'ces_supply') &&
      identical(fit@metadata$stan_source_hash, kernel_hash) && length(data_fields) > 0L &&
      identical(fit@metadata$stan_data_hash, digest::digest(fit@stan_data[data_fields], algo = 'sha256')) && observed_verified
    manifest$attempts[[active]]$identity <- identity
    manifest$attempts[[active]]$identity_verified <- verified
    manifest$attempts[[active]]$observed_arguments_verified <- observed_verified
    manifest$attempts[[active]]$status <- 'fit_complete'; save()
    getFromNamespace('validate_fit_spec', 'antitrustBayes')(fit)
    if (!verified) stop('Fresh fit failed source/data/backend identity checks: ', backend)
    fits[[backend]] <- fit
    p <- keep(helper$nutpie_validation_parameters(fit), backend, 'parameter_summary')
    d <- keep(helper$nutpie_validation_diagnostics(fit), backend, 'diagnostics')
    saveRDS(p, file.path(directory, 'parameters.rds')); saveRDS(d, file.path(directory, 'diagnostics.rds'))
    if (inherits(p, 'error')) manifest$attempts[[active]]$parameter_error <- conditionMessage(p)
    if (inherits(d, 'error')) manifest$attempts[[active]]$diagnostic_error <- conditionMessage(d)
    if (!inherits(p, 'error')) parameters[[backend]] <- p$summary
    if (!inherits(d, 'error')) diagnostics[[backend]] <- d$summary
    if (Sys.time() <= manifest$deadline) {
      cat('BRIDGE finite-observation-error', active, format(Sys.time()), '\n'); flush.console()
      set.seed(unname(bridge_seeds[backend])); started <- proc.time()[['elapsed']]
      ev <- keep(antitrustBayes::bayes_compare(candidate = fit, bridge = TRUE,
        bridge_repetitions = 3L, bridge_method = 'normal'), backend, 'bridge')
      manifest$attempts[[active]]$bridge_elapsed_seconds <- proc.time()[['elapsed']] - started
    } else ev <- simpleError('Execution budget exhausted before bridge sampling; no evidence fabricated.')
    saveRDS(ev, file.path(directory, 'evidence.rds')); evidence[[backend]] <- ev
    if (inherits(ev, 'error')) manifest$attempts[[active]]$evidence_error <- conditionMessage(ev)
    parsed <- helper$.relative_collect_evidence(ev, case, backend, requested = 3L)
    fields <- c('backend', 'share_basis', 'share_observation_measure', 'inference', 'stan_kernel',
      'stan_source_hash', 'stan_data_hash', 'model_spec_version', 'prior_manifest', 'resolved_priors')
    parsed$summary$identity_verified <- verified && isTRUE(parsed$summary$saved_summary_consistent) &&
      isTRUE(parsed$summary$saved_backend_verified) &&
      identical(parsed$summary$context_kernel, 'ces_supply') &&
      identical(parsed$summary$context_saved_hash, kernel_hash) &&
      identical(parsed$summary$context_current_hash, kernel_hash) &&
      all(vapply(fields, function(field) !is.null(parsed$metadata[[field]]) &&
        identical(parsed$metadata[[field]], fit@metadata[[field]]), logical(1)))
    evidence_summary[[backend]] <- parsed$summary; evidence_attempts[[backend]] <- parsed$attempts
    manifest$attempts[[active]]$status <- 'attempt_complete'; save()
  }
  tables <- list(parameters = helper$.relative_collect_bind(lapply(names(parameters), function(b)
      cbind(case = case, backend = b, parameters[[b]]))),
    diagnostics = helper$.relative_collect_bind(lapply(names(diagnostics), function(b)
      cbind(case = case, backend = b, diagnostics[[b]]))),
    evidence = helper$.relative_collect_bind(evidence_summary),
    bridge_attempts = helper$.relative_collect_bind(evidence_attempts))
  paired <- all(c('cmdstanr', 'nutpie') %in% names(fits))
  identity_equal <- paired && identical(helper$nutpie_validation_identity(fits$cmdstanr),
    helper$nutpie_validation_identity(fits$nutpie))
  expected_sources <- c(manifest$runtime$source_hashes, source_hashes)
  current_source_verified <- identical(unname(vapply(file.path(repository, names(expected_sources)),
    digest::digest, character(1), file = TRUE, algo = 'sha256')), unname(expected_sources))
  identity_equal <- isTRUE(identity_equal) && current_source_verified
  tables$paired_identity <- data.frame(case, both_fits_present = paired,
    current_source_verified, source_data_prior_identity_equal = identity_equal)
  tables$attempts <- do.call(rbind, lapply(names(calls), function(backend) {
    attempt <- manifest$attempts[[paste(case, backend, sep = '_')]]
    data.frame(backend, status = helper$.relative_collect_scalar(attempt$status, 'not_attempted'),
      fit_elapsed_seconds = helper$.relative_collect_scalar(attempt$fit_elapsed_seconds),
      bridge_elapsed_seconds = helper$.relative_collect_scalar(attempt$bridge_elapsed_seconds),
      reason = helper$.relative_collect_scalar(attempt$reason, ''),
      parameter_error = helper$.relative_collect_scalar(attempt$parameter_error, ''),
      diagnostic_error = helper$.relative_collect_scalar(attempt$diagnostic_error, ''),
      evidence_error = helper$.relative_collect_scalar(attempt$evidence_error, ''))
  }))
  if (all(c('cmdstanr', 'nutpie') %in% names(parameters)))
    tables$parameter_comparisons <- helper$.relative_collect_compare(parameters$cmdstanr,
      parameters$nutpie, identity_equal, helper$nutpie_validation_compare)
  if (all(c('cmdstanr', 'nutpie') %in% names(evidence_summary)))
    tables$evidence_comparison <- helper$.relative_collect_evidence_compare(evidence_summary$cmdstanr,
      evidence_summary$nutpie, identity_equal && all(vapply(evidence_summary,
        function(x) isTRUE(x$identity_verified), logical(1))))
  if (all(c('cmdstanr', 'nutpie') %in% names(diagnostics))) {
    fields <- c('divergences', 'depth_hits', 'max_rhat', 'min_bulk_ess', 'min_tail_ess')
    aggregate_diagnostic <- function(x, field) {
      value <- x[[field]]
      if (!all(is.finite(value))) return(NA_real_)
      if (field %in% c('divergences', 'depth_hits')) sum(value) else
        if (field == 'max_rhat') max(value) else min(value)
    }
    tables$diagnostic_comparison <- do.call(rbind, lapply(fields, function(field) {
      reference <- aggregate_diagnostic(diagnostics$cmdstanr, field)
      alternative <- aggregate_diagnostic(diagnostics$nutpie, field)
      data.frame(metric = field, cmdstanr = reference, nutpie = alternative,
        nutpie_minus_cmdstanr = alternative - reference)
    }))
  }
  for (name in names(tables)) write.csv(tables[[name]], file.path(output_root,
    paste0('noisy_ces_panel_', name, '.csv')), row.names = FALSE)
  result <- list(created = Sys.time(), experiment_id = manifest$experiment_id, tables = tables,
    uncertainty_definition = paste('Uses the reviewed relative-share collector helpers.',
      'Parameter tolerances are four combined chain-aware MCSEs.',
      'Evidence uncertainty combines within-repetition MCSE and repetition SD without dividing by sqrt(3).',
      'Diagnostics are reported without permissions or sampler-quality gates.'),
    original_noiseless_result = manifest$original_noiseless$status_note)
  saveRDS(result, file.path(output_root, 'noisy-ces-panel-results.rds'))
  manifest$current_source_verified <- current_source_verified
  manifest$status <- if (length(manifest$attempts) == 2L &&
    all(vapply(manifest$attempts, function(x) identical(x$status, 'attempt_complete'), logical(1))))
      'completed_attempts' else if (identical(manifest$status, 'budget_exhausted')) 'budget_exhausted' else 'incomplete_attempts'
  manifest$ended <- Sys.time(); save(); ended <- TRUE; active <- NULL
  invisible(result)
}
