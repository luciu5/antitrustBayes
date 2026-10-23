#!/usr/bin/env Rscript
# Read-only postprocessing of relative_share_backend_validation.R artifacts.
# Sourcing defines functions only. With compatible libraries selected:
# pkgload::load_all('.')
# source('doc/validation/relative_share_validation_collect.R')
# relative_validation_collect(sampling_root, report_dir)
# No sampler, compiler, bridge estimator, or policy solver is called here.

.relative_collect_bind <- function(rows) {
  rows <- Filter(function(x) is.data.frame(x) && nrow(x), rows)
  if (!length(rows)) return(data.frame())
  fields <- unique(unlist(lapply(rows, names)))
  do.call(rbind, lapply(rows, function(x) {
    x[setdiff(fields, names(x))] <- NA
    x[fields]
  }))
}

.relative_collect_scalar <- function(x, default = NA_real_) {
  if (length(x) == 1L && !is.na(x)) x else default
}

.relative_collect_evidence <- function(value, case, backend, requested = 3L) {
  failure <- if (inherits(value, 'error')) conditionMessage(value) else ''
  context <- if (is.data.frame(value)) attr(value, 'bridge_context_info') else NULL
  context <- if (length(context) == 1L) context[[1L]] else NULL
  metadata <- if (is.data.frame(value)) attr(value, 'backend_metadata') else NULL
  metadata <- if (length(metadata) == 1L) metadata[[1L]] else NULL
  attempts <- context$repetitions
  if (is.null(attempts)) attempts <- data.frame()
  if (nrow(attempts) && (anyDuplicated(attempts$repetition) ||
      any(!is.finite(attempts$repetition)))) stop('Invalid saved bridge repetition indices.')
  ids <- sort(unique(c(seq_len(requested), attempts$repetition)))
  rows <- lapply(ids, function(i) {
    row <- attempts[attempts$repetition == i, , drop = FALSE]
    if (!nrow(row)) row <- data.frame(repetition = i,
      status = if (nzchar(failure)) 'saved_error' else 'missing_attempt_record',
      logml = NA_real_, mcse_logml = NA_real_, error = failure, warnings = '')
    cbind(case = case, backend = backend, row)
  })
  attempts <- .relative_collect_bind(rows)
  finite <- is.finite(attempts$logml)
  available <- finite & attempts$status == 'available'
  estimate <- if (any(available)) stats::median(attempts$logml[available]) else NA_real_
  repetition_sd <- if (sum(available) > 1L) stats::sd(attempts$logml[available]) else NA_real_
  within <- attempts$mcse_logml[available]
  # Repetitions reuse posterior draws: do not divide their MCSE by sqrt(3).
  within <- if (length(within) && all(is.finite(within) & within >= 0)) max(within) else NA_real_
  reported_logml <- if (is.data.frame(value)) .relative_collect_scalar(value$logml) else NA_real_
  reported_sd <- if (is.data.frame(value)) .relative_collect_scalar(value$logml_sd) else NA_real_
  same_number <- function(x, y) (is.na(x) && is.na(y)) ||
    (is.finite(x) && is.finite(y) && abs(x - y) <= 1e-12 * max(1, abs(x), abs(y)))
  consistent <- is.data.frame(value) && nrow(value) == 1L &&
    same_number(estimate, reported_logml) && same_number(repetition_sd, reported_sd) &&
    identical(as.integer(value$completed_repetitions), as.integer(sum(available))) &&
    identical(as.integer(value$requested_repetitions), as.integer(requested)) &&
    nrow(attempts) == requested && !any(finite & !available)
  summary <- data.frame(case, backend, logml = estimate, repetition_sd,
    within_mcse = within, uncertainty = sqrt(within^2 + repetition_sd^2),
    requested_repetitions = requested, completed_repetitions = sum(available),
    numerical_status = if (all(available)) 'complete' else if (any(available)) 'partial' else 'unavailable',
    saved_summary_consistent = consistent,
    context_source = .relative_collect_scalar(context$source, NA_character_),
    context_kernel = .relative_collect_scalar(context$kernel, NA_character_),
    context_saved_hash = .relative_collect_scalar(context$saved_hash, NA_character_),
    context_current_hash = .relative_collect_scalar(context$current_hash, NA_character_),
    saved_backend_verified = is.data.frame(value) && identical(as.character(value$backend), backend) &&
      identical(context$backend, backend),
    reason = failure, stringsAsFactors = FALSE)
  list(summary = summary, attempts = attempts, context = context, metadata = metadata)
}

.relative_collect_compare <- function(reference, alternative, identity_equal, compare) {
  result <- compare(reference, alternative, identity_equal)
  result$mcse_available <- is.finite(result$combined_mcse) & result$combined_mcse > 0
  result$deterministic_match <- FALSE
  result$reason <- ''
  for (i in seq_len(nrow(result))) {
    r <- reference[match(result$variable[i], reference$variable), , drop = FALSE]
    a <- alternative[match(result$variable[i], alternative$variable), , drop = FALSE]
    constant <- isTRUE(r$sd == 0) && isTRUE(a$sd == 0)
    indicator <- result$metric[i] %in% c('prob_gt_0', 'prob_gt_5') &&
      isTRUE(result$reference[i] %in% c(0, 1)) && isTRUE(result$alternative[i] %in% c(0, 1))
    equal <- isTRUE(result$difference[i] == 0)
    if (!isTRUE(identity_equal)) next
    if (!is.finite(result$reference[i]) || !is.finite(result$alternative[i])) {
      result$status[i] <- 'unavailable_estimate'
    } else if (constant && equal) {
      result$status[i] <- 'deterministic_match'
      result$deterministic_match[i] <- TRUE
      result$reason[i] <- 'Both saved draw distributions are constant and identical.'
    } else if (constant && !equal) {
      result$status[i] <- 'deterministic_mismatch'
      result$reason[i] <- 'Both saved draw distributions are constant and differ.'
    } else if (indicator && equal && !result$mcse_available[i]) {
      result$status[i] <- 'matching_constant_indicators_mcse_unavailable'
      result$reason[i] <- 'Matching empirical constant indicators do not establish population tail equivalence.'
    } else if (!result$mcse_available[i]) {
      result$status[i] <- 'unavailable_mcse'
      result$reason[i] <- 'No finite positive combined chain-aware MCSE; no equivalence claim.'
    }
  }
  result
}

.relative_collect_evidence_compare <- function(reference, alternative, identity_equal) {
  difference <- alternative$logml - reference$logml
  combined <- sqrt(reference$uncertainty^2 + alternative$uncertainty^2)
  data.frame(reference = reference$logml, alternative = alternative$logml,
    difference, reference_minus_alternative = -difference, combined_mcse = combined,
    tolerance = 4 * combined,
    reference_minus_alternative_combined_mcse = if (is.finite(combined) && combined > 0) -difference / combined else NA_real_,
    identity_equal,
    reference_completed_repetitions = reference$completed_repetitions,
    alternative_completed_repetitions = alternative$completed_repetitions,
    status = if (!isTRUE(identity_equal)) 'identity_mismatch' else
      if (!is.finite(difference)) 'unavailable_evidence' else
      if (!is.finite(combined) || combined <= 0) 'unavailable_mcse' else
      if (abs(difference) <= 4 * combined) 'within_mcse_tolerance' else 'outside_mcse_tolerance')
}

.relative_collect_source_status <- function(result, current_source_verified) {
  result$historical_identity_equal <- result$identity_equal
  result$historical_status <- result$status
  result$current_source_verified <- current_source_verified
  result$identity_equal <- result$historical_identity_equal & current_source_verified
  if (!current_source_verified) result$status[result$historical_identity_equal] <- 'historical_only_source_mismatch'
  result
}

.relative_collect_policy_identity <- function(fit, policy) {
  antitrustBayes:::.lq_validate_plan(fit, policy$plan)
  draws <- policy$draws; plan <- policy$plan
  n <- length(plan$draw_indices)
  stopifnot(is.data.frame(draws), is.data.frame(policy$summary),
    !anyDuplicated(draws[c('regime', 'draw_position')]),
    setequal(unique(draws$regime), c('merger', 'tariff25')),
    nrow(draws) == 2L * n,
    all(draws$draw_position %in% seq_len(n)),
    identical(as.integer(draws$draw_index), as.integer(plan$draw_indices[draws$draw_position])))
  paired <- intersect(draws$draw_position[draws$regime == 'merger' & draws$status == 'ok'],
    draws$draw_position[draws$regime == 'tariff25' & draws$status == 'ok'])
  ok <- draws$status == 'ok'
  stopifnot(all(is.finite(draws$price_pct[ok])),
    all(is.finite(draws$cost_identity_error[ok]) & draws$cost_identity_error[ok] < 1e-8),
    all(draws$plan_fingerprint[ok] == plan$fit_signature$latent_state_fingerprint))
  quality_hashes <- vapply(seq_len(n), function(i)
    digest::digest(as.numeric(plan$quality_draws[i, ]), algo = 'sha256'), character(1))
  stopifnot(identical(as.character(draws$quality_hash[ok]),
    unname(quality_hashes[draws$draw_position[ok]])))
  for (regime in c('merger', 'tariff25')) {
    summary <- policy$summary[policy$summary$regime == regime, , drop = FALSE]
    values <- draws$price_pct[draws$regime == regime & draws$draw_position %in% paired]
    quantiles <- if (length(values)) as.numeric(stats::quantile(values, c(.025, .5, .975))) else rep(NA_real_, 3L)
    expected <- c(mean = if (length(values)) mean(values) else NA_real_,
      sd = if (length(values) > 1L) stats::sd(values) else NA_real_,
      q025 = quantiles[1L], median = quantiles[2L], q975 = quantiles[3L],
      prob_gt_0 = if (length(values)) mean(values > 0) else NA_real_,
      prob_gt_5 = if (length(values)) mean(values > 5) else NA_real_)
    stopifnot(nrow(summary) == 1L, summary$n_planned == n, summary$n_paired == length(paired),
      summary$n_valid == sum(draws$regime == regime & ok),
      summary$n_policy_failures == sum(draws$regime == regime & !ok),
      isTRUE(all.equal(as.numeric(summary[names(expected)]), as.numeric(expected), tolerance = 1e-12)))
    if (length(paired) != n) stopifnot(all(is.na(unlist(summary[grep('^mcse_', names(summary))]))))
  }
  data.frame(plan_fingerprint = plan$fit_signature$latent_state_fingerprint,
    posterior_draw_count = n, n_paired = length(paired),
    attrition_fraction = 1 - length(paired) / n,
    max_cost_identity_error = if (any(ok)) max(draws$cost_identity_error[ok]) else NA_real_,
    plan_and_saved_summary_verified = TRUE)
}

.relative_collect_policy_conditional_mcse <- function(policy, draw_shape) {
  n <- prod(draw_shape)
  stopifnot(length(draw_shape) == 2L, all(draw_shape > 0),
    setequal(policy$plan$draw_indices, seq_len(n)))
  draws <- policy$draws
  paired <- intersect(draws$draw_position[draws$regime == 'merger' & draws$status == 'ok'],
    draws$draw_position[draws$regime == 'tariff25' & draws$status == 'ok'])
  p <- length(paired) / n
  rows <- lapply(c('merger', 'tariff25'), function(regime) {
    selected <- draws$regime == regime & draws$draw_position %in% paired
    index <- draws$draw_index[selected]; values <- draws$price_pct[selected]
    stopifnot(length(index) == length(paired), !anyDuplicated(index),
      all(index %in% seq_len(n)), all(is.finite(values)))
    indicator <- numeric(n); indicator[index] <- 1
    chain_mcse <- function(x) suppressWarnings(posterior::mcse_mean(
      matrix(x, nrow = draw_shape[1L], ncol = draw_shape[2L])))
    ratio_mcse <- function(influence) {
      # Invalid outcomes are never imputed: their ratio numerator influence
      # is zero by definition. Every original transition remains in its chain.
      if (length(values) < 2L || any(!is.finite(influence)) ||
          all(influence == influence[1L])) return(NA_real_)
      chain_mcse(influence) / p
    }
    mean_value <- if (length(values)) mean(values) else NA_real_
    mean_influence <- numeric(n); mean_influence[index] <- values - mean_value
    mean_mcse <- ratio_mcse(mean_influence)
    probabilities <- vapply(c(0, 5), function(threshold) {
      if (length(values)) mean(values > threshold) else NA_real_
    }, numeric(1))
    probability_mcse <- vapply(seq_along(probabilities), function(i) {
      influence <- numeric(n)
      influence[index] <- 1 * (values > c(0, 5)[i]) - probabilities[i]
      ratio_mcse(influence)
    }, numeric(1))
    sd_mcse <- NA_real_
    if (length(values) > 1L) {
      variance <- mean((values - mean_value)^2)
      if (is.finite(variance) && variance > 0) {
        correction <- sqrt(length(values) / (length(values) - 1))
        sd_value <- correction * sqrt(variance)
        influence <- -sd_value * (indicator - p) / (2 * (length(values) - 1))
        influence[index] <- influence[index] + correction *
          ((values - mean_value)^2 - variance) / (2 * sqrt(variance))
        sd_mcse <- ratio_mcse(influence)
      }
    }
    ordinary_mean_mcse <- NA_real_; all_valid_agrees <- NA
    if (length(values) == n) {
      original_order <- numeric(n); original_order[index] <- values
      ordinary_mean_mcse <- chain_mcse(original_order)
      all_valid_agrees <- (is.na(mean_mcse) && is.na(ordinary_mean_mcse)) ||
        (is.finite(mean_mcse) && is.finite(ordinary_mean_mcse) &&
          abs(mean_mcse - ordinary_mean_mcse) <= 1e-12 * max(1, abs(ordinary_mean_mcse)))
      stopifnot(all_valid_agrees)
    }
    data.frame(regime, estimand = 'conditional_on_paired_validity', n_planned = n,
      n_paired = length(values), paired_validity_probability = p,
      mcse_paired_validity_probability = chain_mcse(indicator),
      conditional_mean = mean_value, ratio_mcse_mean = mean_mcse, ratio_mcse_sd = sd_mcse,
      conditional_prob_gt_0 = probabilities[1L], conditional_prob_gt_5 = probabilities[2L],
      ratio_mcse_prob_gt_0 = probability_mcse[1L], ratio_mcse_prob_gt_5 = probability_mcse[2L],
      ordinary_all_valid_mcse_mean = ordinary_mean_mcse,
      all_valid_mean_mcse_agrees = all_valid_agrees,
      mcse_method = 'ratio_influence_function_original_chain_order',
      reason = paste('First-order MCMC ratio MCSE, conditional on paired validity; requires usual CLT and finite moments.',
        'Sample-SD delta influence includes the Bessel-correction derivative and requires finite fourth moments.',
        'Constant or insufficient influences have unavailable MCSE. Quantile MCSE is not inferred.'))
  })
  .relative_collect_bind(rows)
}

# Small independent estimator checks; no fits, solvers or bridge estimates.
relative_validation_collect_checks <- function() {
  if (!requireNamespace('posterior', quietly = TRUE)) stop('Checks require posterior.')
  had_seed <- exists('.Random.seed', envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get('.Random.seed', envir = .GlobalEnv)
  on.exit(if (had_seed) assign('.Random.seed', old_seed, envir = .GlobalEnv) else
    if (exists('.Random.seed', envir = .GlobalEnv, inherits = FALSE))
      rm('.Random.seed', envir = .GlobalEnv), add = TRUE)
  set.seed(419L)
  shape <- c(128L, 4L); n <- prod(shape)
  y <- 4 + 2 * as.numeric(apply(matrix(stats::rnorm(n), nrow = shape[1L]), 2L,
    function(z) as.numeric(stats::filter(z, .65, method = 'recursive'))))
  make_policy <- function(valid) {
    permutation <- sample.int(n)
    rows <- expand.grid(regime = c('merger', 'tariff25'), draw_position = seq_len(n), stringsAsFactors = FALSE)
    rows$draw_index <- permutation[rows$draw_position]
    rows$status <- ifelse(valid[rows$draw_index], 'ok', 'failed')
    rows$price_pct <- y[rows$draw_index]; rows$price_pct[rows$status != 'ok'] <- NA_real_
    list(draws = rows, plan = list(draw_indices = permutation))
  }
  complete <- .relative_collect_policy_conditional_mcse(make_policy(rep(TRUE, n)), shape)
  ordinary <- posterior::mcse_mean(matrix(y, nrow = shape[1L]))
  stopifnot(all(complete$all_valid_mean_mcse_agrees),
    all(abs(complete$ratio_mcse_mean - ordinary) < 1e-12))
  valid <- rep(TRUE, n); valid[c(3, 10, 99, 132, 257, 400)] <- FALSE
  p <- mean(valid); mu <- mean(y[valid]); h <- numeric(n); h[valid] <- y[valid] - mu
  partial <- .relative_collect_policy_conditional_mcse(make_policy(valid), shape)
  stopifnot(all(partial$n_paired == sum(valid)), all(partial$estimand == 'conditional_on_paired_validity'),
    all(abs(partial$ratio_mcse_mean - posterior::mcse_mean(matrix(h, nrow = shape[1L])) / p) < 1e-12))
  for (threshold in c(0, 5)) {
    h <- numeric(n); h[valid] <- 1 * (y[valid] > threshold) - mean(y[valid] > threshold)
    expected <- posterior::mcse_mean(matrix(h, nrow = shape[1L])) / p
    stopifnot(all(abs(partial[[paste0('ratio_mcse_prob_gt_', threshold)]] - expected) < 1e-12))
  }
  # Differentiate the statistic in independent raw ratio moments, rather than
  # copying its conditional-variance influence implementation.
  raw <- cbind(I = 1 * valid, IY = ifelse(valid, y, 0), IY2 = ifelse(valid, y^2, 0))
  center <- colMeans(raw)
  statistic <- function(z) sqrt((z[3L] / z[1L] - (z[2L] / z[1L])^2) * n * z[1L] / (n * z[1L] - 1))
  gradient <- vapply(1:3, function(i) {
    up <- down <- center; up[i] <- up[i] + 1e-6; down[i] <- down[i] - 1e-6
    (statistic(up) - statistic(down)) / 2e-6
  }, numeric(1))
  influence <- as.numeric(sweep(raw, 2L, center) %*% gradient)
  expected_sd <- posterior::mcse_mean(matrix(influence, nrow = shape[1L]))
  stopifnot(all(abs(partial$ratio_mcse_sd - expected_sd) < 1e-8))
  historical <- data.frame(identity_equal = TRUE, status = 'within_mcse_tolerance', reference = -10, alternative = -9.99)
  stale <- .relative_collect_source_status(historical, FALSE)
  stopifnot(stale$historical_identity_equal, !stale$identity_equal,
    stale$status == 'historical_only_source_mismatch', stale$historical_status == historical$status,
    stale$reference == historical$reference, stale$alternative == historical$alternative)
  message('Collector checks passed: original-chain ratio MCSE, ordinary all-valid equality, independent SD delta derivative, historical source gate.')
  invisible(TRUE)
}

.relative_collect_weighted_policy <- function(policies, weights, scenario) {
  rows <- lapply(c('merger', 'tariff25'), function(regime) {
    values <- mass <- numeric(); components <- list()
    for (i in seq_along(weights)) {
      policy <- policies[[i]]
      if (is.null(policy) || inherits(policy, 'error')) return(data.frame(
        scenario, regime, status = 'unavailable_component_policy', retained_mass = NA_real_,
        failure_mass = NA_real_, reason = 'Every original positive-prior k policy component is required.'))
      draws <- policy$draws
      if (anyDuplicated(draws[c('regime', 'draw_position')])) stop('Duplicate conditional policy draw records.')
      paired <- intersect(draws$draw_position[draws$regime == 'merger' & draws$status == 'ok'],
        draws$draw_position[draws$regime == 'tariff25' & draws$status == 'ok'])
      keep <- draws$regime == regime & draws$draw_position %in% paired
      x <- draws$price_pct[keep]
      n <- length(policy$plan$draw_indices)
      if (!n || any(!is.finite(x)) || length(x) > n) stop('Invalid saved conditional policy support.')
      values <- c(values, x)
      # Failed mass remains missing; only then condition the empirical measure.
      mass <- c(mass, rep(weights[[i]] / n, length(x)))
      components[[i]] <- data.frame(scenario, regime, k = as.integer(names(weights)[i]),
        posterior_k = weights[[i]], n_planned = n, n_paired = length(x),
        retained_mass = weights[[i]] * length(x) / n,
        failure_mass = weights[[i]] * (1 - length(x) / n))
    }
    retained <- sum(mass)
    mean_value <- sd_value <- NA_real_; quantiles <- rep(NA_real_, 3L)
    prob0 <- prob5 <- NA_real_
    if (retained > 0) {
      mass <- mass / retained
      mean_value <- sum(mass * values)
      sd_value <- sqrt(sum(mass * (values - mean_value)^2))
      order <- order(values); cumulative <- cumsum(mass[order]); cumulative[length(cumulative)] <- 1
      quantiles <- vapply(c(.025, .5, .975), function(p) values[order][which(cumulative >= p)[1L]], numeric(1))
      prob0 <- sum(mass[values > 0]); prob5 <- sum(mass[values > 5])
    }
    result <- data.frame(scenario, regime,
      status = if (retained >= 1 - 1e-12) 'complete_empirical_mixture' else 'conditional_valid_empirical_mixture',
      retained_mass = retained, failure_mass = max(0, 1 - retained),
      n_paired = length(values), mean = mean_value, sd = sd_value, median = quantiles[2L],
      q025 = quantiles[1L], q975 = quantiles[3L], prob_gt_0 = prob0, prob_gt_5 = prob5,
      mcse_mean = NA_real_, mcse_sd = NA_real_, mcse_median = NA_real_,
      mcse_q025 = NA_real_, mcse_q975 = NA_real_, mcse_prob_gt_0 = NA_real_, mcse_prob_gt_5 = NA_real_,
      reason = 'Exact weighted empirical distribution conditional on paired validity; inverse-CDF quantiles. No joint-chain mixture MCSE or bridge-weight uncertainty propagation.')
    attr(result, 'components') <- .relative_collect_bind(components)
    result
  })
  list(summary = .relative_collect_bind(rows),
    attrition = .relative_collect_bind(lapply(rows, attr, 'components')))
}

relative_validation_collect <- function(root, report_dir, repository = '.',
    model_prior = stats::setNames(rep(1 / 7, 7L),
      paste0('logit_', c('MC', 'B', 'C', 'BF', 'CF', 'BS', 'CS')))) {
  root <- normalizePath(root, mustWork = TRUE)
  repository <- normalizePath(repository, mustWork = TRUE)
  hash <- function(path) if (file.exists(path)) digest::digest(file = path, algo = 'sha256') else NA_character_
  manifest_path <- file.path(root, 'sampling-manifest.rds')
  manifest_hash <- hash(manifest_path)
  collector_hash <- hash(file.path(repository, 'doc/validation/relative_share_validation_collect.R'))
  manifest <- readRDS(manifest_path)
  expected_cases <- c(paste0('logit_', c('MC', 'B', 'C', 'BF', 'CF', 'BS', 'CS')),
    'logit_CF_k1', 'logit_CF_k3', 'ces_C', 'ces_CS', 'logit_C_M4T4', 'ces_C_M4T4')
  if (!setequal(names(manifest$fixtures), expected_cases)) stop('Expected the declared thirteen-case sampling manifest.')
  universe <- expected_cases[seq_len(7L)]
  if (!is.numeric(model_prior) || !identical(names(model_prior), universe) ||
      any(!is.finite(model_prior) | model_prior < 0) || abs(sum(model_prior) - 1) > 1e-10)
    stop('model_prior must explicitly name all seven Logit candidates and sum to one.')
  if (!requireNamespace('antitrustBayes', quietly = TRUE)) stop('Load compatible antitrustBayes development source first.')
  helper <- new.env(parent = environment())
  sys.source(file.path(repository, 'doc/validation/nutpie_backend_validation.R'), envir = helper)
  inputs <- data.frame(path = manifest_path, sha256 = manifest_hash, stringsAsFactors = FALSE)
  read_input <- function(path) {
    inputs <<- rbind(inputs, data.frame(path = path, sha256 = hash(path)))
    if (!file.exists(path)) return(simpleError('Missing saved artifact.'))
    tryCatch(readRDS(path), error = identity)
  }
  source_paths <- c(manifest$runtime$source_hashes, manifest$kernel_hashes)
  names(source_paths) <- c(names(manifest$runtime$source_hashes), names(manifest$kernel_hashes))
  expected_sources <- c(paste0('R/', list.files(file.path(repository, 'R'), '[.]R$')),
    'doc/validation/nutpie_backend_validation.R',
    paste0('inst/stan/', list.files(file.path(repository, 'inst/stan'), '[.]stan$')))
  source_names <- unique(c(names(source_paths), expected_sources))
  sources <- data.frame(path = source_names,
    manifest_sha256 = unname(source_paths[match(source_names, names(source_paths))]),
    current_sha256 = vapply(source_names, function(p) hash(file.path(repository, p)), character(1)))
  sources$matches <- sources$manifest_sha256 == sources$current_sha256
  sources$matches[is.na(sources$matches)] <- FALSE
  current_source_verified <- all(sources$matches)
  fits <- policies <- evidences <- list()
  identities <- posterior <- diagnostic <- evidence_summary <- bridge_attempts <- policy_summary <- policy_failures <- list()
  artifact_errors <- policy_identities <- conditional_policy_mcse <- list()
  safe <- function(expr, case, backend, artifact) tryCatch(expr, error = function(e) {
    artifact_errors[[length(artifact_errors) + 1L]] <<- data.frame(case, backend, artifact, reason = conditionMessage(e))
    NULL
  })
  policy_cases <- c('logit_C', 'logit_CF', 'logit_CS', 'logit_CF_k1', 'logit_CF_k3', 'ces_C', 'ces_CS')
  for (case in expected_cases) for (backend in c('cmdstanr', 'nutpie')) {
    key <- paste(case, backend, sep = '_')
    directory <- file.path(root, key)
    attempt <- manifest$attempts[[key]]
    path <- file.path(directory, 'fit.rds')
    fit <- read_input(path)
    valid_fit <- methods::is(fit, 'BayesLogit')
    fit_identity <- if (valid_fit) safe(helper$nutpie_validation_identity(fit), case, backend, 'fit_identity') else NULL
    metadata <- if (valid_fit) fit@metadata else list()
    kernel_path <- paste0('inst/stan/', metadata$stan_kernel, '.stan')
    kernel_index <- match(kernel_path, names(manifest$kernel_hashes))
    metadata_fields <- metadata$stan_data_fields
    data_hash <- if (valid_fit && length(metadata_fields) && all(metadata_fields %in% names(fit@stan_data)))
      digest::digest(fit@stan_data[metadata_fields], algo = 'sha256') else NA_character_
    fit_sha_ok <- valid_fit && identical(hash(path), attempt$fit_sha256)
    source_ok <- valid_fit && !is.na(kernel_index) &&
      identical(metadata$stan_source_hash, unname(manifest$kernel_hashes[kernel_index]))
    data_ok <- valid_fit && identical(metadata$stan_data_hash, data_hash)
    prior_present <- all(vapply(c('prior_manifest', 'resolved_priors', 'prior_normalization'),
      function(field) !is.null(metadata[[field]]), logical(1)))
    identity_ok <- !is.null(fit_identity) && identical(fit_identity, attempt$identity)
    backend_ok <- identical(metadata$backend, backend)
    basis_ok <- identical(metadata$share_basis, 'relative') &&
      identical(metadata$share_observation_measure, 'ilr_quantity_shares_v1')
    verified <- all(fit_sha_ok, source_ok, data_ok, prior_present, identity_ok, backend_ok, basis_ok)
    identities[[key]] <- data.frame(case, backend,
      attempt_status = .relative_collect_scalar(attempt$status, 'not_attempted'),
      fit_present = valid_fit, fit_sha256 = hash(path), manifest_fit_sha256 = .relative_collect_scalar(attempt$fit_sha256, NA_character_),
      fit_sha_verified = fit_sha_ok, kernel = .relative_collect_scalar(metadata$stan_kernel, NA_character_),
      source_sha256 = .relative_collect_scalar(metadata$stan_source_hash, NA_character_), source_matches_manifest = source_ok,
      data_sha256 = .relative_collect_scalar(metadata$stan_data_hash, NA_character_), data_hash_verified = data_ok,
      prior_identity_present = prior_present, manifest_identity_verified = identity_ok,
      backend_verified = backend_ok, relative_basis_verified = basis_ok,
      historical_identity_verified = verified, current_source_verified,
      identity_verified = verified && current_source_verified,
      identity_status = if (!verified) 'unverified_historical_identity' else
        if (current_source_verified) 'current_identity_verified' else 'historical_only_source_mismatch',
      model_spec_version = .relative_collect_scalar(metadata$model_spec_version),
      reason = if (inherits(fit, 'error')) conditionMessage(fit) else '')
    if (valid_fit) {
      fits[[key]] <- fit
      parameters <- safe(helper$nutpie_validation_parameters(fit), case, backend, 'posterior')
      if (!is.null(parameters)) posterior[[key]] <- cbind(case, backend, parameters$summary)
      diagnostics <- safe(helper$nutpie_validation_diagnostics(fit), case, backend, 'diagnostics')
      if (!is.null(diagnostics)) diagnostic[[key]] <- cbind(case, backend, diagnostics$summary)
    }
    evidence <- read_input(file.path(directory, 'evidence.rds'))
    parsed <- safe(.relative_collect_evidence(evidence, case, backend), case, backend, 'evidence')
    if (is.null(parsed)) parsed <- .relative_collect_evidence(simpleError('Invalid saved evidence artifact.'), case, backend)
    parsed$summary$context_source_verified <- verified &&
      identical(parsed$summary$context_kernel, metadata$stan_kernel) &&
      identical(parsed$summary$context_saved_hash, metadata$stan_source_hash) &&
      identical(parsed$summary$context_current_hash, metadata$stan_source_hash)
    evidence_fields <- c('backend', 'share_basis', 'share_observation_measure', 'inference',
      'stan_kernel', 'stan_source_hash', 'stan_data_hash', 'model_spec_version', 'prior_manifest', 'resolved_priors')
    parsed$summary$saved_fit_metadata_verified <- verified && all(vapply(evidence_fields,
      function(field) !is.null(parsed$metadata[[field]]) && identical(parsed$metadata[[field]], metadata[[field]]), logical(1)))
    parsed$summary$historical_identity_verified <- verified && isTRUE(parsed$summary$saved_summary_consistent) &&
      isTRUE(parsed$summary$context_source_verified) && isTRUE(parsed$summary$saved_backend_verified) &&
      isTRUE(parsed$summary$saved_fit_metadata_verified)
    parsed$summary$current_source_verified <- current_source_verified
    parsed$summary$identity_verified <- parsed$summary$historical_identity_verified && current_source_verified
    evidences[[key]] <- parsed
    evidence_summary[[key]] <- parsed$summary
    bridge_attempts[[key]] <- parsed$attempts
    if (case %in% policy_cases) {
      policy <- read_input(file.path(directory, 'policy.rds'))
      policies[[key]] <- policy
      if (is.list(policy) && !inherits(policy, 'error') && is.data.frame(policy$summary)) {
        policy_identity <- if (valid_fit) safe(.relative_collect_policy_identity(fit, policy), case, backend, 'policy_identity') else NULL
        policy_verified <- verified && !is.null(policy_identity)
        if (!is.null(policy_identity)) policy_identities[[key]] <- cbind(case, backend,
          historical_identity_verified = policy_verified, current_source_verified,
          identity_verified = policy_verified && current_source_verified, policy_identity)
        if (policy_verified) {
          shape <- dim(antitrustBayes:::.bayes_backend_draws(fit@stanfit, 'lp__', 'draws_array'))[1:2]
          ratio <- safe(.relative_collect_policy_conditional_mcse(policy, shape), case, backend, 'policy_conditional_mcse')
          if (!is.null(ratio)) {
            conditional_policy_mcse[[key]] <- cbind(case, backend, ratio)
            ratio <- ratio[match(policy$summary$regime, ratio$regime), ]
            policy$summary$estimand <- ratio$estimand
            for (metric in c('mean', 'sd', 'prob_gt_0', 'prob_gt_5')) {
              field <- paste0('mcse_', metric)
              policy$summary[[paste0('saved_', field)]] <- policy$summary[[field]]
              replace <- !is.finite(policy$summary[[field]]) & is.finite(ratio[[paste0('ratio_', field)]])
              policy$summary[[field]][replace] <- ratio[[paste0('ratio_', field)]][replace]
              policy$summary[[paste0(field, '_method')]] <- ifelse(replace,
                'ratio_influence_function_original_chain_order',
                ifelse(is.finite(policy$summary[[field]]), 'saved_intact_chain_mcse', 'unavailable'))
            }
            policies[[key]] <- policy
          }
        }
        policy_summary[[key]] <- cbind(case, backend, historical_identity_verified = policy_verified,
          current_source_verified, identity_verified = policy_verified && current_source_verified, policy$summary)
        if (!policy_verified) policies[[key]] <- simpleError('Saved policy identity verification failed.')
        if (is.data.frame(policy$failures) && nrow(policy$failures)) policy_failures[[key]] <- cbind(case, backend, policy$failures)
        failed <- policy$draws[policy$draws$status != 'ok', , drop = FALSE]
        if (nrow(failed)) policy_failures[[paste0(key, '_regimes')]] <- cbind(case, backend, stage = 'counterfactual', failed)
      } else artifact_errors[[length(artifact_errors) + 1L]] <- data.frame(case, backend,
        artifact = 'policy', reason = if (inherits(policy, 'error')) conditionMessage(policy) else 'Invalid saved policy artifact.')
    }
  }
  identity_table <- .relative_collect_bind(identities)
  posterior_table <- .relative_collect_bind(posterior)
  evidence_table <- .relative_collect_bind(evidence_summary)
  policy_table <- .relative_collect_bind(policy_summary)
  paired_identity <- posterior_comparisons <- evidence_comparisons <- policy_comparisons <- list()
  for (case in expected_cases) {
    keys <- paste(case, c('cmdstanr', 'nutpie'), sep = '_')
    pair <- fits[keys]
    present <- length(pair) == 2L && all(vapply(pair, methods::is, logical(1), 'BayesLogit'))
    science_equal <- present && identical(helper$nutpie_validation_identity(pair[[1L]]), helper$nutpie_validation_identity(pair[[2L]]))
    verified <- all(identity_table$historical_identity_verified[match(keys, paste(identity_table$case, identity_table$backend, sep = '_'))])
    historical_equal <- isTRUE(science_equal) && isTRUE(verified)
    equal <- historical_equal && current_source_verified
    paired_identity[[case]] <- data.frame(case, source_data_prior_equal = science_equal,
      both_manifest_identities_verified = verified, historical_identity_equal = historical_equal,
      current_source_verified, identity_equal = equal)
    if (all(keys %in% names(posterior))) posterior_comparisons[[case]] <- safe(cbind(case,
      .relative_collect_source_status(.relative_collect_compare(posterior[[keys[1L]]], posterior[[keys[2L]]],
        historical_equal, helper$nutpie_validation_compare), current_source_verified)),
      case, 'paired', 'posterior_comparison')
    ev <- lapply(keys, function(key) evidences[[key]]$summary)
    evidence_comparisons[[case]] <- cbind(case, .relative_collect_source_status(
      .relative_collect_evidence_compare(ev[[1L]], ev[[2L]], historical_equal &&
        all(vapply(ev, function(x) isTRUE(x$historical_identity_verified), logical(1)))), current_source_verified))
    if (all(keys %in% names(policy_summary))) {
      ref <- policy_summary[[keys[1L]]]; ref$variable <- ref$regime
      alt <- policy_summary[[keys[2L]]]; alt$variable <- alt$regime
      policy_comparisons[[case]] <- safe(cbind(case,
        .relative_collect_source_status(.relative_collect_compare(ref, alt, historical_equal &&
          all(ref$historical_identity_verified) && all(alt$historical_identity_verified),
          helper$nutpie_validation_compare), current_source_verified)), case, 'paired', 'policy_comparison')
    }
  }
  # Complete seven-model universe only. No successful-subset normalization.
  universe_rows <- list()
  scenarios <- list(cmdstanr = rep('cmdstanr', 7L), nutpie = rep('nutpie', 7L),
    hybrid = rep(c('cmdstanr', 'nutpie'), length.out = 7L))
  for (scenario in names(scenarios)) {
    backends <- scenarios[[scenario]]; keys <- paste(universe, backends, sep = '_')
    ev <- evidence_table[match(keys, paste(evidence_table$case, evidence_table$backend, sep = '_')), ]
    available_fits <- fits[keys]
    required_fits <- available_fits[model_prior > 0]
    complete <- all(vapply(required_fits, methods::is, logical(1), 'BayesLogit'))
    observed_equal <- complete && all(vapply(required_fits, function(f)
      identical(f@data, required_fits[[1L]]@data), logical(1)))
    contract_error <- ''
    contract <- if (complete) tryCatch({
      antitrustBayes:::validate_comparison_contract(required_fits); TRUE
    }, error = function(e) {contract_error <<- conditionMessage(e); FALSE}) else FALSE
    ready <- observed_equal && contract && all(ev$identity_verified[model_prior > 0]) &&
      all(is.finite(ev$logml[model_prior > 0]))
    probability <- rep(NA_real_, 7L)
    if (isTRUE(ready)) {
      required <- model_prior > 0
      lw <- log(model_prior[required]) + ev$logml[required]
      probability[] <- 0; probability[required] <- exp(lw - max(lw)) / sum(exp(lw - max(lw)))
    }
    universe_rows[[scenario]] <- data.frame(scenario, case = universe, backend = backends,
      hypothesis = c('MC', 'B', 'C', 'BF(k=2)', 'CF(k=2)', 'BS(k=2)', 'CS(k=2)'),
      model_prior = as.numeric(model_prior), logml = ev$logml, posterior_probability = probability,
      observations_identical = observed_equal, scientific_contract_verified = contract,
      current_source_verified,
      status = if (!current_source_verified) 'historical_only_source_mismatch_no_normalization' else
        if (isTRUE(ready)) 'complete_positive_prior_universe' else 'incomplete_universe_no_normalization',
      reason = contract_error)
  }
  universe_comparisons <- lapply(c('nutpie', 'hybrid'), function(scenario) {
    ref <- universe_rows$cmdstanr; alt <- universe_rows[[scenario]]
    data.frame(reference_scenario = 'cmdstanr', alternative_scenario = scenario,
      case = universe, reference_probability = ref$posterior_probability,
      alternative_probability = alt$posterior_probability,
      difference = alt$posterior_probability - ref$posterior_probability,
      status = if (all(is.finite(c(ref$posterior_probability, alt$posterior_probability))))
        'descriptive_probability_difference' else 'incomplete_universe_no_comparison')
  })
  # Original CF support and prior are fixed before inspecting any evidence.
  support <- 1:3; prior <- stats::setNames(c(.2, .5, .3), as.character(support))
  component_cases <- c('logit_CF_k1', 'logit_CF', 'logit_CF_k3')
  conditional <- conditional_summary <- conditional_weights <- conditional_policy <- conditional_attrition <- list()
  conditional_scenarios <- list(cmdstanr = rep('cmdstanr', 3L), nutpie = rep('nutpie', 3L), mixed = c('nutpie', 'cmdstanr', 'nutpie'))
  for (scenario in names(conditional_scenarios)) {
    backends <- conditional_scenarios[[scenario]]; keys <- paste(component_cases, backends, sep = '_')
    ev <- evidence_table[match(keys, paste(evidence_table$case, evidence_table$backend, sep = '_')), ]
    logml <- stats::setNames(ev$logml, names(prior)); components <- stats::setNames(fits[keys], names(prior))
    reps <- vapply(keys, function(key) {
      x <- evidences[[key]]$attempts
      x$logml[match(1:3, x$repetition)]
    }, numeric(3L)); colnames(reps) <- names(prior)
    ready <- all(ev$identity_verified) && all(is.finite(logml)) &&
      all(vapply(components, methods::is, logical(1), 'BayesLogit'))
    combined <- if (ready) safe(antitrustBayes::bayes_combine_conditional_k(logml,
      k_support = support, k_prior = prior, repetitions = reps, component_fits = components),
      'conditional_CF', scenario, 'conditional_k_combination') else NULL
    if (is.null(combined)) {
      conditional_summary[[scenario]] <- data.frame(scenario, logml = NA_real_, repetition_sd = NA_real_,
        within_mcse = NA_real_, uncertainty = NA_real_, requested_repetitions = 3L,
        completed_repetitions = 0L, prior_restored_once = NA,
        status = if (current_source_verified) 'incomplete_conditional_support_no_combination' else
          'historical_only_source_mismatch_no_combination')
      conditional_weights[[scenario]] <- data.frame(scenario, k = support, component_backend = backends,
        prior_k = as.numeric(prior), probability = NA_real_, component_logml = as.numeric(logml))
      next
    }
    conditional[[scenario]] <- combined
    lw <- logml + log(prior); max_lw <- max(lw)
    exact_total <- max_lw + log(sum(exp(lw - max_lw)))
    exact_probability <- exp(lw - exact_total)
    prior_once <- identical(combined$k_prior, prior) && abs(combined$logml - exact_total) < 1e-12 &&
      max(abs(combined$posterior_k$probability - exact_probability)) < 1e-12
    # Conservative local propagation of within-repeat errors; the actual saved
    # family repetitions retain their cross-component common-random dependence.
    within <- if (all(is.finite(ev$within_mcse))) sum(exact_probability * ev$within_mcse) else NA_real_
    conditional_summary[[scenario]] <- data.frame(scenario, logml = combined$logml,
      repetition_sd = combined$logml_sd, within_mcse = within,
      uncertainty = sqrt(within^2 + combined$logml_sd^2), requested_repetitions = combined$requested_repetitions,
      completed_repetitions = combined$completed_repetitions, prior_restored_once = prior_once,
      status = combined$numerical_status)
    weight_reps <- t(apply(reps, 1L, function(x) {
      if (any(!is.finite(x))) return(rep(NA_real_, 3L))
      z <- x + log(prior); exp(z - max(z)) / sum(exp(z - max(z)))
    }))
    conditional_weights[[scenario]] <- data.frame(scenario, k = support, component_backend = backends,
      prior_k = as.numeric(prior), probability = exact_probability, component_logml = as.numeric(logml),
      repetition_sd = apply(weight_reps, 2L, function(x) if (sum(is.finite(x)) > 1L) sd(x, na.rm = TRUE) else NA_real_),
      repetition_min = apply(weight_reps, 2L, function(x) if (any(is.finite(x))) min(x, na.rm = TRUE) else NA_real_),
      repetition_max = apply(weight_reps, 2L, function(x) if (any(is.finite(x))) max(x, na.rm = TRUE) else NA_real_))
    mixture <- safe(.relative_collect_weighted_policy(policies[keys], stats::setNames(exact_probability, names(prior)), scenario),
      'conditional_CF', scenario, 'conditional_policy')
    if (!is.null(mixture)) { conditional_policy[[scenario]] <- mixture$summary; conditional_attrition[[scenario]] <- mixture$attrition }
  }
  conditional_comparisons <- conditional_weight_comparisons <- conditional_policy_comparisons <- list()
  for (scenario in c('nutpie', 'mixed')) {
    conditional_comparisons[[scenario]] <- cbind(reference_scenario = 'cmdstanr', alternative_scenario = scenario,
      .relative_collect_evidence_compare(conditional_summary$cmdstanr, conditional_summary[[scenario]],
        all(c('cmdstanr', scenario) %in% names(conditional))))
    ref <- conditional_weights$cmdstanr; alt <- conditional_weights[[scenario]]
    conditional_weight_comparisons[[scenario]] <- data.frame(reference_scenario = 'cmdstanr', alternative_scenario = scenario,
      k = support, reference_probability = ref$probability, alternative_probability = alt$probability,
      difference = alt$probability - ref$probability,
      status = if (all(is.finite(c(ref$probability, alt$probability)))) 'descriptive_probability_difference' else 'unavailable_conditional_evidence',
      reason = 'Three bridge repetitions share posterior draws; no rigorously propagated posterior-k MCSE.')
    if (all(c('cmdstanr', scenario) %in% names(conditional_policy))) {
      ref <- conditional_policy$cmdstanr; alt <- conditional_policy[[scenario]]
      if (all(c('mean', 'sd') %in% names(ref)) && all(c('mean', 'sd') %in% names(alt))) {
        ref$variable <- ref$regime; alt$variable <- alt$regime
        conditional_policy_comparisons[[scenario]] <- cbind(reference_scenario = 'cmdstanr', alternative_scenario = scenario,
          .relative_collect_compare(ref, alt, TRUE, helper$nutpie_validation_compare))
        names(conditional_policy_comparisons[[scenario]]) <- sub('cmdstan_minus_nutpie',
          'reference_minus_alternative', names(conditional_policy_comparisons[[scenario]]), fixed = TRUE)
      }
    }
  }
  # Recheck original inputs immediately before writing any new output.
  inputs$current_sha256 <- vapply(inputs$path, hash, character(1))
  inputs$unchanged <- (is.na(inputs$sha256) & is.na(inputs$current_sha256)) |
    (!is.na(inputs$sha256) & !is.na(inputs$current_sha256) & inputs$sha256 == inputs$current_sha256)
  if (!all(inputs$unchanged)) stop('An original validation artifact changed during collection; no results written.')
  if (!identical(collector_hash, hash(file.path(repository, 'doc/validation/relative_share_validation_collect.R'))))
    stop('Collector source changed during execution; no results written.')
  historical_sampling_complete <- identical(manifest$status, 'completed_attempts') && all(identity_table$historical_identity_verified)
  historical_evidence_complete <- all(evidence_table$historical_identity_verified) && all(evidence_table$numerical_status == 'complete')
  historical_policy_complete <- length(policy_summary) == 2L * length(policy_cases) && all(policy_table$historical_identity_verified)
  sampling_complete <- historical_sampling_complete && current_source_verified
  evidence_complete <- historical_evidence_complete && current_source_verified
  policy_complete <- historical_policy_complete && current_source_verified
  posterior_complete <- length(posterior) == 26L && length(posterior_comparisons) == 13L &&
    all(vapply(posterior_comparisons, is.data.frame, logical(1)))
  diagnostics_complete <- length(diagnostic) == 26L
  collector_operations_complete <- !length(artifact_errors)
  artifact_status <- if (historical_sampling_complete && historical_evidence_complete && historical_policy_complete &&
      posterior_complete && diagnostics_complete && collector_operations_complete)
    'completed_saved_validation' else 'incomplete_saved_validation'
  status <- data.frame(manifest_status = manifest$status, expected_fit_count = 26L,
    saved_fit_count = sum(identity_table$fit_present), sampling_complete, evidence_complete,
    historical_sampling_complete, historical_evidence_complete,
    posterior_extraction_complete = posterior_complete, diagnostics_complete,
    collector_operations_complete,
    required_policy_artifacts_complete = policy_complete,
    historical_policy_artifacts_complete = historical_policy_complete,
    current_source_matches_manifest = current_source_verified, artifact_status,
    status = if (current_source_verified) artifact_status else 'historical_only_source_mismatch')
  tables <- list(status = status, identity = identity_table, paired_identity = .relative_collect_bind(paired_identity),
    posterior = posterior_table, posterior_comparisons = .relative_collect_bind(posterior_comparisons),
    diagnostics = .relative_collect_bind(diagnostic), evidence = evidence_table,
    evidence_attempts = .relative_collect_bind(bridge_attempts), evidence_comparisons = .relative_collect_bind(evidence_comparisons),
    policy = policy_table, policy_identity = .relative_collect_bind(policy_identities),
    policy_conditional_mcse = .relative_collect_bind(conditional_policy_mcse),
    policy_failures = .relative_collect_bind(policy_failures), policy_comparisons = .relative_collect_bind(policy_comparisons),
    model_probabilities = .relative_collect_bind(universe_rows),
    model_probability_comparisons = .relative_collect_bind(universe_comparisons), conditional_k = .relative_collect_bind(conditional_summary),
    conditional_k_weights = .relative_collect_bind(conditional_weights), conditional_k_comparisons = .relative_collect_bind(conditional_comparisons),
    conditional_k_weight_comparisons = .relative_collect_bind(conditional_weight_comparisons),
    conditional_policy = .relative_collect_bind(conditional_policy), conditional_policy_attrition = .relative_collect_bind(conditional_attrition),
    conditional_policy_comparisons = .relative_collect_bind(conditional_policy_comparisons),
    artifact_errors = .relative_collect_bind(artifact_errors), source_identity = sources, input_identity = inputs)
  result <- list(collected = Sys.time(), root = root, manifest = manifest,
    collector_sha256 = collector_hash,
    model_prior = model_prior, universe_backends = scenarios, conditional_backends = conditional_scenarios,
    tables = tables, conditional_k = conditional,
    bridge_contexts = lapply(evidences, `[[`, 'context'),
    uncertainty_definition = paste('Within MCSE is the maximum when every retained finite repetition supplies a finite MCSE, never divided by sqrt(repetitions).',
      'Backend tolerance is four times sqrt(within_ref^2 + repetition_SD_ref^2 + within_alt^2 + repetition_SD_alt^2).',
      'Attrition-conditional policy mean/threshold MCSE uses ratio influence functions in original chain order; SD uses a delta influence with Bessel correction.',
      'These are conditional on paired validity; failure mass is retained and quantile MCSE remains unavailable under attrition.',
      'Conditional-family within MCSE is a conservative local weighted L1 propagation; k weights and empirical policy mixtures have no propagated MCSE.'))
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  for (name in names(tables)) utils::write.csv(tables[[name]], file.path(report_dir,
    paste0('relative_share_', name, '.csv')), row.names = FALSE, na = '')
  raw_dir <- file.path(root, 'collected')
  dir.create(raw_dir, showWarnings = FALSE)
  raw_path <- tempfile(pattern = paste0(format(Sys.time(), '%Y%m%dT%H%M%S'), '-'), tmpdir = raw_dir, fileext = '.rds')
  saveRDS(result, raw_path)
  result$raw_result_path <- raw_path
  message(status$status, ': ', status$saved_fit_count, '/26 fits; raw result: ', raw_path)
  invisible(result)
}

# Explicit CLI only: sourcing never loads package source or collects artifacts.
# Rscript --vanilla doc/validation/relative_share_validation_collect.R \
#   --collect SAMPLING_ROOT REPORT_DIR [REPOSITORY]
if (sys.nframe() == 0L) {
  arguments <- commandArgs(trailingOnly = TRUE)
  if (!length(arguments)) {
    cat('Usage: relative_share_validation_collect.R --collect SAMPLING_ROOT REPORT_DIR [REPOSITORY]\n       relative_share_validation_collect.R --self-check\n')
  } else if (identical(arguments, '--self-check')) {
    relative_validation_collect_checks()
  } else {
    if (!arguments[1L] %in% '--collect' || !length(arguments) %in% c(3L, 4L))
      stop('Use --collect SAMPLING_ROOT REPORT_DIR [REPOSITORY].')
    repository <- if (length(arguments) == 4L) arguments[4L] else '.'
    preparation <- readRDS(file.path(arguments[2L], 'sampling-manifest.rds'))
    libraries <- unique(vapply(preparation$runtime$packages, function(x)
      # Development package roots are not installed-library directories.
      if (isTRUE(x$available) && length(x$path) == 1L &&
          file.exists(file.path(x$path, 'Meta', 'package.rds'))) dirname(x$path) else '', character(1)))
    .libPaths(c(libraries[nzchar(libraries)], .libPaths()))
    pkgload::load_all(repository, quiet = TRUE)
    relative_validation_collect(arguments[2L], arguments[3L], repository = repository)
  }
}
