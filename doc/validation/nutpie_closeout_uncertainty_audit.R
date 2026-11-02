#!/usr/bin/env Rscript
# Saved policy rows only. Never fits, evaluates a target, or solves a policy.
# Keep the frozen closeout runner and its original uncertainty estimates intact.

nutpie_uncertainty_bound <- function(component_state_se, evidence_se) {
  if (!is.matrix(component_state_se) || nrow(component_state_se) != 5L ||
      is.null(colnames(component_state_se)) ||
      !identical(colnames(component_state_se), names(evidence_se)))
    stop("Require every original k component and matching named statistics.")
  if (any(!is.finite(component_state_se)) || any(component_state_se < 0) ||
      any(!is.finite(evidence_se)) || any(evidence_se < 0))
    return(stats::setNames(rep(NA_real_, ncol(component_state_se)), colnames(component_state_se)))
  colSums(component_state_se) + evidence_se
}

nutpie_uncertainty_component <- function(values, weights, component, statistic,
                                        repetitions = 200L, seed) {
  if (length(values) != 5L || !component %in% seq_along(values) ||
      any(!vapply(values, function(x) is.matrix(x) && identical(dim(x), c(1000L, 4L)), logical(1))))
    stop("Require complete original chain axes for all five components.")
  point <- statistic(values, weights)
  boot <- matrix(NA_real_, repetitions, length(point), dimnames = list(NULL, names(point)))
  set.seed(seed)
  for (rep in seq_len(repetitions)) {
    changed <- values
    changed[[component]] <- apply(values[[component]], 2L, function(chain) {
      starts <- sample.int(length(chain), ceiling(length(chain) / 32L), replace = TRUE)
      index <- unlist(lapply(starts, function(start) ((start - 1L + 0:31) %% length(chain)) + 1L))
      chain[index[seq_along(chain)]]
    })
    boot[rep, ] <- tryCatch(statistic(changed, weights), error = function(e) rep(NA_real_, length(point)))
  }
  list(point = point, bootstrap = boot,
    se = if (all(is.finite(boot))) apply(boot, 2L, stats::sd) else
      stats::setNames(rep(NA_real_, length(point)), names(point)),
    seed = seed, block_length = 32L, repetitions = repetitions)
}

nutpie_uncertainty_difference <- function(reference, alternative, caps) {
  metrics <- names(reference$estimate)
  if (!identical(metrics, names(alternative$estimate))) stop("Mixture statistics do not match.")
  difference <- alternative$estimate - reference$estimate
  bound <- reference$bound + alternative$bound
  both_constant <- reference$constant & alternative$constant
  exact <- both_constant & difference == 0
  data.frame(metric = metrics, reference = reference$estimate, alternative = alternative$estimate,
    difference, combined_mcse_bound = bound, four_mcse_bound = 4 * bound,
    precision_cap = caps, exact_constant_identity = exact,
    agreement_passed = ifelse(both_constant, difference == 0, abs(difference) <= pmax(4 * bound, 1e-10)),
    precision_passed = exact | (is.finite(bound) & 4 * bound <= caps),
    precision_status = ifelse(exact, "exempt_exact_identity",
      ifelse(is.finite(bound) & 4 * bound <= caps, "within_approximate_bound_cap", "unavailable_or_imprecise")))
}

nutpie_uncertainty_values <- function(policies, regime) {
  lapply(policies, function(policy) {
    rows <- policy$draws
    if (!is.data.frame(rows) || !all(c("regime", "draw_index", "status", "price_pct") %in% names(rows)) ||
        anyNA(rows$regime) || !setequal(unique(rows$regime), c("merger", "tariff25")) ||
        anyNA(rows$status) || any(!rows$status %in% c("ok", "failed")))
      stop("Invalid saved policy regimes or statuses.")
    if (!is.numeric(rows$price_pct) || any(rows$status == "ok" & !is.finite(rows$price_pct)))
      stop("A successful saved policy row has no finite effect.")
    for (expected_regime in c("merger", "tariff25")) {
      index <- rows$draw_index[rows$regime == expected_regime]
      if (!is.numeric(index) || length(index) != 4000L || any(!is.finite(index)) ||
          any(index != floor(index)) || anyDuplicated(index) || !setequal(index, 1:4000))
        stop("Incomplete or duplicate saved policy draw axes in ", expected_regime, ".")
    }
    valid <- intersect(rows$draw_index[rows$regime == "merger" & rows$status == "ok"],
      rows$draw_index[rows$regime == "tariff25" & rows$status == "ok"])
    selected <- rows[rows$regime == regime & rows$draw_index %in% valid, ]
    out <- rep(NA_real_, 4000L); out[selected$draw_index] <- selected$price_pct
    matrix(out, 1000L, 4L)
  })
}

nutpie_uncertainty_input_contract <- function(inputs, manifest_hash) {
  if (!length(inputs) || !all(vapply(inputs, function(x) identical(x$status, "completed"), logical(1))))
    stop("Every original component policy and family result must be available.")
  if (!all(vapply(inputs, function(x) identical(x$manifest_hash, manifest_hash), logical(1))))
    stop("A saved policy/evidence result belongs to a different producer-run manifest.")
  invisible(TRUE)
}

nutpie_uncertainty_main <- function() {
  repository <- normalizePath(".", mustWork = TRUE)
  root <- normalizePath(Sys.getenv("AB_NUTPIE_CLOSEOUT_OUTPUT", ".scratch/nutpie-closeout"), mustWork = TRUE)
  manifest <- readRDS(file.path(root, "manifest.rds"))
  .libPaths(manifest$libraries)
  helper <- new.env(parent = globalenv())
  sys.source(file.path(repository, "doc/validation/nutpie_closeout_validation.R"), envir = helper)
  hash <- function(path) digest::digest(file = path, algo = "sha256")
  helper_path <- file.path(repository, "doc/validation/nutpie_closeout_validation.R")
  if (!identical(hash(helper_path), unname(manifest$runtime$scripts[[helper_path]])))
    stop("Frozen closeout helper changed; preserve the original run.")
  ids <- c("conditional_evidence", "conditional_policy", paste("policy", rep(0:4, each = 2L),
    rep(c("cmdstanr", "nutpie"), 5L), sep = "_"))
  paths <- file.path(root, ids, "result.rds")
  if (!all(file.exists(paths))) stop("Wait for every original evidence/policy result; none may be omitted.")
  inputs <- stats::setNames(lapply(paths, readRDS), ids)
  nutpie_uncertainty_input_contract(inputs, hash(file.path(root, "manifest.rds")))
  directory <- file.path(root, "uncertainty-audit")
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  seal_path <- file.path(directory, "manifest.rds")
  seal <- list(input_hashes = stats::setNames(vapply(paths, hash, character(1)), paths),
    frozen_manifest_hash = hash(file.path(root, "manifest.rds")), helper_hash = hash(helper_path),
    audit_hash = hash(file.path(repository, "doc/validation/nutpie_closeout_uncertainty_audit.R")),
    r_version = R.version.string, repetitions = 200L, block_length = 32L,
    seed_base = 421199L,
    method = "one component at a time; sum state SEs plus evidence SE; sum bounds across compared routes")
  if (file.exists(seal_path) && !identical(seal, readRDS(seal_path)))
    stop("Uncertainty audit source or saved inputs changed; retain this audit and choose a reviewed new location.")
  if (!file.exists(seal_path)) saveRDS(seal, seal_path)
  original <- inputs$conditional_policy$value
  family <- inputs$conditional_evidence$value
  result <- component_rows <- list()
  for (route in c("cmdstanr", "nutpie", "mixed")) for (regime in c("merger", "tariff25")) {
    key <- paste(route, regime, sep = "_")
    source <- if (route == "mixed") ifelse(0:4 %% 2L == 0L, "nutpie", "cmdstanr") else rep(route, 5L)
    policies <- lapply(paste("policy", 0:4, source, sep = "_"), function(id) inputs[[id]]$value)
    values <- nutpie_uncertainty_values(policies, regime)
    weights <- family$families[[route]]$posterior_k$probability
    point <- helper$nutpie_closeout_weighted_summary(values, weights)
    if (!identical(point, original[[key]]$estimate)) stop("A frozen policy point estimate changed.")
    component_se <- matrix(NA_real_, 5L, length(point), dimnames = list(as.character(0:4), names(point)))
    for (component in 1:5) {
      checkpoint <- file.path(directory, paste0(key, "_k", component - 1L, ".rds"))
      if (file.exists(checkpoint)) calculation <- readRDS(checkpoint) else {
        cat("Auditing", key, "k", component - 1L, "\n"); flush.console()
        calculation <- nutpie_uncertainty_component(values, weights, component,
          helper$nutpie_closeout_weighted_summary, repetitions = seal$repetitions,
          seed = seal$seed_base + 1000L * match(route, c("cmdstanr", "nutpie", "mixed")) +
            100L * match(regime, c("merger", "tariff25")) + component)
        saveRDS(calculation, checkpoint)
      }
      expected_seed <- seal$seed_base + 1000L * match(route, c("cmdstanr", "nutpie", "mixed")) +
        100L * match(regime, c("merger", "tariff25")) + component
      if (!identical(calculation$point, point) || !identical(calculation$seed, expected_seed) ||
          !identical(calculation$block_length, seal$block_length) || !identical(calculation$repetitions, seal$repetitions))
        stop("A component audit checkpoint does not match its declared source/seed.")
      component_se[component, ] <- calculation$se
      component_rows[[length(component_rows) + 1L]] <- data.frame(route, regime, k = component - 1L,
        metric = names(point), mcse = calculation$se, bootstrap_seed = calculation$seed,
        valid_states = sum(is.finite(values[[component]])), planned_states = 4000L)
    }
    historical <- original[[key]]
    evidence_se <- stats::setNames(as.numeric(historical$evidence_mcse), names(point))
    bound <- nutpie_uncertainty_bound(component_se, evidence_se)
    bound[historical$constant] <- NA_real_
    result[[key]] <- list(route = route, regime = regime, estimate = point, bound = bound,
      component_state_se = component_se, summed_state_se = colSums(component_se),
      evidence_se = evidence_se, historical_mcse = historical$mcse,
      constant = historical$constant, failure_mass = historical$failure_mass, status = historical$status)
  }
  rows <- lapply(result, function(x) data.frame(route = x$route, regime = x$regime,
    metric = names(x$estimate), estimate = x$estimate, historical_approximate_mcse = x$historical_mcse,
    summed_component_state_se = x$summed_state_se, evidence_se = x$evidence_se,
    revised_approximate_mcse_bound = x$bound, exact_constant = x$constant,
    failure_mass = x$failure_mass, status = x$status))
  checks <- list()
  for (route in c("nutpie", "mixed")) for (regime in c("merger", "tariff25")) {
    reference <- result[[paste("cmdstanr", regime, sep = "_")]]
    alternative <- result[[paste(route, regime, sep = "_")]]
    caps <- stats::setNames(unname(manifest$design$policy_four_mcse_caps[names(reference$estimate)]), names(reference$estimate))
    caps["sd"] <- manifest$design$policy_four_mcse_caps[["sd_relative"]] * reference$estimate[["sd"]]
    comparison <- nutpie_uncertainty_difference(reference, alternative, caps)
    comparison$route <- route; comparison$regime <- regime
    comparison$any_attrition <- reference$failure_mass > 0 || alternative$failure_mass > 0
    checks[[length(checks) + 1L]] <- comparison
  }
  report <- file.path(repository, "doc/validation")
  utils::write.csv(do.call(rbind, component_rows), file.path(report, "nutpie_closeout_uncertainty_components.csv"), row.names = FALSE)
  utils::write.csv(do.call(rbind, rows), file.path(report, "nutpie_closeout_uncertainty_bounds.csv"), row.names = FALSE)
  comparisons <- do.call(rbind, checks)
  utils::write.csv(comparisons, file.path(report, "nutpie_closeout_uncertainty_comparisons.csv"), row.names = FALSE)
  saveRDS(list(seal = seal, results = result, comparisons = comparisons), file.path(directory, "result.rds"))
  for (path in names(seal$input_hashes)) if (!identical(hash(path), seal$input_hashes[[path]]))
    stop("An original saved policy/evidence result changed during the audit.")
  lines <- c("# Conditional-policy uncertainty audit", "", paste("Collected", format(Sys.time(), tz = "UTC"), "UTC."), "",
    "This saved-policy-only audit preserves every frozen point estimate, k weight, draw plan, failure and original approximate MCSE.",
    "The frozen helper shares its latent-quality seed across components/backends, so independence of policy Monte Carlo errors is not assumed here.", "",
    "For each route/regime, one k component at a time is resampled within each retained chain using circular blocks32 and200 fixed-seed repetitions.",
    "The sum of the five component statistic SEs is a conservative first-order unknown-covariance bound; the retained evidence-weight SE is then added.",
    "Comparison bounds add both routes' total SE bounds, including homogeneous comparisons that share latent-quality noise.",
    "The bootstrap and normal bridge-weight error estimates remain approximate. Minkowski's inequality bounds the linearized sum without assuming independent component errors; it does not turn finite-sample nonlinear quantile bootstraps into exact confidence intervals.",
    "In particular, atoms, low tail ESS, chain nonstationarity and insufficient block length can limit the approximation. Actual divergences/ESS/attrition remain in the original reports.", "",
    "Matching exact constants retain unavailable MCSE and are exempt identities; different constants fail agreement. Existing .25price-point/.10relative-SD/.03probability precision caps are unchanged.",
    "Within-component policy omissions retain original k weights; failure mass remains explicit. No claim of an unconditional policy distribution is made when any component has attrition.", "",
    paste("Revised comparison rows:", nrow(comparisons), "; agreement failures:", sum(!comparisons$agreement_passed, na.rm = TRUE),
      "; unavailable/imprecise:", sum(!comparisons$precision_passed, na.rm = TRUE), "."), "",
    "Tracked tables: nutpie_closeout_uncertainty_components.csv, nutpie_closeout_uncertainty_bounds.csv and nutpie_closeout_uncertainty_comparisons.csv.",
    "Raw component bootstraps, exact source/input seals and the audit result remain under .scratch/nutpie-closeout/uncertainty-audit.",
    "No HMC, native target evaluation, policy reconstruction or new equilibrium solve was performed for this audit.")
  writeLines(lines, file.path(report, "nutpie_closeout_uncertainty_results.md"))
  invisible(comparisons)
}

if (sys.nframe() == 0L) nutpie_uncertainty_main()
