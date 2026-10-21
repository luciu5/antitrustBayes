#!/usr/bin/env Rscript
# Opt-in public-API closeout. Sourcing defines functions; never starts inference.
# Raw attempts are immutable. Original validation fits/manifests are read-only.

nutpie_closeout_hash <- function(path) digest::digest(file = path, algo = "sha256")

nutpie_closeout_design <- function() {
  tasks <- list()
  add <- function(id, stage, ...) tasks[[length(tasks) + 1L]] <<- c(list(id = id, stage = stage), list(...))
  for (universe in c(FALSE, TRUE)) for (route in c("cmdstanr", "nutpie", "mixed"))
    add(paste0(if (universe) "universe_" else "public_", route), "public", route = route, universe = universe)
  for (model in c("CF", "CS")) for (backend in c("cmdstanr", "nutpie"))
    add(paste("fixed", model, backend, sep = "_"), "fixed", model = model, backend = backend)
  for (k in c(0L, 1L, 3L, 4L)) for (backend in c("cmdstanr", "nutpie"))
    add(paste("sample", k, backend, sep = "_"), "sample", k = k, backend = backend)
  for (k in c(0L, 1L, 3L, 4L)) for (backend in c("cmdstanr", "nutpie"))
    add(paste("component", k, backend, sep = "_"), "component", k = k, backend = backend)
  add("conditional_evidence", "conditional")
  for (k in 0:4) for (backend in c("cmdstanr", "nutpie"))
    add(paste("policy", k, backend, sep = "_"), "policy", k = k, backend = backend)
  add("conditional_policy", "mixture")
  for (kernel in c("logit_supply", "logit_stackelberg", "ces_supply", "ces_stackelberg"))
    add(paste0("portable_", kernel), "portable", kernel = kernel)
  list(tasks = tasks, k_support = 0:4, k_prior = c("0" = .10, "1" = .15, "2" = .25, "3" = .30, "4" = .20),
    controls = list(chains = 4L, iter = 2000L, warmup = 1000L, cores = 2L, init = 0,
      refresh = 500L, control = list(adapt_delta = .95, max_treedepth = 13L)),
    sample_seed_base = 420701L, bridge_seed = 420811L, bridge_repetitions = 3L,
    policy_seed = 491871L, policy_draws = 4000L, policy_bootstrap_seed = 420977L,
    policy_bootstrap_repetitions = 200L, bridge_maximum_four_mcse = .25,
    maximum_probability_difference = .05,
    policy_four_mcse_caps = c(mean = .25, median = .25, q025 = .25, q975 = .25,
      sd_relative = .10, prob_gt_0 = .03, prob_gt_5 = .03),
    scope = paste("Finite normalized evidence through public APIs, independent of sampler-quality verdicts;",
      "all original k states and all failures retained. No default/scientific specification changes."))
}

nutpie_closeout_original_paths <- function(old) {
  out <- list()
  for (backend in c("cmdstanr", "nutpie")) {
    variant <- if (backend == "nutpie") "nutpie_diag" else "cmdstanr"
    for (model in c("B", "MC", "C")) out[[paste(model, backend, sep = "_")]] <-
      if (model == "B") file.path(old, paste0(model, "__", variant), "fit.rds") else
        file.path(old, "evidence", paste0("sample_", model, "_", backend), "fit.rds")
    for (model in c("CF", "CS")) out[[paste(model, backend, sep = "_")]] <-
      file.path(old, paste0(model, "__", variant), "fit.rds")
  }
  out$CES_B_nutpie <- file.path(old, "CES_B__nutpie_diag", "fit.rds")
  out$CES_CS_nutpie <- file.path(old, "CES_CS__nutpie_diag", "fit.rds")
  unlist(out)
}

nutpie_closeout_runtime <- function(repository, old) {
  list(backend = nutpie_validation_runtime(repository),
    cmdstan = list(version = as.character(cmdstanr::cmdstan_version()),
      path = normalizePath(cmdstanr::cmdstan_path(), mustWork = TRUE)),
    caches = list(cmdstan = normalizePath(file.path(old, "cmdstan-cache"), mustWork = FALSE),
      nutpie = normalizePath(file.path(old, "nutpie-cache"), mustWork = FALSE)),
    scripts = vapply(file.path(repository, "doc/validation", c("nutpie_closeout_validation.R",
      "run_nutpie_closeout.sh", "nutpie_evidence_validation.R")), nutpie_closeout_hash, character(1)))
}

nutpie_closeout_assert <- function(manifest) {
  stopifnot(identical(manifest$design, nutpie_closeout_design()),
    identical(manifest$runtime, nutpie_closeout_runtime(manifest$repository, manifest$old)))
  current_stan <- vapply(names(manifest$stan_source_hash), function(kernel)
    nutpie_closeout_hash(antitrustBayes:::.cmdstan_kernel_file(kernel)), character(1))
  if (!identical(current_stan, manifest$stan_source_hash)) stop("A maintained Stan source changed after preparation.")
  for (path in names(manifest$original_hashes)) {
    if (!identical(nutpie_closeout_hash(path), manifest$original_hashes[[path]]))
      stop("An original sealed artifact changed: ", path)
  }
  invisible(TRUE)
}

nutpie_closeout_prepare <- function(repository, old, root) {
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(root, "manifest.rds")
  if (file.exists(path)) { value <- readRDS(path); nutpie_closeout_assert(value); return(value) }
  paths <- nutpie_closeout_original_paths(old)
  if (!all(file.exists(paths))) stop("Missing original fit: ", paste(paths[!file.exists(paths)], collapse = ", "))
  original_manifest <- file.path(old, "manifest.rds")
  original <- readRDS(original_manifest)
  nutpie_evidence_manifest_contract(original)
  starting_git_commit <- trimws(system2("git", c("-C", shQuote(repository), "rev-parse", "HEAD"), stdout = TRUE))
  starting_git_branch <- trimws(system2("git", c("-C", shQuote(repository), "branch", "--show-current"), stdout = TRUE))
  if (length(starting_git_commit) != 1L || !grepl("^[0-9a-f]{40}$", starting_git_commit))
    stop("Cannot record the starting Git revision.")
  if (!length(starting_git_branch) || !nzchar(starting_git_branch)) starting_git_branch <- "(detached)"
  stan_source_hash <- vapply(names(original$source_hash), function(kernel)
    nutpie_closeout_hash(antitrustBayes:::.cmdstan_kernel_file(kernel)), character(1))
  all_paths <- unique(c(original_manifest, file.path(old, "evidence", "manifest.rds"), paths))
  hashes <- stats::setNames(vapply(all_paths, nutpie_closeout_hash, character(1)), all_paths)
  value <- list(repository = repository, old = old, root = root, prepared_at = Sys.time(),
    starting_git_commit = starting_git_commit, starting_git_branch = starting_git_branch,
    stan_source_hash = stan_source_hash,
    design = nutpie_closeout_design(), original_paths = paths, original_hashes = hashes,
    cf_arguments = original$fixtures$CF$arguments, runtime = nutpie_closeout_runtime(repository, old),
    libraries = .libPaths())
  # Full feasible support includes both endpoints for the four ranked firms.
  stopifnot(length(unique(value$cf_arguments$group)) == 4L,
    identical(antitrustBayes:::.stack_validate_k_support(0:4, "fringe", 4L), 0:4))
  for (model in c("B", "MC", "C", "CF", "CS")) {
    fits <- lapply(c("cmdstanr", "nutpie"), function(backend) readRDS(paths[[paste(model, backend, sep = "_")]]))
    stopifnot(identical(nutpie_validation_identity(fits[[1L]]), nutpie_validation_identity(fits[[2L]])))
  }
  saveRDS(value, path)
  value
}

nutpie_closeout_result <- function(manifest, id) {
  path <- file.path(manifest$root, id, "result.rds")
  if (!file.exists(path)) return(NULL)
  readRDS(path)
}

nutpie_closeout_require <- function(manifest, ids) {
  out <- lapply(ids, function(id) nutpie_closeout_result(manifest, id))
  if (any(!vapply(out, function(x) identical(x$status, "completed"), logical(1))))
    stop("Incomplete required tasks: ", paste(ids[!vapply(out, function(x) identical(x$status, "completed"), logical(1))], collapse = ", "))
  stats::setNames(out, ids)
}

nutpie_closeout_component_path <- function(manifest, k, backend) {
  if (k == 2L) manifest$original_paths[[paste("CF", backend, sep = "_")]] else
    file.path(manifest$root, paste("sample", k, backend, sep = "_"), "fit.rds")
}

nutpie_closeout_component_fit <- function(manifest, k, backend) {
  path <- nutpie_closeout_component_path(manifest, k, backend)
  if (k != 2L) {
    sample <- nutpie_closeout_require(manifest, paste("sample", k, backend, sep = "_"))[[1L]]$value
    if (!identical(nutpie_closeout_hash(path), sample$fit_hash)) stop("A saved conditional fit changed: ", path)
  }
  fit <- readRDS(path)
  if (k != 2L && !identical(nutpie_validation_identity(fit), sample$identity)) stop("A component target identity changed.")
  fit
}

nutpie_closeout_bridge_summary <- function(repetitions) {
  values <- vapply(repetitions, function(b) if (is.list(b) && length(b$logml) == 1L &&
    is.finite(b$logml)) b$logml else NA_real_, numeric(1))
  mcse <- vapply(repetitions, function(b) {
    if (is.null(b)) return(NA_real_)
    supplied <- b$mcse_logml
    cv <- tryCatch(bridgesampling::error_measures(b)$cv, error = function(e) NA_real_)
    candidates <- c(supplied, cv)
    candidates <- candidates[is.finite(candidates) & candidates >= 0]
    if (length(candidates)) max(candidates) else NA_real_
  }, numeric(1))
  completed <- sum(is.finite(values))
  spread <- if (completed >= 2L) stats::sd(values[is.finite(values)]) else NA_real_
  uncertainty <- if (completed >= 2L && all(is.finite(mcse[is.finite(values)])))
    max(spread, mcse[is.finite(values)]) else NA_real_
  list(logml = if (completed) stats::median(values[is.finite(values)]) else NA_real_,
    logml_sd = spread, uncertainty = uncertainty, mcse_available = is.finite(uncertainty),
    requested = length(repetitions), completed = completed, repetitions = values, mcse = mcse)
}

nutpie_closeout_public <- function(fits, directory, universe = FALSE, seed) {
  set.seed(seed)
  arguments <- c(fits, list(bridge = TRUE, bridge_repetitions = 3L, bridge_method = "normal", bridge_maxiter = 5000L))
  if (universe) {
    descriptor <- antitrustBayes::bayes_model_universe(as.list(stats::setNames(names(fits), names(fits))),
      stats::setNames(rep(1 / length(fits), length(fits)), names(fits)))
    arguments <- list(fits, universe = descriptor, bridge = TRUE, bridge_repetitions = 3L,
      bridge_method = "normal", bridge_maxiter = 5000L)
  }
  warnings <- character()
  comparison <- withCallingHandlers(do.call(antitrustBayes::bayes_compare, arguments), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  })
  saveRDS(comparison, file.path(directory, "comparison.rds"))
  ordinary <- if (universe) attr(comparison, "source_comparison", exact = TRUE) else comparison
  bridges <- attr(ordinary, "bridge_list", exact = TRUE)
  summaries <- lapply(bridges, nutpie_closeout_bridge_summary)
  weights <- antitrustBayes::bayes_weights(comparison, "bridge")
  stopifnot(identical(names(weights), names(fits)), all(is.finite(weights)), abs(sum(weights) - 1) < 1e-12,
    identical(as.character(ordinary$backend), unname(vapply(fits,
      function(f) antitrustBayes:::.bayes_backend_name(f@stanfit), character(1)))))
  list(comparison = comparison, summaries = summaries, weights = weights, bridge_seed = seed,
    identity = lapply(fits, nutpie_validation_identity), warnings = warnings,
    diagnostics = lapply(fits, nutpie_validation_diagnostics))
}

nutpie_closeout_compare <- function(reference, alternative, label, shared_sources = FALSE) {
  difference <- alternative$logml - reference$logml
  uncertainty <- if (shared_sources) reference$uncertainty + alternative$uncertainty else
    sqrt(reference$uncertainty^2 + alternative$uncertainty^2)
  tolerance <- 4 * uncertainty
  complete <- reference$completed == reference$requested && alternative$completed == alternative$requested
  data.frame(label, reference_logml = reference$logml, alternative_logml = alternative$logml,
    difference, combined_mcse = uncertainty, four_mcse = tolerance,
    uncertainty_combination = if (shared_sources) "sum_shared_sources" else "root_sum_squares_independent",
    precision_passed = complete && is.finite(tolerance) && tolerance <= .25,
    agreement_passed = is.finite(tolerance) && abs(difference) <= max(tolerance, 1e-10),
    requested_reference = reference$requested, completed_reference = reference$completed,
    requested_alternative = alternative$requested, completed_alternative = alternative$completed)
}

nutpie_closeout_family <- function(manifest) {
  components <- list()
  for (k in 0:4) {
    left <- nutpie_closeout_component_fit(manifest, k, "cmdstanr")
    right <- nutpie_closeout_component_fit(manifest, k, "nutpie")
    if (!identical(nutpie_validation_identity(left), nutpie_validation_identity(right)))
      stop("Conditional-k backends have different source/data/prior identities at k=", k)
  }
  for (backend in c("cmdstanr", "nutpie")) {
    ids <- vapply(0:4, function(k) if (k == 2L) paste("fixed_CF", backend, sep = "_") else
      paste("component", k, backend, sep = "_"), character(1))
    results <- nutpie_closeout_require(manifest, ids)
    components[[backend]] <- stats::setNames(lapply(results, function(x) x$value$summaries[[1L]]), as.character(0:4))
  }
  components$mixed <- stats::setNames(lapply(0:4, function(k)
    components[[if (k %% 2L == 0L) "nutpie" else "cmdstanr"]][[as.character(k)]]), as.character(0:4))
  families <- lapply(components, function(parts) {
    z <- vapply(parts, `[[`, numeric(1), "logml")
    repetitions <- do.call(cbind, lapply(parts, `[[`, "repetitions")); colnames(repetitions) <- names(parts)
    value <- antitrustBayes::bayes_combine_conditional_k(z, k_support = 0:4,
      k_prior = manifest$design$k_prior, repetitions = repetitions)
    w <- value$posterior_k$probability
    s <- vapply(parts, `[[`, numeric(1), "uncertainty")
    propagated <- if (all(is.finite(s))) sqrt(sum((w * s)^2)) else NA_real_
    value$uncertainty <- if (is.finite(propagated) && is.finite(value$logml_sd)) max(propagated, value$logml_sd) else NA_real_
    value$requested <- 3L; value$completed <- sum(is.finite(value$repetitions))
    value$weight_mcse <- if (all(is.finite(s))) vapply(seq_along(w), function(i)
      sqrt(sum((w[i] * ((seq_along(w) == i) - w) * s)^2)), numeric(1)) else rep(NA_real_, length(w))
    # Independent arithmetic proves the original prior is restored once.
    logw <- log(manifest$design$k_prior) + z; top <- max(logw)
    independent <- top + log(sum(exp(logw - top)))
    stopifnot(abs(independent - value$logml) < 1e-12,
      max(abs(exp(logw - independent) - w)) < 1e-12)
    value
  })
  list(components = components, families = families,
    comparison = nutpie_closeout_compare(families$cmdstanr, families$nutpie, "CF_full_0_to_4"))
}

nutpie_closeout_weighted_summary <- function(values, weights) {
  stopifnot(length(values) == length(weights), all(is.finite(weights)), all(weights > 0), abs(sum(weights) - 1) < 1e-8)
  valid <- lapply(values, function(x) x[is.finite(x)])
  if (any(!lengths(valid))) stop("A positive-weight k component has no valid policy states; mixture unavailable.")
  y <- unlist(valid, use.names = FALSE)
  mass <- unlist(Map(function(x, w) rep(w / length(x), length(x)), valid, weights), use.names = FALSE)
  ordering <- order(y); cumulative <- cumsum(mass[ordering]); sorted <- y[ordering]
  q <- vapply(c(.025, .5, .975), function(p) sorted[which(cumulative >= p)[1L]], numeric(1))
  constant <- length(unique(y)) == 1L
  mean <- if (constant) y[1L] else sum(mass * y)
  threshold_probability <- function(threshold) {
    hit <- y > threshold
    if (all(hit)) 1 else if (!any(hit)) 0 else sum(mass * hit)
  }
  c(mean = mean, sd = if (constant) 0 else sqrt(sum(mass * (y - mean)^2)),
    median = q[2L], q025 = q[1L], q975 = q[3L],
    prob_gt_0 = threshold_probability(0), prob_gt_5 = threshold_probability(5))
}

nutpie_closeout_policy_mixture <- function(manifest) {
  evidence <- nutpie_closeout_require(manifest, "conditional_evidence")[[1L]]$value
  results <- list()
  for (backend in c("cmdstanr", "nutpie", "mixed")) {
    sources <- if (backend == "mixed") ifelse(0:4 %% 2L == 0L, "nutpie", "cmdstanr") else rep(backend, 5L)
    ids <- paste("policy", 0:4, sources, sep = "_")
    policies <- lapply(nutpie_closeout_require(manifest, ids), function(x) x$value)
    family <- evidence$families[[backend]]
    weights <- family$posterior_k$probability
    parts <- evidence$components[[backend]]
    for (regime in c("merger", "tariff25")) {
      values <- lapply(policies, function(p) {
        d <- p$draws[p$draws$regime == regime, ]
        # Use the same paired-valid states for both policy regimes.
        good <- intersect(p$draws$draw_index[p$draws$regime == "merger" & p$draws$status == "ok"],
          p$draws$draw_index[p$draws$regime == "tariff25" & p$draws$status == "ok"])
        out <- rep(NA_real_, 4000L); out[d$draw_index[d$draw_index %in% good]] <- d$price_pct[d$draw_index %in% good]
        matrix(out, 1000L, 4L)
      })
      estimate <- nutpie_closeout_weighted_summary(values, weights)
      n_valid <- vapply(values, function(x) sum(is.finite(x)), integer(1))
      availability <- sum(weights * n_valid / 4000)
      numerical_sd <- vapply(parts, `[[`, numeric(1), "uncertainty")
      boot <- matrix(NA_real_, manifest$design$policy_bootstrap_repetitions, length(estimate),
        dimnames = list(NULL, names(estimate)))
      evidence_boot <- boot
      if (all(is.finite(numerical_sd))) {
        set.seed(manifest$design$policy_bootstrap_seed + match(backend, c("cmdstanr", "nutpie", "mixed")) +
          10L * match(regime, c("merger", "tariff25")))
        for (b in seq_len(nrow(boot))) {
          # Circular moving blocks preserve each chain and NA failure states.
          resampled <- lapply(values, function(x) {
            apply(x, 2L, function(chain) {
              starts <- sample.int(length(chain), ceiling(length(chain) / 32L), replace = TRUE)
              index <- unlist(lapply(starts, function(start) ((start - 1L + 0:31) %% length(chain)) + 1L))
              chain[index[seq_along(chain)]]
            })
          })
          logz <- stats::rnorm(5L, vapply(parts, `[[`, numeric(1), "logml"), numerical_sd)
          logw <- log(manifest$design$k_prior) + logz
          w <- exp(logw - max(logw)); w <- w / sum(w)
          boot[b, ] <- tryCatch(nutpie_closeout_weighted_summary(resampled, weights), error = function(e) rep(NA_real_, length(estimate)))
          evidence_boot[b, ] <- tryCatch(nutpie_closeout_weighted_summary(values, w), error = function(e) rep(NA_real_, length(estimate)))
        }
      }
      state_mcse <- if (all(is.finite(boot))) apply(boot, 2L, stats::sd) else rep(NA_real_, length(estimate))
      evidence_mcse <- if (all(is.finite(evidence_boot))) apply(evidence_boot, 2L, stats::sd) else rep(NA_real_, length(estimate))
      # Same posterior draws support evidence and policies. Their unknown
      # dependence is not assumed zero: use the sum of the separate SEs.
      mcse <- state_mcse + evidence_mcse
      all_values <- unlist(lapply(values, function(x) x[is.finite(x)]), use.names = FALSE)
      constant <- stats::setNames(rep(length(unique(all_values)) == 1L, length(estimate)), names(estimate))
      constant[["prob_gt_0"]] <- length(unique(all_values > 0)) == 1L
      constant[["prob_gt_5"]] <- length(unique(all_values > 5)) == 1L
      mcse[constant] <- NA_real_
      results[[paste(backend, regime, sep = "_")]] <- list(backend = backend, regime = regime,
        estimate = estimate, mcse = mcse, state_mcse = state_mcse, evidence_mcse = evidence_mcse,
        constant = constant, bootstrap = boot, evidence_bootstrap = evidence_boot, weights = weights,
        valid_by_k = n_valid, planned_by_k = rep(4000L, 5L), failure_mass = 1 - availability,
        status = if (availability == 1) "all_states_valid" else "conditional_within_k_original_k_weights_retained",
        uncertainty = paste("Approximate circular moving-block Monte Carlo SE, block32 separately within each",
          "retained chain;200 repetitions;normal log-evidence perturbations using retained bridge uncertainty separately.",
          "Total MCSE bound is state SE plus evidence SE, retaining unknown dependence conservatively.",
          "NA failure indicators are resampled;no failed draw or k component is dropped from provenance."))
    }
  }
  results
}

nutpie_closeout_portable <- function(task, manifest, directory) {
  model <- c(logit_supply = "B", logit_stackelberg = "CS", ces_supply = "CES_B", ces_stackelberg = "CES_CS")[[task$kernel]]
  path <- manifest$original_paths[[paste(model, "nutpie", sep = "_")]]
  fit <- readRDS(path); original <- antitrustBayes:::.bayes_nutpie_bridge_context(fit)
  on.exit(original$cleanup(), add = TRUE)
  count <- min(20L, nrow(original$samples)); index <- unique(as.integer(round(seq(1, nrow(original$samples), length.out = count))))
  points <- original$samples[index, , drop = FALSE]
  center <- colMeans(original$samples); spread <- apply(original$samples, 2L, stats::sd)
  points <- rbind(points, center, center + .1 * spread, center - .1 * spread)
  evaluate <- function(context) vapply(seq_len(nrow(points)), function(i)
    do.call(context$bridge_args$log_posterior, c(list(points[i, ]), list(data = context$bridge_args$data))), numeric(1))
  expected <- evaluate(original)
  moved <- fit
  moved@stanfit$metadata$lib_path <- file.path(directory, "missing-native-library.so")
  cache <- file.path(directory, "empty-native-cache")
  if (dir.exists(cache)) stop("Portability cache must be new for this attempt.")
  old <- options(antitrustBayes.nutpie_cache = cache)
  on.exit(options(old), add = TRUE)
  restored <- antitrustBayes:::.bayes_nutpie_bridge_context(moved)
  on.exit(restored$cleanup(), add = TRUE)
  restoration <- restored$restoration
  stopifnot(identical(restoration$method, "reconstructed"),
    identical(restoration$reason, "saved_library_unavailable"),
    startsWith(normalizePath(restoration$restored_artifact$path, mustWork = TRUE),
      paste0(normalizePath(cache, mustWork = TRUE), .Platform$file.sep)),
    identical(restoration$restored_artifact$source_hash, fit@metadata$stan_source_hash),
    is.list(restored$sampler_density_verification),
    restored$sampler_density_verification$points == 4000L)
  actual <- evaluate(restored)
  tolerance <- 1e-7 + 1e-9 * pmax(1, abs(actual), abs(expected))
  stopifnot(identical(original$samples, restored$samples), all(is.finite(actual)),
    all(abs(actual - expected) <= tolerance), identical(nutpie_closeout_hash(path), manifest$original_hashes[[path]]))
  list(kernel = task$kernel, points = points, expected = expected, actual = actual,
    maximum_difference = max(abs(actual - expected)), restoration = restored$restoration,
    sampler_density_verification = restored$sampler_density_verification,
    source_hash = fit@metadata$stan_source_hash, data_hash = fit@metadata$stan_data_hash,
    original_fit_unchanged = TRUE)
}

nutpie_closeout_task <- function(task, manifest, directory) {
  task_index <- match(task$id, vapply(manifest$design$tasks, `[[`, character(1), "id"))
  if (is.na(task_index)) stop("Task is absent from the frozen design.")
  bridge_seed <- manifest$design$bridge_seed + 1000L * task_index
  if (task$stage == "public") {
    backends <- if (task$route == "mixed") c("cmdstanr", "nutpie", "cmdstanr") else rep(task$route, 3L)
    fits <- stats::setNames(lapply(seq_along(backends), function(i)
      readRDS(manifest$original_paths[[paste(c("B", "MC", "C")[i], backends[i], sep = "_")]])), c("B", "MC", "C"))
    return(nutpie_closeout_public(fits, directory, universe = task$universe, seed = bridge_seed))
  }
  if (task$stage %in% c("fixed", "component")) {
    if (task$stage == "fixed") path <- manifest$original_paths[[paste(task$model, task$backend, sep = "_")]] else {
      nutpie_closeout_require(manifest, paste("sample", task$k, task$backend, sep = "_"))
      path <- nutpie_closeout_component_path(manifest, task$k, task$backend)
    }
    fit <- if (task$stage == "fixed") readRDS(path) else
      nutpie_closeout_component_fit(manifest, task$k, task$backend)
    return(nutpie_closeout_public(stats::setNames(list(fit), if (task$stage == "fixed") task$model else paste0("CF_k", task$k)),
      directory, seed = bridge_seed))
  }
  if (task$stage == "sample") {
    if (!identical(Sys.getenv("AB_NUTPIE_CLOSEOUT_SAMPLE", "false"), "true")) stop("Additional HMC requires AB_NUTPIE_CLOSEOUT_SAMPLE=true.")
    path <- file.path(manifest$root, task$id, "fit.rds")
    seal_path <- file.path(manifest$root, task$id, "fit-seal.rds")
    if (!file.exists(path)) {
      args <- manifest$cf_arguments; args$k_support <- args$k_prior <- NULL
      args$backend <- task$backend; args$stan_options <- manifest$design$controls
      args$stan_options$seed <- manifest$design$sample_seed_base + task$k
      if (task$backend == "nutpie") args$nutpie_options <- list(adaptation = "diag") else
        args$stan_options$output_dir <- file.path(directory, "cmdstan")
      calls <- do.call(antitrustBayes::bayes_conditional_k_arguments,
        c(args, list(k_support = manifest$design$k_support, k_prior = manifest$design$k_prior)))
      args <- calls[[as.character(task$k)]]
      if (task$backend == "cmdstanr") dir.create(args$stan_options$output_dir, recursive = TRUE, showWarnings = FALSE)
      saveRDS(args, file.path(directory, "arguments.rds"))
      fit <- do.call(antitrustBayes::bayes_est, args)
      saveRDS(fit, path)
      saveRDS(list(fit_hash = nutpie_closeout_hash(path), identity = nutpie_validation_identity(fit),
        arguments = args, manifest_hash = nutpie_closeout_hash(file.path(manifest$root, "manifest.rds"))), seal_path)
    }
    fit <- readRDS(path)
    if (!file.exists(seal_path)) stop("Saved component has no completed fit seal; retain it for explicit inspection.")
    seal <- readRDS(seal_path)
    stopifnot(identical(nutpie_closeout_hash(path), seal$fit_hash),
      identical(nutpie_validation_identity(fit), seal$identity),
      identical(seal$manifest_hash, nutpie_closeout_hash(file.path(manifest$root, "manifest.rds"))))
    stopifnot(identical(as.integer(fit@stan_data$partition_component_k), as.integer(task$k)),
      all(fit@stan_data$partition_prior_probability == 1))
    return(list(fit_path = path, fit_hash = nutpie_closeout_hash(path), identity = nutpie_validation_identity(fit),
      diagnostics = nutpie_validation_diagnostics(fit)))
  }
  if (task$stage == "conditional") return(nutpie_closeout_family(manifest))
  if (task$stage == "policy") {
    if (task$k != 2L) nutpie_closeout_require(manifest, paste("sample", task$k, task$backend, sep = "_"))
    fit <- nutpie_closeout_component_fit(manifest, task$k, task$backend)
    policy <- nutpie_validation_policy(fit)
    stopifnot(length(policy$plan$draw_indices) == manifest$design$policy_draws)
    return(policy)
  }
  if (task$stage == "mixture") return(nutpie_closeout_policy_mixture(manifest))
  if (task$stage == "portable") return(nutpie_closeout_portable(task, manifest, directory))
  stop("Unknown closeout task.")
}

nutpie_closeout_worker <- function(request_path) {
  request <- readRDS(request_path); .libPaths(request$libraries)
  pkgload::load_all(request$repository, quiet = TRUE)
  source(file.path(request$repository, "doc/validation/nutpie_backend_validation.R"))
  source(file.path(request$repository, "doc/validation/nutpie_evidence_validation.R"))
  manifest <- readRDS(request$manifest); nutpie_closeout_assert(manifest)
  options(antitrustBayes.cmdstan_cache = manifest$runtime$caches$cmdstan,
    antitrustBayes.nutpie_cache = manifest$runtime$caches$nutpie)
  warnings <- character(); started <- Sys.time()
  value <- tryCatch(withCallingHandlers(nutpie_closeout_task(request$task, manifest, request$directory),
    warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }), error = identity)
  outcome <- list(status = if (inherits(value, "error")) "failed" else "completed",
    reason = if (inherits(value, "error")) conditionMessage(value) else "", value = if (inherits(value, "error")) NULL else value,
    warnings = warnings, started_at = started, finished_at = Sys.time(), task = request$task,
    manifest_hash = nutpie_closeout_hash(request$manifest))
  saveRDS(outcome, file.path(request$directory, "result.rds"))
  saveRDS(outcome, file.path(manifest$root, request$task$id, "result.rds"))
  nutpie_closeout_assert(manifest)
  if (outcome$status == "failed") quit(status = 1L)
  invisible(outcome)
}

nutpie_closeout_collect <- function(manifest, report_dir) {
  results <- stats::setNames(lapply(manifest$design$tasks, function(task) {
    value <- nutpie_closeout_result(manifest, task$id)
    if (is.null(value)) list(status = "not_run", reason = "No completed execution artifact.") else value
  }), vapply(manifest$design$tasks, `[[`, character(1), "id"))
  status <- do.call(rbind, lapply(seq_along(results), function(i)
    data.frame(task = names(results)[i], stage = manifest$design$tasks[[i]]$stage,
      status = results[[i]]$status, reason = results[[i]]$reason)))
  utils::write.csv(status, file.path(report_dir, "nutpie_closeout_status.csv"), row.names = FALSE)
  completed <- results[vapply(results, function(x) identical(x$status, "completed"), logical(1))]
  evidence <- comparisons <- diagnostics <- probability_checks <- list()
  for (id in names(completed)) {
    value <- completed[[id]]$value
    if (is.null(value$summaries)) next
    for (model in names(value$summaries)) {
      s <- value$summaries[[model]]
      evidence[[length(evidence) + 1L]] <- data.frame(task = id, model, logml = s$logml,
        repetition_sd = s$logml_sd, uncertainty = s$uncertainty, requested = s$requested, completed = s$completed,
        posterior_weight = value$weights[[model]], bridge_seed = value$bridge_seed)
      d <- value$diagnostics[[model]]$summary; d$task <- id; d$model <- model
      diagnostics[[length(diagnostics) + 1L]] <- d
    }
  }
  if (length(evidence)) utils::write.csv(do.call(rbind, evidence), file.path(report_dir, "nutpie_closeout_evidence.csv"), row.names = FALSE)
  if (length(diagnostics)) utils::write.csv(do.call(rbind, diagnostics), file.path(report_dir, "nutpie_closeout_diagnostics.csv"), row.names = FALSE)
  for (route in c("nutpie", "mixed")) for (prefix in c("public_", "universe_")) {
    left <- completed[[paste0(prefix, "cmdstanr")]]$value
    right <- completed[[paste0(prefix, route)]]$value
    if (is.null(left) || is.null(right)) next
    for (model in names(left$summaries)) comparisons[[length(comparisons) + 1L]] <-
      nutpie_closeout_compare(left$summaries[[model]], right$summaries[[model]], paste0(prefix, route, "_", model),
        shared_sources = route == "mixed")
    probability_checks[[length(probability_checks) + 1L]] <- data.frame(route = paste0(prefix, route),
      model = names(left$weights), reference_probability = as.numeric(left$weights),
      alternative_probability = as.numeric(right$weights), difference = as.numeric(right$weights - left$weights),
      tolerance = manifest$design$maximum_probability_difference,
      agreement_passed = abs(right$weights - left$weights) <= manifest$design$maximum_probability_difference)
  }
  for (model in c("CF", "CS")) {
    a <- completed[[paste0("fixed_", model, "_cmdstanr")]]$value$summaries[[1L]]
    b <- completed[[paste0("fixed_", model, "_nutpie")]]$value$summaries[[1L]]
    if (!is.null(a) && !is.null(b)) comparisons[[length(comparisons) + 1L]] <- nutpie_closeout_compare(a, b, paste0("fixed_", model))
  }
  family <- completed$conditional_evidence$value
  if (!is.null(family)) {
    comparisons[[length(comparisons) + 1L]] <- family$comparison
    for (k in 0:4) comparisons[[length(comparisons) + 1L]] <- nutpie_closeout_compare(
      family$components$cmdstanr[[as.character(k)]], family$components$nutpie[[as.character(k)]], paste0("CF_component_k", k))
    comparisons[[length(comparisons) + 1L]] <- nutpie_closeout_compare(
      family$families$cmdstanr, family$families$mixed, "CF_full_0_to_4_mixed", shared_sources = TRUE)
    rows <- lapply(names(family$families), function(backend) {
      f <- family$families[[backend]]
      data.frame(backend, k = 0:4, prior = manifest$design$k_prior,
        component_logml = f$component_logml, posterior_k = f$posterior_k$probability,
        posterior_k_mcse = f$weight_mcse, family_logml = f$logml, family_mcse = f$uncertainty)
    })
    utils::write.csv(do.call(rbind, rows), file.path(report_dir, "nutpie_closeout_conditional_k.csv"), row.names = FALSE)
    for (route in c("nutpie", "mixed")) {
      a <- family$families$cmdstanr; b <- family$families[[route]]
      difference <- b$posterior_k$probability - a$posterior_k$probability
      combined <- if (route == "mixed") a$weight_mcse + b$weight_mcse else sqrt(a$weight_mcse^2 + b$weight_mcse^2)
      probability_checks[[length(probability_checks) + 1L]] <- data.frame(route = paste0("conditional_", route),
        model = paste0("k", 0:4), reference_probability = a$posterior_k$probability,
        alternative_probability = b$posterior_k$probability, difference,
        tolerance = pmin(4 * combined, manifest$design$maximum_probability_difference),
        agreement_passed = is.finite(combined) & abs(difference) <= pmin(4 * combined, manifest$design$maximum_probability_difference))
    }
  }
  if (length(probability_checks)) utils::write.csv(do.call(rbind, probability_checks),
    file.path(report_dir, "nutpie_closeout_probabilities.csv"), row.names = FALSE)
  if (length(comparisons)) utils::write.csv(do.call(rbind, comparisons), file.path(report_dir, "nutpie_closeout_comparisons.csv"), row.names = FALSE)
  policy_tasks <- completed[vapply(completed, function(x) identical(x$task$stage, "policy"), logical(1))]
  if (length(policy_tasks)) {
    rows <- lapply(policy_tasks, function(x) {
      table <- x$value$summary; table$k <- x$task$k; table$backend <- x$task$backend; table
    })
    utils::write.csv(do.call(rbind, rows), file.path(report_dir, "nutpie_closeout_component_policy.csv"), row.names = FALSE)
    failures <- lapply(policy_tasks, function(x) {
      rows <- x$value$draws[x$value$draws$status != "ok", c("regime", "reason"), drop = FALSE]
      if (!nrow(rows)) return(NULL)
      rows$k <- x$task$k; rows$backend <- x$task$backend; rows
    })
    failures <- failures[!vapply(failures, is.null, logical(1))]
    if (length(failures)) {
      failure_rows <- do.call(rbind, failures)
      failure_rows$count <- 1L
      totals <- aggregate(count ~ regime + reason + k + backend, data = failure_rows, sum)
      utils::write.csv(totals, file.path(report_dir, "nutpie_closeout_policy_failures.csv"), row.names = FALSE)
    }
  }
  mixture <- completed$conditional_policy$value
  if (!is.null(mixture)) {
    rows <- lapply(mixture, function(x) data.frame(backend = x$backend, regime = x$regime,
      metric = names(x$estimate), estimate = as.numeric(x$estimate), mcse = as.numeric(x$mcse),
      state_mcse = as.numeric(x$state_mcse), evidence_mcse = as.numeric(x$evidence_mcse),
      exact_constant = x$constant, failure_mass = x$failure_mass, status = x$status))
    utils::write.csv(do.call(rbind, rows), file.path(report_dir, "nutpie_closeout_policy.csv"), row.names = FALSE)
    checks <- list()
    for (route in c("nutpie", "mixed")) for (regime in c("merger", "tariff25")) {
      a <- mixture[[paste("cmdstanr", regime, sep = "_")]]; b <- mixture[[paste(route, regime, sep = "_")]]
      combined <- if (route == "mixed") a$mcse + b$mcse else sqrt(a$mcse^2 + b$mcse^2)
      caps <- stats::setNames(unname(manifest$design$policy_four_mcse_caps[names(a$estimate)]), names(a$estimate))
      caps["sd"] <- manifest$design$policy_four_mcse_caps[["sd_relative"]] * a$estimate[["sd"]]
      difference <- b$estimate - a$estimate
      both_constant <- a$constant & b$constant
      exact_identity <- both_constant & difference == 0
      checks[[length(checks) + 1L]] <- data.frame(route, regime, metric = names(difference), difference,
        combined_mcse = combined, four_mcse = 4 * combined, precision_cap = caps,
        uncertainty_combination = if (route == "mixed") "sum_shared_sources" else "root_sum_squares_independent",
        exact_constant_identity = exact_identity,
        agreement_passed = ifelse(both_constant, difference == 0, abs(difference) <= pmax(4 * combined, 1e-10)),
        precision_passed = exact_identity | (is.finite(combined) & 4 * combined <= caps),
        precision_status = ifelse(exact_identity, "exempt_exact_identity",
          ifelse(is.finite(combined) & 4 * combined <= caps, "within_precision_cap", "unavailable_or_imprecise")),
        any_attrition = a$failure_mass > 0 || b$failure_mass > 0)
    }
    utils::write.csv(do.call(rbind, checks), file.path(report_dir, "nutpie_closeout_policy_comparisons.csv"), row.names = FALSE)
  }
  portable <- completed[startsWith(names(completed), "portable_")]
  if (length(portable)) utils::write.csv(do.call(rbind, lapply(portable, function(x)
    data.frame(kernel = x$value$kernel, points = nrow(x$value$points), maximum_difference = x$value$maximum_difference,
      original_fit_unchanged = x$value$original_fit_unchanged,
      source_hash = x$value$source_hash, data_hash = x$value$data_hash,
      compile_lookup_seconds = x$value$restoration$restored_artifact$compile_lookup_seconds,
      sampler_points = x$value$sampler_density_verification$points,
      sampler_maximum_difference = x$value$sampler_density_verification$max_absolute_difference))),
    file.path(report_dir, "nutpie_closeout_portability.csv"), row.names = FALSE)
  lines <- c("# Nutpie public-API closeout", "", paste("Collected", format(Sys.time(), tz = "UTC"), "UTC."), "",
    paste(sum(status$status == "completed"), "of", nrow(status), "tasks completed;", sum(status$status == "failed"), "failed."), "",
    "Execution completion, finite evidence, numerical precision, sampler diagnostics and policy support are separate outcomes.",
    "No divergence threshold gates calculations. Original targets, priors, controls and failed states remain recorded.", "",
    "Protocol: J=4 CF full support 0:4, prior (.10,.15,.25,.30,.20), eight new singleton fits; original k=2 fits reused.",
    "Each public comparison requests three normal bridge repetitions. Agreement uses four combined numerical uncertainties, capped at .25 log units for precision claims.",
    "Policy calculations retain all 4,000 draws/component with fixed plans; conditional mixtures retain original evidence-based k weights after within-component omissions.",
    "Mixture MCSE bounds add approximate state block-bootstrap and bridge-weight SEs without assuming their dependence is zero; failures are retained and resampled, not replaced.", "",
    "See nutpie_closeout_status.csv, evidence/diagnostics/comparisons/probabilities, conditional_k, component_policy/policy_failures, policy/policy_comparisons and portability CSVs beside this report.",
    "Raw public return objects, repetition objects, plans, draw failures, source/fit seals and every attempt remain in .scratch/nutpie-closeout.")
  writeLines(lines, file.path(report_dir, "nutpie_closeout_results.md"))
  invisible(status)
}

nutpie_closeout_main <- function(action = "--prepare") {
  repository <- normalizePath(".", mustWork = TRUE)
  pkgload::load_all(repository, quiet = TRUE)
  source(file.path(repository, "doc/validation/nutpie_backend_validation.R"))
  source(file.path(repository, "doc/validation/nutpie_evidence_validation.R"))
  old <- normalizePath(Sys.getenv("AB_NUTPIE_CLOSEOUT_ORIGINAL", ".scratch/nutpie-validation"), mustWork = TRUE)
  root <- Sys.getenv("AB_NUTPIE_CLOSEOUT_OUTPUT", file.path(repository, ".scratch/nutpie-closeout"))
  manifest <- nutpie_closeout_prepare(repository, old, root)
  if (action == "--run") {
    if (!identical(Sys.getenv("AB_NUTPIE_CLOSEOUT_EXECUTE", "false"), "true")) stop("Explicit execution opt-in required.")
    seconds <- suppressWarnings(as.integer(Sys.getenv("AB_NUTPIE_CLOSEOUT_SECONDS", "14400")))
    if (length(seconds) != 1L || is.na(seconds) || seconds < 30L || seconds > 14400L) stop("Invocation seconds must be between30 and14400.")
    deadline <- Sys.time() + seconds
    invocation <- tempfile("invocation-", tmpdir = manifest$root, fileext = ".rds")
    saveRDS(list(started_at = Sys.time(), deadline = deadline, action = action), invocation)
    for (task in manifest$design$tasks) {
      current <- nutpie_closeout_result(manifest, task$id)
      if (!is.null(current) && current$status == "completed") next
      if (!is.null(current) && current$status == "failed" && !identical(Sys.getenv("AB_NUTPIE_CLOSEOUT_RETRY_FAILED", "false"), "true")) next
      remaining <- floor(as.numeric(difftime(deadline, Sys.time(), units = "secs")) - 5)
      if (remaining < 10L) break
      task_root <- file.path(manifest$root, task$id); dir.create(task_root, recursive = TRUE, showWarnings = FALSE)
      directory <- tempfile("attempt-", tmpdir = task_root); dir.create(directory)
      request <- file.path(directory, "request.rds")
      saveRDS(list(manifest = file.path(manifest$root, "manifest.rds"), repository = repository,
        libraries = manifest$libraries, task = task, directory = directory), request)
      cat("Starting", task$id, "\n"); flush.console()
      log <- file.path(directory, "worker.log")
      code <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
        c("--vanilla", shQuote(file.path(repository, "doc/validation/nutpie_closeout_validation.R")), "--worker", shQuote(request)),
        stdout = log, stderr = log, timeout = remaining))
      if (!file.exists(file.path(directory, "result.rds"))) {
        outcome <- list(status = "interrupted", reason = paste("Workerexit", code, "without final result; attempt retained."), task = task)
        saveRDS(outcome, file.path(directory, "result.rds")); saveRDS(outcome, file.path(task_root, "result.rds"))
      }
      nutpie_closeout_collect(manifest, file.path(repository, "doc/validation"))
      nutpie_closeout_assert(manifest)
    }
  }
  status <- nutpie_closeout_collect(manifest, file.path(repository, "doc/validation"))
  print(table(status$stage, status$status)); invisible(status)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 2L && args[[1L]] == "--worker") nutpie_closeout_worker(args[[2L]]) else
    nutpie_closeout_main(if (length(args)) args[[1L]] else "--prepare")
}
