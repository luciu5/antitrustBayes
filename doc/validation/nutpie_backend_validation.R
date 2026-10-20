#!/usr/bin/env Rscript
# Fixed paired backend validation. Sourcing defines functions; no inference.
# Prepare: Rscript --vanilla doc/validation/nutpie_backend_validation.R
# Execute: bash doc/validation/run_nutpie_validation.sh --run
# Raw fits/logs persist in .scratch/nutpie-validation; never uses paper results.

nutpie_validation_controls <- function() {
  list(chains = 4L, iter = 2000L, warmup = 1000L, cores = 2L,
       seed = 419771L, refresh = 100L, init = 0,
       control = list(adapt_delta = .95, max_treedepth = 13L))
}

nutpie_validation_budget <- function() {
  value <- suppressWarnings(as.numeric(Sys.getenv("AB_NUTPIE_MAX_SECONDS", "14400")))
  if (length(value) != 1L || !is.finite(value) || value < 30 || value != floor(value))
    stop("AB_NUTPIE_MAX_SECONDS must be an integer >= 30.")
  value
}

nutpie_validation_runtime <- function(repository) {
  source <- sort(list.files(file.path(repository, "R"), pattern = "[.]R$", full.names = TRUE))
  source <- c(source, file.path(repository, "doc", "validation", "nutpie_backend_validation.R"))
  hashes <- vapply(source, digest::digest, character(1), file = TRUE, algo = "sha256")
  names(hashes) <- substring(source, nchar(repository) + 2L)
  packages <- c("antitrustBayes", "antitrust", "coordination", "cmdstanr", "nutpieR", "bridgestan", "posterior", "trade")
  descriptions <- lapply(packages, function(package) {
    if (!requireNamespace(package, quietly = TRUE)) return(list(package = package, available = FALSE))
    description <- utils::packageDescription(package)
    list(package = package, available = TRUE, version = as.character(utils::packageVersion(package)),
      path = find.package(package), remote_sha = description$RemoteSha,
      github_sha = description$GithubSHA1)
  })
  names(descriptions) <- packages
  list(r_version = R.version.string, platform = R.version$platform,
    source_hashes = hashes, packages = descriptions)
}

nutpie_validation_bind <- function(rows) {
  if (!length(rows)) return(data.frame())
  fields <- unique(unlist(lapply(rows, names)))
  rows <- lapply(rows, function(x) {
    x[setdiff(fields, names(x))] <- NA
    x[fields]
  })
  do.call(rbind, rows)
}

nutpie_validation_fixture <- function(conduct, pooled = FALSE) {
  markets <- if (pooled) 3L else 1L
  periods <- if (pooled) 2L else 1L
  rows <- truths <- gates <- list()
  first <- NULL
  for (m in seq_len(markets)) {
    args <- list(conduct = conduct, J = 4L, mode = "parameter",
      alpha = 2.4 * exp(.12 * (m - 2)), costs = c(.8, .95, 1.1, 1.25),
      conditional_shares = c(.38, .29, .20, .13), seed = 419701L + m,
      inside_size = 1000)
    if (conduct %in% c("BF", "CF", "BS", "CS")) args$k <- 2L
    base <- suppressWarnings(do.call(antitrustBayes::bayes_calibrated_market, args))
    if (is.null(first)) first <- base
    for (tt in seq_len(periods)) {
      truth <- if (tt == 1L) base else suppressWarnings(antitrustBayes::bayes_resolve_market(
        base, costs = base$truth$costs * c(1.035, .975, 1.02, .99)))
      d <- truth$truth$data
      d$market <- paste0("m", m); d$year <- paste0("t", tt)
      # Deterministic small observation errors shared by all three fits.
      id <- seq_len(nrow(d))
      q <- d$quantities * exp(.025 * sin(id + 2 * m + tt))
      d$quantities <- 1000 * q / sum(q)
      d$margins <- d$margins * exp(.025 * cos(id + m + 2 * tt))
      stopifnot(all(d$margins > 0 & d$margins < 1))
      rows[[length(rows) + 1L]] <- d[c("prices", "quantities", "margins", "market", "year", "group")]
      td <- truth$truth$data; td$market <- d$market; td$year <- d$year
      truths[[length(truths) + 1L]] <- td
      gates[[length(gates) + 1L]] <- data.frame(conduct, pooled, market = m, period = tt,
        max_price_error = truth$fresh_oracle$max_relative_error,
        max_cost_error = truth$fresh_oracle$max_cost_error,
        max_share_error = truth$fresh_oracle$max_share_error,
        passed = truth$fresh_oracle$passed)
    }
  }
  data <- do.call(rbind, rows)
  common <- first$estimation$inputs
  common[names(data)] <- as.list(data)
  common$margin_observed <- rep(TRUE, nrow(data))
  common$use_rho <- TRUE
  common$share_likelihood <- common$margin_likelihood <- "gaussian"
  common$quality_parameterization <- "collapsed"
  common$quality_scope <- "market"
  common$market_intercept_parameterization <- if (pooled) "collapsed" else "sampled"
  common$use_alpha_market_pool <- pooled
  common$s0_pooling <- if (pooled) "market" else "none"
  common$priors <- antitrustBayes::bayes_prior_manifest()
  common$inference <- "hmc"
  if (conduct %in% c("BF", "CF", "BS", "CS")) {
    # Fixed conditional conduct cases, not a changed prior over uncertain k.
    common$k_support <- 2L; common$k_prior <- 1
  }
  list(id = paste0(conduct, if (pooled) "_pooled" else ""), conduct = conduct,
    pooled = pooled, arguments = common, truth = do.call(rbind, truths),
    equilibrium = do.call(rbind, gates))
}

nutpie_validation_design <- function() {
  cases <- c("MC", "B", "C", "BF", "CF", "BS", "CS")
  fixtures <- c(lapply(cases, nutpie_validation_fixture),
                lapply(c("C", "CF", "CS"), nutpie_validation_fixture, pooled = TRUE),
                lapply(c("B", "CS"), nutpie_validation_ces_fixture))
  names(fixtures) <- vapply(fixtures, `[[`, character(1), "id")
  tasks <- expand.grid(fixture = names(fixtures), variant = c("cmdstanr", "nutpie_diag", "nutpie_low_rank"),
                        stringsAsFactors = FALSE)
  tasks <- tasks[!(grepl("^CES_", tasks$fixture) & tasks$variant == "nutpie_low_rank"), ]
  tasks <- tasks[order(match(tasks$fixture, names(fixtures)),
                       match(tasks$variant, c("cmdstanr", "nutpie_diag", "nutpie_low_rank"))), ]
  rownames(tasks) <- NULL
  tasks$id <- paste(tasks$fixture, tasks$variant, sep = "__")
  list(fixtures = fixtures, tasks = tasks)
}

nutpie_validation_ces_fixture <- function(conduct) {
  # Existing public specified-model route supplies a consistent CES economy.
  fixture <- nutpie_validation_fixture(conduct)
  data <- as.data.frame(fixture$arguments[c("prices", "quantities", "margins", "market", "year", "group")])
  parameters <- list(gamma = 6, s0 = .2, supply_model = if (conduct == "B") 1L else 3L)
  if (conduct == "CS") parameters <- c(parameters, list(structure = "stackelberg",
    partition_prior = "ordered_k", k_support = 2L, k_prior = 1))
  market <- antitrustBayes::bayes_market(prices = data$prices, quantities = data$quantities,
    margins = data$margins, group = data$group, outside_price = 1, insideSize = 1000,
    market_side = "output")
  specified <- antitrustBayes::bayes_specify(do.call(antitrustBayes::bayes_params, parameters), market)
  prepared <- suppressWarnings(antitrustBayes::posterior_policy_models(specified, n_draws = 1L))
  model <- prepared$models[[1L]][[1L]]
  if (methods::is(model, "AntitrustFit")) model <- model@model
  margin <- as.numeric(antitrust::calcMargins(model, preMerger = TRUE))
  cost <- as.numeric(antitrust::calcMC(model, preMerger = TRUE))
  model@priceStart <- model@pricePre * 1.13 + .037
  solved <- as.numeric(antitrust::calcPrices(model, preMerger = TRUE))
  price_error <- max(abs(solved - model@pricePre))
  stopifnot(all(is.finite(margin)), all(margin > 0 & margin < 1), all(cost > 0), price_error < 1e-6)
  fixture$id <- paste0("CES_", conduct)
  fixture$arguments$demand_model <- "ces"
  fixture$arguments$margins <- margin
  fixture$equilibrium <- data.frame(conduct, pooled = FALSE, market = 1L, period = 1L,
    max_price_error = price_error,
    max_cost_error = max(abs(cost - data$prices * (1 - margin))),
    max_share_error = NA_real_, passed = TRUE)
  fixture$truth <- data.frame(prices = data$prices, margins = margin, physical_cost = cost, gamma = 6)
  fixture
}

nutpie_validation_parameters <- function(fit) {
  available <- antitrustBayes:::.bayes_backend_variable_names(fit@stanfit)
  structural <- if (identical(as.character(fit@demand_model), "ces")) "gamma_structural_market" else "alpha_structural_market"
  required <- c("mu_log_a", structural, "s0_prob", "sigma_share", "sigma_margin", "rho")
  missing <- setdiff(required, available)
  if (length(missing)) stop("Required validation parameters absent: ", paste(missing, collapse = ", "))
  variables <- intersect(c(required, "sigma_b_strat_economic", "sigma_b_strat", "sigma_log_a_vec",
    "sigma_s0_pool_vec", "mu_b_market_vec", "sigma_b_market_vec",
    "posterior_component_prob", "prior_component_prob", "firm_leader_prob"), available)
  draws <- antitrustBayes:::.bayes_backend_draws(fit@stanfit, variables, "draws_array")
  rows <- lapply(seq_len(dim(draws)[3L]), function(j) {
    x <- draws[, , j, drop = FALSE]
    dim(x) <- dim(draws)[1:2]
    q <- as.numeric(stats::quantile(x, c(.025, .5, .975)))
    mc <- suppressWarnings(posterior::mcse_quantile(x, probs = c(.025, .5, .975)))
    data.frame(variable = dimnames(draws)[[3L]][[j]], mean = mean(x),
      median = q[2L], q025 = q[1L], q975 = q[3L], sd = stats::sd(as.numeric(x)),
      mcse_mean = suppressWarnings(posterior::mcse_mean(x)),
      mcse_sd = suppressWarnings(posterior::mcse_sd(x)),
      mcse_median = mc[2L], mcse_q025 = mc[1L], mcse_q975 = mc[3L],
      rhat = suppressWarnings(posterior::rhat(x)),
      ess_bulk = suppressWarnings(posterior::ess_bulk(x)),
      ess_tail = suppressWarnings(posterior::ess_tail(x)))
  })
  list(summary = do.call(rbind, rows), draws = draws)
}

nutpie_validation_diagnostics <- function(fit) {
  raw <- antitrustBayes:::.bayes_backend_diagnostics(fit@stanfit)
  fields <- dimnames(raw)[[3L]]
  depth_limit <- fit@metadata$sampler_controls$max_depth
  valid_depth <- is.numeric(depth_limit) && length(depth_limit) == 1L && is.finite(depth_limit)
  common <- antitrustBayes:::.bayes_backend_diagnostic_summary(fit@stanfit,
    max_depth = if (valid_depth) depth_limit else NULL)
  summary <- common$posterior
  finite_extreme <- function(x, fun) if (any(is.finite(x))) fun(x[is.finite(x)]) else NA_real_
  rows <- lapply(seq_len(dim(raw)[2L]), function(chain) {
    values <- function(name) if (name %in% fields) as.numeric(raw[, chain, name]) else numeric()
    data.frame(chain, transitions = dim(raw)[1L],
      divergences = common$chains$divergences[chain],
      max_depth = if (valid_depth) depth_limit else NA_real_,
      depth_hits = common$chains$depth_hits[chain],
      leapfrog_steps = if (length(values("n_leapfrog__"))) sum(values("n_leapfrog__")) else NA_real_,
      bfmi = common$chains$bfmi[chain],
      stepsize_median = common$chains$stepsize_median[chain],
      max_rhat = finite_extreme(summary$rhat, max),
      min_bulk_ess = finite_extreme(summary$ess_bulk, min),
      min_tail_ess = finite_extreme(summary$ess_tail, min),
      available_diagnostics = paste(fields, collapse = ";"))
  })
  list(summary = do.call(rbind, rows), raw = raw, common = common)
}

nutpie_validation_policy_chains <- function(records, shape) {
  if (anyDuplicated(records$draw_index)) stop("Duplicate posterior rows in policy results.")
  if (!setequal(records$draw_index, seq_len(prod(shape)))) return(NULL)
  ordered <- records$price_pct[match(seq_len(prod(shape)), records$draw_index)]
  if (any(!is.finite(ordered))) return(NULL)
  matrix(ordered, nrow = shape[1L], ncol = shape[2L])
}

nutpie_validation_policy <- function(fit) {
  draw_shape <- dim(antitrustBayes:::.bayes_backend_draws(fit@stanfit, "lp__", "draws_array"))[1:2]
  n_draws <- prod(draw_shape)
  plan <- antitrustBayes::posterior_policy_draws(fit, n_draws = n_draws, seed = 491871L)
  prepared <- suppressWarnings(antitrustBayes::posterior_policy_models(fit, draw_plan = plan,
    market = "m1", year = "t1", insideSize = 1000, demand_anchoring = "latent_quality"))
  data <- fit@data[fit@data$market == "m1" & fit@data$year == "t1", ]
  weights <- data$quantities / sum(data$quantities)
  firms <- aggregate(data$quantities, list(firm = data$group), sum)
  pair <- as.character(firms$firm[order(-firms$x)[1:2]])
  owner <- as.character(data$group); owner[owner == pair[2L]] <- pair[1L]
  tariff <- ifelse(data$group %in% pair, 0, .25 / 1.25)
  results <- list()
  for (pos in seq_along(plan$draw_indices)) {
    source <- prepared$models[[1L]][[pos]]
    for (regime in c("merger", "tariff25")) {
      record <- data.frame(regime, draw_position = pos, draw_index = plan$draw_indices[pos],
        status = "failed", price_pct = NA_real_, cost_identity_error = NA_real_,
        quality_hash = NA_character_, plan_fingerprint = NA_character_, reason = "")
      record <- tryCatch({
        if (is.null(source)) stop("Baseline realization failed; see retained failures.")
        identity <- attr(source, "bayes_policy_draw", exact = TRUE)
        stopifnot(identical(identity$draw_index, plan$draw_indices[pos]))
        promoted <- trade::as_trade_fit(source, policy = "tariff", tariffPre = 0,
          cost_basis = "effective", margin_basis = "net_revenue")
        stopifnot(identical(identity, attr(promoted, "bayes_policy_draw", exact = TRUE)))
        owners <- if (regime == "merger") owner else as.character(data$group)
        remap <- function(roles) {
          if (regime == "merger") roles[roles == pair[2L]] <- pair[1L]
          unique(roles)
        }
        native <- !methods::is(source, "AntitrustFit")
        roles <- if (!native) list() else if (methods::is(source, "StackelbergLogit"))
          list(leadersPost = remap(promoted@tariff_state$leadersPre)) else
          list(corePost = remap(promoted@tariff_state$corePre))
        solved <- suppressWarnings(do.call(trade::simulate, c(list(promoted,
          ownerPost = owners, tariffPost = if (regime == "tariff25") tariff else 0), roles)))
        before <- trade::tariff_accounts(solved, preMerger = TRUE)
        after <- trade::tariff_accounts(solved, preMerger = FALSE)
        error <- max(abs(before$physical_cost - after$physical_cost))
        stopifnot(is.finite(error), error < 1e-8, all(is.finite(after$price)))
        record$status <- "ok"
        record$price_pct <- 100 * sum(weights * (after$price - before$price)) / sum(weights * data$prices)
        record$cost_identity_error <- error
        record$quality_hash <- digest::digest(identity$quality, algo = "sha256")
        record$plan_fingerprint <- identity$plan_fingerprint
        record
      }, error = function(e) { record$reason <- conditionMessage(e); record })
      results[[length(results) + 1L]] <- record
    }
  }
  rows <- do.call(rbind, results)
  paired <- intersect(rows$draw_position[rows$status == "ok" & rows$regime == "merger"],
                      rows$draw_position[rows$status == "ok" & rows$regime == "tariff25"])
  for (field in c("draw_index", "quality_hash", "plan_fingerprint")) {
    stopifnot(identical(rows[rows$regime == "merger" & rows$draw_position %in% paired, field],
                        rows[rows$regime == "tariff25" & rows$draw_position %in% paired, field]))
  }
  summary <- do.call(rbind, lapply(c("merger", "tariff25"), function(regime) {
    values <- rows$price_pct[rows$regime == regime & rows$draw_position %in% paired]
    quantiles <- if (length(values)) as.numeric(quantile(values, c(.025, .5, .975))) else rep(NA_real_, 3L)
    mcse_mean <- NA_real_; mcse_quantile <- rep(NA_real_, 3L)
    mcse_sd <- mcse_prob_gt_0 <- mcse_prob_gt_5 <- NA_real_
    complete <- length(paired) == n_draws
    if (complete) {
      records <- rows[rows$regime == regime, ]
      values_by_chain <- nutpie_validation_policy_chains(records, draw_shape)
      stopifnot(!is.null(values_by_chain))
      mcse_mean <- suppressWarnings(posterior::mcse_mean(values_by_chain))
      mcse_sd <- suppressWarnings(posterior::mcse_sd(values_by_chain))
      mcse_quantile <- suppressWarnings(posterior::mcse_quantile(values_by_chain, probs = c(.025, .5, .975)))
      mcse_prob_gt_0 <- suppressWarnings(posterior::mcse_mean(1 * (values_by_chain > 0)))
      mcse_prob_gt_5 <- suppressWarnings(posterior::mcse_mean(1 * (values_by_chain > 5)))
    }
    data.frame(regime, n_planned = n_draws, n_realization_failures = nrow(prepared$failures),
      n_paired = length(paired), n_policy_failures = sum(rows$status != "ok" & rows$regime == regime),
      n_valid = sum(rows$regime == regime & rows$status == "ok"),
      mean = if (length(values)) mean(values) else NA_real_, median = quantiles[2L],
      q025 = quantiles[1L], q975 = quantiles[3L], sd = if (length(values) > 1L) sd(values) else NA_real_,
      prob_gt_0 = if (length(values)) mean(values > 0) else NA_real_,
      prob_gt_5 = if (length(values)) mean(values > 5) else NA_real_,
      mcse_mean, mcse_median = mcse_quantile[2L], mcse_q025 = mcse_quantile[1L], mcse_q975 = mcse_quantile[3L],
      mcse_sd, mcse_prob_gt_0, mcse_prob_gt_5,
      status = if (complete) "complete_chain_order_restored" else "conditional_valid_descriptive",
      reason = if (complete) "" else "Attrition retained; conditional paired-valid distribution has no intact-chain MCSE.")
  }))
  list(summary = summary, draws = rows, plan = plan, failures = prepared$failures)
}

nutpie_validation_identity <- function(fit) {
  fields <- c("stan_kernel", "stan_source_hash", "stan_data_hash", "model_spec_version",
              "prior_manifest", "resolved_priors", "prior_normalization",
              "k_support_resolved", "k_prior_resolved")
  values <- fit@metadata[fields]
  if (any(vapply(values[c("stan_kernel", "stan_source_hash", "stan_data_hash", "model_spec_version")], is.null, logical(1))))
    stop("Fit lacks required source/data/model identity.")
  list(metadata = values, data = fit@stan_data, observed = fit@data)
}

nutpie_validation_compare <- function(reference, alternative, identity_equal) {
  if (anyDuplicated(reference$variable) || anyDuplicated(alternative$variable) ||
      !setequal(as.character(reference$variable), as.character(alternative$variable)))
    stop("Paired fits retained different validation parameter columns.")
  alternative <- alternative[match(reference$variable, alternative$variable), , drop = FALSE]
  metrics <- c("mean", "median", "q025", "q975")
  if (all(c("sd", "mcse_sd") %in% intersect(names(reference), names(alternative)))) metrics <- c(metrics, "sd")
  metrics <- c(metrics, intersect(c("prob_gt_0", "prob_gt_5"), intersect(names(reference), names(alternative))))
  rows <- lapply(metrics, function(metric) {
    difference <- alternative[[metric]] - reference[[metric]]
    combined <- sqrt(alternative[[paste0("mcse_", metric)]]^2 + reference[[paste0("mcse_", metric)]]^2)
    tolerance <- pmax(4 * combined, 1e-10)
    data.frame(variable = reference$variable, metric, reference = reference[[metric]],
      alternative = alternative[[metric]], difference, combined_mcse = combined, tolerance,
      cmdstan_minus_nutpie = -difference,
      cmdstan_minus_nutpie_combined_mcse = ifelse(is.finite(combined) & combined > 0, -difference / combined, NA_real_),
      identity_equal = identity_equal,
      status = if (!identity_equal) rep("identity_mismatch", length(difference)) else ifelse(!is.finite(tolerance),
        "unavailable_mcse", ifelse(abs(difference) <= tolerance, "within_mcse_tolerance", "outside_mcse_tolerance")))
  })
  do.call(rbind, rows)
}

nutpie_validation_worker <- function(task_file) {
  if (!identical(Sys.getenv("AB_NUTPIE_EXECUTE", "false"), "true"))
    stop("Inference workers require AB_NUTPIE_EXECUTE=true.")
  task <- readRDS(task_file)
  .libPaths(task$library_paths)
  pkgload::load_all(task$repository, quiet = TRUE)
  if (!is.null(task$runtime) && !identical(task$runtime, nutpie_validation_runtime(task$repository)))
    stop("Runtime/source provenance changed after task preparation; preserve this run and use a new output root.")
  options(antitrustBayes.cmdstan_cache = file.path(task$root, "cmdstan-cache"),
          antitrustBayes.nutpie_cache = file.path(task$root, "nutpie-cache"))
  cmdstan <- Sys.getenv("AB_NUTPIE_CMDSTAN", "")
  if (nzchar(cmdstan)) cmdstanr::set_cmdstan_path(cmdstan)
  args <- task$fixture$arguments
  backend <- if (task$variant == "cmdstanr") "cmdstanr" else "nutpie"
  kernel <- paste0(task$fixture$arguments$demand_model,
    if (task$fixture$conduct %in% c("BF", "CF", "BS", "CS")) "_stackelberg" else "_supply")
  state_path <- file.path(task$directory, "status.rds")
  fit_path <- file.path(task$directory, "fit.rds")
  reuse_fit <- file.exists(fit_path)
  state <- if (reuse_fit && file.exists(state_path)) readRDS(state_path) else
    list(status = "started", fixture = task$fixture$id, variant = task$variant,
    compile_seconds = NA_real_, estimation_seconds = NA_real_, sampling_seconds = NA_real_,
    warmup_seconds = NA_real_, native_sampling_including_warmup = NA_real_, reason = "")
  state$fit_reused <- reuse_fit
  controls <- if (is.null(task$sampler_controls)) nutpie_validation_controls() else task$sampler_controls
  state$chains <- controls$chains
  state$warmup_per_chain <- controls$warmup
  state$retained_per_chain <- controls$iter - controls$warmup
  state$cores <- controls$cores
  state$target_accept <- controls$control$adapt_delta
  state$maximum_depth <- controls$control$max_treedepth
  state$compile_timing_scope <- "Model preparation includes native cache lookup; cache warmth is not inferred."
  state$versions <- paste(vapply(c("antitrustBayes", "antitrust", "coordination", "cmdstanr", "posterior", "trade", "nutpieR"),
    function(p) paste0(p, "=", if (requireNamespace(p, quietly = TRUE))
      as.character(utils::packageVersion(p)) else "missing"), character(1)), collapse = ";")
  started <- proc.time()[["elapsed"]]
  previous_seconds <- if (reuse_fit && is.numeric(state$total_seconds)) state$total_seconds else 0
  on.exit({ state$invocation_seconds <- proc.time()[["elapsed"]] - started
            state$total_seconds <- previous_seconds + state$invocation_seconds
            saveRDS(state, file.path(task$directory, "status.rds")) }, add = TRUE)
  tryCatch({
    args$backend <- backend
    args$stan_options <- if (is.null(task$sampler_controls)) nutpie_validation_controls() else task$sampler_controls
    if (backend == "cmdstanr" && isTRUE(task$fixture$pooled))
      args$stan_options$control$metric <- "dense_e"
    if (backend == "cmdstanr") args$stan_options$output_dir <- file.path(task$directory, "cmdstan")
    if (backend == "nutpie") args$nutpie_options <- list(adaptation = if (task$variant == "nutpie_diag") "diag" else "low_rank")
    dir.create(file.path(task$directory, "cmdstan"), recursive = TRUE, showWarnings = FALSE)
    argument_path <- file.path(task$directory, "arguments.rds")
    if (reuse_fit) {
      if (!file.exists(argument_path) || !identical(args, readRDS(argument_path)))
        stop("Saved fit arguments differ; retain this run and select a new AB_NUTPIE_OUTPUT directory.")
      fit <- readRDS(fit_path)
      current_source <- digest::digest(file = antitrustBayes:::.cmdstan_kernel_file(kernel), algo = "sha256")
      if (!identical(current_source, fit@metadata$stan_source_hash))
        stop("Saved fit Stan source hash differs from the current unchanged kernel.")
    } else {
      compile_start <- proc.time()[["elapsed"]]
      if (backend == "cmdstanr") antitrustBayes:::.cmdstan_model_for_kernel(kernel) else
        antitrustBayes:::.nutpie_model_for_kernel(kernel)
      state$compile_seconds <- proc.time()[["elapsed"]] - compile_start
      saveRDS(args, argument_path)
      sampling_start <- proc.time()[["elapsed"]]
      fit <- do.call(antitrustBayes::bayes_est, args)
      state$estimation_seconds <- proc.time()[["elapsed"]] - sampling_start
      state$adaptation <- fit@metadata$sampler_controls$adaptation
      saveRDS(fit, fit_path)
    }
    if (backend == "cmdstanr") {
      timing <- fit@stanfit$time()
      saveRDS(timing, file.path(task$directory, "native-timing.rds"))
      if (is.data.frame(timing$chains)) {
        if ("warmup" %in% names(timing$chains)) state$warmup_seconds <- max(timing$chains$warmup)
        if ("sampling" %in% names(timing$chains)) state$sampling_seconds <- max(timing$chains$sampling)
      }
    } else {
      native <- fit@stanfit$metadata$sampling_time_native
      if (is.numeric(native) && length(native) == 1L) state$native_sampling_including_warmup <- native
      saveRDS(fit@stanfit$metadata, file.path(task$directory, "native-timing.rds"))
    }
    # Serialized fits must retain exactly the same posterior states.
    reloaded <- readRDS(fit_path)
    antitrustBayes:::.bayes_backend_sampling_contract(reloaded@stanfit,
      requested_chains = controls$chains, inference = "hmc")
    retained_shape <- dim(antitrustBayes:::.bayes_backend_draws(reloaded@stanfit))[1:2]
    stopifnot(identical(as.integer(retained_shape),
      c(as.integer(controls$iter - controls$warmup), as.integer(controls$chains))))
    stopifnot(identical(antitrustBayes:::.bayes_backend_draws(fit@stanfit),
                        antitrustBayes:::.bayes_backend_draws(reloaded@stanfit)))
    state$status <- "fit_saved"
    saveRDS(state, file.path(task$directory, "status.rds"))
    result <- list(parameters = nutpie_validation_parameters(reloaded),
      diagnostics = nutpie_validation_diagnostics(reloaded), identity = nutpie_validation_identity(reloaded))
    result$parameters$summary$bulk_ess_per_sampling_second <-
      result$parameters$summary$ess_bulk / state$sampling_seconds
    result$parameters$summary$bulk_ess_per_estimation_second <-
      result$parameters$summary$ess_bulk / state$estimation_seconds
    result$diagnostics$summary$bulk_ess_per_sampling_second <-
      result$diagnostics$summary$min_bulk_ess / state$sampling_seconds
    result$diagnostics$summary$bulk_ess_per_native_sampling_and_warmup_second <-
      result$diagnostics$summary$min_bulk_ess / state$native_sampling_including_warmup
    result$diagnostics$summary$bulk_ess_per_estimation_second <-
      result$diagnostics$summary$min_bulk_ess / state$estimation_seconds
    saveRDS(result, file.path(task$directory, "validation.rds"))
    result$policy <- nutpie_validation_policy(reloaded)
    saveRDS(result, file.path(task$directory, "validation.rds"))
    state$status <- "completed"
    state$reason <- ""
  }, error = function(e) {
    state$status <<- if (file.exists(file.path(task$directory, "fit.rds"))) "postfit_failed" else "failed"
    state$reason <<- conditionMessage(e)
    message(state$reason)
  })
  invisible(state)
}

nutpie_validation_collect <- function(manifest, root, report_dir) {
  tables <- list(status = list(), parameters = list(), diagnostics = list(), policy = list(), comparisons = list())
  results <- list()
  for (i in seq_len(nrow(manifest$tasks))) {
    task <- manifest$tasks[i, ]
    directory <- file.path(root, task$id)
    status_path <- file.path(directory, "status.rds")
    state <- if (file.exists(status_path)) readRDS(status_path) else
      list(status = "not_run", reason = "No execution artifact; resume starts a fresh invocation budget.")
    state$fixture <- task$fixture; state$variant <- task$variant; state$run_id <- manifest$run_id
    tables$status[[i]] <- as.data.frame(state, stringsAsFactors = FALSE)
    result_path <- file.path(directory, "validation.rds")
    if (!file.exists(result_path)) next
    result <- readRDS(result_path); results[[task$id]] <- result
    for (field in c("parameters", "diagnostics", "policy")) {
      if (is.null(result[[field]])) next
      value <- result[[field]]$summary
      value$fixture <- task$fixture; value$variant <- task$variant; value$run_id <- manifest$run_id
      tables[[field]][[length(tables[[field]]) + 1L]] <- value
    }
  }
  for (fixture in names(manifest$fixtures)) {
    reference <- results[[paste(fixture, "cmdstanr", sep = "__")]]
    variants <- intersect(c("nutpie_diag", "nutpie_low_rank"), manifest$tasks$variant[manifest$tasks$fixture == fixture])
    for (variant in variants) {
      alternative <- results[[paste(fixture, variant, sep = "__")]]
      if (is.null(reference) || is.null(alternative)) {
        value <- data.frame(status = "unavailable_pair", reason = "Both complete validation artifacts are required.")
      } else {
        equal <- identical(reference$identity, alternative$identity)
        value <- nutpie_validation_compare(reference$parameters$summary, alternative$parameters$summary, equal)
        value$quantity <- "parameter"
        if (!is.null(reference$policy) && !is.null(alternative$policy)) {
          policy_reference <- reference$policy$summary
          policy_alternative <- alternative$policy$summary
          policy_reference$variable <- policy_reference$regime
          policy_alternative$variable <- policy_alternative$regime
          policy_comparison <- nutpie_validation_compare(policy_reference, policy_alternative, equal)
          policy_comparison$quantity <- "policy_conditional_paired_valid"
          availability <- data.frame(variable = policy_reference$regime, metric = "availability",
            reference = policy_reference$n_valid / policy_reference$n_planned,
            alternative = policy_alternative$n_valid / policy_alternative$n_planned,
            status = ifelse(policy_reference$n_valid == policy_alternative$n_valid,
              "same_availability_count", "different_availability_count"), quantity = "policy_unconditional_availability")
          value <- nutpie_validation_bind(list(value, policy_comparison, availability))
        }
      }
      value$fixture <- fixture; value$variant <- variant; value$run_id <- manifest$run_id
      tables$comparisons[[length(tables$comparisons) + 1L]] <- value
    }
  }
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  for (name in names(tables)) {
    values <- nutpie_validation_bind(tables[[name]])
    if (!nrow(values)) values <- data.frame(status = "not_run", reason = "No executed validation artifacts.")
    utils::write.csv(values, file.path(report_dir, paste0("nutpie_backend_", name, ".csv")), row.names = FALSE)
  }
  invisible(tables)
}

nutpie_validation_main <- function() {
  repository <- normalizePath(".", mustWork = TRUE)
  pkgload::load_all(repository, quiet = TRUE)
  root <- Sys.getenv("AB_NUTPIE_OUTPUT", file.path(repository, ".scratch", "nutpie-validation"))
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  root <- normalizePath(root, mustWork = TRUE)
  report_dir <- file.path(repository, "doc", "validation")
  dependencies <- data.frame(package = c("nutpieR", "cmdstanr", "trade", "posterior", "pkgload"),
    minimum = c("1.9.1", "0.9.0", "0", "0", "0"), stringsAsFactors = FALSE)
  dependencies$version <- vapply(dependencies$package, function(p)
    if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else NA_character_, character(1))
  dependencies$available <- vapply(seq_len(nrow(dependencies)), function(i)
    !is.na(dependencies$version[i]) && utils::compareVersion(dependencies$version[i], dependencies$minimum[i]) >= 0L, logical(1))
  utils::write.csv(dependencies, file.path(report_dir, "nutpie_backend_dependencies.csv"), row.names = FALSE)
  manifest_path <- file.path(root, "manifest.rds")
  if (file.exists(manifest_path)) {
    manifest <- readRDS(manifest_path)
    stopifnot(identical(manifest$fixture_hash, digest::digest(manifest$fixtures, algo = "sha256")))
    current_source <- vapply(names(manifest$source_hash), function(kernel)
      digest::digest(file = antitrustBayes:::.cmdstan_kernel_file(kernel), algo = "sha256"), character(1))
    stopifnot(identical(manifest$source_hash, current_source))
    current_controls <- nutpie_validation_controls()
    attempted <- any(file.exists(file.path(root, manifest$tasks$id, "task.rds")))
    if (is.null(manifest$sampler_controls) && !attempted) {
      # The interrupted preparation had never fitted a model. Preserve its
      # scientific fixtures and old runner identity, then record the new budget.
      manifest$prepared_script_hash <- manifest$script_hash
      manifest$sampler_controls <- current_controls
    }
    if (!identical(manifest$sampler_controls, current_controls))
      stop("Saved run sampler controls differ; select a new AB_NUTPIE_OUTPUT directory.")
    manifest$runner_hashes <- unique(c(manifest$runner_hashes, manifest$script_hash,
      digest::digest(file = file.path(repository, "doc", "validation", "nutpie_backend_validation.R"), algo = "sha256")))
    saveRDS(manifest, manifest_path)
  } else {
    design <- nutpie_validation_design()
    manifest <- c(design, list(run_id = format(Sys.time(), "%Y%m%dT%H%M%S%z"),
      script_hash = digest::digest(file = file.path(repository, "doc", "validation", "nutpie_backend_validation.R"), algo = "sha256"),
      fixture_hash = digest::digest(design$fixtures, algo = "sha256"),
      source_hash = vapply(c("logit_supply", "logit_stackelberg", "ces_supply", "ces_stackelberg"),
        function(kernel) digest::digest(file = antitrustBayes:::.cmdstan_kernel_file(kernel), algo = "sha256"), character(1)),
      versions = vapply(c("antitrustBayes", "antitrust", "coordination", "cmdstanr", "posterior", "trade"),
        function(p) if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else NA_character_, character(1)),
      git_sha = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
      sampler_controls = nutpie_validation_controls(),
      started_at = NULL, deadline = NULL, invocations = list(), max_cores = 2L))
    saveRDS(manifest, manifest_path)
    utils::write.csv(do.call(rbind, lapply(design$fixtures, function(fixture) {
      value <- fixture$equilibrium; value$fixture <- fixture$id; value
    })),
      file.path(report_dir, "nutpie_backend_fixtures.csv"), row.names = FALSE)
  }
  execute <- identical(Sys.getenv("AB_NUTPIE_EXECUTE", "false"), "true")
  runtime <- nutpie_validation_runtime(repository)
  attempted <- any(file.exists(file.path(root, manifest$tasks$id, "task.rds")))
  if (!identical(manifest$runtime, runtime)) {
    if (attempted) stop("R source, runner, or dependency provenance changed after execution began. ",
      "Preserve this run and select a new AB_NUTPIE_OUTPUT directory.")
    history <- file.path(root, "preparation-history")
    dir.create(history, showWarnings = FALSE)
    snapshot <- file.path(history, paste0("manifest-", format(Sys.time(), "%Y%m%dT%H%M%OS6"), ".rds"))
    saveRDS(manifest, snapshot)
    manifest$runtime <- runtime
    manifest$script_hash <- digest::digest(file = file.path(repository,
      "doc", "validation", "nutpie_backend_validation.R"), algo = "sha256")
    saveRDS(manifest, manifest_path)
  }
  if (!execute) {
    nutpie_validation_collect(manifest, root, report_dir)
    message("Prepared ", length(manifest$fixtures), " fixed fixtures and ", nrow(manifest$tasks),
      " paired tasks; no compilation or inference. Set AB_NUTPIE_EXECUTE=true to execute.")
    return(invisible(manifest))
  }
  if (!all(dependencies$available)) {
    nutpie_validation_collect(manifest, root, report_dir)
    stop("Paired inference withheld: missing/incompatible dependencies: ",
         paste(dependencies$package[!dependencies$available], collapse = ", "), ".")
  }
  if (!nzchar(Sys.which("timeout"))) stop("GNU timeout is required to enforce the per-invocation deadline.")
  max_seconds <- nutpie_validation_budget()
  manifest$started_at <- Sys.time(); manifest$deadline <- manifest$started_at + max_seconds
  invocation <- list(started_at = manifest$started_at, deadline = manifest$deadline,
    git_sha = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
    max_seconds = max_seconds, script_hash = digest::digest(file = file.path(repository,
      "doc", "validation", "nutpie_backend_validation.R"), algo = "sha256"))
  manifest$invocations <- c(manifest$invocations, list(invocation))
  saveRDS(manifest, manifest_path)
  for (i in seq_len(nrow(manifest$tasks))) {
    task <- manifest$tasks[i, ]; directory <- file.path(root, task$id)
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    task_path <- file.path(directory, "task.rds")
    # Completed fits are reused. Postfit work resumes without resampling.
    # Failed inference needs an explicit retry; retain each attempted artifact.
    if (file.exists(task_path)) {
      state_path <- file.path(directory, "status.rds")
      prior <- if (file.exists(state_path)) readRDS(state_path) else list(status = "interrupted")
      if (identical(prior$status, "completed")) next
      saved_fit <- file.exists(file.path(directory, "fit.rds"))
      if (!saved_fit && !identical(Sys.getenv("AB_NUTPIE_RETRY_FAILED", "false"), "true")) next
      archive <- file.path(directory, paste0("attempt-", format(Sys.time(), "%Y%m%dT%H%M%OS6")))
      dir.create(archive, recursive = TRUE, showWarnings = FALSE)
      previous <- list.files(directory, full.names = TRUE)
      previous <- previous[!dir.exists(previous) & basename(previous) != "fit.rds"]
      if (length(previous) && !all(file.copy(previous, archive))) stop("Cannot preserve prior task attempt.")
    }
    remaining <- as.numeric(difftime(manifest$deadline, Sys.time(), units = "secs"))
    if (remaining <= 5) break
    saveRDS(list(fixture = manifest$fixtures[[task$fixture]], variant = task$variant,
      root = root, directory = directory, library_paths = .libPaths(), repository = repository,
      sampler_controls = manifest$sampler_controls, runtime = manifest$runtime), task_path)
    child_args <- c("--signal=TERM", "--kill-after=5s", as.character(floor(remaining - 5)),
      shQuote(file.path(R.home("bin"), "Rscript")), "--vanilla",
      shQuote(file.path(repository, "doc", "validation", "nutpie_backend_validation.R")), "--worker", shQuote(task_path))
    code <- system2("timeout", child_args, stdout = file.path(directory, "worker.log"),
      stderr = file.path(directory, "worker.log"), env = c("OMP_NUM_THREADS=1", "OPENBLAS_NUM_THREADS=1",
        "MKL_NUM_THREADS=1", "CARGO_BUILD_JOBS=2", "MAKEFLAGS=-j2", "CMAKE_BUILD_PARALLEL_LEVEL=2"))
    status_path <- file.path(directory, "status.rds")
    state <- if (file.exists(status_path)) readRDS(status_path) else list()
    if (code != 0L) {
      state$last_completed_stage <- state$status
      state$status <- if (code %in% c(124L, 137L)) "timed_out" else "worker_failed"
      state$reason <- paste("Worker exit status", code, state$reason)
    }
    state$worker_exit_status <- as.integer(code)
    saveRDS(state, status_path)
    nutpie_validation_collect(manifest, root, report_dir)
  }
  manifest$invocations[[length(manifest$invocations)]]$ended_at <- Sys.time()
  saveRDS(manifest, manifest_path)
  nutpie_validation_collect(manifest, root, report_dir)
  invisible(manifest)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 2L && identical(args[[1L]], "--worker")) nutpie_validation_worker(args[[2L]]) else
    nutpie_validation_main()
}
