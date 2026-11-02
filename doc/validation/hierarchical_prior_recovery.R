#!/usr/bin/env Rscript
# Paired known-truth validation of the package alpha/s0 hierarchy SD defaults.
# No production submission or paper-DGP change. Source defines functions only.
# AB_HIERARCHY_SCOPE=prepare creates and validates fixtures without HMC.
# AB_HIERARCHY_OUTPUT, AB_HIERARCHY_CMDSTAN and AB_HIERARCHY_CACHE select paths.

hierarchy_lse <- function(x) { top <- max(x); top + log(sum(exp(x - top))) }

hierarchy_fixture <- function(case, seed = 581713L, M = 8L, T = 4L, J = 5L) {
  scales <- switch(case, small = c(alpha = .05, s0 = .08),
    moderate = c(alpha = .30, s0 = .40), stop("Unknown fixture: ", case))
  set.seed(seed)
  z_alpha <- rnorm(M); z_s0 <- rnorm(M); z_intercept <- rnorm(M)
  quality <- matrix(rnorm(M * J, 0, .6), M, J)
  quality <- quality - rowMeans(quality)
  relative_price <- array(0, c(M, T, J))
  for (m in seq_len(M)) for (tt in seq_len(T))
    relative_price[m, tt, ] <- .12 * (seq_len(J) - (J + 1)/2) +
      .10 * sin(seq_len(J) * tt + m) + runif(J, -.025, .025)
  z_share <- array(rnorm(M * T * J), c(M, T, J))
  z_margin <- array(rnorm(M * T * J), c(M, T, J))
  alpha <- 2.5 * exp(scales[["alpha"]] * z_alpha)
  eta_s0 <- qlogis(.2) + scales[["s0"]] * z_s0
  s0 <- plogis(eta_s0)
  # Enforce the fitted model's global price-centering convention exactly.
  # The independent market intercept innovations retain known SD .25.
  terms <- weights <- numeric(M * T)
  for (m in seq_len(M)) for (tt in seq_len(T)) {
    index <- (m - 1L) * T + tt
    terms[index] <- (eta_s0[m] + .25 * z_intercept[m] +
      hierarchy_lse(quality[m, ] - alpha[m] * relative_price[m, tt, ]))/alpha[m] +
      mean(relative_price[m, tt, ])
    weights[index] <- 1/alpha[m]
  }
  mu_intercept <- -mean(terms)/mean(weights)
  intercept <- mu_intercept + .25 * z_intercept
  price_center <- 2.2
  meanval <- quality + intercept + alpha * price_center
  rows <- truth_rows <- list()
  for (m in seq_len(M)) for (tt in seq_len(T)) {
    relative <- relative_price[m, tt, ]
    offset <- (eta_s0[m] + intercept[m] +
      hierarchy_lse(quality[m, ] - alpha[m] * relative))/alpha[m]
    price <- price_center + relative + offset
    score <- meanval[m, ] - alpha[m] * price
    shares <- exp(score - hierarchy_lse(c(0, score)))
    conditional <- shares/sum(shares)
    markup <- 1/(alpha[m] * (1 - shares))
    cost <- price - markup
    margin <- markup/price
    observed <- exp(log(conditional) + .025 * z_share[m, tt, ])
    observed <- observed/sum(observed)
    id <- (m - 1L) * T + tt
    rows[[id]] <- data.frame(prices = price, quantities = 1000 * observed,
      margins = margin * exp(.035 * z_margin[m, tt, ]), market = paste0("m", m),
      year = paste0("t", tt), group = paste0("f", seq_len(J)), stringsAsFactors = FALSE)
    truth_rows[[id]] <- data.frame(market = paste0("m", m), year = paste0("t", tt),
      group = paste0("f", seq_len(J)), price = price, conditional_share = conditional,
      share = shares, margin = margin, cost = cost, meanval = meanval[m, ], alpha = alpha[m],
      s0 = s0[m], stringsAsFactors = FALSE)
  }
  data <- do.call(rbind, rows); truth <- do.call(rbind, truth_rows)
  stopifnot(abs(mean(data$prices) - price_center) < 1e-12,
    all(truth$cost > 0), all(data$margins > 0 & data$margins < 1),
    max(abs(truth$price * truth$margin - 1/(truth$alpha * (1 - truth$share)))) < 1e-12)
  for (m in seq_len(M)) for (tt in seq_len(T)) {
    cell <- truth[truth$market == paste0("m", m) & truth$year == paste0("t", tt), ]
    stopifnot(abs(1 - sum(cell$share) - s0[m]) < 1e-12)
  }
  list(case = case, data = data, truth = truth, alpha = alpha, s0 = s0,
    sigma_alpha = scales[["alpha"]], sigma_s0 = scales[["s0"]], sigma_intercept = .25,
    realized_sigma_alpha = sd(log(alpha)), realized_sigma_s0 = sd(qlogis(s0)),
    innovations = list(alpha = z_alpha, s0 = z_s0, intercept = z_intercept,
      quality = quality, relative_price = relative_price, share = z_share, margin = z_margin),
    seed = seed, n_markets = M, n_periods = T, n_products = J)
}

hierarchy_oracle <- function(fixture) {
  errors <- list(); target <- NULL
  for (m in seq_len(fixture$n_markets)) for (tt in seq_len(fixture$n_periods)) {
    cell <- fixture$truth[fixture$truth$market == paste0("m", m) &
      fixture$truth$year == paste0("t", tt), ]
    model <- antitrust::specify(demand = "logit", conduct = "bertrand",
      prices = cell$price, parameters = list(alpha = -cell$alpha[1L], meanval = cell$meanval),
      ownerPre = cell$group, insideSize = 1000, priceOutside = 0,
      priceStart = cell$price * 1.17 + .05, labels = cell$group, output = TRUE)
    solved <- as.numeric(antitrust::calcPrices(model@model, preMerger = TRUE))
    reproduced_shares <- as.numeric(antitrust::calcShares(model@model, preMerger = TRUE))
    recovered_cost <- as.numeric(antitrust::calcMC(model@model, preMerger = TRUE))
    errors[[length(errors) + 1L]] <- data.frame(case = fixture$case, market = m, year = tt,
      price_error = max(abs(solved - cell$price)), share_error = max(abs(reproduced_shares - cell$share)),
      cost_error = max(abs(recovered_cost - cell$cost)))
    if (m == 1L && tt == 1L) target <- model
  }
  errors <- do.call(rbind, errors)
  stopifnot(max(as.matrix(errors[c("price_error", "share_error", "cost_error")])) < 1e-7)
  cell <- fixture$truth[fixture$truth$market == "m1" & fixture$truth$year == "t1", ]
  target_pair <- cell$group[order(-cell$conditional_share)[1:2]]
  owner_post <- cell$group; owner_post[owner_post == target_pair[2L]] <- target_pair[1L]
  observed <- fixture$data[fixture$data$market == "m1" & fixture$data$year == "t1", ]
  weights <- observed$quantities/sum(observed$quantities)
  policy <- list(owner = owner_post, tariff = ifelse(cell$group %in% target_pair, 0, .25/1.25),
    weights = weights, observed_price = observed$prices, pair = target_pair)
  effect <- hierarchy_policy_effects(target, policy)
  stopifnot(all(effect$status == "ok"))
  list(validation = errors, policy = policy, effects = effect)
}

hierarchy_policy_effects <- function(model, policy) {
  trade_fit <- trade::as_trade_fit(model, policy = "tariff", tariffPre = rep(0, length(policy$owner)),
    cost_basis = "effective", margin_basis = "net_revenue")
  results <- lapply(c("merger", "tariff25"), function(regime) {
    tryCatch({
      solved <- trade::simulate(trade_fit,
        ownerPost = if (regime == "merger") policy$owner else NULL,
        tariffPost = if (regime == "tariff25") policy$tariff else rep(0, length(policy$owner)))
      pre <- trade::tariff_accounts(solved, preMerger = TRUE)
      post <- trade::tariff_accounts(solved, preMerger = FALSE)
      stopifnot(all(is.finite(post$price)), all(post$price > 0),
        max(abs(pre$physical_cost - post$physical_cost)) < 1e-7)
      data.frame(regime = regime, status = "ok", price_pct =
        100 * sum(policy$weights * (post$price - pre$price))/sum(policy$weights * policy$observed_price),
        reason = "", stringsAsFactors = FALSE)
    }, error = function(e) data.frame(regime = regime, status = "failed", price_pct = NA_real_,
      reason = conditionMessage(e), stringsAsFactors = FALSE))
  })
  do.call(rbind, results)
}

hierarchy_fit_arguments <- function(fixture, family, directory, seed) {
  # Keep the original comparison's intercept prior fixed across package versions.
  priors <- list(prior_sigma_b_market_family = "half_normal", prior_sigma_b_market_scale = .5)
  if (family == "half_normal") priors <- c(priors, list(prior_sigma_alpha_family = "half_normal",
    prior_sigma_alpha = .5, prior_sigma_s0_pool_family = "half_normal", prior_sigma_s0_pool = .5))
  c(as.list(fixture$data), list(demand_model = "logit", market_side = "output",
    supply_model = "bertrand", use_rho = FALSE, use_alpha_market_pool = TRUE,
    s0_pooling = "market", use_centered_alpha = FALSE, use_centered_s0_pool = FALSE,
    share_likelihood = "gaussian", margin_likelihood = "gaussian", quality_scope = "market",
    quality_parameterization = "collapsed", market_intercept_parameterization = "collapsed",
    priors = priors, stan_options = list(chains = 4L, iter = 4000L, warmup = 2000L,
      cores = 4L, seed = seed, init = 0, refresh = 500L, output_dir = directory,
      control = list(metric = "dense_e", adapt_delta = .995, max_treedepth = 14L))))
}

hierarchy_recovery <- function(fit, fixture, family) {
  extract <- function(name) {
    draws <- as.matrix(fit@stanfit$draws(variables = name, format = "matrix"))
    # CmdStanR may return draws_matrix, whose column indexing retains an array.
    matrix(as.numeric(draws), nrow = nrow(draws), ncol = ncol(draws), dimnames = dimnames(draws))
  }
  alpha <- extract("alpha_structural_market"); outside <- extract("s0_prob")
  market_names <- sort(unique(as.character(fit@data$market)))
  parameter <- function(draws, name, truth, realized_truth = NA_real_) {
    q <- as.numeric(quantile(draws, c(.025, .5, .975)))
    data.frame(case = fixture$case, family = family, parameter = name, truth = truth,
      realized_truth = realized_truth, median = q[2L], q025 = q[1L], q975 = q[3L],
      error = q[2L] - truth, covered_95 = q[1L] <= truth && truth <= q[3L])
  }
  rows <- list()
  for (m in seq_along(market_names)) {
    true_m <- match(market_names[m], paste0("m", seq_along(fixture$alpha)))
    column <- which(fit@stan_data$mky_to_market == m)[1L]
    rows[[length(rows) + 1L]] <- parameter(alpha[, m], paste0("alpha_", market_names[m]), fixture$alpha[true_m])
    rows[[length(rows) + 1L]] <- parameter(outside[, column], paste0("s0_", market_names[m]), fixture$s0[true_m])
    all_years <- which(fit@stan_data$mky_to_market == m)
    stopifnot(max(abs(outside[, all_years, drop = FALSE] - outside[, column])) < 1e-12)
  }
  rows[[length(rows) + 1L]] <- parameter(extract("sigma_log_a_vec")[, 1L], "sigma_log_alpha",
    fixture$sigma_alpha, fixture$realized_sigma_alpha)
  rows[[length(rows) + 1L]] <- parameter(extract("sigma_s0_pool_vec")[, 1L], "sigma_logit_s0",
    fixture$sigma_s0, fixture$realized_sigma_s0)
  rows[[length(rows) + 1L]] <- parameter(extract("sigma_b_market_vec")[, 1L], "sigma_market_intercept",
    fixture$sigma_intercept, sd(.25 * fixture$innovations$intercept))
  do.call(rbind, rows)
}

hierarchy_diagnostics <- function(fit, fixture, family, elapsed) {
  variables <- c("alpha_structural_market", "s0_prob", "sigma_log_a_vec", "sigma_s0_pool_vec",
    "sigma_b_market_vec", "sigma_share", "sigma_margin")
  summary <- fit@stanfit$summary(variables = variables)
  full_summary <- fit@stanfit$summary()
  sampler <- fit@stanfit$sampler_diagnostics(format = "draws_array")
  bfmi <- vapply(seq_len(dim(sampler)[2L]), function(chain) {
    energy <- sampler[, chain, "energy__"]; mean(diff(energy)^2)/var(energy)
  }, numeric(1))
  data.frame(case = fixture$case, family = family, status = "completed", elapsed_seconds = elapsed,
    divergences = sum(sampler[, , "divergent__"]), transitions = prod(dim(sampler)[1:2]),
    max_rhat = max(full_summary$rhat, na.rm = TRUE), min_bulk_ess = min(full_summary$ess_bulk, na.rm = TRUE),
    min_tail_ess = min(full_summary$ess_tail, na.rm = TRUE), min_bfmi = min(bfmi),
    selected_max_rhat = max(summary$rhat, na.rm = TRUE), selected_min_bulk_ess = min(summary$ess_bulk, na.rm = TRUE),
    selected_min_tail_ess = min(summary$ess_tail, na.rm = TRUE), n_summary_variables = nrow(full_summary),
    depth_hits = sum(sampler[, , "treedepth__"] >= 14), reason = "")
}

hierarchy_policies <- function(fit, fixture, oracle, family) {
  plan <- antitrustBayes::posterior_policy_draws(fit, n_draws = 200L, seed = 171951L)
  prepared <- antitrustBayes::posterior_policy_models(fit, draw_plan = plan,
    market = "m1", year = "t1", labels = paste0("f", seq_len(fixture$n_products)),
    insideSize = 1000, demand_anchoring = "latent_quality")
  models <- prepared$models[[1L]]; rows <- list(); identities <- list()
  for (pos in seq_along(models)) {
    if (is.null(models[[pos]])) next
    identity <- attr(models[[pos]], "bayes_policy_draw", exact = TRUE)
    stopifnot(identical(identity$draw_index, plan$draw_indices[pos]))
    identities[[length(identities) + 1L]] <- identity
    result <- tryCatch(hierarchy_policy_effects(models[[pos]], oracle$policy), error = function(e)
      data.frame(regime = c("merger", "tariff25"), status = "failed", price_pct = NA_real_,
        reason = conditionMessage(e)))
    result$draw_position <- pos; result$draw_index <- plan$draw_indices[pos]
    result$plan_fingerprint <- identity$plan_fingerprint
    result$quality_hash <- digest::digest(identity$quality, algo = "sha256")
    rows[[length(rows) + 1L]] <- result
  }
  draws <- if (length(rows)) do.call(rbind, rows) else data.frame()
  paired <- if (nrow(draws)) intersect(draws$draw_position[draws$regime == "merger" & draws$status == "ok"],
    draws$draw_position[draws$regime == "tariff25" & draws$status == "ok"]) else integer()
  for (field in c("draw_index", "plan_fingerprint", "quality_hash"))
    if (length(paired)) stopifnot(identical(
      draws[draws$regime == "merger" & draws$draw_position %in% paired, field],
      draws[draws$regime == "tariff25" & draws$draw_position %in% paired, field]))
  summary <- do.call(rbind, lapply(c("merger", "tariff25"), function(regime) {
    x <- if (length(paired)) draws$price_pct[draws$regime == regime & draws$draw_position %in% paired] else numeric()
    q <- if (length(x)) as.numeric(quantile(x, c(.025, .975))) else c(NA_real_, NA_real_)
    true <- oracle$effects$price_pct[oracle$effects$regime == regime]
    data.frame(case = fixture$case, family = family, regime = regime, oracle = true,
      mean = if (length(x)) mean(x) else NA_real_, q025 = q[1L], q975 = q[2L],
      error = if (length(x)) mean(x) - true else NA_real_, covered_95 = q[1L] <= true && true <= q[2L],
      n_planned = length(plan$draw_indices), n_realization_failures = nrow(prepared$failures),
      n_policy_failures = if (nrow(draws)) sum(draws$status != "ok") else 0L, n_paired = length(paired),
      status = if (length(paired)) "completed" else "no_valid_pairs", reason = "")
  }))
  list(summary = summary, draws = draws, failures = prepared$failures, plan = plan, identities = identities)
}

hierarchy_validation_main <- function() {
  if (!requireNamespace("antitrustBayes", quietly = TRUE) || !requireNamespace("trade", quietly = TRUE))
    stop("Install current antitrustBayes, antitrust and trade before validation.")
  # R removes tempdir() on exit; keep the default output in its persistent parent.
  root <- Sys.getenv("AB_HIERARCHY_OUTPUT", file.path(dirname(tempdir()), "antitrustBayes-hierarchy-recovery"))
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  manifest_path <- file.path(root, "manifest.rds")
  resume <- identical(Sys.getenv("AB_HIERARCHY_RESUME", "false"), "true")
  existing <- file.exists(manifest_path)
  if (existing && !resume) stop("This output contains an inference attempt; use explicit resume without resetting its deadline.")
  if (resume && !existing) stop("No inference manifest exists to resume.")
  scope <- Sys.getenv("AB_HIERARCHY_SCOPE", "full")
  if (!scope %in% c("prepare", "full")) stop("AB_HIERARCHY_SCOPE must be prepare or full.")
  if (existing) {
    if (scope != "full") stop("Resume uses full scope and its original saved fixtures.")
    saved <- readRDS(file.path(root, "fixtures.rds")); fixtures <- saved$fixtures; oracles <- saved$oracles
  } else {
    fixtures <- lapply(c("small", "moderate"), hierarchy_fixture)
    names(fixtures) <- c("small", "moderate")
    stopifnot(identical(fixtures$small$innovations, fixtures$moderate$innovations))
    oracles <- lapply(fixtures, hierarchy_oracle)
    saveRDS(list(fixtures = fixtures, oracles = oracles), file.path(root, "fixtures.rds"))
    write.csv(do.call(rbind, lapply(oracles, `[[`, "validation")), file.path(root, "equilibrium_validation.csv"), row.names = FALSE)
  }
  if (scope == "prepare") { cat("Two common-innovation fixtures and 64 fresh equilibrium solves validated; no HMC.\n"); return(invisible(root)) }
  path <- Sys.getenv("AB_HIERARCHY_CMDSTAN", "")
  if (nzchar(path)) cmdstanr::set_cmdstan_path(path)
  options(antitrustBayes.cmdstan_cache = Sys.getenv("AB_HIERARCHY_CACHE", file.path(root, "stan-cache")))
  started <- Sys.time()
  manifest <- if (existing) readRDS(manifest_path) else list(started = started, deadline = started + 7200, package_version = as.character(packageVersion("antitrustBayes")),
    purpose = "paired_Bertrand_known_truth_prior_validation_not_paper_production",
    fixtures_hash = digest::digest(fixtures, algo = "sha256"), families = c("gamma", "half_normal"),
    scope = scope, status = "started", tasks = expand.grid(case = names(fixtures), family = c("gamma", "half_normal")),
    note = "Package defaults Gamma(2,rate=5) for both SDs; historical HalfNormal(.5). Intercept prior unchanged.")
  if (existing) {
    stopifnot(identical(manifest$package_version, as.character(packageVersion("antitrustBayes"))),
      identical(manifest$fixtures_hash, digest::digest(fixtures, algo = "sha256")))
    if (Sys.time() >= manifest$deadline) stop("The original validation deadline has expired.")
  }
  saveRDS(manifest, manifest_path)
  diagnostics <- recovery <- policy <- evidence <- list(); reference_data <- list()
  flush <- function() {
    for (name in c("diagnostics", "recovery", "policy", "evidence")) {
      value <- get(name, inherits = TRUE)
      if (length(value)) write.csv(do.call(rbind, value), file.path(root, paste0(name, ".csv")), row.names = FALSE)
    }
  }
  save_fit <- function(fit, path) {
    temporary <- tempfile(tmpdir = dirname(path), fileext = ".rds")
    saveRDS(fit, temporary)
    if (!file.rename(temporary, path)) stop("Could not publish saved fit: ", path)
  }
  for (i in seq_len(nrow(manifest$tasks))) {
    task <- manifest$tasks[i, ]; case <- as.character(task$case); family <- as.character(task$family)
    if (Sys.time() >= manifest$deadline) break
    fixture <- fixtures[[case]]; oracle <- oracles[[case]]
    directory <- file.path(root, paste(case, family, sep = "_"))
    argument_path <- file.path(directory, "arguments.rds")
    fit_path <- file.path(directory, "fit.rds")
    args <- hierarchy_fit_arguments(fixture, family, file.path(directory, "cmdstan"), seed = 518921L)
    if (file.exists(argument_path)) {
      stopifnot(identical(readRDS(argument_path), args))
      if (!file.exists(fit_path)) stop("An attempted fit has no saved posterior; preserve its failure instead of silently refitting: ", directory)
    } else {
      dir.create(directory)
      dir.create(args$stan_options$output_dir)
      saveRDS(args, argument_path)
    }
    cat("START", case, family, format(Sys.time(), tz = "UTC"), "\n")
    task_started <- proc.time()[["elapsed"]]
    warnings <- character()
    fit <- if (file.exists(fit_path)) {
      cat("REUSE saved posterior; no HMC repeat:", fit_path, "\n")
      readRDS(fit_path)
    } else withCallingHandlers(tryCatch(do.call(antitrustBayes::bayes_est, args), error = identity),
      warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
    if (length(warnings) || !file.exists(file.path(directory, "warnings.txt")))
      writeLines(warnings, file.path(directory, "warnings.txt"))
    if (inherits(fit, "error")) {
      diagnostics[[i]] <- data.frame(case, family, status = "failed", elapsed_seconds = proc.time()[["elapsed"]] - task_started,
        divergences = NA_integer_, transitions = NA_integer_, max_rhat = NA_real_, min_bulk_ess = NA_real_,
        min_tail_ess = NA_real_, min_bfmi = NA_real_, selected_max_rhat = NA_real_,
        selected_min_bulk_ess = NA_real_, selected_min_tail_ess = NA_real_, n_summary_variables = NA_integer_,
        depth_hits = NA_integer_, reason = conditionMessage(fit))
      flush(); next
    }
    # Verify that the changed priors reach active Stan hierarchy parameters.
    stopifnot(fit@stan_data$use_alpha_market_pool == 1L, fit@stan_data$use_s0_market_pool == 1L,
      fit@stan_data$prior_sigma_alpha_family == if (family == "gamma") 2L else 1L,
      fit@stan_data$prior_sigma_s0_pool_family == if (family == "gamma") 2L else 1L,
      fit@stan_data$prior_sigma_b_market_family == 1L)
    fields <- setdiff(names(fit@stan_data), grep("^prior_", names(fit@stan_data), value = TRUE))
    if (is.null(reference_data[[case]])) reference_data[[case]] <- fit@stan_data[fields] else
      stopifnot(isTRUE(all.equal(reference_data[[case]], fit@stan_data[fields], tolerance = 0)))
    save_fit(fit, file.path(directory, "fit.rds"))
    elapsed <- fit@stanfit$time()$total
    diagnostics[[i]] <- hierarchy_diagnostics(fit, fixture, family, elapsed)
    recovery[[i]] <- hierarchy_recovery(fit, fixture, family)
    flush()
    policy_result <- tryCatch(hierarchy_policies(fit, fixture, oracle, family), error = identity)
    if (inherits(policy_result, "error")) {
      saveRDS(list(status = "failed", reason = conditionMessage(policy_result)), file.path(directory, "policy.rds"))
      policy[[i]] <- data.frame(case, family, regime = c("merger", "tariff25"),
        oracle = oracle$effects$price_pct, mean = NA_real_, q025 = NA_real_, q975 = NA_real_,
        error = NA_real_, covered_95 = NA, n_planned = 200L, n_realization_failures = NA_integer_,
        n_policy_failures = NA_integer_, n_paired = 0L, status = "failed", reason = conditionMessage(policy_result))
    } else {
      policy[[i]] <- policy_result$summary; saveRDS(policy_result, file.path(directory, "policy.rds"))
    }
    set.seed(18815L)
    comparison <- tryCatch(antitrustBayes::bayes_compare(Bertrand = fit, bridge = TRUE,
      bridge_repetitions = 3L, bridge_maxiter = 5000L, bridge_method = "normal", prior_weights = 1), error = identity)
    evidence[[i]] <- if (inherits(comparison, "error")) data.frame(case, family, logml = NA_real_,
      logml_sd = NA_real_, status = "failed", reason = conditionMessage(comparison)) else data.frame(
      case, family, logml = comparison$logml[1L], logml_sd = comparison$logml_sd[1L], status = "completed", reason = "")
    saveRDS(comparison, file.path(directory, "evidence.rds"))
    save_fit(fit, file.path(directory, "fit.rds"))
    flush(); cat("FINISH", case, family, "\n"); print(diagnostics[[i]], row.names = FALSE)
  }
  manifest$status <- if (length(diagnostics) == nrow(manifest$tasks)) "completed" else "budget_exhausted"
  manifest$finished <- Sys.time(); manifest$tasks_finished <- length(diagnostics)
  saveRDS(manifest, manifest_path)
  cat("Results:", normalizePath(root), "\n")
  print(do.call(rbind, diagnostics), row.names = FALSE)
  invisible(root)
}

if (sys.nframe() == 0L) hierarchy_validation_main()
