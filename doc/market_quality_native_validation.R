#!/usr/bin/env Rscript
# Fixed-parameter validation only. No HMC, compilation, or covariance >100 rows.
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(antitrustBayes))
suppressPackageStartupMessages(library(rstan))
baseline_mode <- identical(Sys.getenv("MQ_BASELINE", "0"), "1")
stopifnot(antitrustBayes:::.antitrust_bayes_model_spec_version ==
  if (baseline_mode) 9L else 10L)
repo <- normalizePath(Sys.getenv("MQ_REPO", "/mnt/d/projects/antitrustbayes"))
out <- Sys.getenv("MQ_OUT", file.path(repo, "doc/statistical_toolkit_runs/20260916",
  "performance-audit/recovery/market-quality/native-validation"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
ref <- new.env(parent = globalenv())
sys.source(file.path(repo, "doc/market_quality_reference.R"), ref)

# Import only already-reviewed utility definitions, never the old driver body.
utility_names <- c("cga_as_draw", "cga_init_one", "cga_set_all", "cga_common_init",
  "cga_parameter_labels", "cga_model_for", "cga_fixed", "cga_log_target")
utility <- new.env(parent = globalenv())
for (expr in parse(file.path(repo, "doc/collapsed_gaussian_production_validation.R"))) {
  if (is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
      is.symbol(expr[[2L]]) && as.character(expr[[2L]]) %in% utility_names)
    eval(expr, utility)
}
foc_names <- c("cga_supply_name", "cga_role_conduct_logit", "cga_role_inverse_lerner_ces")
for (expr in parse(file.path(repo, "doc/collapsed_gaussian_benchmark.R"))) {
  if (is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
      is.symbol(expr[[2L]]) && as.character(expr[[2L]]) %in% foc_names)
    eval(expr, utility)
}

records <- list()
record <- function(case, check, error, tolerance = 2e-9, note = "") {
  records[[length(records) + 1L]] <<- data.frame(case = case, check = check,
    error = error, tolerance = tolerance,
    passed = is.finite(error) && error <= tolerance, note = note)
  write.csv(do.call(rbind, records), file.path(out, "validation.csv"), row.names = FALSE)
}
lse <- ref$logsumexp
getdraw <- utility$cga_as_draw
paired <- new.env(parent = emptyenv())
fixture <- function(M = 2L, T = 2L, singleton = FALSE, entry = FALSE) {
  rows <- list()
  for (m in seq_len(M)) for (tt in seq_len(T)) {
    gg <- if (singleton) rep(1L, 2L) else if (m == 1L) c(1L, 1L, 2L, 3L) else
      c(1L, 1L, 2L, 3L, 4L)
    if (entry && tt == 2L && !singleton) gg <- gg[gg != 3L]
    j <- seq_along(gg)
    rows[[length(rows) + 1L]] <- data.frame(prices = 1.03 + .13*j + .02*m + .01*tt,
      quantities = 23 - 2*j + m + tt, margins = .22 + .018*j,
      market = m, year = tt, group = gg, margin_observed = j %% 3L != 2L)
  }
  data <- do.call(rbind, rows)
  data$margins[!data$margin_observed] <- NA_real_
  # Interleave rows so factorization cannot rely on input sorting.
  data <- data[order(seq_len(nrow(data)) %% 3L, seq_len(nrow(data))), ]
  as.list(data)
}
prepare <- function(fx, demand, conduct, structure, scope, mode, side = "output") {
  exact <- structure %in% c("fringe", "stackelberg")
  counts <- tapply(fx$group, interaction(fx$market, fx$year),
                   function(g) length(unique(g)))
  kmax <- min(3L, min(counts))
  ks <- if (!exact) NULL else if (structure == "fringe") 0:kmax else seq_len(kmax)
  kp <- if (is.null(ks)) NULL else seq_along(ks) / sum(seq_along(ks))
  args <- c(fx, list(supply_model = switch(conduct, bertrand = 1L,
    monopolistic = 2L, cournot = 3L, auction = 4L),
    market_side = side, use_cutoff = FALSE, use_rho = TRUE,
    use_alpha_market_pool = FALSE, use_s0_market_pool = FALSE,
    use_s0_market_year_pool = FALSE, use_market_resid_scale = FALSE,
    use_hmt = FALSE, priors = list(), X_s0 = NULL, outside_price = 1.3,
    likelihood_mode = 1L, loo_cluster = fx$market, estimate_df = FALSE,
    resid_df_fixed = 5, structure = structure,
    partition_prior = if (exact) "ordered_k" else "cutoff",
    k_scope = "market", k_support = ks, k_prior = kp,
    quality_scope = scope, share_likelihood = "gaussian",
    quality_parameterization = mode, margin_likelihood = "gaussian"))
  fn <- if (demand == "logit") antitrustBayes:::prepare_stan_data else
    antitrustBayes:::prepare_ces_stan_data
  if (baseline_mode) args$quality_scope <- NULL
  data <- do.call(fn, args)
  if (baseline_mode) {
    # Reference annotations only: the spec-9 native program does NOT consume
    # these fields. Baseline comparisons are restricted to historical M=1.
    stopifnot(data$N_market == 1L)
    data$quality_scope <- 1L
    data$N_quality <- data$N_group
    data$N_quality_block <- 1L
    data$N_quality_coord <- data$N_group-1L
    data$quality_block_start <- 1L
    data$quality_block_n <- data$N_group
    data$obs_to_quality <- data$group
  }
  data
}
bundle <- function(data, demand, conduct, structure, seed = 20260917L, overrides = list()) {
  model <- utility$cga_model_for(demand, structure)
  initial_fit <- utility$cga_fixed(model, data, init = 0, seed = seed)
  init <- utility$cga_common_init(initial_fit, demand, conduct, rho = .31,
                                 tau = .43, s0 = .21)
  init <- utility$cga_set_all(init, "sigma_b_market_vec", .43)
  init <- utility$cga_set_all(init, "mu_b_market_vec", .09)
  init <- utility$cga_set_all(init, "b_market_raw", seq(-.1, .13,
    length.out = length(init$b_market_raw)))
  init <- utility$cga_set_all(init, "b_group_centered", .07)
  for (name in names(overrides)) init <- utility$cga_set_all(init, name, overrides[[name]])
  fit <- utility$cga_fixed(model, data, init = list(init), seed = seed + 1L)
  list(fit = fit, init = init, u = rstan::unconstrain_pars(fit, init), data = data)
}
common_prior <- function(b, demand, mode) {
  d <- b$data; i <- b$init; f <- b$fit
  dn <- function(x, mu, sd) sum(dnorm(as.numeric(x), mu, sd, log = TRUE))
  positive <- function(x, mu, sd) dn(x, mu, sd) - length(x)*pnorm(0, mu, sd,
    lower.tail = FALSE, log.p = TRUE)
  p <- dn(getdraw(f, "mu_log_a"), d$prior_alpha_mean, d$prior_alpha_sd)
  p <- p + dn(i$mu_b_market_vec, 0, 1) + positive(i$sigma_b_market_vec, 0, .5) +
    dn(i$b_market_raw, 0, 1)
  p <- p + dn(i$b_group_centered, 0, i$sigma_b_strat) +
    positive(i$sigma_b_strat, 1, d$prior_sigma_meanval_strat)
  p <- p + dn(i$s0_logit, d$s0_prior_mean, d$s0_prior_sd) +
    dn(i$beta_s0, 0, d$prior_sigma_beta_s0) +
    dn(i$log_sigma_share, log(d$prior_sigma_share), d$prior_log_sigma_share_sd) +
    dn(i$log_sigma_margin, log(d$prior_sigma_margin), d$prior_log_sigma_margin_sd)
  if (length(i$rho_raw)) {
    rho <- getdraw(f, "rho")
    p <- p + (d$prior_lkj-1)*log1p(-rho^2) + log1p(-tanh(i$rho_raw)^2) -
      lbeta(.5, d$prior_lkj) - pbeta(.95^2, .5, d$prior_lkj, log.p = TRUE) + log(.95)
  }
  p
}
predictions <- function(b, demand, conduct, structure, k = NULL) {
  d <- b$data; f <- b$fit
  s0 <- getdraw(f, "s0_prob")
  coef <- getdraw(f, if (demand == "logit") "alpha_structural_market" else
                    "gamma_structural_market")
  ans <- numeric(d$N)
  for (mky in seq_len(d$N_market_year)) {
    ii <- which(d$market_year_idx == mky)
    m <- d$mky_to_market[mky]
    if (structure == "simultaneous") {
      metric <- if (demand == "logit") d$firm_shareIn[ii] else d$firm_revShare[ii]
      if (demand == "logit") ans[ii] <- d$price[ii]*coef[m]*
        switch(conduct, bertrand = (1-(1-s0[mky])*metric),
          cournot = s0[mky]/(s0[mky]+(1-s0[mky])*metric),
          monopolistic = rep(1, length(ii)),
          auction = {S <- (1-s0[mky])*metric; ifelse(S<1e-10, 1, S/(-log1p(-S)))})
      else {
        R0 <- plogis(qlogis(s0[mky])+log(d$outside_price[mky])-log(d$inside_price_index[mky]))
        R <- metric*(1-R0); gamma <- coef[m]
        ans[ii] <- switch(conduct, bertrand = gamma-(gamma-1)*R,
          cournot = 1/(1/gamma+(gamma-1)*R/(gamma*(1+(gamma-1)*R0))),
          monopolistic = rep(gamma, length(ii)))
      }
    } else {
      jj <- d$firm_mky_start[mky]+seq_len(d$firm_mky_n[mky])-1L
      metric <- if (demand == "logit") d$firm_mky_shareIn[jj] else d$firm_mky_revShare[jj]
      rank <- d$firm_mky_rank[jj]
      map <- match(d$obs_to_firm_mky[ii], jj)
      if (demand == "logit") ans[ii] <- d$price[ii]*coef[m]*
        utility$cga_role_conduct_logit(metric, s0[mky], conduct, structure, rank, k)[map]
      else {
        R0 <- plogis(qlogis(s0[mky])+log(d$outside_price[mky])-log(d$inside_price_index[mky]))
        ans[ii] <- utility$cga_role_inverse_lerner_ces(metric, R0, coef[m],
          conduct, structure, rank, k)[map]
      }
    }
  }
  log(ans+1e-9)
}
likelihood_reference <- function(b, demand, conduct, structure, mode) {
  d <- b$data; f <- b$fit; rs <- d$log_shareIn-getdraw(f, "share_mean")
  ss <- getdraw(f, "sigma_share"); sm <- getdraw(f, "sigma_margin")
  tau <- getdraw(f, "sigma_b_strat"); rho <- getdraw(f, "rho")
  idx <- as.integer(d$margin_idx)
  by_market <- numeric(d$N_market); probabilities <- numeric(d$N_partition_component %||% 0)
  if (structure == "simultaneous") {
    mu <- predictions(b, demand, conduct, structure)
    rm <- d$log_marginNormInv[idx]-mu[idx]
    if (mode == "sampled") return(list(total = sum(dnorm(rs, 0, ss, log = TRUE))+
      sum(dnorm(rm-rho*sm/ss*rs[idx], 0, sm*sqrt(1-rho^2), log = TRUE))))
    layout <- ref$make_market_layout(d$market, d$group)
    if (d$quality_scope == 1L) layout <- ref$make_market_layout(rep(1L, d$N), d$group)
    total <- ref$dense_original_joint_log_likelihood(rs, rm, idx, layout, tau, ss, sm, rho)
    return(list(total = total))
  }
  for (m in seq_len(d$N_market)) {
    ii <- which(d$market == m); qq <- which(d$market[idx] == m)
    u <- which(d$partition_unit_market == m)
    cc <- d$partition_unit_start[u]+seq_len(d$partition_unit_count[u])-1L
    local <- numeric(length(cc))
    layout <- ref$make_market_layout(rep(1L, length(ii)), d$group[ii])
    for (kk in seq_along(cc)) {
      mu <- predictions(b, demand, conduct, structure, d$partition_component_k[cc[kk]])
      rm <- d$log_marginNormInv[idx[qq]]-mu[idx[qq]]
      if (mode == "collapsed") ll <- ref$dense_original_joint_log_likelihood(rs[ii], rm,
        match(idx[qq], ii), layout, tau, ss, sm, rho)
      else ll <- sum(dnorm(rs[ii], 0, ss, log = TRUE))+
        sum(dnorm(rm-rho*sm/ss*rs[idx[qq]], 0, sm*sqrt(1-rho^2), log = TRUE))
      local[kk] <- log(d$partition_prior_probability[cc[kk]])+ll
    }
    by_market[m] <- lse(local)
    probabilities[cc] <- exp(local-lse(local))
  }
  list(total = sum(by_market), by_market = by_market, probabilities = probabilities)
}

# Oracle only: small fixtures enumerate joint k states and integrate the full
# dense joint covariance. Production code must never use this enumeration.
joint_k_oracle <- function(b, demand, conduct, structure, mode) {
  d <- b$data; f <- b$fit
  stopifnot(d$N_market <= 3L, d$N <= 100L)
  units <- match(seq_len(d$N_market), d$partition_unit_market)
  components <- lapply(units, function(u)
    d$partition_unit_start[u]+seq_len(d$partition_unit_count[u])-1L)
  states <- expand.grid(lapply(components, seq_along), KEEP.OUT.ATTRS = FALSE)
  stopifnot(nrow(states) <= 125L)
  rs <- d$log_shareIn-getdraw(f, "share_mean")
  idx <- d$margin_idx; ss <- getdraw(f, "sigma_share")
  sm <- getdraw(f, "sigma_margin"); rho <- getdraw(f, "rho")
  layout <- ref$make_market_layout(d$market, d$group)
  scores <- numeric(nrow(states))
  for (s in seq_len(nrow(states))) {
    mu <- numeric(d$N); prior <- 0
    for (m in seq_len(d$N_market)) {
      cc <- components[[m]][states[s, m]]
      local_mu <- predictions(b, demand, conduct, structure, d$partition_component_k[cc])
      ii <- which(d$market == m); mu[ii] <- local_mu[ii]
      prior <- prior+log(d$partition_prior_probability[cc])
    }
    rm <- d$log_marginNormInv[idx]-mu[idx]
    ll <- if (mode == "collapsed") ref$dense_original_joint_log_likelihood(
      rs, rm, idx, layout, getdraw(f, "sigma_b_strat"), ss, sm, rho) else
      sum(dnorm(rs, 0, ss, log = TRUE))+
      sum(dnorm(rm-rho*sm/ss*rs[idx], 0, sm*sqrt(1-rho^2), log = TRUE))
    scores[s] <- prior+ll
  }
  lse(scores)
}

# At arbitrary sampled beta, expanded_target = collapsed_target + the
# conditional beta log density. Exact k requires a mixture of Gaussian quality
# posteriors within each market, rather than a single Gaussian/MAP state.
check_gaussian_completion <- function(sampled, collapsed, case, demand, conduct, structure) {
  d <- collapsed$data; f <- collapsed$fit
  rs <- d$log_shareIn-getdraw(f, "share_mean")
  ss <- getdraw(f, "sigma_share"); sm <- getdraw(f, "sigma_margin")
  tau <- getdraw(f, "sigma_b_strat"); rho <- getdraw(f, "rho")
  post <- likelihood_reference(collapsed, demand, conduct, structure, "collapsed")
  posterior_log_density <- 0; beta_start <- 1L
  for (m in seq_len(d$N_market)) {
    ii <- which(d$market == m); qq <- which(d$market[d$margin_idx] == m)
    local_layout <- ref$make_market_layout(rep(1L, length(ii)), d$group[ii])
    G <- local_layout$n_group[[1L]]; dimension <- G-1L
    beta <- sampled$init$b_group_centered[beta_start+seq_len(dimension)-1L]
    beta_start <- beta_start+dimension
    cc <- if (structure == "simultaneous") NA_integer_ else {
      u <- which(d$partition_unit_market == m)
      d$partition_unit_start[u]+seq_len(d$partition_unit_count[u])-1L
    }
    log_densities <- numeric(length(cc))
    for (j in seq_along(cc)) {
      k <- if (is.na(cc[j])) NULL else d$partition_component_k[cc[j]]
      mu <- predictions(collapsed, demand, conduct, structure, k)
      rm <- d$log_marginNormInv[d$margin_idx[qq]]-mu[d$margin_idx[qq]]
      moments <- ref$posterior_quality_fixed_state_dense(rs[ii], rm,
        match(d$margin_idx[qq], ii), local_layout, tau, ss, sm, rho)[[1L]]
      log_densities[j] <- if (!dimension) 0 else ref$log_mvnorm_dense(beta,
        mean = moments$contrast_mean, covariance = moments$contrast_covariance)
      if (!is.na(cc[j])) log_densities[j] <- log_densities[j]+log(post$probabilities[cc[j]])
    }
    posterior_log_density <- posterior_log_density+lse(log_densities)
  }
  record(case, "native_expanded_collapsed_completion", abs(
    utility$cga_log_target(sampled$fit, sampled$u)-
    utility$cga_log_target(collapsed$fit, collapsed$u)-posterior_log_density), 2e-9)
}
# Locally needed null-coalescing utility, independent of package namespace.
`%||%` <- function(x, y) if (is.null(x)) y else x

check_quality_recovery <- function(b, case, fx, demand, conduct, structure,
                                   scope, side) {
  d <- b$data
  shell <- methods::new("BayesLogit", stanfit = b$fit, stan_data = d,
    data = as.data.frame(fx), price_sd = d$price_sd %||% d$log_price_sd,
    supply_model = as.integer(d$supply_model), demand_model = demand,
    market_side = side, call = quote(fixed_parameter_validation()),
    metadata = list(model_spec_version = 10L, quality_scope = scope,
      quality_parameterization = "collapsed", share_likelihood = "gaussian",
      margin_likelihood = "gaussian", posterior_draws = 1L,
      structure = structure, partition_prior = if (structure == "simultaneous")
        "cutoff" else "ordered_k", k_scope = "market"))
  pq <- posterior_quality(shell, seed = 542L)
  rs <- d$log_shareIn-getdraw(b$fit, "share_mean")
  ss <- getdraw(b$fit, "sigma_share"); sm <- getdraw(b$fit, "sigma_margin")
  rho <- getdraw(b$fit, "rho"); tau <- getdraw(b$fit, "sigma_b_strat")
  for (m in seq_len(d$N_market)) {
    ii <- which(d$market == m); qq <- which(d$market[d$margin_idx] == m)
    k <- if (structure == "simultaneous") NULL else pq$realized_k[1L, m]
    mu <- predictions(b, demand, conduct, structure, k)
    rm <- d$log_marginNormInv[d$margin_idx[qq]]-mu[d$margin_idx[qq]]
    layout <- ref$make_market_layout(rep(1L, length(ii)), d$group[ii])
    moments <- ref$posterior_quality_fixed_state_dense(rs[ii], rm,
      match(d$margin_idx[qq], ii), layout, tau, ss, sm, rho)[[1L]]
    block <- if (scope == "global") 1L else m
    units <- d$quality_block_start[block]+seq_len(d$quality_block_n[block])-1L
    record(case, paste0("quality_conditional_mean_market", m), max(abs(
      pq$conditional_means[1L, units]-moments$group_mean)), 4e-10)
    record(case, paste0("quality_conditional_variance_market", m), max(abs(
      pq$conditional_sds[1L, units]^2-diag(moments$group_covariance))), 4e-10)
    record(case, paste0("conditional_quality_draw_center_market", m),
      abs(sum(pq$draws[1L, units])), 2e-13)
  }
  saveRDS(pq, file.path(out, paste0(case, "-quality-recovery.rds")))
  invisible(pq)
}

check_case <- function(case, fx, demand, conduct, structure, scope, mode,
                       side = "output", gradients = TRUE, overrides = list(),
                       target_tolerance = 2e-9) {
  d <- prepare(fx, demand, conduct, structure, scope, mode, side)
  b <- bundle(d, demand, conduct, structure, overrides = overrides)
  saveRDS(list(data = d, init = b$init,
    model_code_sha256 = digest::digest(b$fit@stanmodel@model_code, algo = "sha256")),
    file.path(out, paste0(case, "-fixed-state.rds")))
  reference <- likelihood_reference(b, demand, conduct, structure, mode)
  target <- utility$cga_log_target(b$fit, b$u)
  saved <- list(target = target,
    bridge_target = rstan::log_prob(b$fit, b$u, adjust_transform = TRUE),
    gradient = rstan::grad_log_prob(b$fit, b$u, adjust_transform = TRUE),
    predictions = lapply(c("share_mean", "pred_shareIn", "pred_marginNormInv", "s0_prob", "a_market"),
      function(n) getdraw(b$fit, n)))
  saveRDS(saved, file.path(out, paste0(case, "-native-output.rds")))
  record(case, "normalized_full_target", abs(target-common_prior(b, demand, mode)-reference$total),
    target_tolerance)
  if (mode == "collapsed") {
    record(case, "generated_collapsed_joint_density", abs(
      getdraw(b$fit, "log_lik_joint_collapsed")-reference$total), target_tolerance)
    if (!baseline_mode && scope == "market") {
      record(case, "generated_market_joint_density_sum", abs(sum(
        getdraw(b$fit, "log_lik_joint_collapsed_market"))-reference$total), target_tolerance)
      record(case, "generated_market_density_valid", abs(
        getdraw(b$fit, "collapsed_market_lik_valid")-1), 0)
    }
  }
  baseline_path <- file.path(dirname(out), "native-baseline", paste0(case, "-native-output.rds"))
  if (!baseline_mode && scope == "global" && file.exists(baseline_path)) {
    old <- readRDS(baseline_path)
    record(case, "prechange_native_target", abs(saved$target-old$target), 2e-10)
    record(case, "prechange_bridge_target", abs(saved$bridge_target-old$bridge_target), 2e-10)
    record(case, "prechange_AD_gradient", max(abs(saved$gradient-old$gradient)), 2e-10)
    record(case, "prechange_GQ_predictions", max(abs(unlist(saved$predictions)-
      unlist(old$predictions))), 2e-10)
  }
  if (structure != "simultaneous") {
    record(case, "component_probabilities", max(abs(getdraw(b$fit, "posterior_component_prob")-
      reference$probabilities)))
    record(case, "market_probability_sums", max(abs(tapply(reference$probabilities,
      d$partition_component_unit, sum)-1)))
    if (!baseline_mode && scope == "market" && d$N_market > 1L)
      record(case, "brute_force_joint_k_dense_oracle", abs(reference$total-
        joint_k_oracle(b, demand, conduct, structure, mode)), 2e-9)
  } else record(case, "ordinary_inverse_margin_prediction",
    max(abs(log(getdraw(b$fit, "pred_marginNormInv")+1e-9)-
            predictions(b, demand, conduct, structure))))
  if (mode == "sampled") {
    quality <- getdraw(b$fit, "b_quality_all")
    for (block in seq_len(d$N_quality_block)) {
      jj <- d$quality_block_start[block]+seq_len(d$quality_block_n[block])-1L
      record(case, paste0("centering_block_", block), abs(sum(quality[jj])), 2e-13)
    }
    record(case, "contrast_dimension", abs(length(b$init$b_group_centered)-d$N_quality_coord), 0)
  }
  if (gradients) {
    ad <- as.numeric(rstan::grad_log_prob(b$fit, b$u, adjust_transform = FALSE))
    labels <- utility$cga_parameter_labels(b$fit)
    for (h in c(2e-5, 4e-6)) {
      fd <- vapply(seq_along(b$u), function(j) {
        plus <- minus <- b$u; plus[j] <- plus[j]+h; minus[j] <- minus[j]-h
        (utility$cga_log_target(b$fit, plus)-utility$cga_log_target(b$fit, minus))/(2*h)
      }, numeric(1))
      err <- abs(ad-fd)/(1+abs(ad))
      record(case, paste0("AD_finite_difference_h", h), max(err), 2e-5,
        paste(labels[which.max(err)], "absolute error", max(abs(ad-fd))))
      write.csv(data.frame(parameter = labels, ad = ad, fd = fd, relative_error = err),
        file.path(out, paste0(case, "-gradient-", h, ".csv")), row.names = FALSE)
    }
  }
  if (!baseline_mode) {
    key <- sub(paste0("_", mode, "_"), "_", case, fixed = TRUE)
    if (mode == "sampled") assign(key, b, paired)
    else if (exists(key, paired, inherits = FALSE))
      check_gaussian_completion(get(key, paired), b, key, demand, conduct, structure)
    if (mode == "collapsed") check_quality_recovery(b, case, fx, demand,
      conduct, structure, scope, side)
  }
  b
}

if (!"--extra-conduct-only" %in% commandArgs(trailingOnly = TRUE)) {
for (demand in c("logit", "ces")) for (conduct in c("bertrand", "cournot"))
  for (structure in c("simultaneous", "fringe", "stackelberg"))
    for (mode in c("sampled", "collapsed")) {
      case <- paste(demand, conduct, structure, mode, sep = "_")
      tryCatch({
        global <- check_case(paste0(case, "_single_global"), fixture(1L), demand,
          conduct, structure, "global", mode)
        if (baseline_mode) next
        market <- check_case(paste0(case, "_single_market"), fixture(1L), demand,
          conduct, structure, "market", mode)
        record(case, "single_market_target_equality", abs(utility$cga_log_target(global$fit, global$u)-
          utility$cga_log_target(market$fit, market$u)), 2e-10)
        record(case, "single_market_bridge_target_equality", abs(rstan::log_prob(global$fit,
          global$u, adjust_transform = TRUE)-rstan::log_prob(market$fit, market$u,
          adjust_transform = TRUE)), 2e-10)
        record(case, "single_market_gradient_equality", max(abs(
          rstan::grad_log_prob(global$fit, global$u, adjust_transform = TRUE)-
          rstan::grad_log_prob(market$fit, market$u, adjust_transform = TRUE))), 2e-10)
        for (name in c("share_mean", "pred_shareIn", "pred_marginNormInv", "s0_prob", "a_market"))
          record(case, paste0("single_market_GQ_", name), max(abs(
            getdraw(global$fit, name)-getdraw(market$fit, name))), 2e-10)
        check_case(paste0(case, "_multiple_market"), fixture(2L), demand, conduct,
          structure, "market", mode)
      }, error = function(e) record(case, "execution_error", Inf, 0, conditionMessage(e)))
    }

if (!baseline_mode) for (singleton in c(FALSE, TRUE)) for (entry in if (singleton) FALSE else c(FALSE, TRUE))
  tryCatch(check_case(paste0("input_singleton", singleton, "_entry", entry),
    fixture(2L, singleton = singleton, entry = entry), "logit", "bertrand",
    "simultaneous", "market", "collapsed", side = "input"),
    error = function(e) record("input_edge", "execution_error", Inf, 0, conditionMessage(e)))

if (!baseline_mode) for (edge in c("rho_zero", "tau_small", "share_scale_small", "missing_last_market")) {
  fx <- fixture(2L)
  overrides <- switch(edge, rho_zero = list(rho_raw = 0),
    tau_small = list(sigma_b_strat = 1e-7),
    share_scale_small = list(log_sigma_share = log(1e-4)), list())
  if (edge == "missing_last_market") {
    fx$margin_observed[fx$market == 2L] <- FALSE
    fx$margins[fx$market == 2L] <- NA_real_
  }
  tryCatch(check_case(paste0("edge_", edge), fx, "logit", "bertrand", "simultaneous",
    "market", "collapsed", overrides = overrides,
    target_tolerance = if (edge == "share_scale_small") 2e-5 else 2e-9),
    error = function(e) record(edge, "execution_error", Inf, 0, conditionMessage(e)))
}

if (!baseline_mode) for (demand in c("logit", "ces"))
  for (conduct in c("bertrand", "cournot"))
    for (structure in c("fringe", "stackelberg")) {
      case <- paste("entry_exit", demand, conduct, structure, sep = "_")
      tryCatch(check_case(case, fixture(2L, entry = TRUE), demand, conduct,
        structure, "market", "collapsed"),
        error = function(e) record(case, "execution_error", Inf, 0, conditionMessage(e)))
    }
}

if (!baseline_mode && "--extra-conduct-only" %in% commandArgs(trailingOnly = TRUE))
  for (demand in c("logit", "ces"))
    for (conduct in if (demand == "logit") c("monopolistic", "auction") else "monopolistic")
      for (scope in c("global", "market")) for (mode in c("sampled", "collapsed")) {
        case <- paste("extra", demand, conduct, mode, "single", scope, sep = "_")
        tryCatch(check_case(case, fixture(1L), demand, conduct, "simultaneous",
          scope, mode), error = function(e) record(case, "execution_error", Inf, 0,
            conditionMessage(e)))
      }

result <- do.call(rbind, records)
print(result)
cat("Checks:", nrow(result), "failed:", sum(!result$passed), "\n")
stopifnot(all(result$passed))
