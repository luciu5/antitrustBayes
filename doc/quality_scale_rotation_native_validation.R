#!/usr/bin/env Rscript

## Native Stan target, generated-quantity, and autodiff validation for the
## contrast-total-scale/share-variance-fraction coordinates.

options(stringsAsFactors = FALSE)
lib <- Sys.getenv("QUALITY_SCALE_LIB", unset = "/tmp/atb-scale-lib")
suppressPackageStartupMessages(library(antitrustBayes, lib.loc = lib))
suppressPackageStartupMessages(library(rstan))

make_fixture <- function() {
  n <- c(3L, 4L, 5L)
  market <- rep(seq_along(n), n)
  N <- length(market)
  list(
    prices = seq(.86, 1.42, length.out = N),
    quantities = unlist(lapply(n, function(k) seq(31, 11, length.out = k))),
    margins = .24 + .025 * sin(seq_len(N)),
    market = market, year = rep(1L, N),
    group = unlist(lapply(n, seq_len), use.names = FALSE),
    margin_observed = !seq_len(N) %in% c(2L, 7L, 12L),
    loo_cluster = market
  )
}

fx <- make_fixture()

model_for <- function(demand, exact) {
  antitrustBayes:::stanmodels[[paste0(
    demand, if (exact) "_stackelberg" else "_supply"
  )]]
}

prepare_data <- function(demand, exact, scale_coordinate) {
  common <- list(
    prices = fx$prices, quantities = fx$quantities, margins = fx$margins,
    market = fx$market, year = fx$year, group = fx$group,
    supply_model = 1L, use_cutoff = FALSE, market_side = "output",
    use_rho = TRUE, use_market_resid_scale = FALSE,
    use_alpha_market_pool = FALSE, use_s0_market_year_pool = FALSE,
    use_s0_market_pool = TRUE, use_centered_alpha = FALSE,
    use_centered_s0_pool = FALSE, market_s0_parameterization = "joint",
    use_hmt = FALSE, priors = list(), X_s0 = NULL,
    margin_observed = fx$margin_observed, outside_price = 1,
    likelihood_mode = 1L, loo_cluster = fx$loo_cluster,
    estimate_df = FALSE, resid_df_fixed = 5,
    structure = if (exact) "stackelberg" else "simultaneous",
    partition_prior = if (exact) "ordered_k" else "cutoff",
    k_support = if (exact) 1:2 else NULL,
    k_prior = if (exact) c(.4, .6) else NULL,
    k_scope = "market", share_likelihood = "gaussian",
    margin_likelihood = "gaussian", quality_parameterization = "collapsed",
    quality_scope = "market",
    quality_scale_parameterization = scale_coordinate
  )
  if (demand == "ces") {
    do.call(antitrustBayes:::prepare_ces_stan_data, common)
  } else {
    do.call(antitrustBayes:::prepare_stan_data, common)
  }
}

init_template <- function(model, data) {
  fit <- rstan::sampling(model, data = data, chains = 1L, iter = 1L,
    warmup = 0L, algorithm = "Fixed_param", init = 0,
    seed = 20260920L, refresh = 0)
  init <- rstan::get_inits(fit)
  if (is.list(init) && length(init) == 1L && is.list(init[[1L]])) {
    init <- init[[1L]]
  }
  list(fit = fit, init = init[intersect(fit@model_pars, names(init))])
}

set_one <- function(init, name, value) {
  if (!name %in% names(init)) return(init)
  old <- init[[name]]
  value <- rep(as.numeric(value), length.out = length(old))
  if (!is.null(dim(old))) dim(value) <- dim(old)
  init[[name]] <- value
  init
}

make_state <- function(template, coordinate, sigma, tau, rho) {
  init <- template
  init <- set_one(init, "mu_b_market_vec", .17)
  init <- set_one(init, "sigma_b_market_vec", .43)
  init <- set_one(init, "b_market_raw", c(-.7, .25, 1.1))
  init <- set_one(init, "mu_s0_pool_vec", -1.12)
  init <- set_one(init, "sigma_s0_pool_vec", .37)
  init <- set_one(init, "c_s0_raw", c(.45, -.8, .2))
  if (coordinate == "legacy") {
    init <- set_one(init, "sigma_b_strat", tau)
    init <- set_one(init, "log_sigma_share", log(sigma))
  } else {
    total <- sqrt(sigma^2 + tau^2)
    fraction_logit <- 2 * (log(sigma) - log(tau))
    init <- set_one(init, "sigma_b_strat", total)
    init <- set_one(init, "log_sigma_share", fraction_logit)
  }
  init <- set_one(init, "log_sigma_margin", log(.14))
  init <- set_one(init, "rho_raw", atanh(rho / .95))
  init <- set_one(init, "mu_log_a_free", log(.84))
  init <- set_one(init, "mu_log_a", log(.84))
  init
}

fixed_state <- function(model, data, init) {
  rstan::sampling(model, data = data, chains = 1L, iter = 1L,
    warmup = 0L, algorithm = "Fixed_param", init = list(init),
    seed = 20260921L, refresh = 0)
}

finite_gradient <- function(fit, u, rel_step = 1e-5) {
  vapply(seq_along(u), function(j) {
    h <- rel_step * max(1, abs(u[[j]]))
    up <- um <- u
    up[[j]] <- up[[j]] + h
    um[[j]] <- um[[j]] - h
    (rstan::log_prob(fit, up, adjust_transform = TRUE) -
       rstan::log_prob(fit, um, adjust_transform = TRUE)) / (2 * h)
  }, numeric(1))
}

rows <- list()
gradients <- list()
state_rows <- list()
scenarios <- list(
  ordinary = c(sigma = .35, tau = .8, rho = .4),
  small_share_positive_rho = c(sigma = .01, tau = 1.3, rho = .94),
  small_share_negative_rho = c(sigma = .01, tau = 1.3, rho = -.94),
  small_quality = c(sigma = .6, tau = .01, rho = .4),
  both_small = c(sigma = .012, tau = .008, rho = -.8)
)

for (demand in c("logit", "ces")) for (exact in c(FALSE, TRUE)) {
  kernel <- paste0(demand, if (exact) "_stackelberg" else "_supply")
  model <- model_for(demand, exact)
  old_data <- prepare_data(demand, exact, "legacy")
  new_data <- prepare_data(demand, exact, "total_fraction")
  old_template <- init_template(model, old_data)
  new_template <- init_template(model, new_data)
  for (scenario in names(scenarios)) {
    value <- scenarios[[scenario]]
    old_init <- make_state(old_template$init, "legacy",
      value[["sigma"]], value[["tau"]], value[["rho"]])
    new_init <- make_state(new_template$init, "total_fraction",
      value[["sigma"]], value[["tau"]], value[["rho"]])
    old_u <- rstan::unconstrain_pars(old_template$fit, old_init)
    new_u <- rstan::unconstrain_pars(new_template$fit, new_init)
    old_lp <- rstan::log_prob(old_template$fit, old_u,
      adjust_transform = TRUE)
    new_lp <- rstan::log_prob(new_template$fit, new_u,
      adjust_transform = TRUE)
    error <- abs(new_lp - old_lp + log(2))
    rows[[length(rows) + 1L]] <- data.frame(
      kernel = kernel, scenario = scenario, quantity = "normalized_target",
      legacy_log_target = old_lp,
      total_fraction_log_target = new_lp,
      max_abs_error = error, tolerance = 3e-8,
      passed = is.finite(error) && error <= 3e-8
    )

    old_draw <- fixed_state(model, old_data, old_init)
    new_draw <- fixed_state(model, new_data, new_init)
    fields <- intersect(c(
      "sigma_b_strat_economic", "sigma_share", "share_variance_fraction",
      "share_mean", "pred_marginNormInv", "s0_prob", "b_market",
      "posterior_component_prob", "log_lik_share_collapsed",
      "log_lik_joint_collapsed", "log_lik_share_collapsed_market",
      "log_lik_joint_collapsed_market"
    ), intersect(old_draw@model_pars, new_draw@model_pars))
    old_gq <- rstan::extract(old_draw, pars = fields, permuted = TRUE)
    new_gq <- rstan::extract(new_draw, pars = fields, permuted = TRUE)
    for (field in fields) {
      old_value <- as.numeric(old_gq[[field]])
      new_value <- as.numeric(new_gq[[field]])
      field_error <- max(abs(as.numeric(old_gq[[field]]) -
                             as.numeric(new_gq[[field]])))
      rows[[length(rows) + 1L]] <- data.frame(
        kernel = kernel, scenario = scenario,
        quantity = paste0("economic_", field),
        legacy_log_target = NA_real_, total_fraction_log_target = NA_real_,
        max_abs_error = field_error, tolerance = 3e-8,
        passed = is.finite(field_error) && field_error <= 3e-8
      )
      state_rows[[length(state_rows) + 1L]] <- data.frame(
        kernel = kernel, scenario = scenario, coordinate = "legacy",
        field = field, element = seq_along(old_value), value = old_value
      )
      state_rows[[length(state_rows) + 1L]] <- data.frame(
        kernel = kernel, scenario = scenario, coordinate = "total_fraction",
        field = field, element = seq_along(new_value), value = new_value
      )
    }

    if (scenario %in% c("ordinary", "small_share_negative_rho")) {
      ad <- as.numeric(rstan::grad_log_prob(new_template$fit, new_u,
        adjust_transform = TRUE))
      fd <- finite_gradient(new_template$fit, new_u)
      grad_error <- max(abs(ad - fd) / (1 + abs(ad)))
      gradients[[length(gradients) + 1L]] <- data.frame(
        kernel = kernel, scenario = scenario,
        max_relative_error = grad_error, tolerance = 3e-4,
        passed = is.finite(grad_error) && grad_error <= 3e-4
      )
    }
  }
}

targets <- do.call(rbind, rows)
gradient_rows <- do.call(rbind, gradients)
write.csv(targets, "doc/quality_scale_rotation_native_target_checks.csv",
          row.names = FALSE)
write.csv(gradient_rows,
          "doc/quality_scale_rotation_native_gradient_checks.csv",
          row.names = FALSE)
state_output <- Sys.getenv(
  "QUALITY_SCALE_STATE_OUTPUT",
  unset = "doc/quality_scale_rotation_native_state_values.csv"
)
write.csv(do.call(rbind, state_rows), state_output, row.names = FALSE)
cat(sprintf("native targets/GQ %d/%d; gradients %d/%d\n",
  sum(targets$passed), nrow(targets),
  sum(gradient_rows$passed), nrow(gradient_rows)))
if (any(!targets$passed) || any(!gradient_rows$passed)) {
  print(subset(targets, !passed))
  print(subset(gradient_rows, !passed))
  stop("native quality-scale validation failed")
}
