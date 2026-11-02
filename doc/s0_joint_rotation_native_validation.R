#!/usr/bin/env Rscript

## Native Stan target, generated-quantity, and autodiff gate for the joint
## market/s0 rotation.  This uses only tiny Fixed_param contexts; it does not
## run HMC.  Set S0_JOINT_LIB to a library containing the freshly built package.

options(stringsAsFactors = FALSE)
lib <- Sys.getenv("S0_JOINT_LIB", unset = "/tmp/atb-joint-lib")
suppressPackageStartupMessages(library(antitrustBayes, lib.loc = lib))
suppressPackageStartupMessages(library(rstan))

make_fixture <- function() {
  M <- 3L
  n <- c(3L, 4L, 5L)
  market <- rep(seq_len(M), n)
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
  name <- paste0(demand, if (exact) "_stackelberg" else "_supply")
  antitrustBayes:::stanmodels[[name]]
}

prepare_data <- function(demand, exact, coordinate) {
  common <- list(
    prices = fx$prices, quantities = fx$quantities, margins = fx$margins,
    market = fx$market, year = fx$year, group = fx$group,
    supply_model = 1L, use_cutoff = FALSE, market_side = "output",
    use_rho = TRUE, use_market_resid_scale = FALSE,
    use_alpha_market_pool = FALSE, use_s0_market_year_pool = FALSE,
    use_s0_market_pool = TRUE, use_centered_alpha = FALSE,
    use_centered_s0_pool = FALSE,
    market_s0_parameterization = coordinate,
    use_hmt = FALSE, priors = list(), X_s0 = NULL,
    margin_observed = fx$margin_observed, outside_price = 1,
    likelihood_mode = 1L, loo_cluster = fx$loo_cluster,
    estimate_df = FALSE, resid_df_fixed = 5,
    structure = if (exact) "stackelberg" else "simultaneous",
    partition_prior = if (exact) "ordered_k" else "cutoff",
    k_support = if (exact) 1:2 else NULL,
    k_prior = if (exact) c(.4, .6) else NULL,
    k_scope = "market", share_likelihood = "gaussian",
    margin_likelihood = "gaussian",
    quality_parameterization = "collapsed", quality_scope = "market"
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
    seed = 190919L, refresh = 0)
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

stable_hypot <- function(x, y) {
  scale <- max(abs(c(x, y)))
  scale * sqrt((x / scale)^2 + (y / scale)^2)
}

make_state <- function(template, demand, coordinate, sigma_b, sigma_s0) {
  init <- template
  u <- c(-.7, .25, 1.1)
  z <- c(.45, -.8, .2)
  S <- stable_hypot(sigma_b, sigma_s0)
  v <- (sigma_b * u + sigma_s0 * z) / S
  h <- (-sigma_s0 * u + sigma_b * z) / S
  init <- set_one(init, "mu_b_market_vec", .17)
  init <- set_one(init, "sigma_b_market_vec", sigma_b)
  init <- set_one(init, "b_market_raw", if (coordinate == "joint") v else u)
  init <- set_one(init, "mu_s0_pool_vec", -1.12)
  init <- set_one(init, "sigma_s0_pool_vec", sigma_s0)
  init <- set_one(init, "s0_market_raw", z)
  init <- set_one(init, "c_s0_raw", h)
  init <- set_one(init, "sigma_b_strat", .21)
  init <- set_one(init, "log_sigma_share", log(.19))
  init <- set_one(init, "log_sigma_margin", log(.14))
  init <- set_one(init, "rho_raw", atanh(.37 / .95))
  init <- set_one(init, "mu_log_a_free", log(.84))
  init <- set_one(init, "mu_log_a", log(.84))
  init
}

fixed_state <- function(model, data, init) {
  rstan::sampling(model, data = data, chains = 1L, iter = 1L,
    warmup = 0L, algorithm = "Fixed_param", init = list(init),
    seed = 190920L, refresh = 0)
}

finite_gradient <- function(fit, u, rel_step = 1e-5) {
  vapply(seq_along(u), function(j) {
    h <- rel_step * max(1, abs(u[j]))
    up <- um <- u
    up[j] <- up[j] + h
    um[j] <- um[j] - h
    (rstan::log_prob(fit, up, adjust_transform = TRUE) -
       rstan::log_prob(fit, um, adjust_transform = TRUE)) / (2 * h)
  }, numeric(1))
}

rows <- list()
grows <- list()
add <- function(kernel, scenario, quantity, error, tolerance, detail) {
  rows[[length(rows) + 1L]] <<- data.frame(
    scenario = paste(kernel, scenario, sep = ":"), quantity = quantity,
    max_abs_error = error, tolerance = tolerance,
    passed = is.finite(error) && error <= tolerance,
    method = "native_stan_fixed_state", detail = detail
  )
}
addg <- function(kernel, scenario, error, tolerance) {
  grows[[length(grows) + 1L]] <<- data.frame(
    scenario = paste(kernel, scenario, sep = ":"),
    quantity = "stan_autodiff_vs_finite_difference",
    max_abs_error = error, tolerance = tolerance,
    passed = is.finite(error) && error <= tolerance,
    method = "native_stan_autodiff",
    detail = "maximum relative error over every unconstrained coordinate"
  )
}

scenarios <- list(
  ordinary = c(.43, .37),
  s0_near_zero = c(.43, 1e-5),
  b_near_zero = c(1e-5, .37),
  both_small = c(2e-5, 1e-5)
)

for (demand in c("logit", "ces")) for (exact in c(FALSE, TRUE)) {
  kernel <- paste0(demand, if (exact) "_stackelberg" else "_supply")
  model <- model_for(demand, exact)
  data_old <- prepare_data(demand, exact, "legacy")
  data_new <- prepare_data(demand, exact, "joint")
  old_template <- init_template(model, data_old)
  new_template <- init_template(model, data_new)

  for (scenario in names(scenarios)) {
    scales <- scenarios[[scenario]]
    old_init <- make_state(old_template$init, demand, "legacy", scales[1], scales[2])
    new_init <- make_state(new_template$init, demand, "joint", scales[1], scales[2])
    old_u <- rstan::unconstrain_pars(old_template$fit, old_init)
    new_u <- rstan::unconstrain_pars(new_template$fit, new_init)
    old_lp <- as.numeric(rstan::log_prob(old_template$fit, old_u,
                                         adjust_transform = TRUE))
    new_lp <- as.numeric(rstan::log_prob(new_template$fit, new_u,
                                         adjust_transform = TRUE))
    add(kernel, scenario, "target", abs(old_lp - new_lp), 2e-8,
        "legacy versus joint native target at the same economic state")

    old_draw <- fixed_state(model, data_old, old_init)
    new_draw <- fixed_state(model, data_new, new_init)
    available <- intersect(c("b_market", "s0_prob", "share_mean",
      "pred_marginNormInv", "alpha_structural_market",
      "gamma_structural_market", "R0_prob", "c_effective_mky",
      "posterior_component_prob", "prior_component_prob"),
      intersect(old_draw@model_pars, new_draw@model_pars))
    old_gq <- rstan::extract(old_draw, pars = available, permuted = TRUE)
    new_gq <- rstan::extract(new_draw, pars = available, permuted = TRUE)
    for (name in available) {
      add(kernel, scenario, paste0("generated_", name),
          max(abs(as.numeric(old_gq[[name]]) - as.numeric(new_gq[[name]]))),
          2e-8, "legacy versus joint transformed/generated quantity")
    }

    ## Full AD checks cover ordinary and boundary states for the ordinary
    ## kernels, plus one exact-k state for each demand system.
    if ((!exact && scenario %in% c("ordinary", "s0_near_zero", "both_small")) ||
        (exact && identical(scenario, "ordinary"))) {
      ad <- as.numeric(rstan::grad_log_prob(new_template$fit, new_u,
                                            adjust_transform = TRUE))
      fd <- finite_gradient(new_template$fit, new_u)
      addg(kernel, scenario, max(abs(ad - fd) / (1 + abs(ad))), 2e-4)
    }
  }
}

target_rows <- do.call(rbind, rows)
gradient_rows <- do.call(rbind, grows)
old_targets <- read.csv("doc/s0_joint_rotation_target_checks.csv",
                        stringsAsFactors = FALSE)
old_gradients <- read.csv("doc/s0_joint_rotation_gradient_checks.csv",
                          stringsAsFactors = FALSE)
write.csv(rbind(old_targets, target_rows),
          "doc/s0_joint_rotation_target_checks.csv", row.names = FALSE)
write.csv(rbind(old_gradients, gradient_rows),
          "doc/s0_joint_rotation_gradient_checks.csv", row.names = FALSE)

if (any(!target_rows$passed) || any(!gradient_rows$passed)) {
  print(subset(target_rows, !passed))
  print(subset(gradient_rows, !passed))
  stop("native joint-rotation validation failed")
}
cat(sprintf("passed %d native target/GQ and %d native gradient checks\n",
            nrow(target_rows), nrow(gradient_rows)))
