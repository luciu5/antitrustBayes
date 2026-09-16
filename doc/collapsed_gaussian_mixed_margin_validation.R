#!/usr/bin/env Rscript
## Optional D-path validation: same collapsed Gaussian shares, fixed-df t margins.
## Uses installed models, no compilation, HMC, changed priors, or economic edits.
rows <- list()
record <- function(demand, structure, check, error) {
  rows[[length(rows) + 1L]] <<- data.frame(demand, structure, check,
    max_abs_error = error, tolerance = 1e-8,
    status = if (is.finite(error) && error <= 1e-8) "passed" else "failed")
  write.csv(do.call(rbind, rows),
    "doc/collapsed_gaussian_mixed_margin_validation.csv", row.names = FALSE)
}
for (demand in c("logit", "ces")) for (structure in c("simultaneous", "stackelberg")) {
  base <- list(prices = seq(1, 2, length.out = 6), quantities = 1:6,
    margins = c(.3, .3, NA, .3, .3, .3),
    margin_observed = c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE), market = rep("M", 6),
    group = c("A", "A", "B", "C", "D", "D"),
    demand_model = demand, supply_model = "bertrand", structure = structure,
    market_side = "output", df = 5, use_rho = FALSE, use_hmt = FALSE,
    use_market_resid_scale = FALSE, use_alpha_market_pool = FALSE,
    use_s0_market_pool = FALSE, use_s0_market_year_pool = FALSE,
    share_likelihood = "gaussian", quality_parameterization = "collapsed",
    stan_options = list(chains = 1, iter = 1, warmup = 0,
      algorithm = "Fixed_param", init = 0, seed = 16092026, refresh = 0))
  if (structure == "stackelberg") {
    base$k_support <- 2L; base$k_prior <- 1; base$partition_prior <- "ordered_k"
  }
  tryCatch({
    b <- do.call(antitrustBayes::bayes_est, c(base, list(margin_likelihood = "gaussian")))
    d <- do.call(antitrustBayes::bayes_est, c(base, list(margin_likelihood = "student_t")))
    draw <- function(fit, name) as.numeric(rstan::extract(fit@stanfit, name)[[name]])
    ii <- which(b@stan_data$margin_obs == 1)
    y <- b@stan_data$log_marginNormInv[ii]
    mu <- log(draw(b, "pred_marginNormInv")[ii] + 1e-9)
    sigma <- draw(b, "sigma_margin")
    expected <- sum(dt((y - mu) / (sigma * sqrt(3/5)), df = 5, log = TRUE) -
      log(sigma * sqrt(3/5)) - dnorm(y, mu, sigma, log = TRUE))
    target <- function(fit) {
      sf <- fit@stanfit
      u <- rstan::unconstrain_pars(sf, rstan::get_inits(sf)[[1]])
      rstan::log_prob(sf, u, adjust_transform = FALSE)
    }
    record(demand, structure, "normalized_target_difference",
      abs(target(d) - target(b) - expected))
    record(demand, structure, "share_density_unchanged",
      abs(draw(d, "log_lik_share_collapsed") - draw(b, "log_lik_share_collapsed")))
    record(demand, structure, "joint_generated_quantity_difference",
      abs(draw(d, "log_lik_joint_collapsed") - draw(b, "log_lik_joint_collapsed") - expected))
  }, error = function(e) {
    record(demand, structure, paste("error:", conditionMessage(e)), Inf)
  })
}
stopifnot(all(do.call(rbind, rows)$status == "passed"))
cat("Optional mixed-margin validation passed.\n")
