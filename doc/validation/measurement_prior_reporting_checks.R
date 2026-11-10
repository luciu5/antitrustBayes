# Explicit extended-tier checks: variance allocation can change economic states
# despite an unchanged observable covariance. No sampling or production edits.
library(testthat)
root <- normalizePath(if (file.exists("DESCRIPTION")) "." else "../..", mustWork = TRUE)
pkgload::load_all(root, quiet = TRUE)
e <- new.env(parent = globalenv())
source(file.path(root, "doc/validation/demand_shock_prototype.R"), local = e)
source(file.path(root, "doc/validation/measurement_prior_prototype.R"), local = e)


test_that("prior upper-tail report is the complement of the lognormal CDF", {
  # A reversed reporting tail would falsely imply most prior mass exceeds .20.
  env <- new.env(parent = e)
  env$ds_private <- function(name) {
    stopifnot(name == ".bayes_backend_extract")
    function(fit, variables) {
      m <- c(.05, .08, .11, .2); omega <- c(.17, .2, .18, .15)
      list(a_market = matrix(1:4), s0_prob = matrix(c(.1, .2, .15, .17)),
        sigma_b_strat_economic = c(.3, .5, .4, .6),
        sigma_share = sqrt(m^2 + omega^2), sigma_margin = c(.2, .3, .25, .4),
        rho = c(.1, .2, .3, .4), sigma_demand_shock = omega,
        rho_demand_shock = c(.2, .3, .4, .5), diagnostic_measurement_sd = m,
        log_sigma_measurement = log(m))
    }
  }
  report <- e$mp_parameters; environment(report) <- env
  output <- tempfile(); dir.create(output); on.exit(unlink(output, recursive = TRUE))
  report(NULL, list(data = list(price_sd = 1)), output,
    list(measurement_prior_log_mean = log(.1), measurement_prior_log_sd = .5))
  p <- read.csv(file.path(output, "measurement_prior_posterior.csv"))
  prior <- p[p$Distribution == "prior", ]
  expect_equal(prior$Pr_Above_20, plnorm(.2, log(.1), .5, lower.tail = FALSE))
  expect_lt(prior$Pr_Above_20, .1)
  expect_equal(prior$Pr_Below_05, prior$Pr_Above_20)
})
