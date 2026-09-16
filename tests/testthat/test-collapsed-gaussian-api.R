test_that("collapsed Gaussian API preserves defaults and records resolved modes", {
  f <- list(
    prices = c(1.0, 1.2, 0.9, 1.1),
    quantities = c(1.1, 0.9, 1.2, 0.8),
    margins = c(0.20, 0.30, 0.25, 0.15),
    market = c(1, 1, 2, 2), year = rep(1, 4), group = c(1, 2, 1, 2),
    supply_model = 1L, market_side = "input"
  )

  d_default <- do.call(prepare_stan_data, f)
  expect_equal(d_default$share_likelihood, 1L)
  expect_equal(d_default$margin_likelihood, 1L)
  expect_equal(d_default$quality_parameterization, 1L)

  d_a <- do.call(prepare_stan_data, c(f, list(
    share_likelihood = "gaussian", quality_parameterization = "sampled",
    margin_likelihood = "gaussian"
  )))
  expect_equal(unname(unlist(d_a[c("share_likelihood", "margin_likelihood",
                                  "quality_parameterization")])), c(2L, 2L, 1L))

  d_b <- do.call(prepare_stan_data, c(f, list(
    share_likelihood = "gaussian", quality_parameterization = "collapsed",
    margin_likelihood = "gaussian", use_rho = FALSE
  )))
  expect_equal(unname(unlist(d_b[c("share_likelihood", "margin_likelihood",
                                  "quality_parameterization")])), c(2L, 2L, 2L))
  expect_equal(attr(d_b, "model_spec_version"), 9L)
  expect_equal(attr(d_b, "likelihood_mode_metadata")$mode,
               "gaussian/collapsed/gaussian")

  incidence <- f
  incidence$group <- c(1, 1, 2, 3)
  incidence$share_likelihood <- "gaussian"
  incidence$quality_parameterization <- "collapsed"
  incidence$margin_likelihood <- "gaussian"
  incidence$use_rho <- FALSE
  d_incidence <- do.call(prepare_stan_data, incidence)
  expect_equal(d_incidence$N_group, 3L)
  expect_equal(tabulate(d_incidence$group, nbins = d_incidence$N_group),
               c(2L, 1L, 1L))

  d_d <- do.call(prepare_stan_data, c(f, list(
    share_likelihood = "gaussian", quality_parameterization = "collapsed",
    margin_likelihood = "student_t", use_rho = FALSE
  )))
  expect_equal(d_d$margin_likelihood, 1L)

  expect_error(do.call(prepare_stan_data, c(f, list(
    share_likelihood = "student_t", margin_likelihood = "gaussian"
  ))), "Student-t shares")
  expect_error(do.call(prepare_stan_data, c(f, list(
    share_likelihood = "student_t", quality_parameterization = "collapsed"
  ))), "collapsed")
  expect_error(do.call(prepare_stan_data, c(f, list(
    share_likelihood = "gaussian", quality_parameterization = "sampled",
    margin_likelihood = "student_t", use_rho = TRUE
  ))), "collapsed quality")
  expect_error(do.call(prepare_stan_data, c(f, list(
    share_likelihood = "gaussian", quality_parameterization = "collapsed",
    use_market_resid_scale = TRUE
  ))), "use_market_resid_scale")
  expect_error(do.call(prepare_stan_data, c(f, list(
    share_likelihood = "gaussian", quality_parameterization = "sampled",
    margin_likelihood = "gaussian", estimate_df = TRUE
  ))), "only available")
})

test_that("missing margins and all four direct model paths encode the same contract", {
  common <- list(
    prices = c(1.0, 1.2, 0.9, 1.1), quantities = c(1.1, 0.9, 1.2, 0.8),
    margins = c(0.20, 0.30, 0.25, 0.15), market = c(1, 1, 2, 2),
    year = rep(1, 4), group = c(1, 2, 1, 2), supply_model = 1L,
    use_cutoff = FALSE, use_hmt = FALSE, priors = list(),
    market_side = "input", share_likelihood = "gaussian",
    quality_parameterization = "collapsed", margin_likelihood = "gaussian",
    use_rho = FALSE, margin_observed = c(TRUE, FALSE, TRUE, FALSE)
  )
  d <- do.call(prepare_stan_data, common)
  expect_equal(d$N_margin, 2L)
  expect_equal(as.integer(d$margin_idx), c(1L, 3L))

  ces <- common
  ces$market_side <- "output"
  ces$margins <- c(0.20, 0.30, 0.25, 0.15)
  ces$priors <- list()
  d_ces <- do.call(prepare_ces_stan_data, ces)
  expect_equal(d_ces$N_margin, 2L)
  expect_equal(d_ces$quality_parameterization, 2L)
  expect_equal(attr(d_ces, "model_spec_version"), 9L)
})

test_that("exact ordered-k rho restriction distinguishes rho zero from active rho", {
  f <- list(
    prices = c(1.0, 1.2, 0.9, 1.1), quantities = c(1.1, 0.9, 1.2, 0.8),
    margins = c(0.20, 0.30, 0.25, 0.15), market = c(1, 1, 2, 2),
    year = rep(1, 4), group = c(1, 2, 1, 2), supply_model = 1L,
    market_side = "input", structure = "fringe", partition_prior = "ordered_k",
    k_support = 0:1, share_likelihood = "gaussian",
    quality_parameterization = "collapsed", margin_likelihood = "gaussian"
  )
  expect_no_error(do.call(prepare_stan_data, c(f, list(use_rho = FALSE))))
  expect_error(do.call(prepare_stan_data, c(f, list(use_rho = TRUE))),
               "supported only for one market")

  f$partition_prior <- "cutoff"
  f$k_support <- NULL
  expect_no_error(do.call(prepare_stan_data, c(f, list(use_rho = TRUE))))
  f$structure <- "simultaneous"
  expect_no_error(do.call(prepare_stan_data, c(f, list(use_rho = TRUE))))
})

test_that("fringe and Stackelberg Logit/CES preparations retain BF/CF/BS/CS data", {
  f <- list(
    prices = c(1.0, 1.2, 0.9, 1.1), quantities = c(1.1, 0.9, 1.2, 0.8),
    margins = c(0.20, 0.30, 0.25, 0.15), market = rep(1, 4), year = rep(1, 4),
    group = c(1, 2, 1, 2), supply_model = 1L, use_cutoff = FALSE,
    use_hmt = FALSE,
    structure = "fringe",
    partition_prior = "ordered_k", k_support = 0:1, use_rho = FALSE,
    share_likelihood = "gaussian", quality_parameterization = "collapsed",
    margin_likelihood = "gaussian", priors = list()
  )
  logit_fringe <- do.call(prepare_stan_data, c(f, list(market_side = "input")))
  expect_true(all(c("N_firm_mky", "N_partition_component", "firm_mky_rank") %in%
                  names(logit_fringe)))
  f$structure <- "stackelberg"
  f$k_support <- 1:2
  logit_stack <- do.call(prepare_stan_data, c(f, list(market_side = "input")))
  expect_true(all(c("N_firm_market", "firm_market_rank") %in%
                  names(logit_stack)))

  f$structure <- "fringe"
  f$k_support <- 0:1
  ces_fringe <- do.call(prepare_ces_stan_data, c(f, list(market_side = "output")))
  expect_true(all(c("N_firm_mky", "firm_mky_revShare") %in%
                  names(ces_fringe)))
  f$structure <- "stackelberg"
  f$k_support <- 1:2
  ces_stack <- do.call(prepare_ces_stan_data, c(f, list(market_side = "output")))
  expect_true(all(c("N_firm_market", "firm_market_revShare") %in%
                  names(ces_stack)))
})
