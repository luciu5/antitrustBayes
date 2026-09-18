## Focused API contract checks for the bounded defaults promotion.  These
## tests prepare data only; they do not compile or sample a Stan model.

test_that("new entry points expose the Gaussian collapsed market defaults", {
  expect_identical(formals(bayes_est)$use_rho, TRUE)
  expect_identical(eval(formals(bayes_est)$share_likelihood),
                   c("gaussian", "student_t"))
  expect_identical(eval(formals(bayes_est)$quality_parameterization),
                   c("collapsed", "sampled"))
  expect_identical(eval(formals(bayes_est)$quality_scope),
                   c("market", "global"))

  expect_identical(eval(formals(prepare_stan_data)$share_likelihood),
                   c("gaussian", "student_t"))
  expect_identical(eval(formals(prepare_stan_data)$quality_parameterization),
                   c("collapsed", "sampled"))
  expect_identical(eval(formals(prepare_stan_data)$quality_scope),
                   c("market", "global"))

  expect_identical(eval(formals(prepare_ces_stan_data)$share_likelihood),
                   c("gaussian", "student_t"))
  expect_identical(eval(formals(prepare_ces_stan_data)$quality_parameterization),
                   c("collapsed", "sampled"))
  expect_identical(eval(formals(prepare_ces_stan_data)$quality_scope),
                   c("market", "global"))
})

test_that("omitted likelihood settings resolve to the new contract", {
  modes <- .resolve_likelihood_modes()
  expect_identical(modes[c("share_likelihood", "margin_likelihood",
                           "quality_parameterization", "quality_scope")],
                   list(share_likelihood = "gaussian",
                        margin_likelihood = "gaussian",
                        quality_parameterization = "collapsed",
                        quality_scope = "market"))
  expect_false(modes$df_active)
  expect_identical(modes$mode, "gaussian/collapsed/gaussian")

  f <- list(
    prices = c(1.0, 1.2, 0.9, 1.1),
    quantities = c(1.1, 0.9, 1.2, 0.8),
    margins = c(0.20, 0.30, 0.25, 0.15),
    market = c(1, 1, 2, 2), year = rep(1, 4), group = c(1, 2, 1, 2),
    supply_model = 1L, market_side = "input"
  )
  d <- do.call(prepare_stan_data, f)
  expect_identical(d$use_rho, 1L)
  expect_identical(d$share_likelihood, 2L)
  expect_identical(d$margin_likelihood, 2L)
  expect_identical(d$quality_parameterization, 2L)
  expect_identical(d$quality_scope, 2L)
  expect_identical(attr(d, "quality_scope"), "market")
  expect_identical(attr(d, "likelihood_mode_metadata")$mode,
                   "gaussian/collapsed/gaussian")
  expect_identical(d$estimate_df, 0L)
  expect_null(attr(d, "likelihood_mode_metadata")$resid_df_fixed)
})

test_that("explicit legacy likelihood and global scope remain available", {
  f <- list(
    prices = c(1.0, 1.2, 0.9, 1.1),
    quantities = c(1.1, 0.9, 1.2, 0.8),
    margins = c(0.20, 0.30, 0.25, 0.15),
    market = c(1, 1, 2, 2), year = rep(1, 4), group = c(1, 2, 1, 2),
    supply_model = 1L, market_side = "input", use_rho = FALSE,
    share_likelihood = "student_t", quality_parameterization = "sampled",
    quality_scope = "global"
  )
  d <- do.call(prepare_stan_data, f)
  expect_identical(d$share_likelihood, 1L)
  expect_identical(d$margin_likelihood, 1L)
  expect_identical(d$quality_parameterization, 1L)
  expect_identical(d$quality_scope, 1L)
  expect_identical(attr(d, "quality_scope"), "global")
  expect_identical(attr(d, "likelihood_mode_metadata")$mode,
                   "student_t/sampled/student_t")

  ## Historical fits without scope metadata continue to decode as global.
  stan_data <- list(group = c(1L, 2L), N_group = 2L)
  fit <- methods::new(
    "BayesLogit", stanfit = NULL,
    data = data.frame(group = c("a", "b")), stan_data = stan_data,
    price_sd = 1, supply_model = 1L, demand_model = "logit",
    market_side = "input", call = quote(fake()), metadata = list()
  )
  expect_identical(.quality_scope_from_fit(fit), "global")
})
