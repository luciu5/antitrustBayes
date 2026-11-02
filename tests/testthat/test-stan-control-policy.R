# Intercept the public estimator immediately before CmdStan sampling. No HMC.
.capture_sampling_controls <- function(options = list(), demand = "logit",
                                       supply = "cournot") {
  captured <- NULL
  testthat::local_mocked_bindings(
    .cmdstan_model_for_kernel = function(kernel) list(
      variables = function() list(data = NULL),
      sample = function(data, ...) {
        captured <<- list(...)
        stop("sampling controls captured", call. = FALSE)
      }
    ), .package = "antitrustBayes")
  grid <- expand.grid(firm = 1:3, market = 1:4)
  expect_error(bayes_est(
    prices = seq(1, 2, length.out = nrow(grid)), quantities = 4 - grid$firm,
    margins = rep(.25, nrow(grid)), market = grid$market,
    group = paste(grid$market, grid$firm, sep = "_"), market_side = "output",
    demand_model = demand, supply_model = supply, use_alpha_market_pool = TRUE,
    use_rho = FALSE, stan_options = options
  ), "sampling controls captured", fixed = TRUE)
  captured
}

test_that("a metric override preserves package adaptation and depth defaults", {
  for (demand in c("logit", "ces")) for (supply in c("bertrand", "cournot")) {
    baseline <- .capture_sampling_controls(demand = demand, supply = supply)
    dense <- .capture_sampling_controls(list(control = list(metric = "dense_e")),
      demand = demand, supply = supply)
    expect_identical(dense$adapt_delta, .95)
    expect_identical(dense$max_treedepth, 13L)
    expect_identical(dense$metric, "dense_e")
    dense$metric <- NULL
    expect_identical(dense, baseline)
  }
})

test_that("partial adaptation and depth overrides preserve other defaults", {
  adapted <- .capture_sampling_controls(list(control = list(adapt_delta = .999)))
  expect_identical(adapted$adapt_delta, .999)
  expect_identical(adapted$max_treedepth, 13L)
  depth <- .capture_sampling_controls(list(control = list(max_treedepth = 16L)))
  expect_identical(depth$adapt_delta, .95)
  expect_identical(depth$max_treedepth, 16L)
  empty <- .capture_sampling_controls(list(control = list()))
  expect_identical(empty, .capture_sampling_controls())
})

test_that("complete controls and top-level aliases retain their explicit values", {
  explicit <- .capture_sampling_controls(list(chains = 2L, iter = 4000L,
    warmup = 2000L, cores = 2L, seed = 123L, init = .1, control = list(
      adapt_delta = .999, max_treedepth = 16L, metric = "dense_e",
      stepsize = .02, adapt_engaged = FALSE)))
  expect_identical(explicit$chains, 2L)
  expect_identical(explicit$parallel_chains, 2L)
  expect_identical(explicit$iter_warmup, 2000L)
  expect_identical(explicit$iter_sampling, 2000L)
  expect_identical(explicit$seed, 123L)
  expect_identical(explicit$init, .1)
  expect_identical(explicit$adapt_delta, .999)
  expect_identical(explicit$max_treedepth, 16L)
  expect_identical(explicit$metric, "dense_e")
  expect_identical(explicit$step_size, .02)
  expect_false(explicit$adapt_engaged)
  alias <- .capture_sampling_controls(list(adapt_delta = .99, max_treedepth = 15L,
    control = list(metric = "dense_e", adapt_delta = .97)))
  expect_identical(alias$adapt_delta, .99)
  expect_identical(alias$max_treedepth, 15L)
  expect_identical(alias$metric, "dense_e")
})

test_that("explicit NULL controls continue to delegate to CmdStan", {
  one <- .capture_sampling_controls(list(control = list(metric = "dense_e",
    adapt_delta = NULL)))
  expect_false("adapt_delta" %in% names(one))
  expect_identical(one$max_treedepth, 13L)
  expect_identical(one$metric, "dense_e")
  all <- .capture_sampling_controls(list(control = NULL))
  expect_false(any(c("adapt_delta", "max_treedepth", "metric") %in% names(all)))
})
