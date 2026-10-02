test_that("default Stan cores follow mc.cores and resolved chain count", {
  had_option <- "mc.cores" %in% names(options())
  old_option <- getOption("mc.cores")
  on.exit(if (had_option) options(mc.cores = old_option)
          else options(mc.cores = NULL), add = TRUE)

  options(mc.cores = NULL)
  expect_identical(.resolve_stan_cores(list(), chains = 4L), 1L)

  options(mc.cores = 4L)
  expect_identical(.resolve_stan_cores(list(), chains = 4L), 4L)
  expect_identical(.resolve_stan_cores(list(), chains = 2L), 2L)
  expect_identical(.resolve_stan_cores(list(), chains = 1L), 1L)
})

test_that("explicit Stan cores override the mc.cores default", {
  had_option <- "mc.cores" %in% names(options())
  old_option <- getOption("mc.cores")
  on.exit(if (had_option) options(mc.cores = old_option)
          else options(mc.cores = NULL), add = TRUE)

  options(mc.cores = 4L)
  expect_identical(.resolve_stan_cores(list(cores = 1L), chains = 4L), 1L)
  expect_identical(.resolve_stan_cores(list(cores = 2L), chains = 1L), 2L)
  expect_identical(.resolve_stan_cores(list(cores = NULL), chains = 4L), 4L)
  expect_identical(getOption("mc.cores"), 4L)
  options(mc.cores = "malformed")
  expect_identical(.resolve_stan_cores(list(cores = 3L), chains = 1L), 3L)
  expect_identical(getOption("mc.cores"), "malformed")
})

test_that("invalid core and chain settings fail clearly", {
  old_option <- options("mc.cores")
  on.exit(options(old_option), add = TRUE)

  bad <- list(0, -1, 1.5, NA_real_, Inf, "4", TRUE)
  for (value in bad) {
    expect_error(
      .resolve_stan_cores(list(), chains = value),
      "stan_options\\$chains.*positive finite integer"
    )
    expect_error(
      .resolve_stan_cores(list(cores = value), chains = 4L),
      "stan_options\\$cores.*positive finite integer"
    )
  }
  for (value in bad) {
    options(mc.cores = value)
    expect_error(
      .resolve_stan_cores(list(), chains = 4L),
      "getOption\\('mc\\.cores'\\).*positive finite integer"
    )
  }
})

test_that("bayes_est passes the resolved chain policy to sampling", {
  old_option <- options("mc.cores")
  on.exit(options(old_option), add = TRUE)
  options(mc.cores = 4L)
  captured <- NULL

  testthat::with_mocked_bindings(
    .cmdstan_model_for_kernel = function(kernel) {
      list(
        variables = function() list(data = NULL),
        sample = function(data, ...) {
          captured <<- list(...)
          stop("sampling interception", call. = FALSE)
        }
      )
    },
    .package = "antitrustBayes",
    {
      estimate <- function(stan_options) {
        bayes_est(
          prices = c(1, 1.1, 1.2), quantities = c(4, 3, 2),
          margins = rep(.25, 3L), market = rep("M1", 3L),
          group = c("A", "B", "C"), market_side = "output",
          use_rho = FALSE, stan_options = stan_options
        )
      }
      expect_error(estimate(list(chains = 2L)), "sampling interception")
      expect_identical(captured$chains, 2L)
      expect_identical(captured$parallel_chains, 2L)

      expect_error(estimate(list(chains = 2L, cores = 1L)),
                   "sampling interception")
      expect_identical(captured$parallel_chains, 1L)

      expect_error(estimate(list(chains = NULL)), "sampling interception")
      expect_identical(captured$chains, 4L)
      expect_identical(captured$parallel_chains, 4L)
    }
  )
})
