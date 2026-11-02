## The existence check is deliberately distinct from convergence diagnostics.
.mock_cmdstan_fit <- function(codes = c(0L, 0L), chains = 2L,
                              draws = array(0, c(3L, 2L, 1L)),
                              type = "CmdStanMCMC") {
  structure(list(
    return_codes = function() codes,
    num_chains = function() chains,
    draws = function(...) draws
  ), class = c(type, "CmdStanFit"))
}

test_that("CmdStan fitting rejects incomplete processes and missing draws", {
  expect_error(.assert_stan_sampling_contract(
    .mock_cmdstan_fit(codes = c(0L, 1L)), 2L),
    "did not complete")
  expect_error(.assert_stan_sampling_contract(
    .mock_cmdstan_fit(chains = 1L), 2L),
    "unexpected chain count")
  expect_error(.assert_stan_sampling_contract(
    .mock_cmdstan_fit(draws = array(numeric(), c(0L, 2L, 1L))), 2L),
    "no usable posterior draws")
  expect_error(.assert_stan_sampling_contract(
    .mock_cmdstan_fit(draws = array(Inf, c(2L, 2L, 1L))), 2L),
    "no usable posterior draws")
  expect_error(.assert_stan_sampling_contract(
    .mock_cmdstan_fit(type = "CmdStanPathfinder"), 2L),
    "unexpected fit type")
})

test_that("CmdStan fitting accepts complete HMC and Pathfinder draws", {
  expect_invisible(.assert_stan_sampling_contract(
    .mock_cmdstan_fit(), 2L))
  expect_invisible(.assert_stan_sampling_contract(
    .mock_cmdstan_fit(codes = 0L, chains = 1L,
                      draws = array(0, c(3L, 1L, 1L)),
                      type = "CmdStanPathfinder"),
    1L, inference = "pathfinder"))
})

test_that("failed CmdStan processes retain their return codes and output", {
  fit <- .mock_cmdstan_fit(codes = c(1L, 0L, NA_integer_), chains = 3L)
  queried <- integer()
  fit$output <- function(id) {
    queried <<- c(queried, id)
    cat(if (id == 1L) "Initialization failed; covariance is not positive definite.\n" else
      "Process interrupted.\n")
  }
  error <- tryCatch(.assert_stan_sampling_contract(fit, 3L), error = identity)
  expect_s3_class(error, "error")
  expect_match(conditionMessage(error), "1=1, 3=NA", fixed = TRUE)
  expect_match(conditionMessage(error), "Initialization failed", fixed = TRUE)
  expect_match(conditionMessage(error), "Process interrupted", fixed = TRUE)
  expect_identical(queried, c(1L, 3L))
})

test_that("failure reporting remains usable when process output is unavailable", {
  fit <- .mock_cmdstan_fit(codes = 1L, chains = 1L)
  fit$output <- function(id) stop("output file unavailable")
  expect_error(.assert_stan_sampling_contract(fit, 1L), "1=1", fixed = TRUE)
  expect_error(.assert_stan_sampling_contract(.mock_cmdstan_fit(codes = integer()), 2L),
    "did not complete", fixed = TRUE)
})

test_that("long process output retains the final failure reason", {
  fit <- .mock_cmdstan_fit(codes = 1L, chains = 1L)
  fit$output <- function(id) {
    cat(strrep("a", 3000L), "\nFinal initialization failure.\n", sep = "")
  }
  error <- tryCatch(.assert_stan_sampling_contract(fit, 1L), error = identity)
  expect_match(conditionMessage(error), "Final initialization failure.", fixed = TRUE)
  expect_lt(nchar(conditionMessage(error)), 2200L)
})
