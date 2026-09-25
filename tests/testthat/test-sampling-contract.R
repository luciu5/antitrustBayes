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
