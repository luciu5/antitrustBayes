test_that("all production Stan programs expose the scoped quality contract", {
  stan_root <- testthat::test_path("..", "..", "inst", "stan")
  models <- c("logit_supply", "ces_supply", "logit_stackelberg",
              "ces_stackelberg")
  for (model in models) {
    code <- paste(readLines(file.path(stan_root, paste0(model, ".stan")),
      warn = FALSE), collapse = "\n")
    expect_match(code, "int<lower=1, upper=2> quality_scope")
    expect_match(code, "N_quality_coord")
    expect_match(code, "obs_to_quality")
    expect_match(code, "quality_block_start")
    expect_match(code, "market_margin_order")
    expect_match(code, "collapsed_joint_from_local_share_stats")
  }
})

test_that("scoped quality Stan sources parse", {
  testthat::skip_if_not_installed("rstan")
  stan_root <- testthat::test_path("..", "..", "inst", "stan")
  models <- c("logit_supply", "ces_supply", "logit_stackelberg",
              "ces_stackelberg")
  for (model in models) {
    expect_silent(rstan::stanc(
      file = file.path(stan_root, paste0(model, ".stan")),
      model_name = model,
      allow_undefined = TRUE,
      obfuscate_model_name = FALSE
    ))
  }
})

test_that("market collapsed helpers use the local observation count", {
  stan_root <- testthat::test_path("..", "..", "inst", "stan")
  models <- c("logit_supply", "ces_supply", "logit_stackelberg",
              "ces_stackelberg")
  for (model in models) {
    code <- paste(readLines(file.path(stan_root, paste0(model, ".stan")),
      warn = FALSE), collapse = "\n")
    expect_match(code,
      "return centered_gaussian_from_stats\\(within, a, c, logv,\\s*share_n \\+ margin_n")
  }
})
