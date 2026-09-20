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

test_that("market collapsed helpers use margin-first joint Gaussian conditioning", {
  stan_root <- testthat::test_path("..", "..", "inst", "stan")
  models <- c("logit_supply", "ces_supply", "logit_stackelberg",
              "ces_stackelberg")
  for (model in models) {
    code <- paste(readLines(file.path(stan_root, paste0(model, ".stan")),
      warn = FALSE), collapse = "\n")
    expect_match(code,
      "real slope = rho \\* sigma_share / sigma_margin;")
    expect_match(code,
      "real paired_variance = square\\(sigma_share\\) \\* \\(1.0 - square\\(rho\\)\\);")
    expect_match(code,
      "return margin_ll \\+ centered_gaussian_density")
    expect_match(code, "market_obs_position")
  }
})
