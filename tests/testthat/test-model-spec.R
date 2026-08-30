test_that("Stan likelihoods implement the current model contract", {
  installed_path <- system.file("stan", package = "antitrustBayes")
  root <- testthat::test_path("..", "..")
  stan_dir <- if (nzchar(installed_path)) installed_path else
    file.path(root, "inst", "stan")
  logit_code <- paste(
    readLines(file.path(stan_dir, "logit_supply.stan"), warn = FALSE),
    collapse = "\n"
  )
  ces_code <- paste(
    readLines(file.path(stan_dir, "ces_supply.stan"), warn = FALSE),
    collapse = "\n"
  )

  for (code in list(logit_code, ces_code)) {
    expect_match(code, "clr_shareIn", fixed = TRUE)
    expect_match(code, "clr_normal_lpdf", fixed = TRUE)
    expect_match(code, "student_t_lpdf", fixed = TRUE)
    expect_match(code, "log_lik_period", fixed = TRUE)
    expect_match(code, "log_lik_market", fixed = TRUE)
    expect_match(code, "log_lik_cluster", fixed = TRUE)
    expect_match(code, "bargaining_power", fixed = TRUE)
    expect_false(grepl("prior_lkj", code, fixed = TRUE))
    expect_false(grepl("rho_raw", code, fixed = TRUE))
  }
  expect_match(logit_code, "firm_shareIn", fixed = TRUE)
  expect_match(ces_code, "firm_revShare", fixed = TRUE)
  expect_match(ces_code, "R0_prob", fixed = TRUE)
})
