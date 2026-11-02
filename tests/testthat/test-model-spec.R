test_that("Stan likelihoods implement the current model contract", {
  stan_dir <- .ab_test_path("inst", "stan")
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
    expect_match(code, "student_t_lpdf", fixed = TRUE)
    expect_match(code, "log_lik_period", fixed = TRUE)
    expect_match(code, "log_lik_market", fixed = TRUE)
    expect_match(code, "log_lik_cluster", fixed = TRUE)
    expect_match(code, "bargaining_power", fixed = TRUE)
    expect_match(code, "prior_lkj", fixed = TRUE)
    expect_match(code, "rho_raw", fixed = TRUE)
    expect_match(code, "use_rho", fixed = TRUE)
    expect_match(code, "resid_df_fixed", fixed = TRUE)
    expect_match(code, "estimate_df", fixed = TRUE)
    expect_match(code, "real nu", fixed = TRUE)
    expect_false(grepl("clr_normal_lpdf", code, fixed = TRUE))
  }
  expect_match(logit_code, "firm_shareIn", fixed = TRUE)
  expect_match(ces_code, "firm_revShare", fixed = TRUE)
  expect_match(ces_code, "R0_prob", fixed = TRUE)
})
