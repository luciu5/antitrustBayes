test_that("prior manifests preserve centers and scale dispersion", {
  base <- bayes_prior_manifest("baseline")
  tight <- bayes_prior_manifest("tight")
  diffuse <- bayes_prior_manifest("diffuse")

  expect_equal(base$prior_alpha_mean, 0)
  expect_equal(tight$prior_alpha_mean, base$prior_alpha_mean)
  expect_equal(diffuse$s0_prior_mean, base$s0_prior_mean)
  expect_equal(tight$prior_alpha_sd, 0.5 * base$prior_alpha_sd)
  expect_equal(diffuse$s0_prior_sd, 2 * base$s0_prior_sd)
  expect_equal(base$margin_df, 5)
  expect_equal(base$bargaining_prior_lower, 0)
  expect_equal(base$bargaining_prior_upper, 0.9999)
})

test_that("bank model priors implement institutional and sensitivity sets", {
  institutional <- bank_model_priors()
  expect_equal(
    institutional,
    c(bertrand = 0.30, cournot = 0.30, auction = 0.10,
      monopolistic = 0.30)
  )
  expect_equal(unname(bank_model_priors(scheme = "equal")), rep(0.25, 4))
  expect_equal(
    unname(bank_model_priors(c("bertrand", "cournot", "mc"), "core_three")),
    rep(1 / 3, 3)
  )
  expect_error(
    bank_model_priors(c("bertrand", "bargaining"), "institutional"),
    "does not include bargaining"
  )
})

test_that("invalid bargaining prior bounds are rejected", {
  bad <- bayes_prior_manifest()
  bad$bargaining_prior_lower <- 0.8
  bad$bargaining_prior_upper <- 0.7
  expect_error(validate_priors(bad), "Bargaining prior bounds")
})
