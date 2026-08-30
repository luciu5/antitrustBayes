test_that("Bayes decision responds to asymmetric error costs", {
  draws <- c(rep(-1, 6), rep(1, 4))
  equal <- merger_decision(draws)
  enforcement <- merger_decision(draws, false_negative_cost = 2)
  expect_false(equal$challenge)
  expect_true(enforcement$challenge)
  expect_equal(enforcement$probability_cutoff, 1 / 3)
})

test_that("decision regret is zero for the oracle", {
  truth <- c(-1, 1, 2, -2)
  oracle <- truth > 0
  result <- decision_regret(oracle, truth)
  expect_equal(result$regret, rep(0, length(truth)))

  wrong <- decision_regret(!oracle, truth, false_negative_cost = 2)
  expect_gt(sum(wrong$regret), 0)
})
