test_that("fast conduct validation returns reproducible theory-linked output", {
  a <- conduct_validation("weak", n_rep = 3, seed = 99)
  b <- conduct_validation("weak", n_rep = 3, seed = 99)
  expect_equal(a, b)
  expect_equal(nrow(a), 3)
  expect_true(all(a$true_model_probability >= 0 & a$true_model_probability <= 1))
  expect_true(all(a$min_information >= 0))
})
