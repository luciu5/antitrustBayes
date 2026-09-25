test_that("one-year exact-k plans retain their market state under a year selector", {
  plan <- list(k_cache = list(`31` = c(market_1 = 2L)))
  selected_year <- list(unit_id = "market_1_1", market_id = "market_1")
  expect_identical(
    antitrustBayes:::.lq_plan_k_for_unit(plan, 1L, selected_year), 2L
  )

  panel_plan <- list(k_cache = list(`31` = c(market_1_1 = 1L,
                                             market_1_2 = 3L)))
  expect_identical(
    antitrustBayes:::.lq_plan_k_for_unit(
      panel_plan, 1L,
      list(unit_id = "market_1_2", market_id = "market_1")
    ), 3L
  )
  expect_error(
    antitrustBayes:::.lq_plan_k_for_unit(
      panel_plan, 1L,
      list(unit_id = "market_2_1", market_id = "market_2")
    ), "no role state"
  )
})
