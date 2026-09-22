market_year_size_fixture <- function() {
  list(
    prices = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
    quantities = c(1, 2, 7, 100, 50, 3, 4, 8, 2),
    margins = rep(0.3, 9),
    market = c("A", "A", "A", "B", "B", "A", "A", "B", "B"),
    year = c(2020, 2020, 2020, 2020, 2020, 2021, 2021, 2021, 2021),
    group = letters[seq_len(9)]
  )
}

test_that("market-year totals preserve original units for unequal panel cells", {
  x <- market_year_size_fixture()
  sd <- do.call(prepare_stan_data, c(x, list(
    supply_model = 1L, use_cutoff = FALSE, use_rho = FALSE,
    use_hmt = FALSE, priors = list()
  )))

  expected_totals <- c(10, 7, 150, 10)
  expected_mky <- c(1L, 1L, 1L, 3L, 3L, 2L, 2L, 4L, 4L)
  expect_equal(as.numeric(attr(sd, "N_mt")), expected_totals)
  expect_equal(as.numeric(attr(sd, "inside_quantity_total")), expected_totals)
  expect_equal(as.integer(attr(sd, "market_year_idx")), expected_mky)
  expect_equal(as.numeric(attr(sd, "N_mt_by_observation")),
               expected_totals[expected_mky])
  expect_equal(as.numeric(sd$shareIn),
               x$quantities / expected_totals[expected_mky])
  expect_equal(attr(sd, "market_year_size")$units,
               "original input quantity units")
  expect_match(attr(sd, "market_year_size")$source,
               "original inside quantities")
  expect_false("N_mt" %in% names(sd))
  expect_false("inside_quantity_total" %in% names(sd))
})

test_that("market-year size mapping follows input row order", {
  x <- market_year_size_fixture()
  order <- c(9L, 6L, 3L, 8L, 1L, 7L, 4L, 2L, 5L)
  shuffled <- lapply(x, `[`, order)
  sd <- do.call(prepare_stan_data, c(shuffled, list(
    supply_model = 1L, use_cutoff = FALSE, use_rho = FALSE,
    use_hmt = FALSE, priors = list()
  )))

  expected_row_totals <- with(shuffled, ifelse(market == "A" & year == 2020, 10,
    ifelse(market == "A" & year == 2021, 7,
      ifelse(market == "B" & year == 2020, 150, 10))))
  expect_equal(as.numeric(attr(sd, "N_mt_by_observation")), expected_row_totals)
  expect_equal(as.numeric(sd$shareIn), shuffled$quantities / expected_row_totals)
  expect_equal(as.numeric(sd$shareIn),
               x$quantities[order] / expected_row_totals)
})

test_that("CES retains quantity totals without changing normalized outcomes", {
  x <- market_year_size_fixture()
  sd <- do.call(prepare_ces_stan_data, c(x, list(
    supply_model = 1L, use_cutoff = FALSE, use_rho = FALSE,
    use_hmt = FALSE, priors = list()
  )))
  expected_totals <- c(10, 7, 150, 10)
  mky <- as.integer(attr(sd, "market_year_idx"))
  expect_equal(as.numeric(attr(sd, "N_mt")), expected_totals)
  expect_equal(as.numeric(sd$shareIn), x$quantities / expected_totals[mky])
  expect_equal(as.numeric(sd$revShare),
               x$prices * x$quantities / vapply(seq_len(4), function(m) {
                 sum(x$prices[mky == m] * x$quantities[mky == m])
               }, numeric(1))[mky])
})

test_that("market-year totals reject invalid quantities", {
  x <- market_year_size_fixture()
  x$quantities[4L] <- 0
  expect_error(
    do.call(prepare_stan_data, c(x, list(
      supply_model = 1L, use_cutoff = FALSE, use_rho = FALSE,
      use_hmt = FALSE, priors = list()
    ))),
    "positive and finite"
  )
})

test_that("specified fits expose the one-cell original quantity total", {
  market <- bayes_market(
    prices = c(1, 2), quantities = c(4, 6), margins = c(.3, .3),
    group = c("a", "b"), outside_price = 1, market_side = "input"
  )
  fit <- bayes_specify(
    bayes_params(alpha = 1, s0 = .2, supply_model = 1L), market
  )
  expect_equal(fit@metadata$N_mt, 10)
  expect_equal(fit@metadata$market_year_idx, c(1L, 1L))
  expect_equal(fit@data$N_mt, c(10, 10))
  expect_equal(fit@observed$N_mt, c(10, 10))
})
