test_that("local conduct slopes match the analytical expansions", {
  expect_equal(
    conduct_local_slopes(0.2),
    c(monopolistic = 0, auction = -0.5, bertrand = -1, cournot = -5)
  )
})

test_that("profiling unknown outside share enlarges equivalence", {
  q <- c(0.03, 0.08, 0.15, 0.25)
  known <- conduct_distance(
    q, "cournot", "bertrand", s0_k = 0.2, profile_s0 = FALSE
  )
  unknown <- conduct_distance(
    q, "cournot", "bertrand", s0_k = 0.2, profile_s0 = TRUE
  )
  expect_lte(unknown$distance, known$distance + 1e-8)
  expect_equal(
    conduct_distance(q, "bertrand", "bertrand", profile_s0 = FALSE)$distance,
    0,
    tolerance = 1e-8
  )
})

test_that("CES bargaining is exactly a Bertrand rescaling", {
  q <- c(0.05, 0.1, 0.2)
  beta <- 0.7
  expect_equal(
    conduct_schedule(q, "bargaining", s0 = 0.3, demand_model = "ces",
                     gamma = 4, bargaining_power = beta),
    conduct_schedule(q, "bertrand", s0 = 0.3, demand_model = "ces",
                     gamma = 4) / (1 - beta)
  )
})

test_that("CES Cournot schedule includes the outside revenue share", {
  q <- c(0.08, 0.16, 0.24)
  s0 <- 0.2
  gamma <- 4
  expected <- 1 / gamma + (gamma - 1) * q /
    (gamma * (1 + (gamma - 1) * s0))
  expect_equal(
    1 / conduct_schedule(q, "cournot", s0 = s0, demand_model = "ces",
                         gamma = gamma),
    expected,
    tolerance = 1e-12
  )
})

test_that("pairwise information matrix is well formed", {
  info <- conduct_information_matrix(
    c(0.03, 0.08, 0.15), models = c("bertrand", "auction", "mc"),
    s0_k = 0.3
  )
  expect_equal(unname(diag(info)), rep(0, 3))
  expect_true(all(info >= 0))
})

test_that("information distance supports CES schedules", {
  d <- conduct_distance(
    c(0.04, 0.10, 0.18), "bertrand", "cournot", s0_k = 0.2,
    demand_model = "ces", gamma = 4
  )
  expect_true(is.finite(d$distance))
  expect_identical(d$demand_model, "ces")
})
