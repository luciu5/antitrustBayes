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

test_that("conduct_schedule rejects jointly infeasible firm shares", {
  ## Each share is individually below 1 - s0 = 0.8, but together they exceed it.
  q <- c(0.5, 0.5)
  for (demand in c("logit", "ces")) {
    expect_error(
      conduct_schedule(q, "bertrand", s0 = 0.2, demand_model = demand),
      "sum\\(q\\) <= 1 - s0"
    )
  }
  expect_error(conduct_schedule(c(0.3, 0.3, 0.3), "cournot", s0 = 0.2),
               "sum\\(q\\) <= 1 - s0")
  expect_error(conduct_distance(q, "bertrand", "cournot", s0_k = 0.2),
               "sum\\(q\\) <= 1 - s0")
  expect_error(
    conduct_distance(c(0.45, 0.40), "bertrand", "cournot", s0_k = 0.2,
                     profile_s0 = FALSE),
    "sum\\(q\\) <= 1 - s0"
  )
})

test_that("conduct_schedule accepts partially enumerated and fully enumerated markets", {
  s0 <- 0.2
  ## Strictly less than the inside share: a subset of firms.
  expect_length(conduct_schedule(c(0.1, 0.2), "bertrand", s0 = s0), 2L)
  ## Exactly the inside share, including floating-point noise.
  full <- c(0.1, 0.2, 0.3, 0.2)
  expect_equal(sum(full), 0.8)
  expect_length(conduct_schedule(full, "cournot", s0 = s0), 4L)
  noisy <- c(0.1, 0.2, 0.3) / 0.6 * (1 - s0)
  expect_length(conduct_schedule(noisy, "auction", s0 = s0), 3L)
  expect_length(
    conduct_schedule(noisy, "bertrand", s0 = s0, demand_model = "ces", gamma = 3),
    3L
  )
  ## A single firm may hold the whole inside share (previously accepted, still is).
  expect_equal(conduct_schedule(0.8, "bertrand", s0 = 0.2), 0.2)
  ## Materially infeasible by more than the tolerance is still rejected.
  expect_error(conduct_schedule(c(0.4, 0.4, 0.0001), "bertrand", s0 = 0.2),
               "sum\\(q\\) <= 1 - s0")
})

test_that("conduct_distance profiles s0 up to 1 - sum(q), not 1 - max(q)", {
  ## Old bound 1 - max(q) = 0.7 would let the comparison model choose
  ## s0 up to 0.7 even though sum(q) = 0.6 requires s0 <= 0.4.
  q <- c(0.30, 0.20, 0.10)
  for (model_l in c("bertrand", "cournot", "auction")) {
    d <- conduct_distance(q, "cournot", model_l, s0_k = 0.3,
                          s0_bounds = c(0.001, 0.95))
    expect_lte(d$s0_l, 1 - sum(q) + 1e-8)
    expect_gt(d$s0_l, 0)
    expect_true(is.finite(d$distance))
  }
  ## Profiling remains valid: it can only reduce the fixed-s0 distance.
  known <- conduct_distance(q, "cournot", "bertrand", s0_k = 0.3,
                            profile_s0 = FALSE)
  profiled <- conduct_distance(q, "cournot", "bertrand", s0_k = 0.3,
                               profile_s0 = TRUE)
  expect_lte(profiled$distance, known$distance + 1e-8)
  expect_equal(known$s0_l, 0.3)
  ## A brute-force grid over the feasible interval agrees with the optimizer.
  grid <- seq(0.001, 1 - sum(q), length.out = 400L)
  target <- conduct_schedule(q, "cournot", 0.3)
  grid_distance <- vapply(grid, function(s0) {
    shape <- conduct_schedule(q, "bertrand", s0)
    alpha <- sum(target * shape) / sum(shape^2)
    0.5 * sum((target - alpha * shape)^2)
  }, numeric(1))
  expect_equal(profiled$distance, min(grid_distance), tolerance = 1e-4)
})

test_that("a partially enumerated market profiles over its whole feasible range", {
  ## Only two of the market's firms are listed; the unlisted ones are absorbed
  ## in 1 - s0 - sum(q) and must not restrict the profiled outside share.
  q <- c(0.05, 0.10)
  wide <- conduct_distance(q, "cournot", "bertrand", s0_k = 0.2)
  expect_lte(wide$s0_l, 1 - sum(q) + 1e-8)
  ## The bound is aggregate: 0.45 + 0.30 leaves at most 0.25 for the outside good,
  ## whereas a max(q) bound would have allowed up to 0.55.
  tight <- conduct_distance(c(0.45, 0.30), "cournot", "bertrand", s0_k = 0.2,
                            s0_bounds = c(0.001, 0.95))
  expect_lte(tight$s0_l, 0.25 + 1e-8)
  ## A comparison interval entirely above 1 - sum(q) is impossible.
  expect_error(
    conduct_distance(c(0.45, 0.30), "cournot", "bertrand", s0_k = 0.2,
                     s0_bounds = c(0.30, 0.95)),
    "incompatible with 'q'"
  )
  ## Fixed s0 ignores the profiling interval.
  fixed <- conduct_distance(c(0.45, 0.30), "cournot", "bertrand", s0_k = 0.2,
                            profile_s0 = FALSE, s0_bounds = c(0.30, 0.95))
  expect_equal(fixed$s0_l, 0.2)
  expect_true(is.finite(fixed$distance))
})

test_that("fully enumerated markets profile without numerical failure", {
  q <- c(0.30, 0.25, 0.15)           # sum 0.70, s0 = 0.30 exactly
  d <- conduct_distance(q, "bertrand", "cournot", s0_k = 0.3)
  expect_equal(d$convergence, 0L)
  expect_lte(d$s0_l, 0.3 + 1e-8)
  expect_equal(conduct_distance(q, "bertrand", "bertrand", s0_k = 0.3)$distance,
               0, tolerance = 1e-8)
})

test_that("shipped conduct validation designs remain feasible", {
  out <- conduct_validation("strong", n_rep = 2, seed = 5)
  expect_equal(nrow(out), 2L)
  expect_true(all(out$min_information >= 0))
})
