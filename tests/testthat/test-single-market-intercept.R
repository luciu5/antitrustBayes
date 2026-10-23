test_that("Logit kernels retain one proper inside-market intercept in a single market", {
  for (kernel in c("logit_supply", "logit_stackelberg")) {
    code <- paste(readLines(.ab_test_path("inst", "stan", paste0(kernel, ".stan")),
                           warn = FALSE), collapse = "\n")
    expect_match(code, "vector[1] mu_b_market_vec;", fixed = TRUE)
    expect_match(code, "market_intercept_location = rep_vector(mu_b_market_vec[1], N_market);",
                 fixed = TRUE)
    expect_false(grepl("market_intercept_location = rep_vector(0.0, N_market)",
                       code, fixed = TRUE))
    expect_length(gregexpr("normal_lpdf(mu_b_market_vec | 0, 1)", code,
                          fixed = TRUE)[[1L]], 1L)
    expect_match(code, "target += normal_lpdf(mu_b_market_vec | 0, 1);\n  if (is_single_market == 0)",
                 fixed = TRUE)
    expect_match(code, "vector<lower=0>[is_single_market == 0 ? 1 : 0] sigma_b_market_vec;",
                 fixed = TRUE)
  }
})

test_that("intercept estimation preserves market flags, quality blocks and outside-share pools", {
  families <- c(MC = 2L, B = 1L, C = 3L, A = 4L,
                BF = 1L, CF = 3L, BS = 1L, CS = 3L)
  for (M in c(1L, 4L)) for (T in c(1L, 4L)) for (family in names(families)) {
    structure <- if (family %in% c("BF", "CF")) "fringe" else
      if (family %in% c("BS", "CS")) "stackelberg" else "simultaneous"
    args <- list(prices = rep(seq(1.2, 2.4, length.out = 7L), M * T),
      quantities = rep(7:1, M * T), margins = rep(.3, 7L * M * T),
      market = rep(paste0("m", seq_len(M)), each = 7L * T),
      year = rep(rep(seq_len(T), each = 7L), M),
      group = rep(paste0("p", 1:7), M * T), supply_model = families[[family]],
      market_side = "output", structure = structure, use_cutoff = FALSE,
      partition_prior = if (structure == "simultaneous") "cutoff" else "ordered_k")
    if (structure != "simultaneous") {
      args$k_support <- 3L
      args$k_prior <- 1
    }
    for (pool in c("none", "market", "market_year")) {
      args$use_s0_market_pool <- identical(pool, "market")
      args$use_s0_market_year_pool <- identical(pool, "market_year")
      sampled <- do.call(prepare_stan_data, args)
      expect_identical(sampled$is_single_market, as.integer(M == 1L))
      expect_identical(sampled$N_market, M)
      expect_identical(sampled$N_market_year, M * T)
      expect_identical(sampled$N_quality_block, M)
      expect_identical(sampled$use_s0_market_pool, as.integer(pool == "market"))
      expect_identical(sampled$use_s0_market_year_pool, as.integer(pool == "market_year"))
      expect_identical(attr(sampled, "model_spec_version"), 13L)
      if (M == 1L) {
        collapsed <- do.call(prepare_stan_data,
          c(args, list(market_intercept_parameterization = "collapsed")))
        expect_identical(lapply(collapsed, identity), lapply(sampled, identity))
        expect_identical(collapsed$market_intercept_parameterization, 1L)
        expect_identical(attr(collapsed, "market_intercept_parameterization_effective"), "sampled")
        expect_identical(attr(collapsed, "market_intercept_coordinate_convention"),
                         "sampled_single_market_coordinates")
      }
    }
  }
})

test_that("seven-firm positive-cost equilibria retain fixed conditional shares at 10, 20 and 40 percent outside", {
  r <- c(.25, .20, .16, .13, .11, .09, .06)
  costs <- c(.70, .82, .96, 1.10, 1.27, 1.45, 1.64)
  for (s0 in c(.10, .20, .40)) for (conduct in c("B", "C", "CF")) {
    arguments <- list(conduct = conduct, J = 7L, mode = "parameter", alpha = 2.4,
      conditional_shares = r, outside_share = s0, costs = costs, seed = 453971L)
    if (conduct == "CF") arguments$k <- 3L
    market <- do.call(bayes_calibrated_market, arguments)
    truth <- market$truth
    p <- truth$data$prices
    b <- unname(truth$meanval - mean(truth$meanval))
    intercept <- mean(truth$meanval) - truth$alpha * mean(p)
    price_std <- (p - mean(p)) / stats::sd(p)
    alpha_std <- truth$alpha * stats::sd(p)
    reconstructed <- stats::qlogis(s0) + intercept + b - alpha_std * price_std
    expect_equal(unname(truth$conditional_shares), r, tolerance = 1e-13)
    expect_equal(truth$s0, s0, tolerance = 1e-13)
    expect_equal(unname(truth$costs), costs, tolerance = 1e-13)
    expect_true(all(truth$costs > 0))
    expect_equal(sum(b), 0, tolerance = 1e-13)
    expect_equal(intercept, mean(log(r)) - stats::qlogis(s0), tolerance = 1e-13)
    expect_equal(reconstructed, log(r), tolerance = 1e-13)
    expect_equal(unname(truth$shares), (1 - s0) * r, tolerance = 1e-13)
    expect_equal(unname(p * (1 - truth$margins)), costs, tolerance = 1e-12)
    expect_true(market$fresh_oracle$passed)
    expect_lt(market$fresh_oracle$max_cost_error, 2e-6)
  }
})
