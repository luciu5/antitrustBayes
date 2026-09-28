rotation_fixture <- function(multi_market = FALSE) {
  if (isTRUE(multi_market)) {
    market <- c("M1", "M1", "M2", "M2")
    year <- rep("Y1", 4L)
  } else {
    market <- rep("M1", 4L)
    year <- rep(c("Y1", "Y2"), each = 2L)
  }
  list(
    prices = c(1.0, 1.1, 0.9, 1.05),
    quantities = c(1.0, 0.8, 1.2, 0.9),
    margins = rep(0.25, 4L),
    market = market,
    year = year,
    group = c("G1", "G2", "G1", "G2")
  )
}

test_that("market s0 rotation resolver distinguishes requested and active modes", {
  resolve <- antitrustBayes:::.resolve_market_s0_parameterization

  legacy <- resolve("legacy", N_market = 2L)
  rotated <- resolve("rotated", N_market = 2L)
  one_market <- resolve("rotated", N_market = 1L)
  joint <- resolve("joint", N_market = 2L,
                   use_s0_market_pool = TRUE,
                   use_centered_s0_pool = FALSE)

  expect_equal(legacy$code, 1L)
  expect_equal(legacy$requested, "legacy")
  expect_equal(legacy$effective, "legacy")
  expect_false(legacy$active)

  expect_equal(rotated$code, 2L)
  expect_equal(rotated$requested, "rotated")
  expect_equal(rotated$effective, "rotated")
  expect_true(rotated$active)
  expect_equal(rotated$coordinate_convention,
               "rotated_market_outside_share")

  expect_equal(one_market$code, 2L)
  expect_equal(one_market$requested, "rotated")
  expect_equal(one_market$effective, "legacy")
  expect_false(one_market$active)
  expect_equal(one_market$coordinate_convention,
               "legacy_market_outside_share")
  expect_equal(joint$code, 3L)
  expect_equal(joint$effective, "joint")
  expect_true(joint$active)
  expect_true(joint$joint_active)
  expect_equal(joint$coordinate_convention,
               "joint_orthogonal_market_outside_share")
  expect_error(resolve("joint", N_market = 1L,
                       use_s0_market_pool = TRUE,
                       use_centered_s0_pool = FALSE),
               "requires more than one market")
  expect_error(resolve("joint", N_market = 2L,
                       use_s0_market_pool = FALSE,
                       use_centered_s0_pool = FALSE),
               "requires market-level")
  expect_error(resolve("joint", N_market = 2L,
                       use_s0_market_pool = TRUE,
                       use_centered_s0_pool = TRUE),
               "requires a noncentered")
  expect_error(resolve("invalid", N_market = 2L), "should be one of")
})

test_that("helpers record rotation mode and keep one-market NC pooling a no-op", {
  d <- rotation_fixture()
  args <- list(
    prices = d$prices, quantities = d$quantities, margins = d$margins,
    market = d$market, year = d$year, group = d$group,
    supply_model = 1L, use_rho = FALSE, use_s0_market_pool = TRUE,
    use_centered_s0_pool = FALSE, market_s0_parameterization = "rotated"
  )

  logit <- do.call(prepare_stan_data, args)
  ces <- do.call(prepare_ces_stan_data, c(args, list(
    use_cutoff = FALSE, use_hmt = FALSE, priors = list(),
    market_side = "output"
  )))

  for (sd in list(logit, ces)) {
    expect_equal(sd$market_s0_parameterization, 2L)
    expect_equal(sd$use_s0_market_pool, 1L)
    expect_equal(sd$use_s0_market_year_pool, 0L)
    expect_equal(attr(sd, "market_s0_parameterization_requested"), "rotated")
    expect_equal(attr(sd, "market_s0_parameterization_effective"), "legacy")
    expect_false(attr(sd, "market_s0_rotation_active"))
    expect_equal(attr(sd, "market_s0_coordinate_convention"),
                 "legacy_market_outside_share")
  }
})

test_that("multi-market helpers activate the rotated contract for both demand models", {
  d <- rotation_fixture(multi_market = TRUE)
  args <- list(
    prices = d$prices, quantities = d$quantities, margins = d$margins,
    market = d$market, year = d$year, group = d$group,
    supply_model = 1L, use_rho = FALSE, use_s0_market_pool = TRUE,
    use_centered_s0_pool = FALSE, market_s0_parameterization = "rotated"
  )

  logit <- do.call(prepare_stan_data, args)
  ces <- do.call(prepare_ces_stan_data, c(args, list(
    use_cutoff = FALSE, use_hmt = FALSE, priors = list(),
    market_side = "output"
  )))
  for (sd in list(logit, ces)) {
    expect_equal(sd$market_s0_parameterization, 2L)
    expect_equal(attr(sd, "market_s0_parameterization_requested"), "rotated")
    expect_equal(attr(sd, "market_s0_parameterization_effective"), "rotated")
    expect_true(attr(sd, "market_s0_rotation_active"))
    expect_equal(sd$N_market, 2L)
  }
})

test_that("multi-market helpers expose the narrow joint orthogonal contract", {
  d <- rotation_fixture(multi_market = TRUE)
  args <- list(
    prices = d$prices, quantities = d$quantities, margins = d$margins,
    market = d$market, year = d$year, group = d$group,
    supply_model = 1L, use_rho = FALSE, use_s0_market_pool = TRUE,
    use_centered_s0_pool = FALSE, market_s0_parameterization = "joint"
  )

  logit <- do.call(prepare_stan_data, args)
  ces <- do.call(prepare_ces_stan_data, c(args, list(
    use_cutoff = FALSE, use_hmt = FALSE, priors = list(),
    market_side = "output"
  )))
  for (sd in list(logit, ces)) {
    expect_equal(sd$market_s0_parameterization, 3L)
    expect_equal(attr(sd, "market_s0_parameterization_requested"), "joint")
    expect_equal(attr(sd, "market_s0_parameterization_effective"), "joint")
    expect_true(attr(sd, "market_s0_rotation_active"))
    expect_true(attr(sd, "market_s0_joint_active"))
    expect_equal(attr(sd, "market_s0_coordinate_convention"),
                 "joint_orthogonal_market_outside_share")
  }

  expect_error(do.call(prepare_stan_data, within(args, {
    use_s0_market_pool <- FALSE
  })), "requires market-level")
})

test_that("all four Stan kernels retain rotation diagnostics and economic equations", {
  stan_dir <- .ab_test_path("inst", "stan")
  models <- c("logit_supply", "logit_stackelberg",
              "ces_supply", "ces_stackelberg")

  for (nm in models) {
    code <- paste(readLines(file.path(stan_dir, paste0(nm, ".stan")),
                            warn = FALSE), collapse = "\n")
    flat <- gsub("[[:space:]]+", " ", code)
    expect_match(flat,
                 "int<lower=1, upper=3> market_s0_parameterization;",
                 fixed = TRUE)
    expect_match(flat, "market_s0_parameterization == 2 && is_single_market == 0",
                 fixed = TRUE)
    expect_match(flat, "vector[N_market_year] c_unbounded_mky",
                 fixed = TRUE)
    expect_match(flat, "vector[N_market_year] c_effective_mky",
                 fixed = TRUE)
    expect_match(flat, "s0_prob[mky] = bounded_inv_logit(s0[mky])",
                 fixed = TRUE)
    expect_match(flat,
                 "c_unbounded_mky[mky] = s0[mky] + market_intercept_location[m]",
                 fixed = TRUE)
    expect_match(flat,
                 "c_effective_mky[mky] = logit(s0_prob[mky]) + market_intercept_location[m]",
                 fixed = TRUE)
    expect_match(flat,
                 "if (is_single_market == 0 && use_s0_market_pool == 1 && use_centered_s0_pool == 0)",
                 fixed = TRUE)
    expect_match(flat,
                 "real w = c_s0_raw[m]",
                 fixed = TRUE)
    expect_match(flat,
                 "real w = c_s0_raw[mky]",
                 fixed = TRUE)
    expect_match(flat,
                 "z_s0 = c_s0_raw - b_market / sigma_s0_pool_vec[1]",
                 fixed = TRUE)
    expect_match(flat,
                 "market_s0_parameterization == 3 && is_single_market == 0",
                 fixed = TRUE)
    expect_match(flat,
                 "market_s0_joint_scale = hypot(sigma_b_market_vec[1], sigma_s0_pool_vec[1])",
                 fixed = TRUE)
    expect_match(flat,
                 "sigma_b_market_vec[1] * b_market_raw - sigma_s0_pool_vec[1] * c_s0_raw",
                 fixed = TRUE)
    expect_match(flat,
                 "sigma_s0_pool_vec[1] * b_market_raw[m] + sigma_b_market_vec[1] * c_s0_raw[m]",
                 fixed = TRUE)
    expect_match(flat, "target += std_normal_lpdf(c_s0_raw)", fixed = TRUE)
    if (grepl("stackelberg", nm, fixed = TRUE)) {
      expect_match(flat, "stack_", fixed = TRUE)
      expect_match(flat, "conduct == 1", fixed = TRUE)
      expect_match(flat, "structure == 1", fixed = TRUE)
    } else {
      expect_match(flat, "supply_model == 3", fixed = TRUE)
    }
    expect_match(flat,
                 "share_mean[i] = logit(s0_prob[mky]) + market_intercept_location[m]",
                 fixed = TRUE)
  }

  logit_code <- paste(readLines(file.path(stan_dir, "logit_supply.stan")),
                      collapse = "\n")
  expect_match(logit_code, "log_kappa_vec", fixed = TRUE)
  expect_match(logit_code,
               "mu_log_a = log_kappa_vec[1] - log_s0_ref_global",
               fixed = TRUE)
  expect_no_match(logit_code, "log_kappa_market", fixed = TRUE)
  expect_match(logit_code,
               "log_a_market_derived = mu_log_a + sigma_log_a * r_market_a_raw",
               fixed = TRUE)

  stackelberg_code <- paste(
    readLines(file.path(stan_dir, "logit_stackelberg.stan")), collapse = "\n"
  )
  expect_match(stackelberg_code, "log_kappa_vec", fixed = TRUE)
  expect_no_match(stackelberg_code, "log_kappa_market", fixed = TRUE)
  expect_match(stackelberg_code,
               "log_a_market_derived = mu_log_a + sigma_log_a * r_market_a_raw",
               fixed = TRUE)

  ces_code <- paste(readLines(file.path(stan_dir, "ces_supply.stan")),
                    collapse = "\n")
  expect_match(ces_code, "ces_inv_margin", fixed = TRUE)
  expect_match(ces_code, "gamma_structural_market", fixed = TRUE)
})
