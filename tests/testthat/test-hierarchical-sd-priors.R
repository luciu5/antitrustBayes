.hierarchy_prior_data <- function(demand = "logit", model = "B", n_market = 2L,
                                  pooling = "market", priors = list(),
                                  centered = FALSE, alpha_pool = TRUE) {
  grid <- expand.grid(firm = 1:7, year = 1:2, market = seq_len(n_market))
  codes <- c(MC = 2L, B = 1L, C = 3L, A = 4L, BF = 1L, CF = 3L, BS = 1L, CS = 3L)
  structure <- if (model %in% c("BF", "CF")) "fringe" else
    if (model %in% c("BS", "CS")) "stackelberg" else "simultaneous"
  args <- list(prices = seq(1.1, 2.4, length.out = nrow(grid)),
    quantities = 8 - grid$firm, margins = rep(.3, nrow(grid)),
    market = grid$market, year = grid$year,
    group = paste(grid$market, grid$firm, sep = "_"),
    supply_model = unname(codes[model]), structure = structure,
    use_cutoff = FALSE, use_rho = TRUE, use_hmt = FALSE,
    use_s0_market_pool = pooling == "market",
    use_s0_market_year_pool = pooling == "market_year",
    use_alpha_market_pool = alpha_pool, use_centered_alpha = centered,
    use_centered_s0_pool = centered, priors = priors)
  if (structure != "simultaneous") {
    args$partition_prior <- "ordered_k"; args$k_support <- 1:6
    args$k_prior <- rep(1/6, 6)
  }
  prepare <- if (demand == "logit") prepare_stan_data else prepare_ces_stan_data
  list(data = do.call(prepare, args), arguments = args)
}

test_that("structural hierarchy Gamma defaults and named alternatives resolve", {
  for (manifest in c("baseline", "tight", "diffuse")) {
    for (prefix in c("prior_sigma_alpha", "prior_sigma_s0_pool")) {
      defaults <- resolve_prior_manifest(manifest, list())
      expect_identical(defaults[[paste0(prefix, "_family")]], 2L)
      expect_equal(defaults[[paste0(prefix, "_shape")]], 2)
      expect_equal(defaults[[paste0(prefix, "_rate")]], 5)
      for (family in c("half_normal", "gamma")) {
        p <- resolve_prior_manifest(manifest, setNames(list(family), paste0(prefix, "_family")))
        expect_identical(p[[paste0(prefix, "_family")]], if (family == "gamma") 2L else 1L)
        expect_silent(validate_priors(p))
      }
    }
  }
  for (prefix in c("prior_sigma_alpha", "prior_sigma_s0_pool")) {
    for (bad in list("unknown", NA_character_, c("gamma", "gamma"), 0, 1.5, Inf, TRUE, NULL))
      expect_error(resolve_prior_manifest("baseline", setNames(list(bad), paste0(prefix, "_family"))), prefix)
    for (name in c(prefix, paste0(prefix, c("_shape", "_rate")))) {
      for (bad in list(0, -1, Inf, NA_real_, c(1, 2), NULL))
        expect_error(.hierarchy_prior_data(priors = setNames(list(bad), name)), "Prior")
    }
  }
})

test_that("both demand systems and all conduct families forward both hierarchy priors", {
  p <- list(prior_sigma_alpha_family = "gamma", prior_sigma_alpha_shape = 3,
    prior_sigma_alpha_rate = 10, prior_sigma_s0_pool_family = "half_normal",
    prior_sigma_s0_pool = .7)
  for (demand in c("logit", "ces")) for (model in c("MC", "B", "C", "A", "BF", "CF", "BS", "CS")) {
    for (pooling in c("market", "market_year")) for (centered in c(FALSE, TRUE)) {
      d <- .hierarchy_prior_data(demand, model, pooling = pooling, centered = centered, priors = p)$data
      expect_identical(d$prior_sigma_alpha_family, 2L)
      expect_equal(d$prior_sigma_alpha_shape, 3)
      expect_equal(d$prior_sigma_alpha_rate, 10)
      expect_identical(d$prior_sigma_s0_pool_family, 1L)
      expect_equal(d$prior_sigma_s0_pool, .7)
      expect_equal(d$use_s0_market_pool, as.integer(pooling == "market"))
      expect_equal(d$use_s0_market_year_pool, as.integer(pooling == "market_year"))
    }
  }
})

test_that("public estimation forwards the new families before HMC", {
  captured <- NULL
  fields <- unlist(lapply(c("prior_sigma_alpha", "prior_sigma_s0_pool"),
    function(p) c(p, paste0(p, c("_family", "_shape", "_rate")))))
  testthat::local_mocked_bindings(.cmdstan_model_for_kernel = function(kernel) list(
    variables = function() list(data = setNames(as.list(rep(1, length(fields))), fields)),
    sample = function(data, ...) { captured <<- data; stop("captured before HMC") }),
    .package = "antitrustBayes")
  for (demand in c("logit", "ces")) {
    args <- .hierarchy_prior_data(demand)$arguments
    args$supply_model <- "bertrand"; args$demand_model <- demand; args$market_side <- "output"
    args$priors <- list(prior_sigma_alpha_family = "half_normal", prior_sigma_alpha = .25,
      prior_sigma_s0_pool_family = "gamma", prior_sigma_s0_pool_shape = 3,
      prior_sigma_s0_pool_rate = 7)
    expect_error(do.call(bayes_est, args), "captured before HMC")
    expect_identical(captured$prior_sigma_alpha_family, 1L)
    expect_equal(captured$prior_sigma_alpha, .25)
    expect_identical(captured$prior_sigma_s0_pool_family, 2L)
    expect_equal(captured$prior_sigma_s0_pool_shape, 3)
    expect_equal(captured$prior_sigma_s0_pool_rate, 7)
  }
})

test_that("active prior identity includes exactly the selected hierarchy densities", {
  make_fit <- function(priors = list(), n_market = 2L, pooling = "market", alpha_pool = TRUE) {
    prepared <- .hierarchy_prior_data(n_market = n_market, pooling = pooling,
      alpha_pool = alpha_pool, priors = priors)
    a <- prepared$arguments; d <- prepared$data
    new("BayesLogit", data = data.frame(prices = a$prices, quantities = a$quantities,
      margins = a$margins, market = a$market, year = a$year, group = a$group),
      stan_data = d, supply_model = 1L, demand_model = "logit", market_side = "output",
      price_sd = 1, call = quote(identity()), metadata = list(structure = "simultaneous",
        resolved_priors = d[names(bayes_prior_manifest())]))
  }
  active <- function(fit) .mu_active_priors(fit, gaussian_only = TRUE)
  base <- active(make_fit())
  for (prefix in c("prior_sigma_alpha", "prior_sigma_s0_pool")) {
    expect_false(identical(base, active(make_fit(setNames(list(7), paste0(prefix, "_rate"))))))
    expect_false(identical(base, active(make_fit(setNames(list(3), paste0(prefix, "_shape"))))))
    expect_identical(base, active(make_fit(setNames(list(9), prefix))))
    hn <- setNames(list("half_normal"), paste0(prefix, "_family"))
    expect_false(identical(base, active(make_fit(hn))))
    changed <- hn; changed[[paste0(prefix, "_rate")]] <- 9
    expect_identical(active(make_fit(hn)), active(make_fit(changed)))
    changed <- hn; changed[[prefix]] <- .7
    expect_false(identical(active(make_fit(hn)), active(make_fit(changed))))
    stale <- make_fit(); stale@metadata$resolved_priors[[paste0(prefix, "_rate")]] <- 7
    expect_error(active(stale), "Active prior metadata")
  }
  unused <- list(prior_sigma_alpha_family = "half_normal", prior_sigma_alpha = 9,
    prior_sigma_s0_pool_family = "half_normal", prior_sigma_s0_pool = 9)
  expect_identical(active(make_fit(n_market = 1L)), active(make_fit(unused, n_market = 1L)))
  expect_identical(active(make_fit(pooling = "none", alpha_pool = FALSE)),
    active(make_fit(unused, pooling = "none", alpha_pool = FALSE)))
  # One-market panels still have an active market-year outside-share SD.
  year_pool <- make_fit(n_market = 1L, pooling = "market_year")
  fields <- names(active(year_pool))
  expect_true("prior_sigma_s0_pool_family" %in% fields)
  expect_false("prior_sigma_alpha_family" %in% fields)
  expect_false(identical(active(year_pool), active(make_fit(
    list(prior_sigma_s0_pool_family = "half_normal"), n_market = 1L, pooling = "market_year"))))
})
