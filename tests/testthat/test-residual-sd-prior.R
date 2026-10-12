.residual_sd_data <- function(demand = "logit", model = "B", enabled = TRUE,
                              priors = list(), n_market = 2L) {
  grid <- expand.grid(firm = 1:7, year = 1:2, market = seq_len(n_market))
  codes <- c(MC = 2L, B = 1L, C = 3L, A = 4L, BF = 1L, CF = 3L, BS = 1L, CS = 3L)
  structure <- if (model %in% c("BF", "CF")) "fringe" else
    if (model %in% c("BS", "CS")) "stackelberg" else "simultaneous"
  args <- list(prices = seq(1.1, 2.4, length.out = nrow(grid)),
    quantities = 8 - grid$firm, margins = rep(.3, nrow(grid)),
    market = grid$market, year = grid$year, group = paste(grid$market, grid$firm, sep = "_"),
    supply_model = unname(codes[model]), structure = structure,
    use_cutoff = FALSE, use_rho = TRUE, use_hmt = FALSE,
    use_market_resid_scale = enabled, quality_parameterization = "sampled", priors = priors)
  if (structure != "simultaneous") {
    args$partition_prior <- "ordered_k"; args$k_support <- 1:6
    args$k_prior <- rep(1/6, 6)
  }
  prepare <- if (demand == "logit") prepare_stan_data else prepare_ces_stan_data
  list(data = do.call(prepare, args), arguments = args)
}

test_that("all hierarchy SD defaults are Gamma without enabling residual heterogeneity", {
  for (manifest in c("baseline", "tight", "diffuse")) {
    p <- resolve_prior_manifest(manifest, list())
    for (prefix in c("prior_sigma_alpha", "prior_sigma_s0_pool", "prior_sigma_b_market",
                     "prior_sigma_resid_market")) {
      expect_identical(p[[paste0(prefix, "_family")]], 2L)
      expect_equal(p[[paste0(prefix, "_shape")]], 2)
      expect_equal(p[[paste0(prefix, "_rate")]], if (prefix == "prior_sigma_resid_market") 7 else 5)
    }
  }
  expect_false(formals(bayes_est)$use_market_resid_scale)
  expect_equal(.residual_sd_data(enabled = FALSE)$data$use_market_resid_scale, 0L)
  for (bad in list("invalid", NA_character_, c("gamma", "gamma"), 0, 1.5, Inf, TRUE, NULL))
    expect_error(resolve_prior_manifest("baseline", list(prior_sigma_resid_market_family = bad)),
                 "prior_sigma_resid_market_family")
  for (name in c("prior_sigma_resid_market", "prior_sigma_resid_market_shape", "prior_sigma_resid_market_rate"))
    for (bad in list(0, -1, Inf, NA_real_, c(1, 2), NULL))
      expect_error(.residual_sd_data(priors = setNames(list(bad), name)), "Prior")
  invalid <- bayes_prior_manifest(); invalid$prior_sigma_resid_market_family <- 3
  expect_error(validate_priors(invalid), "prior_sigma_resid_market_family")
})

test_that("all conduct and demand families forward residual SD prior choices", {
  for (demand in c("logit", "ces")) for (model in c("MC", "B", "C", "A", "BF", "CF", "BS", "CS")) {
    for (enabled in c(FALSE, TRUE)) for (family in c("half_normal", "gamma")) {
      d <- .residual_sd_data(demand, model, enabled, priors = list(
        prior_sigma_resid_market_family = family, prior_sigma_resid_market = .45,
        prior_sigma_resid_market_shape = 3, prior_sigma_resid_market_rate = 9))$data
      expect_identical(d$prior_sigma_resid_market_family, if (family == "gamma") 2L else 1L)
      expect_equal(d$prior_sigma_resid_market, .45)
      expect_equal(d$prior_sigma_resid_market_shape, 3)
      expect_equal(d$prior_sigma_resid_market_rate, 9)
      expect_equal(d$use_market_resid_scale, as.integer(enabled))
    }
  }
})

test_that("public estimation forwards residual prior families before HMC", {
  fields <- c("use_market_resid_scale", "prior_sigma_resid_market",
    "prior_sigma_resid_market_family", "prior_sigma_resid_market_shape", "prior_sigma_resid_market_rate")
  captured <- NULL
  testthat::local_mocked_bindings(.cmdstan_model_for_kernel = function(kernel) list(
    variables = function() list(data = setNames(as.list(rep(1, length(fields))), fields)),
    sample = function(data, ...) { captured <<- data; stop("captured before HMC") }),
    .package = "antitrustBayes")
  for (demand in c("logit", "ces")) for (family in c("half_normal", "gamma")) {
    args <- .residual_sd_data(demand)$arguments
    args$supply_model <- "bertrand"; args$demand_model <- demand; args$market_side <- "output"
    args$priors <- list(prior_sigma_resid_market_family = family,
      prior_sigma_resid_market = .45, prior_sigma_resid_market_shape = 3, prior_sigma_resid_market_rate = 9)
    expect_error(do.call(bayes_est, args), "captured before HMC")
    expect_identical(captured$prior_sigma_resid_market_family, if (family == "gamma") 2L else 1L)
    expect_equal(captured$prior_sigma_resid_market, .45)
    expect_equal(captured$prior_sigma_resid_market_shape, 3)
    expect_equal(captured$prior_sigma_resid_market_rate, 9)
  }
})

test_that("evidence identity includes exactly the active residual SD prior", {
  make_fit <- function(enabled = TRUE, priors = list()) {
    prepared <- .residual_sd_data(enabled = enabled, priors = priors)
    a <- prepared$arguments; d <- prepared$data
    new("BayesLogit", data = data.frame(prices = a$prices, quantities = a$quantities,
      margins = a$margins, market = a$market, year = a$year, group = a$group),
      stan_data = d, supply_model = 1L, demand_model = "logit", market_side = "output",
      price_sd = 1, call = quote(identity()), metadata = list(structure = "simultaneous",
        resolved_priors = d[names(bayes_prior_manifest())]))
  }
  active <- function(fit) .mu_active_priors(fit, gaussian_only = TRUE)
  base <- active(make_fit())
  expect_false(identical(base, active(make_fit(priors = list(prior_sigma_resid_market_rate = 9)))))
  expect_false(identical(base, active(make_fit(priors = list(prior_sigma_resid_market_shape = 3)))))
  expect_identical(base, active(make_fit(priors = list(prior_sigma_resid_market = 9))))
  hn <- list(prior_sigma_resid_market_family = "half_normal")
  expect_false(identical(base, active(make_fit(priors = hn))))
  expect_identical(active(make_fit(priors = hn)), active(make_fit(priors = c(hn,
    list(prior_sigma_resid_market_rate = 9, prior_sigma_resid_market_shape = 9)))))
  expect_false(identical(active(make_fit(priors = hn)), active(make_fit(priors = c(hn,
    list(prior_sigma_resid_market = .7))))))
  expect_identical(active(make_fit(FALSE)), active(make_fit(FALSE, c(hn,
    list(prior_sigma_resid_market = 9, prior_sigma_resid_market_shape = 9, prior_sigma_resid_market_rate = 9)))))
  stale <- make_fit(); stale@metadata$resolved_priors$prior_sigma_resid_market_rate <- 9
  expect_error(active(stale), "Active prior metadata")
})
