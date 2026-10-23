.market_sd_data <- function(demand = "logit", model = "B", n_market = 2L,
                            coordinates = "sampled", priors = list()) {
  family <- c(MC = 2L, B = 1L, C = 3L, A = 4L,
              BF = 1L, CF = 3L, BS = 1L, CS = 3L)
  structure <- if (model %in% c("BF", "CF")) "fringe" else
    if (model %in% c("BS", "CS")) "stackelberg" else "simultaneous"
  args <- list(prices = seq(1.1, 2.4, length.out = 7L * n_market),
    quantities = rep(7:1, n_market), margins = rep(.3, 7L * n_market),
    market = rep(paste0("m", seq_len(n_market)), each = 7L),
    year = rep(2001L, 7L * n_market),
    group = paste0("firm", seq_len(7L * n_market)),
    supply_model = unname(family[model]), use_cutoff = FALSE, use_rho = TRUE, use_hmt = FALSE,
    structure = structure, market_intercept_parameterization = coordinates,
    priors = priors, share_basis = "absolute")
  if (structure != "simultaneous") {
    args$partition_prior <- "ordered_k"
    args$k_support <- 1:6
    args$k_prior <- rep(1/6, 6)
  }
  prepare <- if (demand == "logit") prepare_stan_data else prepare_ces_stan_data
  list(data = do.call(prepare, args), arguments = args)
}

test_that("market SD family names and valid hyperparameters resolve canonically", {
  expect_identical(bayes_prior_manifest()$prior_sigma_b_market_family, 2L)
  for (manifest in c("baseline", "tight", "diffuse")) {
    p <- resolve_prior_manifest(manifest, list(prior_sigma_b_market_family = "gamma",
      prior_sigma_b_market_shape = 3, prior_sigma_b_market_rate = 7))
    expect_identical(p$prior_sigma_b_market_family, 2L)
    expect_equal(p$prior_sigma_b_market_shape, 3)
    expect_equal(p$prior_sigma_b_market_rate, 7)
    expect_silent(validate_priors(p))
  }
  expect_identical(resolve_prior_manifest("baseline",
    list(prior_sigma_b_market_family = "half_normal"))$prior_sigma_b_market_family, 1L)
  for (bad in list("unknown", NA_character_, c("gamma", "gamma"), 0, 1.5, Inf, NULL))
    expect_error(resolve_prior_manifest("baseline", list(prior_sigma_b_market_family = bad)),
      "prior_sigma_b_market_family")
  for (name in c("prior_sigma_b_market_scale", "prior_sigma_b_market_shape", "prior_sigma_b_market_rate")) {
    for (bad in list(0, -1, Inf, NA_real_, c(1, 2), NULL)) {
      priors <- setNames(list(bad), name)
      expect_error(.market_sd_data(priors = priors), "Prior")
    }
  }
})

test_that("all conduct and demand kernels receive the selected market SD prior", {
  p <- list(prior_sigma_b_market_family = "gamma",
    prior_sigma_b_market_shape = 2.5, prior_sigma_b_market_rate = 6)
  for (demand in c("logit", "ces")) for (model in c("MC", "B", "C", "A", "BF", "CF", "BS", "CS")) {
    for (coordinates in c("sampled", "collapsed")) {
      d <- .market_sd_data(demand, model, coordinates = coordinates, priors = p)$data
      expect_identical(d$prior_sigma_b_market_family, 2L)
      expect_equal(d$prior_sigma_b_market_shape, 2.5)
      expect_equal(d$prior_sigma_b_market_rate, 6)
    }
    default <- .market_sd_data(demand, model)$data
    hn <- .market_sd_data(demand, model, priors = list(
      prior_sigma_b_market_family = "half_normal", prior_sigma_b_market_scale = .5))$data
    expect_identical(hn$prior_sigma_b_market_family, 1L)
    expect_equal(hn$prior_sigma_b_market_scale, .5)
    expect_identical(hn[setdiff(names(hn), "prior_sigma_b_market_family")],
                     default[setdiff(names(default), "prior_sigma_b_market_family")])
  }
})

test_that("the public estimation route forwards priors before sampling", {
  captured <- NULL
  fields <- c("prior_sigma_b_market_family", "prior_sigma_b_market_scale",
              "prior_sigma_b_market_shape", "prior_sigma_b_market_rate")
  testthat::local_mocked_bindings(.cmdstan_model_for_kernel = function(kernel) list(
    variables = function() list(data = setNames(as.list(rep(1, 4)), fields)),
    sample = function(data, ...) { captured <<- data; stop("captured before HMC") }),
    .package = "antitrustBayes")
  for (demand in c("logit", "ces")) {
    args <- .market_sd_data(demand)$arguments
    args$supply_model <- "bertrand"; args$demand_model <- demand
    args$market_side <- "output"
    args$priors <- list(prior_sigma_b_market_family = "gamma",
      prior_sigma_b_market_shape = 3, prior_sigma_b_market_rate = 4)
    expect_error(do.call(bayes_est, args), "captured before HMC")
    expect_identical(captured$prior_sigma_b_market_family, 2L)
    expect_equal(captured$prior_sigma_b_market_shape, 3)
    expect_equal(captured$prior_sigma_b_market_rate, 4)
  }
})

test_that("evidence identity tracks active market SD choices only", {
  make_fit <- function(n_market, priors) {
    prepared <- .market_sd_data(n_market = n_market, priors = priors)
    a <- prepared$arguments; d <- prepared$data
    new("BayesLogit", data = data.frame(prices = a$prices, quantities = a$quantities,
      margins = a$margins, market = a$market, year = a$year, group = a$group),
      stan_data = d, supply_model = 1L, demand_model = "logit", market_side = "output",
      price_sd = 1, call = quote(identity()), metadata = list(structure = "simultaneous",
        resolved_priors = d[names(bayes_prior_manifest())]))
  }
  active <- function(fit) .mu_active_priors(fit, gaussian_only = TRUE)
  hn <- make_fit(2L, list(prior_sigma_b_market_family = "half_normal"))
  gamma <- make_fit(2L, list(prior_sigma_b_market_family = "gamma"))
  expect_false(identical(active(hn), active(gamma)))
  expect_false(identical(active(gamma), active(make_fit(2L, list(
    prior_sigma_b_market_family = "gamma", prior_sigma_b_market_rate = 6)))))
  expect_identical(active(gamma), active(make_fit(2L, list(
    prior_sigma_b_market_family = "gamma", prior_sigma_b_market_scale = 9))))
  expect_identical(active(hn), active(make_fit(2L, list(
    prior_sigma_b_market_family = "half_normal", prior_sigma_b_market_shape = 9))))
  expect_false(identical(active(hn), active(make_fit(2L, list(
    prior_sigma_b_market_family = "half_normal", prior_sigma_b_market_scale = .7)))))
  stale <- gamma; stale@metadata$resolved_priors$prior_sigma_b_market_rate <- 6
  expect_error(active(stale), "Active prior metadata")
  expect_identical(active(make_fit(1L, list())), active(make_fit(1L,
    list(prior_sigma_b_market_family = "gamma", prior_sigma_b_market_rate = 9))))
})
