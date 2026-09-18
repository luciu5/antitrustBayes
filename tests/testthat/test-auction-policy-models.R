.auction_policy_fixture <- function() {
  params <- bayes_params(s0 = .2, alpha = rep(1.5, 3), supply_model = 4L)
  market <- bayes_market(prices = c(5, 6, 7), quantities = c(40, 30, 20),
    margins = c(.3, .3, .3), group = c("A", "A", "B"), outside_price = 1,
    market_side = "output", insideSize = 90, labels = c("one", "two", "three"))
  bayes_specify(params, market)
}

test_that("auction realization preserves the public cost-based baseline", {
  fit <- .auction_policy_fixture()
  before <- serialize(fit, NULL)
  plan <- posterior_policy_draws(fit, n_draws = 2L, seed = 12L)
  observed <- suppressWarnings(posterior_policy_models(fit, draw_plan = plan))
  latent <- suppressWarnings(posterior_policy_models(fit, draw_plan = plan,
    demand_anchoring = "latent_quality"))
  expect_equal(nrow(observed$failures), 0L)
  expect_equal(nrow(latent$failures), 0L)
  expect_identical(observed$draw_plan, latent$draw_plan)
  for (position in seq_along(plan$draw_indices)) {
    one <- observed$models$unit_1[[position]]
    two <- latent$models$unit_1[[position]]
    expect_s4_class(one@model, "Auction2ndLogit")
    expect_equal(unname(one@model@mcPre), unname(two@model@mcPre))
    expect_equal(unname(one@model@pricePre), unname(two@model@pricePre), tolerance = 1e-9)
    expect_equal(antitrust::calcShares(one@model), antitrust::calcShares(two@model), tolerance = 1e-9)
    expect_identical(attr(two, "bayes_latent_quality")$demand_representation,
                     "auction_cost_based_value_index")
    expect_identical(attr(two, "bayes_policy_draw")$draw_index, plan$draw_indices[[position]])
  }
  expect_identical(serialize(fit, NULL), before)
})

test_that("latent auction value and costs remain fixed under cost and ownership policies", {
  fit <- .auction_policy_fixture()
  plan <- posterior_policy_draws(fit, n_draws = 1L, seed = 10L)
  prepared <- suppressWarnings(posterior_policy_models(fit, draw_plan = plan))
  source <- prepared$models$unit_1[[1L]]
  original_cost <- as.numeric(source@model@mcPre)
  score <- log(antitrust::calcShares(source@model) / .2) + c(.15, .15, -.15)
  result <- .lq_apply_to_fit(source,
    list(score = score, s0 = .2, quality = c(.15, -.15), draw_index = 1L), plan)
  model <- result$fit@model
  normalized <- exp(score-max(score)); normalized <- .8 * normalized/sum(normalized)
  expect_equal(unname(model@slopes$meanval), unname(log(normalized/.2)), tolerance = 1e-12)
  expect_equal(as.numeric(model@mcPre), original_cost, tolerance = 0)
  shares <- antitrust::calcShares(model)
  firm_share <- as.numeric(model@ownerPre %*% shares)
  markup <- -log1p(-firm_share)/(abs(model@slopes$alpha)*firm_share)
  expect_equal(unname(model@pricePre), original_cost + markup, tolerance = 1e-10)
  value <- model@slopes$meanval - model@slopes$alpha * original_cost
  expect_equal(result$metadata$auction_value, unname(value), tolerance = 1e-12)
  changed <- suppressWarnings(antitrust::simulate(result$fit,
    ownerPost = c("merged", "merged", "merged"),
    mcDelta = original_cost*c(.25, .25, 0)))
  expect_equal(as.numeric(changed@mcPost), original_cost*c(1.25, 1.25, 1), tolerance = 1e-12)
  index <- changed@slopes$meanval + changed@slopes$alpha*changed@mcDelta
  expect_equal(index - changed@slopes$alpha*changed@mcPost, value, tolerance = 1e-12)
  expected <- exp(index); expected <- expected/(1+sum(expected))
  expect_equal(unname(antitrust::calcShares(changed, preMerger = FALSE)), unname(expected), tolerance = 1e-12)
  expect_equal(unname(changed@pricePost - changed@mcPost),
    unname(antitrust::calcMargins(changed, preMerger = FALSE, level = TRUE)), tolerance = 1e-12)
})

test_that("auction policy realization does not broaden unsupported demand or side", {
  fit <- .auction_policy_fixture()
  fit@market_side <- "input"
  expect_error(posterior_policy_models(fit), "output Logit/CES")
  params <- bayes_params(s0 = .2, gamma = rep(5, 3), supply_model = 4L)
  market <- bayes_market(prices = c(5, 6, 7), quantities = c(40, 30, 20),
    margins = rep(.3, 3), group = c("A", "A", "B"), outside_price = 1)
  expect_error(posterior_policy_models(bayes_specify(params, market)), "output Logit/CES")
})

test_that("collapsed joint quality is realized once in an auction policy plan", {
  data <- data.frame(prices = c(5, 6, 7), quantities = c(40, 30, 20),
    margins = rep(.3, 3), group = c("A", "A", "B"),
    market = "market_1", year = "1", outside_price = 1,
    margin_observed = TRUE, stringsAsFactors = FALSE)
  incidence <- .build_quality_incidence(rep(1L, 3), c(1L, 1L, 2L),
    "market_1", c("A", "B"), "market", margin_idx = 1:3)
  stan_data <- c(list(N = 3L, N_group = 2L, group = c(1L, 1L, 2L),
    N_market = 1L, N_market_year = 1L, market = rep(1L, 3),
    market_year_idx = rep(1L, 3), margin_idx = 1:3,
    log_shareIn = log(data$quantities/sum(data$quantities)),
    log_marginNormInv = log(1/data$margins),
    quality_scope = 2L, quality_parameterization = 2L,
    share_likelihood = 2L, margin_likelihood = 2L,
    likelihood_mode = 1L, use_rho = 1L, use_cutoff = 0L,
    use_market_resid_scale = 0L),
    incidence[setdiff(names(incidence), c("quality_scope", "quality_scope_code", "quality_lookup"))])
  attr(stan_data, "quality_incidence") <- incidence
  attr(stan_data, "quality_scope") <- "market"
  attr(stan_data, "group_levels") <- c("A", "B")
  attr(stan_data, "market_levels") <- "market_1"
  attr(stan_data, "model_spec_version") <- 10L
  fit <- methods::new("BayesLogit", stanfit = NULL, data = data,
    stan_data = stan_data, price_sd = 1, supply_model = 4L,
    demand_model = "logit", market_side = "output", call = quote(fixed_draw_fixture()),
    metadata = list(model_spec_version = 10L, quality_scope = "market",
      quality_parameterization = "collapsed", share_likelihood = "gaussian",
      margin_likelihood = "gaussian", structure = "simultaneous",
      partition_prior = "cutoff", posterior_draws = 2L))
  firm_s <- c(rep(.8*70/90, 2), .8*20/90)
  markup <- -log1p(-firm_s)/(1.5*firm_s)
  draws <- list(a_market = matrix(rep(1.5, 2), ncol = 1L),
    s0_prob = matrix(rep(.2, 2), ncol = 1L),
    share_mean = matrix(rep(qlogis(.2)-1.5*c(-1, 0, 1), 2), nrow = 2, byrow = TRUE),
    sigma_b_strat = rep(.8, 2), sigma_share = rep(.2, 2),
    sigma_margin = rep(.3, 2), rho = rep(.25, 2),
    pred_marginNormInv = matrix(rep(data$prices/markup, 2), nrow = 2, byrow = TRUE))
  testthat::local_mocked_bindings(.posterior_extract = function(...) draws,
    .package = "antitrustBayes")
  before <- serialize(fit, NULL)
  plan <- posterior_policy_draws(fit, n_draws = 2L, seed = 47L)
  expect_equal(unname(rowSums(plan$quality_draws)), c(0, 0), tolerance = 1e-12)
  expect_identical(plan$quality_conditioning, "observed log shares and matched Gaussian margins")
  latent <- suppressWarnings(posterior_policy_models(fit, draw_plan = plan,
    demand_anchoring = "latent_quality", labels = c("one", "two", "three"), insideSize = 90))
  observed <- suppressWarnings(posterior_policy_models(fit, draw_plan = plan,
    labels = c("one", "two", "three"), insideSize = 90))
  expect_equal(nrow(latent$failures), 0L)
  expect_equal(nrow(observed$failures), 0L)
  for (position in seq_along(plan$draw_indices)) {
    draw <- plan$draw_indices[[position]]
    model <- latent$models[[1L]][[position]]@model
    expected <- draws$share_mean[draw, ] + plan$quality_draws[position, c(1, 1, 2)]
    expected <- exp(expected-max(expected)); expected <- .8*expected/sum(expected)
    expect_equal(unname(model@slopes$meanval), unname(log(expected/.2)), tolerance = 1e-12)
    expect_equal(as.numeric(model@mcPre), as.numeric(observed$models[[1L]][[position]]@model@mcPre))
    expect_equal(attr(latent$models[[1L]][[position]], "bayes_policy_draw")$quality,
      as.numeric(plan$quality_draws[position, ]))
  }
  expect_identical(serialize(fit, NULL), before)
})
