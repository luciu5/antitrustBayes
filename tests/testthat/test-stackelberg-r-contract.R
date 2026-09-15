test_that("Stackelberg preparation preserves unique firm-MY and partition maps", {
  prices <- c(2, 2.1, 3, 3.1, 2.2, 2.3, 3.2, 3.3)
  quantities <- c(10, 8, 12, 9, 11, 7, 13, 10)
  margins <- rep(c(.30, .28, .27, .25), 2)
  market <- rep(c("m1", "m2", "m1", "m2"), 2)
  year <- rep(c("y1", "y1", "y1", "y1", "y2", "y2", "y2", "y2"), 1)
  group <- c("A", "C", "B", "D", "A", "C", "B", "D")
  ## Product rows are interleaved by market while each cell has two firms.
  d <- prepare_stan_data(
    prices, quantities, margins, market, year, group,
    supply_model = 1L, structure = "stackelberg",
    k_support = c(1L, 2L), k_prior = c(.6, .4),
    use_rho = FALSE, use_hmt = FALSE
  )
  expect_equal(d$N_firm_mky, 8L)
  expect_equal(d$mky_to_partition_unit, array(c(1L, 1L, 2L, 2L)))
  expect_equal(d$N_partition_unit, 2L)
  expect_equal(d$N_partition_component, 4L)
  expect_equal(d$partition_component_k, array(c(1L, 2L, 1L, 2L)))
  common_fields <- c("shareIn", "firm_shareIn", "marginNormInv",
                     "log_shareIn", "clr_shareIn", "log_marginNormInv",
                     "price", "price_std")
  expect_true(all(vapply(d[common_fields], function(x) length(dim(x)) == 1L,
                         logical(1))))
  new_fields <- c(
    "firm_mky_shareIn", "firm_mky_revShare", "firm_mky_market",
    "firm_mky_market_year", "firm_mky_group", "firm_mky_rank",
    "firm_mky_start", "firm_mky_n", "obs_to_firm_mky",
    "partition_unit_market", "partition_unit_market_year",
    "partition_unit_start", "partition_unit_count", "partition_unit_k_min",
    "partition_unit_k_max", "mky_to_partition_unit",
    "partition_component_unit", "partition_component_k",
    "partition_prior_logprob", "partition_prior_probability"
  )
  expect_true(all(vapply(d[new_fields], function(x) length(dim(x)) == 1L,
                         logical(1))))

  expect_error(
    prepare_stan_data(
      prices, quantities, margins, market, year, group,
      supply_model = 1L, structure = "stackelberg",
      k_support = c(1L, 2L), k_prior = c(.6, .4),
      k_scope = "period", use_rho = FALSE, use_hmt = FALSE
    ),
    "one partition per market"
  )

  expect_error(
    prepare_stan_data(
      prices, quantities, margins, market, year, group,
      supply_model = 1L, structure = "stackelberg",
      k_support = c(1L, 2L), cutoff_scale = .037,
      use_rho = FALSE, use_hmt = FALSE
    ),
    "not defined"
  )
})

test_that("fringe cutoff scales remain prespecified and are carried into data", {
  d <- prepare_stan_data(
    prices = c(1.5, 1.8), quantities = c(10, 8), margins = c(.3, .25),
    market = c("m", "m"), year = c("y", "y"), group = c("a", "b"),
    supply_model = 1L, structure = "fringe", partition_prior = "cutoff",
    cutoff_scale = .37,
    use_rho = FALSE, use_hmt = FALSE
  )
  expect_equal(d$cutoff_scale, .37)
  expect_equal(attr(d, "cutoff_scale"), .37)
  expect_equal(attr(d, "model_spec_version"), 8L)
})

test_that("singleton CES common vectors retain one-dimensional Stan shapes", {
  d <- prepare_ces_stan_data(
    prices = 2, quantities = 1, margins = .3,
    market = "m", year = "y", group = factor("z", levels = c("z", "a")),
    supply_model = 1L, structure = "stackelberg", k_support = 1L,
    market_side = "output",
    use_cutoff = FALSE, use_rho = FALSE, use_hmt = FALSE, priors = list()
  )
  common_fields <- c("shareIn", "revShare", "firm_revShare", "marginNormInv",
                     "log_shareIn", "clr_shareIn", "log_marginNormInv",
                     "price", "log_price_std")
  expect_true(all(vapply(d[common_fields], function(x) length(dim(x)) == 1L,
                         logical(1))))
})

test_that("Stack summaries distinguish cutoff mixtures and direct ordered-k", {
  market <- bayes_market(
    prices = c(2, 2, 2), quantities = c(1, 1, 1),
    margins = c(.3, .3, .3), group = c("a", "b", "c"),
    outside_price = 1, market_side = "input", insideSize = 10,
    labels = c("p1", "p2", "p3")
  )
  cutoff_fit <- bayes_specify(
    bayes_params(alpha = c(.7, .7), s0 = .2,
                 cutoff_share = c(.1, .5), cutoff_scale = .03,
                 supply_model = 1L, structure = "fringe",
                 partition_prior = "cutoff"), market
  )
  cutoff <- bayes_cutoff_summary(cutoff_fit)$by_market_year
  expect_equal(length(cutoff$firm_ids[[1L]]), 3L)
  expect_null(cutoff$firm_leader_probability[[1L]])
  expect_true(length(cutoff$firm_strategic_probability[[1L]]) == 3L)
  expect_equal(cutoff$firm_core_probability[[1L]],
               cutoff$firm_strategic_probability[[1L]])
  expect_true(cutoff$raw_k_pmf[[1L]][["0"]] > 0)
  expect_true(cutoff$hard_k_pmf[[1L]][["3"]] > 0)
  expect_true(cutoff$core_share_conditional[[1L]] < 1)
  expect_identical(cutoff_fit@metadata$prior_normalization$mode,
                   "supplied_draws")
  expect_false("external_offset" %in% names(cutoff_fit@metadata$prior_normalization))

  ordered_fit <- bayes_specify(
    bayes_params(alpha = .7, s0 = .2, supply_model = 1L,
                 structure = "stackelberg", k_support = 1:3,
                 k_prior = c(.2, .3, .5)), market
  )
  ordered <- bayes_cutoff_summary(ordered_fit)
  expect_null(ordered$cutoff_posterior)
  expect_true(all(is.na(ordered$by_market_year$cutoff_mean)))
  expect_equal(sum(ordered$by_market_year$realized_k_posterior[[1L]]), 1,
               tolerance = 1e-12)

  unsorted <- bayes_specify(
    bayes_params(alpha = .7, s0 = .2, supply_model = 1L,
                 structure = "fringe", k_support = c(2L, 0L, 1L),
                 k_prior = c(.4, .3, .3)), market
  )
  uq <- bayes_cutoff_summary(unsorted)$by_market_year
  expect_equal(uq$realized_k_median, 1)
  expect_equal(uq$realized_k_lower, 0)
  expect_equal(uq$realized_k_upper, 2)
})

test_that("component role probabilities include firms absent from a period", {
  got <- .stack_component_role_probability(
    k = c(1L, 2L), probability = c(.25, .75), n_firm = 3L
  )
  expect_equal(got, c(1, .75, 0))
})

test_that("fringe smooth core summaries retain cutoff and outside-share dependence", {
  market <- bayes_market(
    prices = c(2, 2.2, 2.4), quantities = c(8, 2, 1),
    margins = c(.3, .3, .3), group = c("a", "b", "c"),
    outside_price = 1, market_side = "input", insideSize = 10,
    labels = c("p1", "p2", "p3")
  )
  fit <- bayes_specify(
    bayes_params(alpha = c(.7, .7), s0 = c(.1, .8), cutoff_share = c(.2, .9),
                 cutoff_scale = c(.1, .2), supply_model = 1L,
                 structure = "fringe", partition_prior = "cutoff"), market
  )
  got <- bayes_cutoff_summary(fit)$by_market_year
  expect_null(got$firm_leader_probability[[1L]])
  expect_true(length(got$firm_strategic_probability[[1L]]) == 3L)
  expect_equal(got$firm_core_probability[[1L]],
               got$firm_strategic_probability[[1L]])
  sizes <- c(8, 2, 1) / 11
  tau <- c(.1, .2); cuts <- c(.2, .9); s0 <- c(.1, .8)
  smooth <- vapply(seq_along(cuts), function(i) {
    sum(sizes * plogis((sizes - cuts[[i]]) / tau[[i]]))
  }, numeric(1))
  expect_equal(got$smooth_core_share_conditional[[1L]], mean(smooth), tolerance = 1e-12)
  expect_equal(got$smooth_core_share_unconditional[[1L]],
               mean(smooth * (1 - s0)), tolerance = 1e-12)
})

test_that("role offsets reject selected tied boundaries only", {
  component <- list(k = 1:3, probability = rep(1 / 3, 3), cutoff = .4)
  expect_error(
    .stack_bridge_select_k(component, list(partition_prior = "cutoff"),
                           c(.5, .5, .2), "threshold", NULL, -1L),
    "tied firm-size boundary"
  )
  ## The cutoff equals a size, but the selected offset boundary is distinct.
  expect_equal(
    .stack_bridge_select_k(component, list(partition_prior = "cutoff"),
                           c(.6, .4, .2), "threshold", NULL, 1L),
    2L
  )
  direct_tie <- list(k = 1:3, probability = c(.5, .3, .2), cutoff = NULL)
  expect_equal(
    .stack_bridge_select_k(direct_tie,
                           list(partition_prior = "ordered_k"),
                           c(.5, .5, .2), "map", NULL, 1L),
    2L
  )
  direct_component <- list(k = 1:3, probability = c(.2, .5, .3), cutoff = NULL)
  expect_error(
    .stack_bridge_select_k(direct_component,
                           list(partition_prior = "ordered_k"),
                           c(.4, .3, .2, .1), "map", 4L, 0L),
    "supplied support"
  )
  expect_error(
    .stack_bridge_select_k(direct_component,
                           list(partition_prior = "ordered_k"),
                           c(.4, .3, .2, .1), "map", NULL, 2L),
    "supplied support"
  )
})

test_that("bridge and canonical preparation agree on tied firm ranks", {
  ## Factor levels deliberately put B before A; the declared IDs must still
  ## resolve the exact tie as A then B in both paths.
  g <- factor(c("A", "B", "A", "B"), levels = c("B", "A"))
  prices <- c(2, 2.1, 3, 3.1)
  quantities <- c(1, 1, 2, 2)
  d <- prepare_stan_data(
    prices, quantities, c(.3, .3, .3, .3),
    market = c("m", "m", "m", "m"), year = rep("y", 4), group = g,
    supply_model = 1L, structure = "stackelberg", k_support = 1:2,
    k_prior = c(.5, .5), use_rho = FALSE,
    use_hmt = FALSE
  )
  expected <- .stack_bridge_sizes(data.frame(
    prices = prices, quantities = quantities, group = as.character(g)
  ), "logit")
  expect_equal(names(expected), c("A", "B"))
  expect_equal(as.numeric(expected), as.numeric(d$firm_mky_shareIn), tolerance = 1e-15)
  expect_equal(as.integer(d$firm_mky_rank), c(1L, 2L))
})

test_that("ordered-k reconstruction preserves the stored market ranking", {
  df <- data.frame(
    prices = c(2, 2, 2, 2),
    quantities = c(1, 7, 6, 2),
    group = c("a", "b", "c", "d")
  )
  fit <- methods::new(
    "BayesLogit", data = data.frame(), supply_model = 1L,
    demand_model = "logit", market_side = "output", call = quote(fake()),
    metadata = list(
      partition_prior = "ordered_k",
      firm_market_ranking = list(
        firm_market_market = rep(1L, 5),
        firm_market_group_label = c("a", "b", "c", "d", "e"),
        firm_market_rank = 1:5
      )
    ),
    stanfit = list(), stan_data = list(), price_sd = 1
  )
  unit <- list(market_index = 1L)
  got <- .stack_bridge_ranked_sizes(fit, unit, df, "logit")
  expect_identical(names(got), letters[1:5])
  expect_equal(unname(got), c(1, 7, 6, 2, 0) / 16)
  ## Period-specific size order is b, c, d, a, whereas the fixed market
  ## ordering remains a, b, c, d, e.
  expect_identical(names(.stack_bridge_sizes(df, "logit")),
                   c("b", "c", "d", "a"))
})

test_that("estimated Stack bridge fails closed when component GQ is absent", {
  skip_if_not_installed("rstan")
  fake <- methods::new(
    "BayesLogit", data = data.frame(), supply_model = 1L,
    demand_model = "logit", market_side = "input", call = quote(fake()),
    metadata = list(structure = "stackelberg"), stanfit = list(),
    stan_data = list(
      mky_to_partition_unit = array(1L, dim = 1L),
      partition_unit_start = array(1L, dim = 1L),
      partition_unit_count = array(1L, dim = 1L),
      partition_component_k = array(1L, dim = 1L)
    ), price_sd = 1
  )
  unit <- list(market_year_index = 1L, market = "m1")
  expect_error(
    .stack_bridge_component_prob(
      fake, unit, list(k_scope = "market", partition_prior = "ordered_k"),
      sizes = 1, draw = 1L,
      draws = list(a_market = 0, posterior_component_prob = NULL)
    ),
    "missing generated quantity"
  )
})

test_that("shared-market cache uses full component support for selected years", {
  df <- data.frame(
    prices = c(2, 2.1, 2.2, 2.3, 2, 2.1, 2.2),
    quantities = c(4, 3, 2, 1, 4, 3, 2), margins = rep(.3, 7),
    market = rep("m", 7), year = c(rep("y1", 4), rep("y2", 3)),
    group = c("a", "b", "c", "d", "a", "b", "c"),
    outside_price = rep(1, 7)
  )
  stan_data <- list(
    mky_to_partition_unit = array(c(1L, 1L), dim = 2L),
    partition_unit_start = array(1L, dim = 1L),
    partition_unit_count = array(3L, dim = 1L),
    partition_component_k = array(1:3, dim = 3L), use_cutoff = 0L
  )
  attr(stan_data, "model_spec_version") <- 8L
  fit <- methods::new(
    "BayesLogit", data = df, supply_model = 1L, demand_model = "logit",
    market_side = "input", call = quote(fake()),
    metadata = list(
      structure = "stackelberg", partition_prior = "ordered_k",
      k_scope = "market", k_prior = NULL, cutoff_scale = NULL,
      firm_market_ranking = list(
        firm_market_market = rep(1L, 4),
        firm_market_group_label = c("a", "b", "c", "d"),
        firm_market_rank = 1:4
      )
    ),
    stanfit = list(), stan_data = stan_data, price_sd = 1
  )
  setup <- list(
    data = df,
    units = data.frame(unit_id = "m_y1", market = "m", year = "y1"),
    info = list(list(unit_id = "m_y1", market_id = "m", year_id = "y1",
                     idx = 1:4, market_index = 1L, market_year_index = 1L))
  )
  draws <- list(
    a_market = matrix(.7, nrow = 1L, ncol = 1L),
    s0_prob = matrix(.2, nrow = 1L, ncol = 2L),
    posterior_component_prob = matrix(c(.1, .2, .7), nrow = 1L)
  )
  expect_error(
    .stack_k_cache(fit, draws, setup, 1L, "map", fixed_k = 4L,
                   role_offset = 0L),
    "shared market support"
  )
})

test_that("specified Stackelberg simulation uses the native one-step route", {
  skip_if_not_installed("coordination")
  market <- bayes_market(
    prices = c(2, 2.2, 2.4, 1.8), quantities = c(10, 8, 6, 9),
    margins = c(.3, .28, .25, .32), group = c("a", "b", "c", "d"),
    outside_price = 1, market_side = "input", insideSize = 100,
    labels = paste0("p", 1:4)
  )
  fit <- bayes_specify(
    bayes_params(alpha = c(.7, .8), s0 = .2, supply_model = 1L,
                 structure = "stackelberg", k_support = 1:3,
                 k_prior = rep(1 / 3, 3)), market
  )
  native <- bayes_stackelberg_model(fit, draw = 1L, fixed_k = 1L,
                                    role_mapping = "map")
  expect_identical(native@diagnostics$bayes_reconstruction$normalization$mode,
                   "supplied_draws")
  expect_equal(native@diagnostics$bayes_reconstruction$normalization$model_spec_version,
               8L)
  out <- bayes_sim(
    fit, ownerPost = antitrust::counterfactual(
      costs = stats::setNames(.02, "p1")
    ), n_draws = 1L, mc.cores = 1L
  )
  expect_identical(out@metadata$simulation_route,
                   "coordination::stackelberg_simulate")
  expect_equal(out@n_valid[[1L]], 1L)
})

test_that("exact fringe cutoff summary reports ordered strategic-core states", {
  market <- bayes_market(
    prices = c(2, 2.2, 2.4), quantities = c(10, 8, 6),
    margins = c(.3, .28, .25), group = c("a", "b", "c"),
    outside_price = 1, market_side = "output", insideSize = 100
  )
  fit <- bayes_specify(
    bayes_params(
      alpha = 1, s0 = .2, supply_model = 1L, structure = "fringe",
      k_support = 0:2, k_prior = c(.2, .3, .5)
    ),
    market
  )
  got <- bayes_cutoff_summary(fit)$by_market_year
  expect_null(got$firm_leader_probability[[1L]])
  expect_equal(
    unname(got$firm_strategic_probability[[1L]]),
    c(.8, .5, 0)
  )
  expect_equal(
    unname(got$realized_k_posterior[[1L]]),
    c(.2, .3, .5)
  )
  expect_true(is.na(got$cutoff_mean))
})

test_that("exact fringe reconstructs and simulates through coordination", {
  skip_if_not_installed("coordination", minimum_version = "0.0.0.9003")
  market <- bayes_market(
    prices = c(10, 11, 12), quantities = c(10, 8, 6),
    margins = c(.3, .28, .25), group = c("a", "b", "c"),
    outside_price = 1, market_side = "output", insideSize = 100,
    labels = c("p1", "p2", "p3")
  )
  fit <- bayes_specify(
    bayes_params(
      alpha = 1.4, s0 = .25, supply_model = 1L, structure = "fringe",
      k_support = 0:3, k_prior = c(.1, .2, .3, .4)
    ),
    market
  )
  none <- bayes_core_fringe_model(fit, fixed_k = 0L, role_mapping = "map")
  expect_s4_class(none, "CoreFringeLogit")
  expect_length(none@corePre, 0L)
  expect_equal(unname(none@pricePre), market$prices, tolerance = 1e-7)

  two <- bayes_core_fringe_model(fit, fixed_k = 2L, role_mapping = "map")
  expect_identical(two@corePre, c("a", "b"))
  expect_equal(two@diagnostics$bayes_reconstruction$selected_k, 2L)

  sparse_support <- bayes_specify(
    bayes_params(
      alpha = 1.4, s0 = .25, supply_model = 1L, structure = "fringe",
      k_support = c(0L, 2L), k_prior = c(.4, .6)
    ),
    market
  )
  expect_error(
    bayes_core_fringe_model(sparse_support, fixed_k = 1L),
    "supplied support"
  )

  sim <- bayes_sim(
    fit, n_draws = 1L, fixed_k = 1L, mcDelta = c(.02, 0, 0),
    mc.cores = 1L
  )
  expect_equal(sim@n_valid[[1L]], 1L)
  expect_identical(sim@metadata$simulation_route,
                   "coordination::core_fringe_simulate")
  expect_identical(sim@metadata$k_cache[[1L]][[1L]], 1L)
  expect_true(all(is.finite(sim@price_post[[1L]][1, ])))

  fit_ces <- bayes_specify(
    bayes_params(
      gamma = 2.4, s0 = .25, supply_model = 1L, structure = "fringe",
      k_support = 0:3, k_prior = rep(.25, 4)
    ),
    market
  )
  ces <- bayes_core_fringe_model(fit_ces, fixed_k = 2L)
  expect_s4_class(ces, "CoreFringeCES")
  expect_identical(ces@corePre, c("a", "b"))
  expect_equal(unname(ces@pricePre), market$prices, tolerance = 1e-7)
})

test_that("exact fringe maps merged core ownership and permits k zero", {
  skip_if_not_installed("coordination", minimum_version = "0.0.0.9003")
  market <- bayes_market(
    prices = c(10, 11, 12), quantities = c(12, 8, 5),
    margins = c(.3, .28, .25), group = c("a", "b", "c"),
    outside_price = 1, market_side = "output", insideSize = 100,
    labels = c("p1", "p2", "p3")
  )
  fit <- bayes_specify(
    bayes_params(
      alpha = 1.3, s0 = .2, supply_model = 3L, structure = "fringe",
      k_support = 0:3, k_prior = rep(.25, 4)
    ),
    market
  )
  merged <- bayes_sim(
    fit, ownerPost = c("ab", "ab", "c"), n_draws = 1L,
    fixed_k = 1L, mc.cores = 1L
  )
  expect_identical(merged@sims[[1L]][[1L]]@corePost, "ab")

  all_moncom <- bayes_sim(
    fit, ownerPost = c("ab", "ab", "c"), n_draws = 1L,
    fixed_k = 0L, corePost = character(), mc.cores = 1L
  )
  expect_length(all_moncom@sims[[1L]][[1L]]@corePost, 0L)
  expect_identical(all_moncom@metadata$k_cache[[1L]][[1L]], 0L)
})
