## Independent public-API regression fixture.  The posterior object is a
## deliberately small estimated-fit shell; the CmdStanR draw adapter is mocked so this
## checks the public draw-plan/simulation path without running Stan.

library(antitrustBayes)

.nonzero_input_fit <- function(supply_model = 1L) {
  prices <- c(1.1, 1.2, 1.3)
  quantities <- c(40, 30, 20)
  data <- data.frame(
    prices = prices,
    quantities = quantities,
    margins = c(.50, .40, .30),
    market = rep("m1", 3L),
    year = rep("1", 3L),
    group = c("a", "b", "c"),
    loo_cluster = rep("m1", 3L),
    margin_observed = rep(TRUE, 3L),
    outside_price = rep(1, 3L),
    stringsAsFactors = FALSE
  )
  stan_data <- list(
    quality_parameterization = 1L,
    N_group = 3L,
    group = 1:3
  )
  attr(stan_data, "model_spec_version") <- 13L
  attr(stan_data, "quality_scope") <- "global"
  attr(stan_data, "group_levels") <- c("a", "b", "c")
  methods::new(
    "BayesLogit", stanfit = NULL, stan_data = stan_data,
    data = data, supply_model = as.integer(supply_model), demand_model = "logit",
    market_side = "input", price_sd = 1,
    call = quote(fake_nonzero_quality_fit()),
    metadata = list(
      quality_parameterization = "sampled",
      posterior_draws = 2L,
      structure = "simultaneous",
      partition_prior = "cutoff",
      model_spec_version = 13L
    )
  )
}

.nonzero_input_draws <- function() {
  q <- rbind(
    c(.20, -.10, -.10),
    c(-.15, .05, .10)
  )
  base_score <- log(c(40, 30, 20) / 90)
  list(
    a_market = matrix(c(.8, .8), nrow = 2L, ncol = 1L),
    s0_prob = matrix(c(.20, .20), nrow = 2L, ncol = 1L),
    b_quality_all = q,
    ## Quality is part of the latent score here, so the fixture exercises a
    ## nonzero q state rather than merely carrying an unused generated field.
    share_mean = matrix(rep(base_score, 2L), nrow = 2L, byrow = TRUE) + q
  )
}

test_that("public input Logit keeps nonzero quality fixed through a value shock", {
  draws <- .nonzero_input_draws()
  owner <- c("a", "b", "c")
  value_shock <- antitrust::counterfactual(
    ownership = owner, costs = rep(.10, 3L)
  )

  testthat::with_mocked_bindings(
    .cmdstan_extract = function(fit, variables = NULL, ...) draws,
    .package = "antitrustBayes",
    {
      for (supply_model in c(1L, 3L)) {
      fit <- .nonzero_input_fit(supply_model)
      params <- bayes_params(alpha = .8, s0 = .2, supply_model = supply_model)
      market <- bayes_market(
        prices = c(1.1, 1.2, 1.3), quantities = c(40, 30, 20),
        margins = c(.50, .40, .30), group = owner, outside_price = 1,
        market_side = "input", insideSize = 90,
        labels = paste0("p", seq_along(owner))
      )
      baseline <- suppressWarnings(antitrustBayes:::build_conduct_object(
        params, 1L, market, ownerPre = owner, ownerPost = owner
      ))
      value_pre <- unname(antitrustBayes:::as_antitrust_fit(baseline)@model@mcPre)
      plan <- posterior_policy_draws(fit, n_draws = 2L, seed = 19L)
      pre <- suppressWarnings(simulate(
        fit, owner, draw_plan = plan, mc.cores = 1L,
        demand_anchoring = "latent_quality"
      ))
      post <- suppressWarnings(simulate(
        fit, value_shock, draw_plan = plan, mc.cores = 1L,
        demand_anchoring = "latent_quality"
      ))

      expect_equal(as.integer(pre@n_valid), 2L)
      expect_equal(as.integer(post@n_valid), 2L)
      expect_equal(pre@metadata$draw_indices, post@metadata$draw_indices)
      expect_equal(pre@metadata$latent_diagnostics$m1[[1L]]$quality,
                   post@metadata$latent_diagnostics$m1[[1L]]$quality)
      expect_true(any(abs(plan$quality_draws) > 1e-8))

      for (i in seq_len(2L)) {
        pre_model <- pre@sims$m1[[i]]
        post_model <- post@sims$m1[[i]]
        pre_diag <- pre@metadata$latent_diagnostics$m1[[i]]
        post_diag <- post@metadata$latent_diagnostics$m1[[i]]

        ## The observed value primitive is recovered before either policy.
        expect_equal(unname(pre_model@mcPre), value_pre, tolerance = 1e-10)
        expect_equal(unname(post_model@mcPre), unname(pre_model@mcPre),
                     tolerance = 1e-10)
        expect_equal(unname(pre_model@mcPost), unname(pre_model@mcPre),
                     tolerance = 1e-10)

        ## A 10 percent value/cost shock is applied after the pre-policy
        ## primitive is recovered, and q remains fixed across both scenarios.
        expect_equal(unname(post_model@mcPost),
                     1.1 * unname(post_model@mcPre), tolerance = 1e-10)
        expect_equal(pre_diag$quality, post_diag$quality, tolerance = 1e-12)
        expect_true(any(abs(pre_diag$quality) > 1e-8))
        expect_equal(sum(pre_diag$reference_quantity_shares), .8,
                     tolerance = 1e-10)
        expect_equal(sum(post_diag$reference_quantity_shares), .8,
                     tolerance = 1e-10)

        ## q changes the reconstructed baseline state, while the solved input
        ## FOC still reports a positive dollar markdown.
        expect_true(pre_diag$max_abs_pre_price_discrepancy > 1e-8)
        expect_true(post_diag$max_abs_pre_price_discrepancy > 1e-8)
        expect_lt(pre_diag$max_abs_foc_at_latent_pre, 1e-6)
        expect_lt(post_diag$max_abs_foc_at_latent_pre, 1e-6)
        pre_resid <- unname(pre_model@mcPre) - unname(pre_model@pricePre) -
          antitrust::calcMargins(pre_model, preMerger = TRUE, level = TRUE)
        post_resid <- unname(post_model@mcPost) - unname(post_model@pricePost) -
          antitrust::calcMargins(post_model, preMerger = FALSE, level = TRUE)
        expect_lt(max(abs(pre_resid)), 1e-6)
        expect_lt(max(abs(post_resid)), 1e-6)
      }
      }
    }
  )
})
