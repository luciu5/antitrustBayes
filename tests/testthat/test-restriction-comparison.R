make_restriction_fit <- function(markets = c("m1", "m2"),
                                 support = c(1L, 2L),
                                 posterior_draws = 2L) {
  n_market <- length(markets)
  n_component <- n_market * length(support)
  stan_data <- list(
    N_partition_unit = as.integer(n_market),
    partition_unit_market = seq_len(n_market),
    partition_unit_start = seq.int(1L, n_component, by = length(support)),
    partition_unit_count = rep.int(length(support), n_market),
    partition_component_unit = rep(seq_len(n_market), each = length(support)),
    partition_component_k = rep(support, n_market)
  )
  attr(stan_data, "market_levels") <- markets
  methods::new(
    "BayesLogit", stanfit = NULL,
    data = data.frame(market = markets, stringsAsFactors = FALSE),
    stan_data = stan_data, price_sd = 1,
    supply_model = 1L, demand_model = "logit", market_side = "input",
    call = quote(fake()),
    metadata = list(
      structure = "stackelberg", partition_prior = "ordered_k",
      posterior_draws = as.integer(posterior_draws)
    )
  )
}

test_that("one-market composite restrictions preserve support order", {
  fit <- make_restriction_fit("m1", support = c(3L, 1L, 2L))
  draws <- list(
    posterior_component_prob = rbind(c(.2, .5, .3), c(.6, .1, .3)),
    prior_component_prob = rbind(c(.2, .5, .3), c(.2, .5, .3))
  )
  got <- bayes_restriction_compare(
    fit,
    list(singleton = c(2L), composite = c(2L, 3L, 2L, 1L)),
    reference = "composite", draws = draws
  )
  expect_equal(got$prior_mass, c(.3, 1))
  expect_equal(got$posterior_mass, c(.3, 1))
  expect_equal(got$log_evidence_vs_full, c(0, 0))
  expect_equal(got$log_bf_vs_reference, c(0, 0))
  expect_equal(got$bf_vs_reference, c(1, 1))
  expect_identical(attr(got, "method"), "encompassing_prior")
})

test_that("joint restrictions average products within draws", {
  fit <- make_restriction_fit()
  ## The market probabilities move together.  The posterior joint mass is
  ## mean(c(1, 0)) = .5, whereas the product of marginal means is .25.
  draws <- list(
    posterior_component_prob = rbind(
      c(1, 0, 1, 0), c(0, 1, 0, 1)
    ),
    prior_component_prob = matrix(.5, nrow = 2L, ncol = 4L)
  )
  got <- bayes_restriction_compare(
    fit,
    list(joint = list(m1 = 1L, m2 = 1L),
         partial = list(m1 = c(1L, 2L))),
    draws = draws
  )
  expect_equal(got$prior_mass, c(.25, 1))
  expect_equal(got$posterior_mass, c(.5, 1))
  expect_equal(got$log_evidence_vs_full, c(log(2), 0))
})

test_that("share-only posterior probabilities give unit evidence", {
  fit <- make_restriction_fit("m1", support = c(1L, 2L))
  p <- rbind(c(.8, .2), c(.3, .7))
  got <- bayes_restriction_compare(
    fit, list(first = 1L, all = c(1L, 2L)),
    draws = list(posterior_component_prob = p, prior_component_prob = p)
  )
  expect_equal(got$log_evidence_vs_full, c(0, 0))
  expect_equal(got$posterior_mass, got$prior_mass)
})

test_that("zero posterior event mass is reported as negative infinity", {
  fit <- make_restriction_fit("m1", support = c(1L, 2L))
  got <- bayes_restriction_compare(
    fit, list(impossible = 2L),
    draws = list(
      posterior_component_prob = rbind(c(1, 0), c(1, 0)),
      prior_component_prob = matrix(c(.5, .5, .5, .5), nrow = 2L)
    )
  )
  expect_equal(got$prior_mass, .5)
  expect_equal(got$posterior_mass, 0)
  expect_identical(got$log_evidence_vs_full, -Inf)
})

test_that("zero prior restriction fails clearly", {
  fit <- make_restriction_fit("m1", support = c(1L, 2L))
  expect_error(
    bayes_restriction_compare(
      fit, list(outside_support = 3L),
      draws = list(
        posterior_component_prob = matrix(c(.5, .5, .5, .5), nrow = 2L),
        prior_component_prob = rbind(c(1, 0), c(1, 0))
      )
    ),
    "zero prior mass"
  )
})

test_that("specified exact fits use their coherent prior as posterior", {
  fit <- methods::new(
    "BayesSpecifiedFit", params = list(n_draws = 4L), market = list(),
    data = data.frame(market = "unit_1"), supply_model = 1L,
    demand_model = "logit", market_side = "input", call = quote(fake()),
    metadata = list(
      structure = "fringe", partition_prior = "ordered_k",
      k_support_resolved = list(c(0L, 2L)),
      k_prior_resolved = list(c(.25, .75))
    )
  )
  got <- bayes_restriction_compare(fit, list(core = 2L))
  expect_equal(got$prior_mass, .75)
  expect_equal(got$posterior_mass, .75)
  expect_equal(got$log_evidence_vs_full, 0)
  expect_false("bf_vs_reference" %in% names(got))
})

test_that("malformed maps and probabilities fail through exact-k validation", {
  fit <- make_restriction_fit("m1", support = c(1L, 2L))
  bad <- list(
    posterior_component_prob = matrix(c(.5, NA, .5, .5), nrow = 2L),
    prior_component_prob = matrix(.5, nrow = 2L, ncol = 2L)
  )
  expect_error(
    bayes_restriction_compare(fit, list(one = 1L), draws = bad),
    "contains invalid probabilities"
  )
  expect_error(
    bayes_restriction_compare(fit, list(one = 1.5), draws =
      list(posterior_component_prob = matrix(.5, 2L, 2L),
           prior_component_prob = matrix(.5, 2L, 2L))),
    "finite integer k"
  )
})
