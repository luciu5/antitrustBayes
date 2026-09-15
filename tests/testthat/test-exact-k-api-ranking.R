test_that("exact ordered-k validates positional support and prior", {
  expect_error(
    .stack_validate_k_support(c(0, 1), "stackelberg"),
    "strictly positive"
  )
  expect_error(
    .stack_validate_k_support(c(0, 0), "fringe"),
    "unique"
  )
  expect_error(
    .stack_validate_k_support(c(0, 3), "fringe", c(2L, 4L)),
    "feasible in every selected market"
  )
  expect_error(
    .stack_validate_public_k(1:3, c(.2, .8), "stackelberg"),
    "same length"
  )
  expect_error(
    .stack_validate_public_k(1:2, c(.2, .7), "stackelberg"),
    "sum to one"
  )

  firm <- list(firm_market_n = c(3L, 4L), firm_mky_n = c(3L, 3L, 4L))
  got <- .stack_partition_arrays(
    firm, mky_to_market = c(1L, 1L, 2L), market_levels = c("m1", "m2"),
    k_support = c(0L, 2L), k_prior = NULL, k_scope = "market",
    partition_prior = "ordered_k", structure = "fringe"
  )
  expect_identical(got$partition_supports, list(c(0L, 2L), c(0L, 2L)))
  expect_equal(got$k_prior_by_unit, list(c(.5, .5), c(.5, .5)))
  expect_identical(got$partition_component_k,
                   c(0L, 2L, 0L, 2L))
})

test_that("one fixed firm-market rank applies across the complete panel", {
  # In m1, A wins period one and B wins period two by more, so B has the
  # market-panel rank one. C appears only in period one and its absent second
  # period contributes zero. m2 has an exact tie resolved by the firm label.
  market <- c(1, 1, 1, 1, 1, 2, 2)
  mky <- c(1, 1, 1, 2, 2, 3, 3)
  group <- c(1, 1, 2, 1, 2, 1, 2)
  share <- c(.30, .30, .35, .10, .90, .50, .50)
  firm <- .stack_unique_firm_mky(
    mky, group, share, share, market, c("A", "B"), "logit"
  )

  expect_identical(firm$firm_market_group_label, c("B", "A", "A", "B"))
  expect_equal(firm$firm_market_shareIn, c(.625, .35, .5, .5))
  expect_identical(firm$firm_market_rank, c(1L, 2L, 1L, 2L))
  expect_identical(firm$obs_to_firm_market, c(2L, 2L, 1L, 2L, 1L, 3L, 4L))

  inherited <- firm$firm_mky_rank[firm$obs_to_firm_mky]
  expect_identical(inherited, c(2L, 2L, 1L, 2L, 1L, 1L, 2L))
  expect_equal(firm$firm_mky_shareIn[firm$obs_to_firm_mky],
               c(.6, .6, .35, .1, .9, .5, .5))
})

test_that("ranking and observation mappings survive a row permutation", {
  market <- c(1, 1, 1, 1, 2, 2)
  mky <- c(1, 1, 2, 2, 3, 3)
  group <- c(1, 2, 1, 2, 1, 2)
  share <- c(.7, .3, .2, .8, .4, .6)
  labels <- c("A", "B")
  base <- .stack_unique_firm_mky(mky, group, share, share, market, labels)
  perm <- c(6, 2, 4, 1, 5, 3)
  moved <- .stack_unique_firm_mky(
    mky[perm], group[perm], share[perm], share[perm], market[perm], labels
  )
  expect_identical(moved$firm_market_group_label,
                   base$firm_market_group_label)
  expect_equal(moved$firm_market_shareIn, base$firm_market_shareIn)
  expect_identical(moved$firm_mky_group_label, base$firm_mky_group_label)
  expect_equal(moved$firm_mky_shareIn, base$firm_mky_shareIn)
  expect_identical(moved$firm_market_group_label[moved$obs_to_firm_market],
                   labels[group[perm]])
  expect_identical(moved$firm_mky_group_label[moved$obs_to_firm_mky],
                   labels[group[perm]])
})

test_that("Stack data keeps caller rows and uses separate traversal maps", {
  market <- c(2L, 1L, 2L, 1L)
  mky <- c(2L, 1L, 2L, 1L)
  group <- c(2L, 1L, 1L, 2L)
  share <- c(.4, .6, .6, .4)
  shell <- list(N_market = 2L, mky_to_market = c(1L, 2L))
  attr(shell, "market_levels") <- c("m1", "m2")
  got <- .stack_append_data(
    shell, mky, market, group, c("A", "B"), share, share,
    structure = list(name = "fringe"), partition_prior = "ordered_k",
    k_support = c(0L, 1L), k_prior = c(.25, .75)
  )
  expect_identical(got$structure, 1L)
  expect_identical(as.integer(got$market_obs_order), c(2L, 4L, 1L, 3L))
  expect_identical(as.integer(got$market_obs_start), c(1L, 3L))
  expect_identical(as.integer(got$market_obs_n), c(2L, 2L))
  ranking <- attr(got, "firm_market_ranking", exact = TRUE)
  expect_identical(ranking$firm_market_group_label[got$obs_to_firm_market],
                   c("B", "A", "A", "B"))
  expect_equal(got$partition_prior_probability,
               array(c(.25, .75, .25, .75), dim = 4L))
})

test_that("preparers route explicit exact structures and ignore k for simultaneous", {
  args <- list(
    prices = c(2, 2.2, 2.4), quantities = c(5, 3, 2),
    margins = c(.3, .28, .25), market = rep("m", 3),
    year = rep("y", 3), group = c("A", "B", "C"),
    supply_model = 1L, use_rho = FALSE, use_hmt = FALSE
  )
  simultaneous <- do.call(prepare_stan_data, c(
    args,
    list(structure = "simultaneous", k_support = c(NA, -4),
         k_prior = "deliberately ignored")
  ))
  expect_identical(attr(simultaneous, "structure"), "simultaneous")
  expect_null(attr(simultaneous, "k_support", exact = TRUE))
  expect_null(simultaneous$N_partition_component)

  fringe <- do.call(prepare_stan_data, c(
    args,
    list(structure = "fringe", k_support = 0:2,
         k_prior = c(.2, .3, .5))
  ))
  expect_identical(attr(fringe, "partition_prior"), "ordered_k")
  expect_identical(as.integer(fringe$partition_component_k), 0:2)
  expect_identical(fringe$structure, 1L)
  expect_null(attr(fringe, "cutoff_scale", exact = TRUE))
  expect_equal(
    unname(attr(fringe, "partition_metadata")$k_prior_by_unit[[1L]]),
    c(.2, .3, .5)
  )

  expect_error(
    do.call(prepare_stan_data, c(args, list(
      structure = "stackelberg", k_support = 0:2,
      k_prior = c(.2, .3, .5)
    ))),
    "strictly positive"
  )
})

test_that("specified exact fits retain requested and resolved k provenance", {
  params <- bayes_params(
    alpha = .7, s0 = .2, supply_model = 1L,
    structure = "stackelberg", k_support = c(1L, 3L)
  )
  market <- bayes_market(
    prices = c(2, 2.2, 2.4), quantities = c(5, 3, 2),
    margins = c(.3, .28, .25), group = c("A", "B", "C"),
    outside_price = 1, market_side = "output", insideSize = 10
  )
  fit <- bayes_specify(params, market)
  expect_identical(fit@metadata$k_support, c(1L, 3L))
  expect_null(fit@metadata$k_prior)
  expect_identical(fit@metadata$k_support_resolved[[1L]], c(1L, 3L))
  expect_equal(fit@metadata$k_prior_resolved[[1L]], c(.5, .5))
  expect_identical(
    fit@metadata$firm_market_ranking$tie_break,
    "descending_share_then_ascending_firm_label"
  )
})

test_that("CES exact ranking uses conditional inside revenue shares", {
  d <- prepare_ces_stan_data(
    prices = c(1, 4, 2, 2), quantities = c(8, 2, 3, 7),
    margins = rep(.25, 4), market = rep("m", 4),
    year = c("y1", "y1", "y2", "y2"), group = c("A", "B", "A", "B"),
    supply_model = 3L, structure = "fringe", k_support = 0:2,
    k_prior = c(.2, .3, .5), market_side = "output",
    use_cutoff = FALSE, use_rho = FALSE, use_hmt = FALSE, priors = list()
  )
  ranking <- attr(d, "firm_market_ranking", exact = TRUE)
  # Period revenue shares are tied in y1 and are 0.3/0.7 in y2; B ranks first.
  expect_identical(ranking$metric, "conditional_inside_revenue_share")
  expect_identical(ranking$firm_market_group_label, c("B", "A"))
  expect_equal(ranking$firm_market_revShare, c(.6, .4))
  expect_identical(d$structure, 1L)
  expect_identical(as.integer(d$partition_component_k), 0:2)
})
