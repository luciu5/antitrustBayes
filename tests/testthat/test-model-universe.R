make_model_universe_exact_fit <- function(markets = c("m1", "m2"),
                                           n_firms = c(3L, 2L),
                                           support = list(0:2, 0:1),
                                           structure = "fringe",
                                           supply_model = 1L,
                                           share_likelihood = "gaussian",
                                           margin_likelihood = "gaussian",
                                           quality_parameterization = NULL) {
  n_unit <- length(markets)
  stopifnot(length(n_firms) == n_unit, length(support) == n_unit)
  n_component <- sum(lengths(support))
  starts <- cumsum(c(1L, head(lengths(support), -1L)))
  ranking <- list(
    firm_market_market = rep(seq_len(n_unit), n_firms),
    firm_market_group_label = unlist(lapply(seq_len(n_unit), function(u) {
      paste0("m", u, "_firm", seq_len(n_firms[[u]]))
    }), use.names = FALSE),
    firm_market_rank = unlist(lapply(n_firms, seq_len), use.names = FALSE)
  )
  stan_data <- list(
    N_partition_unit = as.integer(n_unit),
    partition_unit_market = seq_len(n_unit),
    partition_unit_start = starts,
    partition_unit_count = lengths(support),
    partition_component_unit = rep(seq_len(n_unit), lengths(support)),
    partition_component_k = unlist(support, use.names = FALSE),
    partition_prior_probability = unlist(lapply(support, function(k) rep(1 / length(k), length(k))),
                                        use.names = FALSE),
    partition_prior_logprob = unlist(lapply(support, function(k) rep(-log(length(k)), length(k))),
                                    use.names = FALSE),
    likelihood_mode = 1L,
    shareIn = rep(.4, 4L), marginNormInv = rep(.5, 4L),
    margin_obs = rep(1L, 4L), market_year_idx = rep(seq_len(n_unit), each = 2L),
    use_rho = 0L, use_s0_market_pool = 0L, use_s0_market_year_pool = 0L,
    use_alpha_market_pool = 0L, use_market_resid_scale = 0L, use_hmt = 0L
  )
  attr(stan_data, "market_levels") <- markets
  attr(stan_data, "k_scope") <- "market"
  attr(stan_data, "structure") <- "fringe"
  attr(stan_data, "partition_prior") <- "ordered_k"
  attr(stan_data, "share_likelihood") <- "gaussian"
  attr(stan_data, "margin_likelihood") <- "gaussian"
  observed_groups <- unlist(lapply(seq_len(n_unit), function(u) {
    paste0("m", u, "_firm", seq_len(n_firms[[u]]))
  }), use.names = FALSE)
  stan_data$N <- length(observed_groups)
  stan_data$market <- rep(seq_len(n_unit), n_firms)
  stan_data$group <- match(observed_groups, sort(unique(observed_groups)))
  stan_data$N_group <- length(unique(observed_groups))
  attr(stan_data, "group_levels") <- sort(unique(observed_groups))
  attr(stan_data, "quality_scope") <- "global"
  attr(stan_data, "quality_incidence") <- antitrustBayes:::.build_quality_incidence(
    stan_data$market, stan_data$group, markets,
    attr(stan_data, "group_levels"), quality_scope = "global")
  data <- data.frame(
    prices = rep(1, length(observed_groups)),
    quantities = rep(1, length(observed_groups)),
    margins = rep(.5, length(observed_groups)),
    market = rep(markets, n_firms), year = rep("1", length(observed_groups)),
    group = observed_groups, margin_observed = TRUE, outside_price = 1,
    stringsAsFactors = FALSE
  )
  methods::new(
    "BayesLogit", stanfit = NULL, data = data, stan_data = stan_data,
    price_sd = 1, supply_model = as.integer(supply_model), demand_model = "logit",
    market_side = "input", call = quote(fake()), metadata = list(
      structure = structure, partition_prior = "ordered_k", k_scope = "market",
      supply_model = as.integer(supply_model), demand_model = "logit", prior_manifest = "baseline",
      resolved_priors = list(prior_alpha_sd = 1), firm_market_ranking = ranking,
      share_likelihood = share_likelihood, margin_likelihood = margin_likelihood,
      quality_parameterization = quality_parameterization, quality_scope = "global",
      likelihood = "joint", posterior_draws = 2L,
      prior_normalization = list(mode = "native_stan_normalized", model_spec_version = 12L)
    )
  )
}

test_that("the public universe route accepts named collections and named fits", {
  fit <- make_model_universe_exact_fit()
  universe <- bayes_model_universe(
    list(BF = list(fit = "source", k = "interior")), c(BF = 1))
  testthat::local_mocked_bindings(
    .bayes_compare_universe = function(args, ...) args,
    .package = "antitrustBayes")
  expect_identical(bayes_compare(list(source = fit), universe = universe),
                   list(source = fit))
  expect_identical(bayes_compare(source = fit, universe = universe),
                   list(source = fit))
})

test_that("prior identity uses active Stan values rather than provenance labels", {
  fit <- make_model_universe_exact_fit()
  fit@stan_data$prior_alpha_sd <- 1
  original <- .mu_fingerprint(fit)
  changed <- fit
  changed@stan_data$prior_alpha_sd <- 2
  expect_error(.mu_fingerprint(changed), "Active prior metadata")
  changed@metadata$resolved_priors$prior_alpha_sd <- 2
  expect_false(identical(.mu_fingerprint(changed), original))
  changed <- fit
  changed@metadata$prior_manifest <- "diffuse"
  expect_identical(.mu_fingerprint(changed), original)
  changed@metadata$prior_normalization$resolved_priors <- list(prior_alpha_sd = 2)
  expect_error(.mu_fingerprint(changed), "Active prior metadata")
  changed <- fit
  changed@stan_data$prior_cutoff_alpha <- 3
  changed@metadata$resolved_priors$prior_cutoff_alpha <- 5
  expect_identical(.mu_fingerprint(changed), original)
})

test_that("the descriptor requires explicit named scientific prior weights", {
  expect_error(bayes_model_universe(list(B = "b")), "prior_weights")
  expect_error(
    bayes_model_universe(list(B = "b", C = "c"), c(.5, .5, .0)),
    "one entry"
  )
  u <- bayes_model_universe(
    list(B = "b", BF = list(fit = "bf", k = "interior")),
    c(B = .6, BF = .4)
  )
  expect_s3_class(u, "bayes_model_universe")
  expect_identical(u$prior_weight_provenance$rule,
                   "explicit_named_scientific_prior")
})

test_that("strict interior uses the ranked firm union in every panel market", {
  fit <- make_model_universe_exact_fit()
  descriptor <- bayes_model_universe(
    list(BF = list(fit = "bf", k = "interior")), c(BF = 1)
  )
  rows <- .mu_fit_descriptor_rows(list(bf = fit), descriptor)
  expect_equal(rows[[1L]]$resolved$allowed, list(m1 = c(1L, 2L), m2 = 1L))

  multiproduct <- fit
  multiproduct@data <- rbind(
    multiproduct@data,
    multiproduct@data[multiproduct@data$market == "m1" &
                        multiproduct@data$group == "m1_firm1", , drop = FALSE]
  )
  multi_rows <- .mu_fit_descriptor_rows(list(bf = multiproduct), descriptor)
  expect_equal(multi_rows[[1L]]$resolved$n_union, c(m1 = 3L, m2 = 2L))

  endpoint <- bayes_model_universe(
    list(BF0 = list(fit = "bf", k = 0L)), c(BF0 = 1)
  )
  expect_error(.mu_fit_descriptor_rows(list(bf = fit), endpoint), "market-to-k")

  inconsistent <- fit
  inconsistent@data$group[inconsistent@data$market == "m1"] <-
    c("m1_firm1", "m1_firm2", "m1_phantom")
  expect_error(
    .mu_fit_descriptor_rows(list(bf = inconsistent), descriptor),
    "ranked firm union"
  )
})

test_that("overlapping exact family rectangles and duplicate atomics are rejected", {
  fit <- make_model_universe_exact_fit()
  one <- bayes_model_universe(
    list(BF1 = list(fit = "bf", k = list(m1 = 1L, m2 = 1L)),
         BF2 = list(fit = "bf", k = "interior")),
    c(BF1 = .5, BF2 = .5)
  )
  expect_error(.mu_fit_descriptor_rows(list(bf = fit), one), "overlap")
})

test_that("multimarket evidence uses a finite-mixture identity", {
  rows <- list(
    list(hypothesis = "A", source = "F", kind = "family", supply = "bertrand",
         structure = "fringe", fingerprint = "same", resolved = list()),
    list(hypothesis = "B", source = "G", kind = "atomic", supply = "cournot",
         structure = "simultaneous", fingerprint = "other", resolved = NULL)
  )
  ## log BF(A,F) = log(.5/.25) = log(2), while a source log evidence of 3
  ## yields log Z_A = 3 + log(2).  The actual restriction helper supplies
  ## the within-draw product before averaging; this fixture checks the final
  ## finite-mixture normalization and explicit prior weights.
  rows[[1L]]$evidence_override <- list(
    prior_mass = .25, posterior_mass = .5, log_evidence_vs_full = log(2),
    method = "encompassing_prior"
  )
  got <- .mu_result_from_evidence(
    rows, data.frame(model = c("F", "G"), logml = c(3, 2), logml_sd = c(.1, .1)),
    c(A = .25, B = .75)
  )$table
  expect_equal(got$logml_restricted, c(3 + log(2), 2))
  expect_equal(sum(got$posterior_probability), 1)
  expect_gt(got$posterior_probability[[1L]], got$posterior_probability[[2L]])
})

test_that("the existing restriction helper supplies the shared-theta correction", {
  fit <- make_model_universe_exact_fit(
    n_firms = c(3L, 3L), support = list(0:1, 0:1)
  )
  descriptor <- bayes_model_universe(
    list(BF_joint = list(fit = "bf", k = list(m1 = 1L, m2 = 1L))),
    c(BF_joint = 1)
  )
  row <- .mu_fit_descriptor_rows(list(bf = fit), descriptor)[[1L]]
  draws <- list(
    posterior_component_prob = rbind(
      c(.2, .8, .2, .8), c(.8, .2, .8, .2)
    ),
    prior_component_prob = matrix(.5, nrow = 2L, ncol = 4L)
  )
  restriction <- .mu_restriction_result(row, draws = draws)
  expect_equal(restriction$prior_mass, .25)
  expect_equal(restriction$posterior_mass, .34)
  expect_equal(restriction$log_evidence_vs_full, log(.34 / .25))

  result <- .mu_result_from_evidence(
    list(row), data.frame(model = "bf", logml = log(6.25), logml_sd = NA_real_),
    c(BF_joint = 1), restriction_results = list(as.list(restriction[1L, ]))
  )
  expect_equal(exp(result$table$logml_restricted), 8.5)
  expect_equal(result$table$posterior_probability, 1)
  ## The draw-level posterior mass is E[.8*.8, .2*.2] = .34; multiplying
  ## marginal draw means would incorrectly give .25.
  expect_equal(.34, mean(c(.8 * .8, .2 * .2)))
  expect_equal(c(16 / 17, 1 / 17), c(.8 * .8, .2 * .2) / .68)
})

test_that("incomplete positive-prior evidence withholds all universe probabilities", {
  rows <- list(list(hypothesis = "A", source = "F", kind = "atomic",
                    supply = "bertrand", structure = "simultaneous",
                    fingerprint = "same", resolved = NULL))
  got <- .mu_result_from_evidence(
    rows, data.frame(model = "F", logml = NA_real_), c(A = 1)
  )
  expect_true(all(is.na(got$table$posterior_probability)))
  expect_identical(got$status, "incomplete_evidence")

  underflow_rows <- list(list(hypothesis = "A", source = "F", kind = "family",
                              supply = "bertrand", structure = "fringe",
                              fingerprint = "same", resolved = list()))
  underflow_rows[[1L]]$evidence_override <- list(
    prior_mass = .25, posterior_mass = 0, log_evidence_vs_full = log(2),
    method = "encompassing_prior"
  )
  underflow <- .mu_result_from_evidence(
    underflow_rows, data.frame(model = "F", logml = 1), c(A = 1)
  )
  expect_equal(underflow$table$logml_restricted, 1 + log(2))
  expect_equal(underflow$table$posterior_probability, 1)
  expect_identical(underflow$status, "complete")
  expect_identical(underflow$table$evidence_status, "restriction_mass_underflow")

  inf <- .mu_result_from_evidence(
    rows, data.frame(model = "F", logml = Inf), c(A = 1)
  )
  expect_identical(inf$status, "incomplete_evidence")
})

test_that("partial source repetitions remain visible without suppressing finite universe probabilities", {
  rows <- lapply(c("F", "G"), function(source) list(hypothesis = source, source = source,
    kind = "atomic", supply_model = "bertrand", supply = "bertrand",
    structure = "simultaneous", fingerprint = source, resolved = NULL))
  evidence <- data.frame(model = c("F", "G"), logml = c(3, 2), logml_sd = c(.1, NA_real_))
  context <- list(F = list(status = "complete", requested_repetitions = 3L,
    completed_repetitions = 3L), G = list(status = "partial", requested_repetitions = 3L,
    completed_repetitions = 1L, repetitions = data.frame(rep = 1:3,
      status = c("complete", "failed", "failed"), logml = c(2, NA, NA))))
  attr(evidence, "bridge_context_info") <- context
  result <- .mu_result_from_evidence(rows, evidence, c(.3, .7))
  expected <- exp(c(3, 2)) * c(.3, .7); expected <- expected / sum(expected)
  expect_equal(result$table$posterior_probability, expected)
  expect_identical(result$status, "complete")
  expect_identical(result$numerical_status, "partial")
  expect_identical(result$table$point_estimate_available, c(TRUE, TRUE))
  expect_identical(result$table$numerical_status, c("complete", "partial"))
  expect_identical(result$table$requested_repetitions, c(3L, 3L))
  expect_identical(result$table$completed_repetitions, c(3L, 1L))
  expect_true(is.na(result$table$logml_sd[2L]))
  evidence$logml[2L] <- NA_real_
  result <- .mu_result_from_evidence(rows, evidence, c(.3, .7))
  expect_true(all(is.na(result$table$posterior_probability)))
  expect_identical(result$status, "incomplete_evidence")
  expect_identical(result$numerical_status, "unavailable")
  no_metadata <- data.frame(model = c("F", "G"), logml = c(3, 2))
  expect_identical(.mu_result_from_evidence(rows, no_metadata, c(.3, .7))$numerical_status, "unreported")
})

test_that("public universe result preserves mixed-backend repetition and sampler metadata", {
  left <- make_model_universe_exact_fit(structure = "simultaneous")
  right <- make_model_universe_exact_fit(structure = "simultaneous", supply_model = 3L)
  attr(left@stan_data, "model_spec_version") <- 12L
  attr(right@stan_data, "model_spec_version") <- 12L
  descriptor <- bayes_model_universe(list(B = "left", C = "right"), c(B = .5, C = .5))
  source <- data.frame(model = c("left", "right"), logml = c(-3, -4),
    logml_sd = c(.1, NA_real_), backend = c("cmdstanr", "nutpie"))
  context <- list(left = list(status = "complete", requested_repetitions = 3L,
    completed_repetitions = 3L), right = list(status = "partial", requested_repetitions = 3L,
    completed_repetitions = 1L, repetitions = data.frame(rep = 1:3,
      status = c("complete", "failed", "failed"), logml = c(-4, NA, NA))))
  sampler <- list(left = list(backend = "cmdstanr", divergences = 2L),
    right = list(backend = "nutpie", divergences = 4L))
  attr(source, "bridge_context_info") <- context
  attr(source, "backend_metadata") <- sampler
  testthat::local_mocked_bindings(bayes_compare = function(...) source, .package = "antitrustBayes")
  result <- .bayes_compare_universe(list(left = left, right = right), descriptor,
    bridge = TRUE, bridge_maxiter = 5000L, bridge_repetitions = 3L,
    bridge_method = "normal", loo_unit = "market")
  expect_identical(attr(result, "status"), "complete")
  expect_identical(attr(result, "numerical_status"), "partial")
  expect_identical(result$backend, c("cmdstanr", "nutpie"))
  expect_identical(attr(result, "bridge_context_info"), context)
  expect_identical(attr(result, "backend_metadata"), sampler)
  expect_true(all(is.finite(result$post_prob)))
})

test_that("sampler provenance cannot create duplicate scientific universe hypotheses", {
  left <- make_model_universe_exact_fit(structure = "simultaneous")
  right <- left
  left@metadata$backend <- "cmdstanr"
  right@metadata$backend <- "nutpie"
  right@metadata$seed <- 8271L
  right@metadata$sampler <- list(adaptation = "low_rank", lib_hash = "other-binary")
  expect_identical(.mu_fingerprint(left), .mu_fingerprint(right))
  descriptor <- bayes_model_universe(list(first = "left", second = "right"), c(first = .5, second = .5))
  expect_error(.mu_fit_descriptor_rows(list(left = left, right = right), descriptor), "duplicate")
})

test_that("matching statistical fingerprints govern duplicate detection", {
  a <- make_model_universe_exact_fit(structure = "simultaneous")
  b <- a
  b@metadata$quality_parameterization <- "collapsed"
  attr(a@stan_data, "likelihood_mode_metadata") <- list(
    share_likelihood = "gaussian", margin_likelihood = "gaussian",
    quality_parameterization = "sampled", quality_parameterization_code = 1L,
    mode = "gaussian/sampled/gaussian"
  )
  attr(b@stan_data, "likelihood_mode_metadata") <- list(
    share_likelihood = "gaussian", margin_likelihood = "gaussian",
    quality_parameterization = "collapsed", quality_parameterization_code = 2L,
    mode = "gaussian/collapsed/gaussian"
  )
  ## Gaussian sampled and collapsed quality are the same joint density.
  expect_identical(.mu_fingerprint(a), .mu_fingerprint(b))
  duplicate <- bayes_model_universe(list(B1 = "a", B2 = "b"), c(.5, .5))
  expect_error(.mu_fit_descriptor_rows(list(a = a, b = b), duplicate), "duplicate")

  centered <- a
  centered@metadata$s0_pooling_requested <- "auto"
  uncentered <- a
  uncentered@metadata$s0_pooling_requested <- "manual"
  expect_identical(.mu_fingerprint(centered), .mu_fingerprint(uncentered))

  different_likelihood <- make_model_universe_exact_fit(
    structure = "simultaneous", margin_likelihood = "student_t",
    quality_parameterization = "sampled"
  )
  different_demand <- different_likelihood
  different_demand@demand_model <- "ces"
  different_demand@metadata$demand_model <- "ces"
  expect_false(identical(.mu_fingerprint(a), .mu_fingerprint(different_likelihood)))
  expect_false(identical(.mu_fingerprint(a), .mu_fingerprint(different_demand)))

  hybrid_sample <- make_model_universe_exact_fit(
    structure = "simultaneous", margin_likelihood = "student_t",
    quality_parameterization = "collapsed"
  )
  hybrid_collapse <- hybrid_sample
  hybrid_sample@stan_data$use_rho <- 0L
  hybrid_collapse@stan_data$use_rho <- 0L
  hybrid_sample@metadata$quality_parameterization <- "sampled"
  attr(hybrid_sample@stan_data, "likelihood_mode_metadata") <- list(
    share_likelihood = "gaussian", margin_likelihood = "student_t",
    quality_parameterization = "sampled", mode = "gaussian/sampled/student_t"
  )
  attr(hybrid_collapse@stan_data, "likelihood_mode_metadata") <- list(
    share_likelihood = "gaussian", margin_likelihood = "student_t",
    quality_parameterization = "collapsed", mode = "gaussian/collapsed/student_t"
  )
  ## For the valid rho=0 Gaussian-share/fixed-t-margin hybrid, collapsing the
  ## Gaussian quality is the same integrated joint density.
  expect_identical(.mu_fingerprint(hybrid_sample), .mu_fingerprint(hybrid_collapse))
  kept <- bayes_model_universe(list(B = "a", B_t = "b"), c(.5, .5))
  expect_silent(.mu_fit_descriptor_rows(list(a = a, b = different_likelihood), kept))
})

test_that("ordered-k fingerprints retain effective full and conditional priors", {
  set_prior <- function(fit, prior) {
    fit@stan_data$partition_prior_probability <- unlist(prior, use.names = FALSE)
    fit@stan_data$partition_prior_logprob <- log(fit@stan_data$partition_prior_probability)
    fit@metadata$k_support_resolved <- list(m1 = 0:2, m2 = 0:1)
    fit@metadata$k_prior_resolved <- stats::setNames(prior, c("m1", "m2"))
    fit
  }
  a <- set_prior(make_model_universe_exact_fit(), list(c(.4, .2, .4), c(.5, .5)))
  b <- set_prior(a, list(c(.7, .1, .2), c(.9, .1)))
  ra <- .mu_resolve_restriction(list(k = "interior"), a, "A")
  rb <- .mu_resolve_restriction(list(k = "interior"), b, "B")
  expect_false(identical(.mu_fingerprint(a), .mu_fingerprint(b)))
  expect_identical(.mu_fingerprint(a, restriction = ra),
                   .mu_fingerprint(b, restriction = rb))
  descriptor <- bayes_model_universe(
    list(A = list(fit = "a", k = "interior"),
         B = list(fit = "b", k = "interior")), c(A = .5, B = .5))
  expect_error(.mu_fit_descriptor_rows(list(a = a, b = b), descriptor),
               "same conditional ordered-k prior")
  c <- set_prior(a, list(c(.3, .3, .4), c(.5, .5)))
  rc <- .mu_resolve_restriction(list(k = "interior"), c, "C")
  expect_false(identical(.mu_fingerprint(a, restriction = ra),
                         .mu_fingerprint(c, restriction = rc)))
  expect_length(.mu_fit_descriptor_rows(list(a = a, b = c), descriptor), 2L)
  reordered <- b
  reordered@metadata$k_support_resolved <- list(m2 = 1:0, m1 = 2:0)
  reordered@metadata$k_prior_resolved <- list(m2 = c(.1, .9), m1 = c(.2, .1, .7))
  expect_identical(.mu_fingerprint(b), .mu_fingerprint(reordered))
  normalized_location <- a
  normalized_location@metadata$prior_normalization$k_support_resolved <-
    normalized_location@metadata$k_support_resolved
  normalized_location@metadata$prior_normalization$k_prior_resolved <-
    normalized_location@metadata$k_prior_resolved
  normalized_location@metadata$k_support_resolved <- NULL
  normalized_location@metadata$k_prior_resolved <- NULL
  expect_identical(.mu_fingerprint(a), .mu_fingerprint(normalized_location))
  dual <- a
  dual@metadata$prior_normalization$k_support_resolved <- list(m1 = 0:2, m2 = 0:1)
  dual@metadata$prior_normalization$k_prior_resolved <- list(m1 = c(.2, .4, .4), m2 = c(.5, .5))
  expect_error(.mu_fingerprint(dual), "prior metadata disagrees")
  dual@metadata$prior_normalization$k_support_resolved <- list(m1 = 1:3, m2 = 0:1)
  expect_error(.mu_fingerprint(dual), "support metadata disagrees")
  stale <- a
  stale@metadata$k_prior_resolved[[1L]] <- c(.2, .4, .4)
  expect_error(.mu_fingerprint(stale), "prior metadata disagrees")
  stale <- a
  stale@metadata$k_support_resolved[[1L]] <- 1:3
  expect_error(.mu_fingerprint(stale), "support metadata disagrees")
  invalid <- a
  invalid@stan_data$partition_prior_logprob[[1L]] <- 0
  expect_error(.mu_fingerprint(invalid), "log probabilities disagree")
  invalid <- a
  invalid@stan_data$partition_prior_probability[[1L]] <- .5
  invalid@stan_data$partition_prior_logprob <- log(invalid@stan_data$partition_prior_probability)
  expect_error(.mu_fingerprint(invalid), "sum to one")
  zero <- set_prior(a, list(c(.4, 0, .6), c(.5, .5)))
  rz <- .mu_resolve_restriction(list(k = list(m1 = 1L, m2 = 1L)), zero, "Z")
  expect_true(.mu_exact_prior_state(zero, rz)$m1$zero_prior_event)
})

test_that("interior restrictions retain supplied support without enlarging it", {
  fit <- make_model_universe_exact_fit(markets = "m1", n_firms = 24L,
                                       support = list(0:5))
  resolved <- .mu_resolve_restriction(list(k = "interior"), fit, "BF")
  expect_identical(resolved$allowed$m1, 1:5)
  expect_error(.mu_resolve_restriction(list(k = 1:6), fit, "BF"), "unsupported")
})

test_that("dimension aggregation sums atomic probabilities without relabeling MC or auction", {
  tab <- data.frame(
    hypothesis = c("MC", "B", "BF", "A"), model = c("MC", "B", "BF", "A"),
    source = c("mc", "b", "bf", "a"), kind = c("atomic", "atomic", "family", "atomic"),
    supply_model = c("monopolistic", "bertrand", "bertrand", "auction"),
    structure = c("simultaneous", "simultaneous", "fringe", "simultaneous"),
    prior_weight = c(.2, .3, .3, .2),
    posterior_probability = c(.2, .3, .3, .2), stringsAsFactors = FALSE
  )
  class(tab) <- c("bayes_universe_comparison", "data.frame")
  attr(tab, "status") <- "complete"
  got <- bayes_dimensions(tab, c("conduct", "strategic_core"))
  expect_equal(got$posterior_probability[got$dimension == "conduct" & got$level == "price"], .6)
  expect_equal(got$posterior_probability[got$dimension == "conduct" & got$level == "monopolistic"], .2)
  expect_equal(got$posterior_probability[got$dimension == "conduct" & got$level == "auction"], .2)
  expect_equal(got$direction[got$dimension == "conduct" & got$level == "monopolistic"], "moncom")
  expect_equal(got$posterior_probability[got$dimension == "strategic_core" & got$level == "interior"], .3)
})

test_that("the universe table retains the ordinary bayes_weights contract", {
  tab <- data.frame(
    model = c("B", "C"), model_prior = c(.5, .5),
    post_prob = c(.25, .75), stringsAsFactors = FALSE
  )
  expect_equal(bayes_weights(tab), c(B = .25, C = .75))
})
test_that("inactive prior settings do not encode additional canonical hypotheses", {
  a <- make_model_universe_exact_fit(structure = "simultaneous")
  a@metadata$resolved_priors <- bayes_prior_manifest("baseline")
  a@metadata$resid_df_fixed <- 5
  a@stan_data$resid_df_fixed <- 5
  changes <- list(
    prior_cutoff_alpha = 987, prior_cutoff_beta = 987,
    prior_sigma_meanval_fringe = 987, prior_sigma_alpha = 987,
    prior_sigma_s0_pool = 987, prior_sigma_beta_s0 = 987,
    prior_sigma_resid_market = 987, prior_nu_shape = 987,
    prior_nu_rate = 987, prior_lkj = 987, bargaining_prior_alpha = 987
  )
  for (name in names(changes)) {
    b <- a
    b@metadata$resolved_priors[[name]] <- changes[[name]]
    expect_identical(.mu_fingerprint(a), .mu_fingerprint(b), info = name)
  }
  b <- a
  b@metadata$resid_df_fixed <- 9
  b@stan_data$resid_df_fixed <- 9
  b@metadata$resolved_priors$resid_df_fixed <- 9
  expect_identical(.mu_fingerprint(a), .mu_fingerprint(b))
  b <- a
  b@metadata$resolved_priors$prior_alpha_sd <- 2
  expect_false(identical(.mu_fingerprint(a), .mu_fingerprint(b)))
  # The same cutoff hyperparameters are active in a smooth fringe model.
  smooth <- a
  smooth@metadata$structure <- "fringe"
  smooth@metadata$partition_prior <- "cutoff"
  smooth@stan_data$use_cutoff <- 1L
  changed <- smooth
  changed@metadata$resolved_priors$prior_cutoff_alpha <- 987
  expect_false(identical(.mu_fingerprint(smooth), .mu_fingerprint(changed)))
})
