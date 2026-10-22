.backend_fixture <- function(columns, iterations = 4L, chains = 2L) {
  draws <- array(seq_len(iterations * chains * length(columns)),
                 c(iterations, chains, length(columns)),
                 dimnames = list(NULL, NULL, columns))
  structure(list(
    draws_array = draws,
    diagnostics_array = array(0, c(iterations, chains, 1L),
      dimnames = list(NULL, NULL, "divergent__")),
    metadata = list(stan_variables = unique(sub("\\[.*$", "", columns)),
                    provenance = list(source_hash = "frozen-source"))
  ), class = "NutpieFit")
}

test_that("backend extraction retains chain order and Stan matrix indices", {
  fit <- .backend_fixture(c("x[2,1]", "x[1,2]", "theta", "x[1,1]", "x[2,2]"))
  got <- .bayes_backend_extract(fit, c("theta", "x"))
  expect_identical(names(got), c("theta", "x"))
  expect_equal(got$theta, as.numeric(17:24))
  expect_equal(dim(got$x), c(8L, 2L, 2L))
  expect_equal(as.numeric(got$x[, 2L, 1L]), as.numeric(1:8))
  expect_equal(as.numeric(got$x[, 1L, 2L]), as.numeric(9:16))
  expect_equal(as.numeric(got$x[, 1L, 1L]), as.numeric(25:32))
  expect_equal(as.numeric(got$x[, 2L, 2L]), as.numeric(33:40))
  expect_equal(dim(.bayes_backend_draws(fit, "x", "draws_array")),
               c(4L, 2L, 4L))
  expect_identical(.bayes_backend_name(fit), "nutpie")
  with_lp <- .backend_fixture(c("theta", "lp__", "x[1]", "x[2]"))
  expect_identical(.bayes_backend_variable_names(with_lp), c("theta", "x"))
  expect_identical(names(.bayes_backend_extract(with_lp)), c("theta", "x"))
  expect_equal(.bayes_backend_extract(with_lp, "lp__")$lp__, as.numeric(9:16))
})

test_that("indexed extraction rejects missing or malformed Stan coordinates", {
  expect_error(.bayes_backend_extract(.backend_fixture(c("x[1,1]", "x[2,2]")), "x"),
               "Incomplete indexed")
  expect_error(.bayes_backend_extract(.backend_fixture(c("x[0]", "x[1]")), "x"),
               "Invalid indexed")
  expect_error(.bayes_backend_extract(.backend_fixture(c("x[1]", "x[1,2]")), "x"),
               "Invalid indexed")
  expect_error(.bayes_backend_extract(.backend_fixture("theta"), "log_lik"),
               "missing posterior field")
  expect_error(.bayes_extract_indexed_draws(
    matrix(1:4, 2L, dimnames = list(NULL, c("x[1]", "x[1]"))), "x"),
    "unique Stan names")
})

test_that("backend extraction reconstructs vectors and higher Stan dimensions", {
  fields <- c("v[2]", "a[1,2,2]", "a[1,1,1]", "v[1]",
              "a[1,2,1]", "a[1,1,2]")
  fit <- .backend_fixture(fields, iterations = 2L, chains = 1L)
  got <- .bayes_backend_extract(fit, c("v", "a"))
  expect_equal(dim(got$v), c(2L, 2L))
  expect_equal(got$v, cbind(as.numeric(7:8), as.numeric(1:2)))
  expect_equal(dim(got$a), c(2L, 1L, 2L, 2L))
  expect_equal(as.numeric(got$a[, 1L, 2L, 2L]), as.numeric(3:4))
  expect_equal(as.numeric(got$a[, 1L, 1L, 1L]), as.numeric(5:6))
  expect_equal(as.numeric(got$a[, 1L, 2L, 1L]), as.numeric(9:10))
  expect_equal(as.numeric(got$a[, 1L, 1L, 2L]), as.numeric(11:12))
})

test_that("saved Nutpie fits summarize and expose diagnostics without native objects", {
  fit <- .backend_fixture(c("theta", "sigma"), iterations = 16L)
  fit <- unserialize(serialize(fit, NULL))
  summary <- suppressWarnings(.bayes_backend_summary(fit, "theta"))
  expect_identical(as.character(summary$variable), "theta")
  expect_equal(summary$mean, mean(1:32))
  expect_true(all(c("rhat", "ess_bulk", "q5", "q95") %in% names(summary)))
  expect_equal(dim(.bayes_backend_diagnostics(fit)), c(16L, 2L, 1L))
  expect_equal(.bayes_backend_provenance(fit)$provenance$source_hash, "frozen-source")
  fit$diagnostics_array <- NULL
  expect_error(.bayes_backend_diagnostics(fit), "no retained sampler diagnostics")
})

test_that("common diagnostics report native depth hits and unavailable metrics", {
  fit <- .backend_fixture("theta", iterations = 4L)
  fit$metadata$num_warmup <- 20L
  fit$metadata$resolved_options <- list(adaptation = "low_rank", max_treedepth = 13L)
  fit$metadata$sampler_config <- '{"maxdepth":13}'
  fit$metadata$sampling_time_native <- 2.5
  fit$diagnostics_array <- array(0, c(4, 2, 5), dimnames = list(NULL, NULL,
    c("divergent__", "treedepth__", "maxdepth_reached__", "energy__", "stepsize__")))
  fit$diagnostics_array[, , "divergent__"] <- c(1, 0, 0, 0, 0, 0, 1, 1)
  fit$diagnostics_array[, , "treedepth__"] <- 12
  fit$diagnostics_array[, , "maxdepth_reached__"] <- c(0, 1, 0, 0, 1, 0, 1, 0)
  fit$diagnostics_array[, , "energy__"] <- c(1, 2, 4, 8, 2, 1, 3, 2)
  fit$diagnostics_array[, , "stepsize__"] <- rep(c(.1, .2), each = 4L)
  result <- suppressWarnings(.bayes_backend_diagnostic_summary(fit))
  expect_equal(result$chains$divergences, c(1, 2))
  expect_equal(result$chains$depth_hits, c(1, 2))
  expect_equal(result$chains$bfmi, c(63 / 115, 9 / 4))
  expect_equal(result$chains$stepsize_median, c(.1, .2))
  expect_length(result$unavailable, 0L)
  expect_identical(result$adaptation, "low_rank")
  expect_true(all(c("rhat", "ess_bulk", "ess_tail") %in% names(result$posterior)))
  expect_identical(result$sampling$warmup_iterations, 20L)
  expect_identical(result$sampling$native_config, fit$metadata$sampler_config)
  expect_equal(result$sampling$timing$total_seconds, 2.5)
  expect_true(is.na(result$sampling$timing$warmup_seconds))
  expect_true(is.na(result$sampling$timing$sampling_seconds))
  expect_match(result$sampling$timing$unavailable$phase_times, "separate warmup")
  fit$diagnostics_array <- fit$diagnostics_array[, , "divergent__", drop = FALSE]
  missing <- suppressWarnings(.bayes_backend_diagnostic_summary(fit))
  expect_true(all(is.na(missing$chains$bfmi)))
  expect_match(missing$unavailable$bfmi, "Hamiltonian energy")
  expect_match(missing$unavailable$depth_hits, "Depth-hit flags")
  expect_match(missing$unavailable$stepsize_median, "step sizes")
})

test_that("legacy CmdStan objects retain their extraction and summary paths", {
  columns <- c("theta", "x[2]", "x[1]")
  native <- .backend_fixture(columns)
  cmdstan <- structure(list(
    metadata = function() list(stan_variables = c("theta", "x"), iter_warmup = 100L),
    draws = function(variables = NULL, format = "draws_array")
      .bayes_backend_draws(native, variables, format),
    summary = function(variables = NULL) data.frame(variable = variables, mean = 42),
    time = function() list(total = 4, chains = data.frame(chain_id = 1:2,
      warmup = c(2, 1.5), sampling = c(1, .9)))
  ), class = "CmdStanFit")
  got <- .bayes_backend_extract(cmdstan, c("theta", "x"))
  expect_equal(got$theta, as.numeric(1:8))
  expect_equal(got$x[, 1L], as.numeric(17:24))
  expect_equal(got$x[, 2L], as.numeric(9:16))
  expect_identical(.bayes_backend_name(cmdstan), "cmdstanr")
  expect_equal(.bayes_backend_summary(cmdstan, "theta")$mean, 42)
  info <- .bayes_backend_sampling_info(cmdstan)
  expect_identical(info$backend, "cmdstanr")
  expect_identical(info$chains, 2L)
  expect_identical(info$retained_draws_per_chain, 4L)
  expect_identical(info$warmup_iterations, 100L)
  expect_equal(info$timing$total_seconds, 4)
  expect_equal(info$timing$warmup_seconds, c(2, 1.5))
  expect_equal(info$timing$sampling_seconds, c(1, .9))
  expect_true(is.na(info$timing$compile_seconds))
})

test_that("collapsed quality can reconstruct directly from a retained Nutpie fit", {
  columns <- c("share_mean[3]", "sigma_share", "share_mean[1]",
               "sigma_b_strat", "share_mean[2]", "share_mean[4]")
  native <- .backend_fixture(columns)
  native$draws_array[] <- 0
  native$draws_array[, , "sigma_share"] <- .4
  native$draws_array[, , "sigma_b_strat"] <- .6
  stan_data <- list(N = 4L, N_group = 2L, group = c(1L, 2L, 1L, 2L),
                    log_shareIn = c(.2, -.1, .3, -.2), likelihood_mode = 2L)
  attr(stan_data, "group_levels") <- c("a", "b")
  fit <- methods::new("BayesLogit", stanfit = native, stan_data = stan_data,
    data = data.frame(group = c("a", "b", "a", "b")), price_sd = 1,
    supply_model = 1L, demand_model = "logit", market_side = "output",
    metadata = list(quality_parameterization = "collapsed", share_likelihood = "gaussian"))
  got <- posterior_quality(fit, seed = 17L)
  ## With centered b=(q,-q), posterior precision for q is
  ## 2/tau^2 + N/sigma^2; its sufficient statistic is sum(z*y)/sigma^2.
  precision <- 2 / .6^2 + 4 / .4^2
  expected <- sum(c(1, -1, 1, -1) * stan_data$log_shareIn) / .4^2 / precision
  expect_equal(as.numeric(got$conditional_means[, 1L]), rep(expected, 8L))
  expect_equal(as.numeric(got$conditional_means[, 2L]), rep(-expected, 8L))
  expect_equal(as.numeric(got$conditional_sds[, 1L]), rep(sqrt(1 / precision), 8L))
  expect_equal(rowSums(got$draws), rep(0, 8L), ignore_attr = TRUE)
})

test_that("saved Nutpie policy plans retain posterior rows and reject changed provenance", {
  columns <- c("share_mean[3]", "b_quality_all[2]", "a_market[1]",
               "s0_prob[1]", "b_quality_all[1]", "share_mean[1]",
               "b_quality_all[3]", "share_mean[2]")
  native <- .backend_fixture(columns, iterations = 2L, chains = 1L)
  quality <- rbind(c(.2, -.1, -.1), c(-.15, .05, .1))
  native$draws_array[, , "a_market[1]"] <- .8
  native$draws_array[, , "s0_prob[1]"] <- .2
  for (j in 1:3) {
    native$draws_array[, , paste0("b_quality_all[", j, "]")] <- quality[, j]
    native$draws_array[, , paste0("share_mean[", j, "]")] <-
      log(c(40, 30, 20)[[j]] / 90) + quality[, j]
  }
  stan_data <- list(quality_parameterization = 1L, N_group = 3L, group = 1:3)
  attr(stan_data, "model_spec_version") <- 12L
  attr(stan_data, "group_levels") <- c("a", "b", "c")
  attr(stan_data, "quality_scope") <- "global"
  fit <- methods::new("BayesLogit", stanfit = native, stan_data = stan_data,
    data = data.frame(prices = c(1.1, 1.2, 1.3), quantities = c(40, 30, 20),
      margins = c(.5, .4, .3), group = c("a", "b", "c"), market = "m1",
      year = "1", outside_price = 1), price_sd = 1, supply_model = 1L,
    demand_model = "logit", market_side = "input",
    metadata = list(backend = "nutpie", structure = "simultaneous",
                    model_spec_version = 12L,
                    quality_parameterization = "sampled"))
  plan <- posterior_policy_draws(fit, n_draws = 2L, seed = 19L)
  expect_equal(unname(plan$quality_draws), quality[plan$draw_indices, , drop = FALSE])
  saved <- unserialize(serialize(list(fit = fit, plan = plan), NULL))
  expect_identical(.lq_validate_plan(saved$fit, saved$plan), saved$plan)
  simulated <- suppressWarnings(simulate(saved$fit, c("a", "b", "c"),
    draw_plan = saved$plan, mc.cores = 1L, demand_anchoring = "latent_quality"))
  expect_equal(as.integer(simulated@n_valid), 2L)
  expect_identical(simulated@metadata$draw_indices, saved$plan$draw_indices)
  saved$fit@metadata$backend <- "cmdstanr"
  expect_error(.lq_validate_plan(saved$fit, saved$plan), "different fit")
})
