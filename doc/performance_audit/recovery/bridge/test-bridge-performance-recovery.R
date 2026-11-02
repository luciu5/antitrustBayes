test_that("bridge RNG helpers preserve seed presence and validate seed types", {
  before <- .bayes_bridge_rng_state()
  on.exit(.bayes_bridge_restore_rng(before), add = TRUE)

  value <- .bayes_bridge_preserve_rng({
    set.seed(731)
    stats::rnorm(3L)
  })
  expect_length(value, 3L)
  expect_identical(.bayes_bridge_rng_state(), before)

  expect_identical(.bayes_bridge_validate_seed(731), 731L)
  expect_null(.bayes_bridge_validate_seed(NULL))
  expect_error(.bayes_bridge_validate_seed("731"), "finite integer")
  expect_error(.bayes_bridge_validate_seed(7.5), "finite integer")
  expect_error(.bayes_bridge_validate_seed(NA_real_), "finite integer")
  expect_error(.bayes_bridge_validate_cores(1.5), "positive integer")
  expect_identical(.bayes_bridge_validate_cores(2), 2L)
})

test_that("explicit bridge streams are independent and deterministic", {
  before <- .bayes_bridge_rng_state()
  on.exit(.bayes_bridge_restore_rng(before), add = TRUE)

  first <- .bayes_bridge_preserve_rng(
    .bayes_bridge_allocate_streams(3L, bridge_seed = 1729L)
  )
  second <- .bayes_bridge_preserve_rng(
    .bayes_bridge_allocate_streams(3L, bridge_seed = 1729L)
  )
  expect_identical(first, second)
  expect_length(first$streams, 3L)
  expect_true(all(vapply(first$streams, is.integer, logical(1))))
  expect_identical(first$seed_source, "explicit")
  expect_identical(first$stream_kind[[1L]], "L'Ecuyer-CMRG")
})

test_that("bridge dispatch records worker errors and restores caller RNG", {
  contexts <- list(
    list(samples = list(), stanfit_model = NULL),
    list(samples = list(), stanfit_model = NULL)
  )
  before <- .bayes_bridge_rng_state()
  on.exit(.bayes_bridge_restore_rng(before), add = TRUE)

  expect_warning(
    dispatched <- .bayes_bridge_dispatch(
      contexts, rep(NA_character_, 2L), c("a", "b"),
      bridge_repetitions = 2L, bridge_maxiter = 17L,
      bridge_cores = 1L, bridge_seed = 1729L,
      bridge_fun = function(...) stop("synthetic bridge failure")
    ),
    "Bridge sampling failed"
  )
  expect_identical(.bayes_bridge_rng_state(), before)
  expect_equal(nrow(dispatched$failures), 4L)
  expect_identical(dispatched$failures$model,
                   rep(c("a", "b"), each = 2L))
  expect_identical(dispatched$failures$repetition, rep(1:2, 2L))
  expect_true(all(vapply(dispatched$bridge_list, function(x) {
    length(x) == 2L && all(vapply(x, is.null, logical(1)))
  }, logical(1))))
})

test_that("parallel explicit streams restore caller RNG on worker errors", {
  skip_on_os("windows")
  contexts <- list(
    list(samples = list(), stanfit_model = NULL),
    list(samples = list(), stanfit_model = NULL)
  )
  before <- .bayes_bridge_rng_state()
  on.exit(.bayes_bridge_restore_rng(before), add = TRUE)

  expect_warning(
    dispatched <- .bayes_bridge_dispatch(
      contexts, rep(NA_character_, 2L), c("a", "b"),
      bridge_repetitions = 2L, bridge_maxiter = 17L,
      bridge_cores = 2L, bridge_seed = 1729L,
      bridge_fun = function(...) stop("synthetic parallel bridge failure")
    ),
    "Bridge sampling failed"
  )
  expect_identical(.bayes_bridge_rng_state(), before)
  expect_identical(dispatched$execution$dispatch, "fork")
  expect_identical(dispatched$execution$requested_cores, 2L)
  expect_identical(dispatched$execution$effective_cores, 2L)
  expect_identical(dispatched$execution$seed_source, "explicit")
  expect_identical(dispatched$execution$stream_kind[[1L]], "L'Ecuyer-CMRG")
  expect_equal(nrow(dispatched$failures), 4L)
})

test_that("serial NULL-seed dispatch retains the legacy execution path", {
  contexts <- list(list(samples = list(), stanfit_model = NULL))
  dispatched <- .bayes_bridge_dispatch(
    contexts, NA_character_, "a", bridge_repetitions = 1L,
    bridge_maxiter = 17L, bridge_cores = 1L, bridge_seed = NULL,
    bridge_fun = function(...) list(logml = 1)
  )
  expect_identical(dispatched$execution$dispatch, "serial")
  expect_identical(dispatched$execution$requested_cores, 1L)
  expect_identical(dispatched$execution$effective_cores, 1L)
  expect_identical(dispatched$execution$seed_source, "none")
  expect_null(dispatched$execution$stream_kind)
  expect_true(is.finite(dispatched$bridge_list[[1L]][[1L]]$logml))
})

.make_bridge_recovery_fit <- function(label, fail = FALSE, prior = 1) {
  stan_data <- list(
    shareIn = c(.4, .6), clr_shareIn = c(-.2, .2),
    marginNormInv = c(2, 2), log_marginNormInv = log(c(2, 2)),
    margin_obs = c(1L, 1L), market_year_idx = c(1L, 1L),
    cluster = c(1L, 1L), likelihood_mode = 1L,
    use_s0_market_pool = 0L, use_s0_market_year_pool = 0L, use_hmt = 0L
  )
  attr(stan_data, "model_spec_version") <- 9L
  new("BayesLogit", stanfit = NULL,
      data = data.frame(
        prices = c(1, 2), quantities = c(4, 6), margins = c(.25, .2),
        market = c("m", "m"), year = c("1", "1"), group = c("a", "b"),
        loo_cluster = c("m", "m"), margin_observed = c(TRUE, TRUE),
        outside_price = c(1, 1)
      ),
      stan_data = stan_data, price_sd = 1, supply_model = 1L,
      demand_model = "logit", market_side = "output", call = call("fake"),
      metadata = list(
        prior_manifest = "baseline",
        resolved_priors = list(prior_alpha_sd = 1),
        prior_normalization = list(mode = "native_stan_normalized",
                                   model_spec_version = 9L),
        quality_parameterization = "collapsed",
        bayes_is_collapsed_quality = TRUE,
        model_prior = prior, bridge_id = label, bridge_fail = fail
      ))
}

.mock_bridge_dispatch_for_recovery <- function(contexts, context_errors,
                                                model_names,
                                                bridge_repetitions, ...) {
  testthat::expect_length(contexts, 3L)
  testthat::expect_length(context_errors, 3L)
  testthat::expect_length(model_names, 3L)
  bridge_list <- lapply(seq_along(contexts), function(i) {
    if (is.null(contexts[[i]])) return(rep(list(NULL), bridge_repetitions))
    lapply(seq_len(bridge_repetitions), function(rep_id) {
      list(logml = 10 + i + rep_id / 10)
    })
  })
  failed <- which(vapply(contexts, is.null, logical(1)))
  rows <- lapply(failed, function(i) data.frame(
    model_index = i, model = model_names[[i]],
    repetition = seq_len(bridge_repetitions), stage = "context",
    error = context_errors[[i]], stringsAsFactors = FALSE
  ))
  failures <- if (length(rows)) do.call(rbind, rows) else
    data.frame(model_index = integer(0), model = character(0),
               repetition = integer(0), stage = character(0),
               error = character(0), stringsAsFactors = FALSE)
  for (i in failed) {
    warning("Bridge context restoration failed for model ", model_names[[i]],
            "; evidence marked incomplete: ", context_errors[[i]],
            call. = FALSE)
  }
  list(
    bridge_list = bridge_list, failures = failures,
    execution = list(requested_cores = 1L, effective_cores = 1L,
                     dispatch = "serial", fallback = NULL,
                     seed_source = "none", realized_seed = NULL,
                     stream_kind = NULL, n_jobs = 3L * bridge_repetitions)
  )
}

.run_public_bridge_recovery <- function(failed, prior_weights = rep(1 / 3, 3L)) {
  fits <- lapply(seq_len(3L), function(i) {
    .make_bridge_recovery_fit(
      c("first", "middle", "last")[[i]],
      fail = identical(i, failed), prior = prior_weights[[i]]
    )
  })
  names(fits) <- c("first", "middle", "last")
  testthat::local_mocked_bindings(
    .bayes_bridge_restore_context = function(fit) {
      if (isTRUE(fit@metadata$bridge_fail)) {
        stop("synthetic restoration failure for ", fit@metadata$bridge_id)
      }
      list(samples = list(), stanfit_model = NULL,
           source = "test", kernel = "logit_supply",
           saved_hash = "saved", current_hash = "current",
           migration = FALSE, zero_width_fields_added = character(0L),
           target_check = TRUE)
    },
    .bayes_bridge_dispatch = .mock_bridge_dispatch_for_recovery,
    .package = "antitrustBayes"
  )
  result <- NULL
  expect_warning(
    result <- do.call(bayes_compare, c(fits, list(
      bridge = TRUE, bridge_repetitions = 2L, bridge_cores = 1L
    ))),
    "Bridge context restoration failed"
  )
  result
}

test_that("public bridge comparison retains failed model slots and repetitions", {
  for (failed in c(1L, 2L, 3L)) {
    comp <- .run_public_bridge_recovery(failed)
    expect_identical(comp$model, c("first", "middle", "last"))
    expect_equal(comp$logml[failed], NA_real_)
    expect_true(all(is.finite(comp$logml[-failed])))
    bridge_list <- attr(comp, "bridge_list")
    expect_length(bridge_list, 3L)
    expect_true(all(vapply(bridge_list, length, integer(1)) == 2L))
    expect_true(all(vapply(bridge_list[[failed]], is.null, logical(1))))
    expect_true(all(vapply(bridge_list[-failed], function(reps) {
      all(vapply(reps, function(x) is.list(x) && is.finite(x$logml), logical(1)))
    }, logical(1))))
    failures <- attr(comp, "bridge_failure_records")
    expect_equal(failures$model, rep(c("first", "middle", "last")[[failed]], 2L))
    expect_identical(failures$repetition, 1:2)
    expect_true(all(is.na(comp$post_prob)))
  }
})

test_that("zero-prior missing bridge evidence remains a zero posterior", {
  comp <- .run_public_bridge_recovery(2L, prior_weights = c(.5, 0, .5))
  expect_true(all(is.finite(comp$post_prob[c(1L, 3L)])))
  expect_identical(comp$post_prob[[2L]], 0)
  expect_equal(sum(comp$post_prob), 1)
})

test_that("public universe comparison forwards bridge execution controls", {
  fit <- .make_bridge_recovery_fit("source")
  universe <- bayes_model_universe(list(H = "source"), c(H = 1))
  testthat::local_mocked_bindings(
    .bayes_compare_universe = function(args, universe, bridge, bridge_maxiter,
                                        bridge_repetitions, bridge_cores,
                                        bridge_seed, loo_unit) {
      list(bridge_repetitions = bridge_repetitions,
           bridge_cores = bridge_cores, bridge_seed = bridge_seed,
           loo_unit = loo_unit)
    },
    .package = "antitrustBayes"
  )
  got <- bayes_compare(
    source = fit, universe = universe, bridge_repetitions = 4L,
    bridge_cores = 2L, bridge_seed = 23L
  )
  expect_identical(got$bridge_repetitions, 4L)
  expect_identical(got$bridge_cores, 2L)
  expect_identical(got$bridge_seed, 23L)
  expect_identical(got$loo_unit, "cluster")
})
