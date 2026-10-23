.nutpie_test_options <- function() list(chains = 2L, cores = 1L, iter = 30L,
  warmup = 10L, refresh = 0L, seed = 12L,
  control = list(adapt_delta = .95, max_treedepth = 13L))

test_that("Nutpie preserves explicit common controls and separates adaptation", {
  options <- .nutpie_test_options()
  options$init <- .25
  options$save_warmup <- TRUE
  translated <- .nutpie_sample_options(options)
  expect_identical(translated$num_draws, 20L)
  expect_identical(translated$num_warmup, 10L)
  expect_identical(translated$num_chains, 2L)
  expect_identical(translated$cores, 1L)
  expect_identical(translated$target_accept, .95)
  expect_identical(translated$max_treedepth, 13L)
  expect_identical(translated$seed, 12L)
  expect_identical(translated$init, .25)
  expect_true(translated$save_warmup)
  expect_identical(translated$adaptation, "diag")
  low_rank <- .nutpie_sample_options(options,
    list(adaptation = "low_rank", mass_matrix_gamma = 1e-4))
  expect_identical(low_rank$adaptation, "low_rank")
  expect_identical(low_rank$mass_matrix_gamma, 1e-4)
  low_rank$adaptation <- "diag"
  low_rank$mass_matrix_gamma <- NULL
  expect_identical(low_rank, translated)
  options$seed <- NULL
  set.seed(97)
  rng_before <- .Random.seed
  expect_false("seed" %in% names(.nutpie_sample_options(options)))
  expect_identical(.Random.seed, rng_before)
})

test_that("Nutpie rejects unsupported and ignored sampler requests", {
  options <- .nutpie_test_options()
  for (control in c("metric", "stepsize", "adapt_engaged")) {
    unsupported <- options
    unsupported$control[[control]] <- 1
    expect_error(.nutpie_sample_options(unsupported), paste0("Unsupported Nutpie control.*", control))
  }
  for (name in c("sig_figs", "output_dir", "surprise")) {
    unsupported <- options
    unsupported[[name]] <- 1
    expect_error(.nutpie_sample_options(unsupported), paste0("Unsupported stan_options.*", name))
  }
  expect_error(.nutpie_sample_options(options, list(target_accept = .99)), "Unsupported nutpie_options")
  expect_error(.nutpie_sample_options(options, list(mass_matrix_gamma = 1e-3)), "require adaptation")
  for (name in c("thin", "threads_per_chain")) {
    unsupported <- options
    unsupported[[name]] <- 2L
    expect_error(.nutpie_sample_options(unsupported), paste0(name, " = 1"), fixed = TRUE)
  }
  options$warmup <- 0L
  expect_error(.nutpie_sample_options(options), "warmup.*integer >= 1")
})

test_that("Nutpie validates integer and probability controls before native work", {
  for (name in c("iter", "warmup", "chains", "cores")) {
    for (bad in list(NA, Inf, 1.1, numeric(), c(1, 2))) {
      options <- .nutpie_test_options()
      options[[name]] <- bad
      expect_error(.nutpie_sample_options(options), "finite integer")
    }
  }
  options <- .nutpie_test_options()
  options$control$adapt_delta <- 1
  expect_error(.nutpie_sample_options(options), "in \\(0, 1\\)")
  expect_error(.nutpie_sample_options(.nutpie_test_options(), list(store_gradient = NA)), "TRUE or FALSE")
  expect_error(.nutpie_sample_options(.nutpie_test_options(), list(progress = "maybe")), "progress must")
  expect_error(.nutpie_sample_options(.nutpie_test_options(), list(mindepth = 14L)),
    "mindepth cannot exceed")
  unnamed <- stats::setNames(list(1), NA_character_)
  expect_error(.nutpie_sample_options(unnamed), "uniquely named list")
  expect_error(.nutpie_sample_options(.nutpie_test_options(), unnamed), "uniquely named list")
})

test_that("Nutpie diagnostic arrays use explicit chain and draw identities", {
  raw <- list(chain = c(2L, 1L, 2L, 1L), draw = c(2L, 1L, 1L, 2L),
    logp = c(-22, -11, -21, -12), energy = c(22, 11, 21, 12),
    diverging = c(TRUE, FALSE, FALSE, FALSE), depth = c(4L, 1L, 3L, 2L),
    n_steps = c(15L, 1L, 7L, 3L))
  result <- .nutpie_diagnostics_array(raw, 2L, 2L)
  expect_equal(unname(result[, , "lp__"]), matrix(c(-11, -12, -21, -22), 2L))
  expect_equal(unname(result[, , "divergent__"]), matrix(c(0, 0, 0, 1), 2L))
  expect_equal(unname(result[, , "treedepth__"]), matrix(1:4, 2L))
  expect_equal(unname(result[, , "n_leapfrog__"]), matrix(c(1, 3, 7, 15), 2L))
  expect_false("accept_stat__" %in% dimnames(result)[[3L]])
  raw$draw[[1L]] <- 1L
  expect_error(.nutpie_diagnostics_array(raw, 2L, 2L), "duplicate diagnostic draws")
  raw$chain[[1L]] <- 0L
  expect_error(.nutpie_diagnostics_array(raw, 2L, 2L), "invalid diagnostic chain/draw")
  raw$chain[[1L]] <- Inf
  expect_error(.nutpie_diagnostics_array(raw, 2L, 2L), "invalid diagnostic chain/draw")
  raw$chain <- as.character(raw$chain)
  expect_error(.nutpie_diagnostics_array(raw, 2L, 2L), "invalid diagnostic chain/draw")
})

test_that("Nutpie completion requires all requested usable posterior draws", {
  draws <- array(c(rep(-1, 4), 1:4), c(2, 2, 2),
    dimnames = list(iteration = 1:2, chain = 1:2, variable = c("lp__", "theta")))
  fit <- structure(list(draws_array = draws,
    metadata = list(num_samples = 2L, num_chains = 2L)), class = "NutpieFit")
  expect_invisible(.nutpie_assert_fit(fit, "hmc", 2L))
  expect_invisible(.nutpie_assert_fit(fit, "hmc"))
  expect_error(.nutpie_assert_fit(fit, "hmc", 3L), "every requested chain")
  expect_error(.nutpie_assert_fit(fit, "pathfinder", 2L), "HMC fits only")
  fit$draws_array[1, 1, "lp__"] <- -Inf
  expect_error(.nutpie_assert_fit(fit, "hmc", 2L), "usable posterior log densities")
  fit$draws_array[1, 1, "lp__"] <- -1
  fit$draws_array[1, 1, "theta"] <- NaN
  expect_error(.nutpie_assert_fit(fit, "hmc", 2L), "generated-quantity draws")
  fit$draws_array[1, 1, "theta"] <- -Inf
  expect_error(.nutpie_assert_fit(fit, "hmc", 2L), "nonfinite posterior field: theta")
  fit$draws_array[1, 1, "theta"] <- Inf
  expect_error(.nutpie_assert_fit(fit, "hmc", 2L), "nonfinite posterior field: theta")
  fit$draws_array[1, 1, "theta"] <- 1
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(fit, path)
  expect_identical(readRDS(path), fit)
  dimnames(fit$draws_array)[[3L]] <- c("lp__", "lp__")
  expect_error(.nutpie_assert_fit(fit, "hmc"), "duplicate posterior field names")
  dimnames(fit$draws_array) <- NULL
  expect_error(.nutpie_assert_fit(fit, "hmc"), "missing or duplicate posterior field names")
})

test_that("only documented zero-prior partition scores may be negative infinite", {
  draws <- array(c(rep(-1, 4), rep(-Inf, 4), rep(0, 4)), c(2, 2, 3),
    dimnames = list(iteration = 1:2, chain = 1:2,
      variable = c("lp__", "component_log_score[1]", "prior_component_prob[1]")))
  fit <- structure(list(draws_array = draws,
    metadata = list(kernel = "logit_stackelberg", num_samples = 2L)), class = "NutpieFit")
  expect_invisible(.nutpie_assert_fit(fit, "hmc", 2L))
  fit$draws_array[1, 1, "prior_component_prob[1]"] <- .1
  expect_error(.nutpie_assert_fit(fit, "hmc", 2L), "nonfinite posterior field")
  fit$draws_array[1, 1, "prior_component_prob[1]"] <- 0
  fit$draws_array[1, 1, "component_log_score[1]"] <- Inf
  expect_error(.nutpie_assert_fit(fit, "hmc", 2L), "nonfinite posterior field")
  fit$draws_array[1, 1, "component_log_score[1]"] <- -Inf
  fit$metadata$kernel <- "logit_supply"
  expect_error(.nutpie_assert_fit(fit, "hmc", 2L), "nonfinite posterior field")
})

test_that("Nutpie enumerates each frozen kernel data universe without compilation", {
  for (kernel in c("logit_supply", "ces_supply", "logit_stackelberg", "ces_stackelberg")) {
    fields <- .nutpie_stan_data_fields(kernel)
    expect_true(all(c("N", "share_basis") %in% fields))
    expect_true(all(c("price", "market", "N_market", "prior_lkj") %in% fields))
    expect_identical(anyDuplicated(fields), 0L)
    expect_false(any(c("b_strat", "log_lik", "lp__") %in% fields))
  }
  expect_true("revShare" %in% .nutpie_stan_data_fields("ces_supply"))
  expect_false("revShare" %in% .nutpie_stan_data_fields("logit_supply"))
})

test_that("native Nutpie cache follows the option and restores the environment", {
  cache <- tempfile("nutpie-cache-")
  dir.create(cache)
  on.exit(unlink(cache, recursive = TRUE), add = TRUE)
  old_options <- options(antitrustBayes.nutpie_cache = cache)
  on.exit(options(old_options), add = TRUE)
  before <- Sys.getenv("R_USER_CACHE_DIR", unset = NA_character_)
  during <- .nutpie_with_cache(Sys.getenv("R_USER_CACHE_DIR"))
  expect_identical(during, normalizePath(file.path(cache, "nutpie-native")))
  expect_identical(Sys.getenv("R_USER_CACHE_DIR", unset = NA_character_), before)
  expect_error(.nutpie_with_cache(stop("cache error")), "cache error")
  expect_identical(Sys.getenv("R_USER_CACHE_DIR", unset = NA_character_), before)
})

test_that("native Nutpie draws retain generated quantities, log density, and provenance", {
  cache <- tempfile("nutpie-native-fixture-")
  dir.create(cache)
  on.exit(unlink(cache, recursive = TRUE), add = TRUE)
  old_options <- options(antitrustBayes.nutpie_cache = cache)
  on.exit(options(old_options), add = TRUE)
  captured <- NULL
  native_draws <- posterior::as_draws_array(array(1:16, c(2, 2, 4),
    dimnames = list(iteration = 1:2, chain = 1:2,
      variable = c("theta", "b[1]", "log_lik[1]", "replicate[1]"))))
  attr(native_draws, "sampler_config") <- '{"num_warmup":1}'
  attr(native_draws, "sampling_time") <- .25
  diagnostics <- list(chain = c(1, 1, 2, 2), draw = c(1, 2, 1, 2),
    logp = c(-1, -2, -3, -4), diverging = rep(FALSE, 4))
  local_mocked_bindings(.nutpie_require = function() invisible(TRUE),
    .nutpie_sample_native = function(...) {
      captured <<- list(...)
      native_draws
    }, .nutpie_diagnostics_native = function(draws) diagnostics,
    .package = "antitrustBayes")
  model <- structure(list(native_model = list(lib_path = "unused"),
    metadata = list(nutpie_version = "1.9.1", source_hash = "source")), class = "NutpieModel")
  options <- .nutpie_test_options()
  options$iter <- 3L
  options$warmup <- 1L
  data <- list(N = 1L, y = array(2, 1L))
  fit <- .nutpie_sample(model, data, options)
  expect_s3_class(fit, "NutpieFit")
  expect_equal(dim(fit$draws_array), c(2, 2, 5))
  expect_equal(as.numeric(fit$draws_array[, , "lp__"]), diagnostics$logp)
  expect_equal(as.numeric(fit$draws_array[, , "log_lik[1]"]), 9:12)
  expect_equal(as.numeric(fit$draws_array[, , "replicate[1]"]), 13:16)
  expect_identical(fit$metadata$stan_data_hash, digest::digest(data, algo = "sha256"))
  expect_identical(fit$metadata$data_json, captured$data)
  expect_identical(fit$metadata$sampler_config, attr(native_draws, "sampler_config"))
  expect_identical(fit$metadata$sampling_time_native, .25)
  expect_identical(captured$seed, 12L)
  expect_identical(captured$adaptation, "diag")
  expect_match(captured$data, '"y": \\[2\\]')
  path <- file.path(cache, "fit.rds")
  saveRDS(fit, path)
  expect_identical(readRDS(path), fit)
  native_draws <- fit$draws_array
  expect_identical(.nutpie_sample(model, data, options)$draws_array, fit$draws_array)
  native_draws[1L, 1L, "lp__"] <- native_draws[1L, 1L, "lp__"] - 1
  expect_error(.nutpie_sample(model, data, options), "log densities disagree")
  native_draws <- native_draws[1L, , , drop = FALSE]
  expect_error(.nutpie_sample(model, data, options), "every requested chain")
})

test_that("Nutpie compiles exactly the shared kernel in an external cache", {
  skip_if_not_installed("nutpieR", "1.9.1")
  cache <- tempfile("nutpie-source-fixture-")
  dir.create(cache)
  on.exit(unlink(cache, recursive = TRUE), add = TRUE)
  old_options <- options(antitrustBayes.nutpie_cache = cache)
  on.exit(options(old_options), add = TRUE)
  source <- .cmdstan_kernel_file("logit_supply")
  original_hash <- digest::digest(file = source, algo = "sha256")
  captured <- NULL
  local_mocked_bindings(.nutpie_compile_native = function(...) {
    captured <<- list(...)
    lib <- file.path(cache, "fixture.so")
    writeBin(as.raw(1:10), lib)
    list(lib_path = lib)
  }, .package = "antitrustBayes")
  model <- .nutpie_model_for_kernel("logit_supply")
  expect_identical(model$metadata$source_hash, original_hash)
  expect_identical(digest::digest(file = captured$stan_file, algo = "sha256"), original_hash)
  expect_true(startsWith(captured$stan_file, normalizePath(cache)))
  expect_false(identical(captured$stan_file, source))
  expect_identical(captured$compile_args, "STAN_THREADS=true")
  expect_identical(model$metadata$compile_args, captured$compile_args)
  expect_identical(digest::digest(file = source, algo = "sha256"), original_hash)
})

test_that("public backend defaults and explicit CmdStanR requests agree", {
  captured <- list()
  local_mocked_bindings(.cmdstan_model_for_kernel = function(kernel) {
    captured[[length(captured) + 1L]] <<- kernel
    stop("CmdStan route selected")
  }, .nutpie_require = function() stop("Nutpie dependency unexpectedly requested"),
  .package = "antitrustBayes")
  args <- list(prices = c(1, 2, 3), quantities = c(3, 2, 1), margins = rep(.25, 3),
    market = rep(1, 3), group = 1:3, market_side = "output")
  expect_error(do.call(bayes_est, args), "CmdStan route selected")
  expect_error(do.call(bayes_est, c(args, list(backend = "cmdstanr"))), "CmdStan route selected")
  expect_identical(captured[[1L]], captured[[2L]])
})

test_that("public Nutpie capability and dependency errors occur before sampling", {
  expect_error(bayes_est(backend = "nutpie", inference = "pathfinder"), "HMC|hmc")
  expect_error(bayes_est(nutpie_options = list(adaptation = "low_rank")), "require backend")
  local_mocked_bindings(.nutpie_require = function() stop("optional nutpieR missing"),
    .package = "antitrustBayes")
  expect_error(bayes_est(backend = "nutpie"), "optional nutpieR missing")
})
