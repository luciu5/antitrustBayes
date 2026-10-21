.nutpie_portability_fixture <- function(directory) {
  kernel <- "logit_supply"
  source <- .cmdstan_kernel_file(kernel)
  data <- list(N = 2L)
  data_hash <- digest::digest(data, algo = "sha256")
  json_file <- file.path(directory, "data.json")
  cmdstanr::write_stan_json(data, json_file)
  json <- paste(readLines(json_file, warn = FALSE), collapse = "\n")
  q <- seq(-1, 1, length.out = 8L)
  draws <- array(c(exp(q), rep(999, 8L), -exp(q) + q), c(4L, 2L, 3L),
    dimnames = list(NULL, NULL, c("scale", "unrelated_generated", "lp__")))
  library <- file.path(directory, "saved.so")
  writeBin(charToRaw("unit-test artifact; never load as native code"), library)
  hash <- digest::digest(file = source, algo = "sha256")
  native <- structure(list(draws_array = draws, metadata = list(
    source_code = paste(readLines(source, warn = FALSE), collapse = "\n"),
    source_hash = hash, inference = "hmc", kernel = kernel,
    num_chains = 2L, num_samples = 4L, lib_path = library,
    lib_hash = digest::digest(file = library, algo = "sha256"),
    stan_data_hash = data_hash, data_json = json,
    data_json_hash = digest::digest(json, algo = "sha256"), seed = 119L,
    platform = R.version$platform, nutpie_version = "1.9.1",
    compile_args = "STAN_THREADS=true", stanc_args = character())), class = "NutpieFit")
  methods::new("BayesLogit", stanfit = native, stan_data = data, metadata = list(
    backend = "nutpie", inference = "hmc", stan_kernel = kernel,
    stan_source_hash = hash, stan_data_fields = "N", stan_data_hash = data_hash))
}

.nutpie_portability_model <- function(version = "2.8.0", offset = 0) {
  pieces <- as.integer(strsplit(version, ".", fixed = TRUE)[[1L]])
  list(model_version = function() list(major = pieces[1L], minor = pieces[2L], patch = pieces[3L]),
    model_info = function() "unit-test evaluator; no native code loaded",
    param_names = function(include_tp, include_gq) "scale",
    param_unc_names = function() "scale", param_unc_num = function() 1L,
    param_unconstrain = log,
    param_constrain = function(q, include_tp, include_gq) exp(q),
    log_density = function(q, propto, jacobian) {
      stopifnot(isTRUE(jacobian))
      -exp(q) + q + offset
    })
}

.nutpie_portability_artifact <- function(fit, directory) {
  library <- file.path(directory, "restored.so")
  writeBin(charToRaw("different verified unit-test artifact"), library)
  list(metadata = list(lib_path = library,
    lib_hash = digest::digest(file = library, algo = "sha256"),
    source_hash = fit@metadata$stan_source_hash, platform = R.version$platform,
    nutpie_version = "1.9.1", compile_args = "STAN_THREADS=true",
    stanc_args = character(), compile_time = .1))
}

test_that("verified saved Nutpie libraries restore without compilation or fit mutation", {
  directory <- withr::local_tempdir()
  withr::local_options(antitrustBayes.nutpie_cache = directory)
  fit <- .nutpie_portability_fixture(directory)
  before <- digest::digest(fit, algo = "sha256")
  local_mocked_bindings(
    .bayes_nutpie_bridge_client_version = function() "2.8.0",
    .bayes_nutpie_bridge_model = function(lib, data_file, seed) {
      expect_identical(basename(lib), paste0("antitrustBayes_", fit@stanfit$metadata$lib_hash, ".so"))
      expect_false(identical(lib, fit@stanfit$metadata$lib_path))
      expect_true(file.exists(data_file))
      expect_identical(seed, 119L)
      .nutpie_portability_model()
    },
    .bayes_nutpie_bridge_rebuild = function(...) stop("Unexpected compilation"))
  context <- .bayes_nutpie_bridge_context(fit)
  expect_false(context$migration)
  expect_identical(context$restoration$method, "saved_library")
  expect_identical(context$restoration$restored_artifact$hash, fit@stanfit$metadata$lib_hash)
  expect_identical(context$restoration$restored_artifact$bridgestan_version, "2.8.0")
  expect_identical(context$sampler_density_verification$points, 8L)
  expect_equal(as.numeric(context$samples), seq(-1, 1, length.out = 8L))
  expect_identical(digest::digest(fit, algo = "sha256"), before)
})

test_that("missing and changed Nutpie libraries reconstruct without loading rejected bytes", {
  directory <- withr::local_tempdir()
  withr::local_options(antitrustBayes.nutpie_cache = directory)
  fit <- .nutpie_portability_fixture(directory)
  artifact <- .nutpie_portability_artifact(fit, directory)
  calls <- list()
  local_mocked_bindings(
    .bayes_nutpie_bridge_client_version = function() "2.8.0",
    .bayes_nutpie_bridge_model = function(lib, data_file, seed) {
      expect_identical(basename(lib), paste0("antitrustBayes_", artifact$metadata$lib_hash, ".so"))
      .nutpie_portability_model()
    },
    .bayes_nutpie_bridge_rebuild = function(kernel, fresh) {
      calls[[length(calls) + 1L]] <<- list(kernel = kernel, fresh = fresh)
      artifact
    })
  missing <- fit
  missing@stanfit$metadata$lib_path <- file.path(directory, "absent.so")
  missing_before <- digest::digest(missing, algo = "sha256")
  context <- .bayes_nutpie_bridge_context(missing)
  expect_true(context$migration)
  expect_identical(context$restoration$reason, "saved_library_unavailable")
  expect_identical(context$restoration$original_artifact$expected_hash, fit@stanfit$metadata$lib_hash)
  expect_identical(context$restoration$restored_artifact$hash, artifact$metadata$lib_hash)
  expect_equal(context$sampler_density_verification$max_absolute_difference, 0, tolerance = 1e-14)
  expect_identical(context$sampler_density_verification$points, 8L)
  expect_false(calls[[1L]]$fresh)
  expect_identical(digest::digest(missing, algo = "sha256"), missing_before)
  writeBin(charToRaw("changed original bytes"), fit@stanfit$metadata$lib_path)
  changed_hash <- digest::digest(file = fit@stanfit$metadata$lib_path, algo = "sha256")
  context <- .bayes_nutpie_bridge_context(fit)
  expect_identical(context$restoration$reason, "saved_library_hash_mismatch")
  expect_true(calls[[2L]]$fresh)
  expect_identical(digest::digest(file = fit@stanfit$metadata$lib_path, algo = "sha256"), changed_hash)
})

test_that("source data JSON and seed failures precede any native restoration", {
  directory <- withr::local_tempdir()
  fit <- .nutpie_portability_fixture(directory)
  local_mocked_bindings(
    .bayes_nutpie_bridge_client_version = function() "2.8.0",
    .bayes_nutpie_bridge_model = function(...) stop("Unexpected native load"),
    .bayes_nutpie_bridge_rebuild = function(...) stop("Unexpected compilation"))
  changed <- fit
  changed@stanfit$metadata$source_hash <- strrep("0", 64L)
  expect_error(.bayes_nutpie_bridge_context(changed), "source hashes differ")
  changed <- fit
  changed@stan_data$N <- 3L
  expect_error(.bayes_nutpie_bridge_context(changed), "Stan data differ")
  changed <- fit
  changed@stanfit$metadata$stan_data_hash <- strrep("0", 64L)
  expect_error(.bayes_nutpie_bridge_context(changed), "Stan data hashes differ")
  changed <- fit
  changed@stanfit$metadata$data_json <- "{}"
  expect_error(.bayes_nutpie_bridge_context(changed), "serialized Stan data differ")
  changed <- fit
  changed@stanfit$metadata$seed <- NA_real_
  expect_error(.bayes_nutpie_bridge_context(changed), "valid saved model-construction seed")
  changed <- fit
  changed@stanfit$metadata$num_chains <- 3L
  expect_error(.bayes_nutpie_bridge_context(changed), "every requested chain")
  expect_error(.bayes_nutpie_bridge_context(fit, list(kernel = "ces_supply")), "contract differs")
})

test_that("reconstructed source or binary mismatch fails before its native load", {
  directory <- withr::local_tempdir()
  fit <- .nutpie_portability_fixture(directory)
  fit@stanfit$metadata$lib_path <- file.path(directory, "absent.so")
  artifact <- .nutpie_portability_artifact(fit, directory)
  local_mocked_bindings(
    .bayes_nutpie_bridge_client_version = function() "2.8.0",
    .bayes_nutpie_bridge_model = function(...) stop("Unexpected native load"),
    .bayes_nutpie_bridge_rebuild = function(...) artifact)
  artifact$metadata$source_hash <- strrep("0", 64L)
  expect_error(.bayes_nutpie_bridge_context(fit), "source/library integrity")
  artifact$metadata$source_hash <- fit@metadata$stan_source_hash
  artifact$metadata$lib_hash <- strrep("0", 64L)
  expect_error(.bayes_nutpie_bridge_context(fit), "source/library integrity")
})

test_that("platform and loader failures reconstruct but incompatible ABIs remain unavailable", {
  directory <- withr::local_tempdir()
  withr::local_options(antitrustBayes.nutpie_cache = directory)
  fit <- .nutpie_portability_fixture(directory)
  artifact <- .nutpie_portability_artifact(fit, directory)
  rebuilt <- loaded <- 0L
  version <- "2.8.0"
  local_mocked_bindings(
    .bayes_nutpie_bridge_client_version = function() "2.8.0",
    .bayes_nutpie_bridge_model = function(lib, ...) {
      loaded <<- loaded + 1L
      if (identical(basename(lib), paste0("antitrustBayes_", fit@stanfit$metadata$lib_hash, ".so")))
        stop("incompatible binary")
      .nutpie_portability_model(version)
    },
    .bayes_nutpie_bridge_rebuild = function(kernel, fresh) {
      expect_true(fresh)
      rebuilt <<- rebuilt + 1L
      artifact
    })
  context <- .bayes_nutpie_bridge_context(fit)
  expect_match(context$restoration$reason, "saved_library_load_failed: incompatible binary")
  expect_identical(loaded, 2L)
  platform <- fit
  platform@stanfit$metadata$platform <- "different-platform"
  loaded <- 0L
  context <- .bayes_nutpie_bridge_context(platform)
  expect_identical(context$restoration$reason, "saved_library_platform_mismatch")
  expect_identical(loaded, 1L)
  version <- "2.9.0"
  expect_error(.bayes_nutpie_bridge_context(platform), "restored model ABI.*2.9.0")
})

test_that("reconstructed targets must reproduce saved sampler densities", {
  directory <- withr::local_tempdir()
  withr::local_options(antitrustBayes.nutpie_cache = directory)
  fit <- .nutpie_portability_fixture(directory)
  artifact <- .nutpie_portability_artifact(fit, directory)
  fit@stanfit$metadata$lib_path <- file.path(directory, "absent.so")
  local_mocked_bindings(
    .bayes_nutpie_bridge_client_version = function() "2.8.0",
    .bayes_nutpie_bridge_model = function(...) .nutpie_portability_model(offset = .01),
    .bayes_nutpie_bridge_rebuild = function(...) artifact)
  expect_error(.bayes_nutpie_bridge_context(fit), "disagrees with saved sampler log densities")
})

test_that("evidence reconstruction restores the user's external cache option after failure", {
  directory <- withr::local_tempdir()
  withr::local_options(antitrustBayes.nutpie_cache = directory)
  local_mocked_bindings(.nutpie_model_for_kernel = function(kernel) stop("compiler unavailable"))
  expect_error(.bayes_nutpie_bridge_rebuild("logit_supply", fresh = TRUE), "compiler unavailable")
  expect_identical(getOption("antitrustBayes.nutpie_cache"), directory)
})

test_that("restored libraries get distinct content names without changing compiler artifacts", {
  directory <- withr::local_tempdir()
  fit <- .nutpie_portability_fixture(directory)
  artifact <- .nutpie_portability_artifact(fit, directory)
  withr::local_options(antitrustBayes.nutpie_cache = directory)
  alias <- .bayes_nutpie_bridge_library_alias(artifact$metadata$lib_path, artifact$metadata$lib_hash)
  expect_false(identical(alias, artifact$metadata$lib_path))
  expect_false(identical(dirname(alias), dirname(artifact$metadata$lib_path)))
  expect_true(startsWith(alias, normalizePath(directory)))
  expect_identical(digest::digest(file = alias, algo = "sha256"), artifact$metadata$lib_hash)
  expect_identical(digest::digest(file = artifact$metadata$lib_path, algo = "sha256"), artifact$metadata$lib_hash)
})

test_that("saved-library density disagreements cannot trigger a compilation rescue", {
  directory <- withr::local_tempdir()
  withr::local_options(antitrustBayes.nutpie_cache = directory)
  fit <- .nutpie_portability_fixture(directory)
  local_mocked_bindings(
    .bayes_nutpie_bridge_client_version = function() "2.8.0",
    .bayes_nutpie_bridge_model = function(...) .nutpie_portability_model(offset = .01),
    .bayes_nutpie_bridge_rebuild = function(...) stop("Unexpected compilation rescue"))
  expect_error(.bayes_nutpie_bridge_context(fit), "disagrees with saved sampler log densities")
})

test_that("same-basename libraries with different contents get different native loader names", {
  directory <- withr::local_tempdir()
  withr::local_options(antitrustBayes.nutpie_cache = file.path(directory, "external-cache"))
  paths <- file.path(directory, c("first", "second"), "model.so")
  for (i in seq_along(paths)) {
    dir.create(dirname(paths[i]))
    writeBin(charToRaw(paste("different artifact", i)), paths[i])
  }
  hashes <- vapply(paths, digest::digest, character(1), file = TRUE, algo = "sha256")
  loaded_paths <- character()
  local_mocked_bindings(.bayes_nutpie_bridge_model = function(lib, ...) {
    loaded_paths <<- c(loaded_paths, lib)
    .nutpie_portability_model()
  })
  for (i in seq_along(paths)) .bayes_nutpie_bridge_load(paths[i], "unused.json", 1L, "2.8.0", unname(hashes[i]))
  expect_identical(basename(paths[1L]), basename(paths[2L]))
  expect_false(identical(basename(loaded_paths[1L]), basename(loaded_paths[2L])))
  expect_equal(unname(vapply(loaded_paths, digest::digest, character(1), file = TRUE, algo = "sha256")), unname(hashes))
})
