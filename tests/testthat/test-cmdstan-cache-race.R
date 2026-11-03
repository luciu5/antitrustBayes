## Concurrent fits that start together can all find the per-source cache
## directory missing. A process that loses the create/copy race must treat
## "already exists" as success rather than a failure. The races are staged with
## real filesystem state, using a compile stub in place of CmdStan.

.cache_race_context <- function(env = parent.frame()) {
  skip_if_not_installed("cmdstanr")
  base <- withr::local_tempdir(.local_envir = env)
  withr::local_options(antitrustBayes.cmdstan_cache = base, .local_envir = env)
  cache_env <- antitrustBayes:::.cmdstan_model_cache
  rm(list = ls(cache_env, all.names = TRUE), envir = cache_env)
  state <- new.env()
  state$compiled <- 0L
  state$version <- "2.39.0"
  state$path <- file.path(base, "cmdstan")
  dir.create(state$path)
  compile <- function(file, ...) {
    state$compiled <- state$compiled + 1L
    exe <- paste0(file, ".exe"); writeLines("x", exe)
    structure(list(exe_file = function() exe, file = file), class = "CmdStanModel")
  }
  testthat::local_mocked_bindings(
    cmdstan_model = compile,
    cmdstan_version = function() state$version,
    cmdstan_path = function() state$path,
    .package = "cmdstanr", .env = env
  )
  list(base = base, state = state)
}

.cache_dir_for <- function(base, expected = 1L) {
  dirs <- list.dirs(file.path(base, "stan"), recursive = FALSE)
  expect_length(dirs, expected)
  dirs
}

test_that("a competitor having already created the cache directory and source is not an error", {
  ctx <- .cache_race_context()
  ## First call stages the hashed directory exactly as a winning competitor would.
  first <- antitrustBayes:::.cmdstan_model_for_kernel("ces_supply")
  cache <- .cache_dir_for(ctx$base)
  expect_true(file.exists(file.path(cache, "ces_supply.stan")))
  ## A second process in a fresh R session finds directory, source and model.rds
  ## present; it must reuse them rather than fail on "already exists".
  rm(list = ls(antitrustBayes:::.cmdstan_model_cache, all.names = TRUE),
     envir = antitrustBayes:::.cmdstan_model_cache)
  second <- antitrustBayes:::.cmdstan_model_for_kernel("ces_supply")
  expect_s3_class(second, "CmdStanModel")
  expect_identical(ctx$state$compiled, 1L)
})

test_that("the disk cache separates CmdStan versions and resolved paths", {
  ctx <- .cache_race_context()
  antitrustBayes:::.cmdstan_model_for_kernel("ces_supply")
  original <- .cache_dir_for(ctx$base)

  ## A different CmdStan version must not reuse an object compiled against the
  ## original toolchain, even though the kernel and cache root are unchanged.
  ctx$state$version <- "2.39.1"
  antitrustBayes:::.cmdstan_model_for_kernel("ces_supply")
  by_version <- .cache_dir_for(ctx$base, expected = 2L)
  expect_true(original %in% by_version)
  expect_length(setdiff(by_version, original), 1L)

  ## Restoring the version while changing the normalized toolchain path is a
  ## distinct identity too; this catches reuse across separate installations.
  ctx$state$version <- "2.39.0"
  ctx$state$path <- file.path(ctx$base, "other-cmdstan")
  dir.create(ctx$state$path)
  antitrustBayes:::.cmdstan_model_for_kernel("ces_supply")
  by_path <- .cache_dir_for(ctx$base, expected = 3L)
  expect_length(setdiff(by_path, by_version), 1L)
  expect_identical(ctx$state$compiled, 3L)
})

test_that("a pre-existing directory without a model compiles once and leaves no staging files", {
  ctx <- .cache_race_context()
  antitrustBayes:::.cmdstan_model_for_kernel("ces_supply")
  cache <- .cache_dir_for(ctx$base)
  ## Simulate the competitor that created the directory and source but has not
  ## finished compiling: remove the model object only.
  unlink(list.files(cache, pattern = "model\\.rds$", full.names = TRUE))
  rm(list = ls(antitrustBayes:::.cmdstan_model_cache, all.names = TRUE),
     envir = antitrustBayes:::.cmdstan_model_cache)
  model <- antitrustBayes:::.cmdstan_model_for_kernel("ces_supply")
  expect_s3_class(model, "CmdStanModel")
  expect_identical(ctx$state$compiled, 2L)
  expect_false(any(grepl("^ces_supply_[^.]+\\.stan$", list.files(cache))))
  expect_identical(readLines(file.path(cache, "ces_supply.stan")),
                   readLines(antitrustBayes:::.cmdstan_kernel_file("ces_supply")))
})

test_that("concurrent first use of the same cache succeeds in every process", {
  skip_on_os("windows")
  ctx <- .cache_race_context()
  ## One trial rarely loses the race, so repeat with fresh caches; the original
  ## implementation failed in roughly 40% of 8-way starts on this filesystem.
  ## Workers wait for a release file so they reach the cache check together.
  failures <- character()
  for (trial in seq_len(12L)) {
    withr::local_options(antitrustBayes.cmdstan_cache = file.path(ctx$base, paste0("t", trial)))
    release <- tempfile("release")
    workers <- lapply(seq_len(8L), function(i) parallel::mcparallel(tryCatch({
      rm(list = ls(antitrustBayes:::.cmdstan_model_cache, all.names = TRUE),
         envir = antitrustBayes:::.cmdstan_model_cache)
      while (!file.exists(release)) Sys.sleep(0.0005)
      antitrustBayes:::.cmdstan_model_for_kernel("ces_supply")
      "ok"
    }, error = function(e) conditionMessage(e))))
    Sys.sleep(0.05)
    file.create(release)
    results <- unlist(parallel::mccollect(workers))
    failures <- c(failures, setdiff(results, "ok"))
  }
  expect_length(failures, 0L)
  expect_identical(readLines(file.path(.cache_dir_for(file.path(ctx$base, "t1")), "ces_supply.stan")),
                   readLines(antitrustBayes:::.cmdstan_kernel_file("ces_supply")))
})

test_that("a cache location that cannot be created still fails clearly", {
  skip_on_os("windows")
  skip_if(identical(Sys.info()[["user"]], "root"), "root can write anywhere")
  ctx <- .cache_race_context()
  blocker <- file.path(ctx$base, "stan")
  writeLines("a file, not a directory", blocker)
  expect_error(antitrustBayes:::.cmdstan_model_for_kernel("ces_supply"),
               "Cannot create CmdStan cache directory")
  expect_identical(ctx$state$compiled, 0L)
})
