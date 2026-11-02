#!/usr/bin/env Rscript

## Runtime validation plan.  This script never compiles a model.  It consumes
## an already compiled stanmodel RDS and a matching data RDS in a fresh R
## process, so the first TBB initialization can be observed without changing
## the caller's process state.

args <- commandArgs(trailingOnly = TRUE)

value_for <- function(flag, default = NA_character_) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^", flag, "="), "", hit[[1L]])
}

if (any(args %in% c("--help", "-h")) || !length(args)) {
  cat(paste(
    "Usage:",
    "validate_runtime.R --model-rds=MODEL.rds --data-rds=DATA.rds",
    "--threads=N [--out=REPORT.txt] [--force-error]",
    "\nN is tested in a fresh process; this script does not compile.\n",
    sep = "\n"
  ))
  quit(save = "no", status = if (length(args)) 0L else 2L)
}

main <- function() {
  model_path <- value_for("--model-rds")
  data_path <- value_for("--data-rds")
  thread_text <- value_for("--threads", "1")
  out_path <- value_for("--out")
  force_error <- "--force-error" %in% args

  if (is.na(model_path) || is.na(data_path))
    stop("--model-rds and --data-rds are required", call. = FALSE)
  threads <- suppressWarnings(as.integer(thread_text))
  if (is.na(threads) || threads < 1L)
    stop("--threads must be a positive integer", call. = FALSE)
  if (!requireNamespace("rstan", quietly = TRUE))
    stop("the rstan package is required", call. = FALSE)

  state <- function(tag) {
    env <- Sys.getenv("STAN_NUM_THREADS", unset = NA_character_)
    list(
      tag = tag,
      rstan_threads_per_chain = as.integer(rstan::rstan_options("threads_per_chain")),
      stan_num_threads_present = !is.na(env),
      stan_num_threads = if (is.na(env)) "<unset>" else env
    )
  }

  old_option <- as.integer(rstan::rstan_options("threads_per_chain"))
  old_env <- Sys.getenv("STAN_NUM_THREADS", unset = NA_character_)
  restore <- function() {
    ## rstan_options() itself writes STAN_NUM_THREADS, so restore the exact
    ## environment state after restoring the RStan option.
    try(rstan::rstan_options(threads_per_chain = old_option), silent = TRUE)
    if (is.na(old_env)) Sys.unsetenv("STAN_NUM_THREADS")
    else Sys.setenv(STAN_NUM_THREADS = old_env)
  }
  on.exit(restore(), add = TRUE)

## This setter must precede model loading/constructor activity.
  rstan::rstan_options(threads_per_chain = threads)
  state_after_setter <- state("after_setter")

  model <- readRDS(model_path)
  data <- readRDS(data_path)
  if (methods::is(model, "stanfit")) model <- rstan::get_stanmodel(model)
  if (!methods::is(model, "stanmodel"))
    stop("MODEL.rds must contain a stanmodel or stanfit", call. = FALSE)

  model_code <- methods::slot(model, "model_code")
  model_name <- methods::slot(model, "model_name")
  translation <- tryCatch({
  ## stanc() is parse-only here; no C++ compilation is requested.
  parsed <- rstan::stanc(model_code = model_code, model_name = model_name,
                         verbose = FALSE)
  cpp <- parsed$cppcode
  list(
    status = "ok",
    has_reduce_sum = grepl("reduce_sum", model_code, fixed = TRUE),
    has_thread_prologue = grepl("RSTAN_THREADING|init_threadpool_tbb", cpp),
    error = NA_character_
  )
  }, error = function(e) {
    list(status = "error", has_reduce_sum = NA,
         has_thread_prologue = NA, error = conditionMessage(e))
  })

  attempt <- tryCatch({
    if (force_error) stop("intentional state-restoration failure probe")
    ## One chain/core isolates within-chain runtime behavior from forked chains.
    fit <- rstan::sampling(
      object = model, data = data, chains = 1L, iter = 1L, warmup = 0L,
      algorithm = "Fixed_param", init = 0, seed = 1L, refresh = 0L, cores = 1L
    )
    list(status = "success", error = NA_character_)
  }, error = function(e) {
    list(status = "error", error = conditionMessage(e))
  })
  state_before_restore <- state("before_restore")
  restore()
  state_after_restore <- state("after_restore")

  report <- list(
    model = normalizePath(model_path, mustWork = FALSE),
    data = normalizePath(data_path, mustWork = FALSE),
    requested_threads = threads,
    forced_error = force_error,
    state_initial = list(
      rstan_threads_per_chain = old_option,
      stan_num_threads_present = !is.na(old_env),
      stan_num_threads = if (is.na(old_env)) "<unset>" else old_env
    ),
    state_after_setter = state_after_setter,
    translation = translation,
    sampling = attempt,
    state_before_restore = state_before_restore,
    state_after_restore = state_after_restore,
    actual_tbb_workers_observed = NA_integer_,
    actual_worker_note = paste(
      "An instrumented reducer/TBB probe is still required;",
      "RStan success does not identify the worker count."
    )
  )

  output <- capture.output(dput(report))
  if (is.na(out_path)) writeLines(output) else writeLines(output, out_path)
  invisible(report)
}

main()
