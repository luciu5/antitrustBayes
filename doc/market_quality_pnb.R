#!/usr/bin/env Rscript
# Frozen PNB comparison; no edits to mergerBayes or its historical outputs.
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(antitrustBayes))
suppressPackageStartupMessages(library(dplyr))
repo <- normalizePath(Sys.getenv("MQ_REPO", "/mnt/d/projects/antitrustbayes"))
merger_repo <- normalizePath(Sys.getenv("MQ_MERGER_REPO", "~/Projects/mergerBayes"))
root <- Sys.getenv("MQ_PNB_ROOT", file.path(repo,
  "doc/statistical_toolkit_runs/20260916/performance-audit/recovery/market-quality/pnb"))
data_source <- Sys.getenv("MQ_PNB_DATA_SOURCE", "mergerBayes")
args <- commandArgs(trailingOnly = TRUE)
models <- data.frame(case = c("B", "MC", "C", "A", "BF", "BS", "CF", "CS"),
  supply = c("bertrand", "monopolistic", "cournot", "auction",
    "bertrand", "bertrand", "cournot", "cournot"),
  structure = c(rep("simultaneous", 4L), "fringe", "stackelberg", "fringe", "stackelberg"))
frozen_options <- list(chains = 4L, iter = 8000L, warmup = 2000L, cores = 4L,
  seed = 1960L, refresh = 0L,
  control = list(adapt_delta = .99, max_treedepth = 12L))
prior <- c(B = 1/12, MC = .25, C = 1/12, A = .25,
  BF = 1/12, BS = 1/12, CF = 1/12, CS = 1/12)
hash <- function(x) digest::digest(x, algo = "sha256")
save_once <- function(x, path) {
  if (file.exists(path)) {
    stopifnot(identical(hash(readRDS(path)), hash(x)))
    return(invisible(path))
  }
  saveRDS(x, paste0(path, ".tmp")); stopifnot(file.rename(paste0(path, ".tmp"), path))
}
prepare <- function() {
  if (identical(data_source, "bundled")) {
    inputs <- new.env(parent = emptyenv())
    load(file.path(repo, "data/pnb.rda"), envir = inputs)
    data <- inputs$pnb
    data$margin_observed <- is.finite(data$margins) & data$margins > 0
    data_source_info <- list(kind = "bundled PNB, separately labeled application",
      path = file.path(repo, "data/pnb.rda"),
      file_sha256 = digest::digest(file.path(repo, "data/pnb.rda"), file = TRUE, algo = "sha256"),
      comparison = "price/share summary agrees with rich_main; margin inputs differ; differences are not likelihood-only")
  } else if (data_source %in% c("mergerBayes", "archived")) {
  # Import only the retrieved preparation function, never execute its runner.
  env <- new.env(parent = globalenv())
  archive <- if (identical(data_source, "archived"))
    normalizePath(Sys.getenv("MQ_PNB_ARCHIVED_FILE"), mustWork = TRUE) else NULL
  if (!is.null(archive) && basename(archive) != "supreme_data.RData")
    stop("The archived runner input must be named supreme_data.RData.")
  env$datadir <- if (is.null(archive)) file.path(merger_repo, "data") else dirname(archive)
  for (expr in parse(file.path(merger_repo, "code/pnb_rich_antitrustbayes.R"))) {
    if (is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
      identical(expr[[2L]], as.name("prepare_pnb_data"))) eval(expr, env)
  }
  data <- env$prepare_pnb_data()
  data_source_info <- if (is.null(archive))
    list(kind = "rich_main original-input preparation",
      path = attr(data, "source"), comparison = "same runner/preparation; preserve input provenance") else
    list(kind = "archived tracked-input reconstruction using rich_main preparation",
      path = archive, file_sha256 = digest::digest(archive, file = TRUE, algo = "sha256"),
      archive_manifest = jsonlite::fromJSON(file.path(dirname(archive), "manifest.json")),
      comparison = "tracked historical snapshot matches rich_main summaries; original process serialized input unavailable")
  } else stop("MQ_PNB_DATA_SOURCE must be mergerBayes, archived or bundled.")
  stopifnot(nrow(data) == 39L, length(unique(data$group)) == 39L)
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  save_once(data, file.path(root, "frozen-data.rds"))
  plan <- list(models = models, prior = prior, k_support = 1:5,
    k_prior = rep(.2, 5L), stan_options = frozen_options, df_active = FALSE,
    share_likelihood = "gaussian", margin_likelihood = "gaussian",
    quality_parameterization = "collapsed", quality_scope = "market",
    use_rho = TRUE, market_s0_parameterization = "legacy", market_side = "input",
    data_source = data_source_info, data_sha256 = hash(data), comparator_git = system2("git",
      c("-C", merger_repo, "rev-parse", "HEAD"), stdout = TRUE),
    comparator_driver_sha256 = digest::digest(file.path(merger_repo,
      "code/pnb_rich_antitrustbayes.R"), file = TRUE, algo = "sha256"))
  save_once(plan, file.path(root, "frozen-plan.rds"))
  write.csv(models, file.path(root, "frozen-models.csv"), row.names = FALSE)
  cat("PNB 39-bank data and settings frozen; no sampling performed.\n")
}
fit_one <- function(case) {
  stopifnot(antitrustBayes:::.antitrust_bayes_model_spec_version == 10L)
  gate <- file.path(dirname(root), "native-validation/approved-gate.rds")
  if (!file.exists(gate) || !isTRUE(readRDS(gate)$approved))
    stop("Native density/gradient, unit tests and independent review gate not approved.")
  p <- readRDS(file.path(root, "frozen-plan.rds"))
  data <- readRDS(file.path(root, "frozen-data.rds"))
  stopifnot(identical(p$data_sha256, hash(data)))
  global_regression <- identical(case, "B_global_regression")
  model <- p$models[p$models$case == if (global_regression) "B" else case, ]
  if (nrow(model) != 1L) stop("Unknown PNB case: ", case)
  files <- file.path(root, paste0(case, c("-started.rds", "-result.rds", "-fit.rds")))
  if (any(file.exists(files)))
    stop("Earlier PNB attempt retained; use a separately labeled operational retry directory.")
  opts <- p$stan_options
  opts$sample_file <- file.path(root, paste0(case, "-stan.csv"))
  call_args <- list(prices = data$prices, quantities = data$quantities,
    margins = data$margins, market = rep("PNB", nrow(data)), group = data$group,
    market_side = "input", demand_model = "logit", supply_model = model$supply,
    structure = model$structure, use_rho = TRUE,
    share_likelihood = "gaussian", margin_likelihood = "gaussian",
    quality_parameterization = "collapsed",
    quality_scope = if (global_regression) "global" else "market",
    prior_manifest = "baseline", market_s0_parameterization = "legacy",
    margin_observed = data$margin_observed, stan_options = opts)
  if (model$structure != "simultaneous") {
    call_args$partition_prior <- "ordered_k"
    call_args$k_support <- p$k_support; call_args$k_prior <- p$k_prior
    call_args$k_scope <- "market"
  }
  driver <- file.path(root, paste0(case, "-driver.R"))
  stopifnot(file.copy(file.path(repo, "doc/market_quality_pnb.R"), driver, overwrite = FALSE))
  save_once(list(case = case, args = call_args, time = Sys.time(), PID = Sys.getpid(),
    driver_sha256 = digest::digest(driver, file = TRUE, algo = "sha256"),
    model_hashes = vapply(antitrustBayes:::stanmodels, function(m)
      digest::digest(m@model_code, serialize = FALSE, algo = "sha256"), character(1))), files[1L])
  warnings <- character(); error <- ""; start <- proc.time()[["elapsed"]]
  fit <- tryCatch(withCallingHandlers(do.call(bayes_est, call_args), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  }), error = function(e) {error <<- conditionMessage(e); NULL})
  elapsed <- proc.time()[["elapsed"]]-start
  if (!is.null(fit)) saveRDS(fit, files[3L])
  saveRDS(list(case = case, elapsed = elapsed, error = error, warnings = warnings,
    native_seconds = if (!is.null(fit)) sum(rstan::get_elapsed_time(fit@stanfit)) else NA_real_,
    completed = Sys.time()), files[2L])
  if (nzchar(error)) stop(error)
  cat(case, "saved; wall seconds:", elapsed, "\n")
}
if ("--prepare" %in% args) prepare()
case <- args[startsWith(args, "--case=")]
if (length(case)) fit_one(sub("^--case=", "", case[[1L]]))
