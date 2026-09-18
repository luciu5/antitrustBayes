#!/usr/bin/env Rscript
# One market count per process so RSS measurements have separate baselines.
# This is a native target/gradient benchmark, not an HMC/recovery claim.
options(stringsAsFactors = FALSE)
repo <- normalizePath(Sys.getenv("MQ_REPO", "/mnt/d/projects/antitrustbayes"))
root <- file.path(repo, "doc/statistical_toolkit_runs/20260916/performance-audit/recovery/market-quality/scaling")
dir.create(root, recursive = TRUE, showWarnings = FALSE)
arguments <- commandArgs(trailingOnly = TRUE)
env <- new.env(parent = globalenv())
env$commandArgs <- function(...) character()
sys.source(file.path(repo, "doc/market_quality_benchmark.R"), env)
prepare_fixture <- function(M) {
  fixture <- env$mq_fixture(M, T = 2L, F = 4L, demand = "logit",
    conduct = "bertrand", structure = "fringe", seed = 161803L)
  env$mq_save_once(fixture, file.path(root, paste0("M", M, "-fixture.rds")))
}
if ("--prepare" %in% arguments) {
  for (M in c(10L, 100L, 500L, 2000L)) prepare_fixture(M)
  cat("Frozen scaling fixtures prepared; no target evaluation or HMC.\n")
}
requested <- arguments[startsWith(arguments, "--M=")]
if (length(requested)) {
  M <- as.integer(sub("^--M=", "", requested[[1L]]))
  stopifnot(M %in% c(10L, 100L, 500L, 2000L))
  suppressPackageStartupMessages(library(antitrustBayes))
  stopifnot(antitrustBayes:::.antitrust_bayes_model_spec_version == 10L)
  gate <- file.path(dirname(root), "native-validation/approved-gate.rds")
  stopifnot(file.exists(gate), isTRUE(readRDS(gate)$approved))
  fx <- readRDS(file.path(root, paste0("M", M, "-fixture.rds")))
  utility <- new.env(parent = globalenv())
  names_needed <- c("cga_as_draw", "cga_init_one", "cga_set_all", "cga_common_init",
    "cga_model_for", "cga_fixed", "cga_log_target")
  for (expr in parse(file.path(repo, "doc/collapsed_gaussian_production_validation.R")))
    if (is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
        is.symbol(expr[[2L]]) && as.character(expr[[2L]]) %in% names_needed) eval(expr, utility)
  args <- c(as.list(fx$data), list(supply_model = 1L, market_side = "output",
    use_rho = TRUE, use_cutoff = FALSE, use_alpha_market_pool = TRUE,
    use_s0_market_pool = TRUE, use_s0_market_year_pool = FALSE,
    use_market_resid_scale = FALSE, use_hmt = FALSE, priors = list(),
    X_s0 = NULL, outside_price = fx$outside_price, likelihood_mode = 1L,
    loo_cluster = fx$data$market, estimate_df = FALSE, resid_df_fixed = 5,
    structure = "fringe", partition_prior = "ordered_k", k_scope = "market",
    k_support = 0:3, k_prior = rep(.25, 4L), quality_scope = "market",
    share_likelihood = "gaussian", margin_likelihood = "gaussian",
    quality_parameterization = "collapsed"))
  start <- proc.time()[["elapsed"]]
  data <- do.call(antitrustBayes:::prepare_stan_data, args)
  setup <- proc.time()[["elapsed"]]-start
  stopifnot(data$N_market == M, data$N_quality == 4L*M,
    data$N_partition_component == 4L*M, !anyDuplicated(names(data)))
  model <- utility$cga_model_for("logit", "fringe")
  start <- proc.time()[["elapsed"]]
  initial <- utility$cga_fixed(model, data, init = 0, seed = 161803L)
  init <- utility$cga_common_init(initial, "logit", "bertrand", rho = .31, tau = .43, s0 = .21)
  fit <- utility$cga_fixed(model, data, init = list(init), seed = 161804L)
  construction <- proc.time()[["elapsed"]]-start
  u <- rstan::unconstrain_pars(fit, init)
  target <- rstan::log_prob(fit, u, adjust_transform = TRUE)
  gradient <- rstan::grad_log_prob(fit, u, adjust_transform = TRUE)
  stopifnot(is.finite(target), all(is.finite(gradient)))
  repetitions <- 25L
  timer <- function(fun) {
    gc()
    start <- proc.time()[["elapsed"]]
    for (rep in seq_len(repetitions)) fun()
    (proc.time()[["elapsed"]]-start)/repetitions
  }
  target_time <- timer(function() rstan::log_prob(fit, u, adjust_transform = TRUE))
  gradient_time <- timer(function() rstan::grad_log_prob(fit, u, adjust_transform = TRUE))
  probabilities <- rstan::extract(fit, pars = "posterior_component_prob")$posterior_component_prob
  stopifnot(max(abs(tapply(as.numeric(probabilities), data$partition_component_unit, sum)-1)) < 1e-12)
  peak <- readLines("/proc/self/status", warn = FALSE)
  peak_kb <- as.numeric(sub("^VmHWM:\\s+([0-9]+).*$", "\\1", peak[grepl("^VmHWM:", peak)]))
  result <- data.frame(M = M, N = data$N, quality_groups = data$N_quality,
    components = data$N_partition_component, sampling_dimensions = length(u),
    preparation_seconds = setup, fixed_state_and_GQ_seconds = construction,
    target_seconds = target_time, gradient_seconds = gradient_time,
    peak_RSS_kb = peak_kb, target = target, HMC_seconds = NA_real_, ESS_second = NA_real_,
    rho = .31, repetitions = repetitions, HMC_status = "not run")
  write.csv(result, file.path(root, paste0("M", M, "-native-timing.csv")), row.names = FALSE)
  saveRDS(list(result = result, data_sha256 = env$mq_hash(fx$data),
    model_code_sha256 = digest::digest(fit@stanmodel@model_code, algo = "sha256", serialize = FALSE),
    init = init, R_version = R.version.string), file.path(root, paste0("M", M, "-timing-provenance.rds")))
  print(result)
}
