#!/usr/bin/env Rscript

## Standalone, frozen hybrid benchmark runner.
##
## This runner sources the baseline benchmark implementation copied beside it,
## retains the seven mode-B primary designs, and changes only the fitted
## margin likelihood to fixed-df=5 Student-t.  It has two explicit phases:
##
##   Rscript run_hybrid_gaussian_t_margin.R --freeze
##   Rscript run_hybrid_gaussian_t_margin.R --run
##
## The second phase refuses to run unless the first phase completed.  All
## artifacts, including failures and warnings, remain below /tmp.

options(warn = 1)

hybrid_args <- commandArgs(trailingOnly = TRUE)
hybrid_out <- normalizePath(dirname(dirname(sub("^--file=", "",
  commandArgs(trailingOnly = FALSE)[startsWith(commandArgs(trailingOnly = FALSE), "--file=")][[1L]]))),
  mustWork = FALSE)
hybrid_repo <- "/mnt/d/projects/antitrustbayes"
hybrid_lib <- "/tmp/cga-Rlib"
hybrid_expected_commit <- "12b26c0"
hybrid_reference_csv <- "/tmp/cga-primary-20260916-configfix/benchmark.csv"
hybrid_baseline_script <- file.path(hybrid_out, "scripts",
  "collapsed_gaussian_benchmark_baseline.R")

if (any(hybrid_args %in% c("--help", "-h"))) {
  cat("Usage: Rscript run_hybrid_gaussian_t_margin.R --freeze|--run [--resume]\n")
  quit(save = "no", status = 0L)
}
if (!length(hybrid_args) || !any(hybrid_args %in% c("--freeze", "--run")))
  stop("Choose exactly one phase: --freeze or --run")
if (sum(hybrid_args %in% c("--freeze", "--run")) != 1L)
  stop("Choose exactly one phase: --freeze or --run")
hybrid_phase <- hybrid_args[hybrid_args %in% c("--freeze", "--run")][[1L]]
hybrid_resume <- "--resume" %in% hybrid_args
if (hybrid_resume && hybrid_phase != "--run")
  stop("--resume requires --run")

if (!file.exists(hybrid_baseline_script))
  stop("Missing frozen baseline script: ", hybrid_baseline_script)
source(hybrid_baseline_script, local = TRUE)

dir.create(hybrid_out, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(hybrid_out, "data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(hybrid_out, "logs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(hybrid_out, "stanfit"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(hybrid_out, "reference"), recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(hybrid_lib)) stop("Missing frozen package library: ", hybrid_lib)
.libPaths(c(hybrid_lib, .libPaths()))
suppressPackageStartupMessages(library(antitrustBayes, lib.loc = hybrid_lib))

hybrid_digest_file <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  digest::digest(file = path, algo = "sha256")
}

hybrid_git <- function(...) {
  out <- tryCatch(system2("git", c(...), stdout = TRUE, stderr = TRUE),
    error = function(e) character())
  if (!length(out)) NA_character_ else paste(out, collapse = "\n")
}

hybrid_commit <- hybrid_git("-C", hybrid_repo, "rev-parse", "HEAD")
if (!identical(substr(hybrid_commit, 1L, nchar(hybrid_expected_commit)),
               hybrid_expected_commit))
  stop("Baseline checkout is not the requested commit; found ", hybrid_commit,
       ", expected prefix ", hybrid_expected_commit)

hybrid_stan_files <- c("logit_supply.stan", "ces_supply.stan",
  "logit_stackelberg.stan", "ces_stackelberg.stan")
hybrid_repo_stan <- file.path(hybrid_repo, "inst", "stan", hybrid_stan_files)
hybrid_inst_stan <- file.path(hybrid_lib, "antitrustBayes", "stan", hybrid_stan_files)
if (any(!file.exists(hybrid_repo_stan)) || any(!file.exists(hybrid_inst_stan)))
  stop("Missing baseline or installed Stan source needed for equivalence check")
hybrid_repo_hashes <- vapply(hybrid_repo_stan, hybrid_digest_file, character(1))
hybrid_inst_hashes <- vapply(hybrid_inst_stan, hybrid_digest_file, character(1))
if (!identical(unname(hybrid_repo_hashes), unname(hybrid_inst_hashes)))
  stop("Installed Stan sources do not match baseline checkout")

hybrid_plan <- function() {
  base <- cga_make_plan()
  out <- base[base$mode == "B", , drop = FALSE]
  if (nrow(out) != 7L) stop("Expected exactly seven mode-B primary designs")
  out$source_mode <- "B"
  out$mode <- "D"
  out$share_likelihood <- "gaussian"
  out$quality_parameterization <- "collapsed"
  out$margin_likelihood <- "student_t"
  out$seed_source_case_id <- paste(out$design, "B", sep = "__")
  out$sampler_seed <- vapply(seq_len(nrow(out)), function(i)
    cga_derived_seed(out$seed[[i]], out$seed_source_case_id[[i]]), integer(1))
  out$case_id <- paste(out$design, "D", sep = "__")
  out$run_profile <- "hybrid_gaussian_share_student_t_margin"
  out$comparison_label <- "D: collapsed Gaussian shares + fixed-df5 Student-t margins"
  out <- out[order(match(out$design,
    c("logit_bertrand_exact_271828", "logit_bertrand_noisy_271828",
      "logit_bertrand_noisy_314159", "logit_cournot_noisy_271828",
      "logit_BF_noisy_271828", "logit_BS_noisy_271828",
      "ces_bertrand_noisy_271828"))), , drop = FALSE]
  rownames(out) <- NULL
  out
}

hybrid_validate_plan <- function(plan) {
  required <- c("case_id", "design", "source_mode", "mode", "demand_model",
    "supply_model", "structure", "share_scenario", "seed", "partition_prior",
    "seed_source_case_id", "sampler_seed",
    "share_likelihood", "margin_likelihood", "quality_parameterization",
    "use_rho", "use_market_resid_scale", "use_alpha_market_pool",
    "use_s0_market_pool", "use_s0_market_year_pool", "use_centered_alpha",
    "use_centered_s0_pool", "use_hmt", "df", "chains", "iter", "warmup",
    "adapt_delta", "max_treedepth", "metric", "k_support", "k_prior")
  if (!all(required %in% names(plan)) || nrow(plan) != 7L)
    stop("Hybrid plan is missing required fields or is not seven rows")
  if (!all(plan$source_mode == "B") || !all(plan$mode == "D"))
    stop("Hybrid labels do not identify the source B plan and D result")
  expected_sampler_seed <- vapply(seq_len(nrow(plan)), function(i)
    cga_derived_seed(plan$seed[[i]], plan$seed_source_case_id[[i]]), integer(1))
  if (!identical(as.integer(plan$sampler_seed), as.integer(expected_sampler_seed)))
    stop("Hybrid sampler seeds do not reproduce the original B case IDs")
  if (!all(plan$share_likelihood == "gaussian") ||
      !all(plan$quality_parameterization == "collapsed") ||
      !all(plan$margin_likelihood == "student_t") || !all(plan$df == 5))
    stop("Hybrid likelihood/API fields are not frozen")
  if (!identical(as.integer(plan$chains), rep(4L, 7L)) ||
      !identical(as.integer(plan$iter), rep(1000L, 7L)) ||
      !identical(as.integer(plan$warmup), rep(500L, 7L)) ||
      !identical(as.integer(plan$max_treedepth), rep(13L, 7L)) ||
      !all(plan$adapt_delta == .95) || !all(plan$metric == "diag_e"))
    stop("Hybrid sampler settings changed")
  flags <- c("use_rho", "use_market_resid_scale", "use_alpha_market_pool",
    "use_s0_market_pool", "use_s0_market_year_pool", "use_centered_alpha",
    "use_centered_s0_pool", "use_hmt")
  if (any(vapply(flags, function(x) any(plan[[x]]), logical(1))))
    stop("Hybrid legacy/rho/HMT flags must all be false")
  invisible(TRUE)
}

hybrid_make_panels <- function(plan, save_files = TRUE) {
  panels <- vector("list", nrow(plan))
  for (i in seq_len(nrow(plan))) {
    case <- plan[i, , drop = FALSE]
    panel <- cga_make_dgp(seed = case$seed,
      demand_model = case$demand_model, supply_model = case$supply_model,
      structure = case$structure, share_scenario = case$share_scenario)
    check <- cga_validate_dgp(panel)
    if (!identical(check$status, "passed"))
      stop("DGP validation failed for ", case$design)
    panels[[i]] <- panel
    if (save_files) {
      csv_path <- file.path(hybrid_out, "data", paste0(case$design, ".csv"))
      rds_path <- file.path(hybrid_out, "data", paste0(case$design, ".rds"))
      utils::write.csv(panel$data, csv_path, row.names = FALSE, na = "")
      saveRDS(panel, rds_path)
    }
  }
  panels
}

hybrid_reference <- function() {
  if (!file.exists(hybrid_reference_csv))
    stop("Missing prior primary CSV: ", hybrid_reference_csv)
  ref <- utils::read.csv(hybrid_reference_csv, stringsAsFactors = FALSE,
    na.strings = c("", "NA"))
  ref <- ref[ref$mode == "B", , drop = FALSE]
  if (nrow(ref) != 7L) stop("Prior primary CSV does not contain seven B rows")
  ref <- ref[order(match(ref$design, hybrid_plan()$design)), , drop = FALSE]
  utils::write.csv(ref, file.path(hybrid_out, "reference",
    "primary_B_metrics.csv"), row.names = FALSE, na = "")
  ref
}

hybrid_freeze <- function(plan) {
  hybrid_validate_plan(plan)
  panels <- hybrid_make_panels(plan, save_files = TRUE)
  ref <- hybrid_reference()
  plan_path <- file.path(hybrid_out, "benchmark_plan.csv")
  utils::write.csv(plan, plan_path, row.names = FALSE, na = "")
  data_hash <- vapply(plan$design, function(design)
    hybrid_digest_file(file.path(hybrid_out, "data", paste0(design, ".csv"))),
    character(1))
  rds_hash <- vapply(plan$design, function(design)
    hybrid_digest_file(file.path(hybrid_out, "data", paste0(design, ".rds"))),
    character(1))
  row_common <- data.frame(
    frozen_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    baseline_commit = hybrid_commit,
    package = "antitrustBayes", package_version = as.character(packageVersion("antitrustBayes")),
    package_library = hybrid_lib,
    baseline_script_sha256 = hybrid_digest_file(hybrid_baseline_script),
    protocol_sha256 = hybrid_digest_file(file.path(hybrid_out, "scripts",
      "collapsed_gaussian_benchmark_protocol_baseline.md")),
    prior_primary_csv = hybrid_reference_csv,
    prior_primary_csv_sha256 = hybrid_digest_file(hybrid_reference_csv),
    source_change = "margin_likelihood only: gaussian -> student_t",
    share_likelihood = "gaussian", quality_parameterization = "collapsed",
    margin_likelihood = "student_t", df = 5, use_rho = FALSE,
    use_market_resid_scale = FALSE, use_alpha_market_pool = FALSE,
    use_s0_market_pool = FALSE, use_s0_market_year_pool = FALSE,
    use_centered_alpha = FALSE, use_centered_s0_pool = FALSE, use_hmt = FALSE,
    chains = 4L, iter = 1000L, warmup = 500L, adapt_delta = .95,
    max_treedepth = 13L, metric = "diag_e", cores = 1L,
    gate = "zero divergences; zero max-depth hits; BFMI >= 0.30; rank Rhat <= 1.01; bulk/tail ESS >= 100",
    policy_status = "not run in this bounded benchmark",
    stringsAsFactors = FALSE)
  common <- row_common[rep(1L, nrow(plan)), , drop = FALSE]
  common$case_id <- plan$case_id
  common$design <- plan$design
  common$source_mode <- plan$source_mode
  common$mode <- plan$mode
  common$seed_source_case_id <- plan$seed_source_case_id
  common$sampler_seed <- plan$sampler_seed
  common$demand_model <- plan$demand_model
  common$supply_model <- plan$supply_model
  common$structure <- plan$structure
  common$share_scenario <- plan$share_scenario
  common$seed <- plan$seed
  common$partition_prior <- plan$partition_prior
  common$k_support <- plan$k_support
  common$k_prior <- plan$k_prior
  common$k_true <- plan$k_true
  common$bridge_requested <- plan$bridge_requested
  common$data_csv <- file.path("data", paste0(plan$design, ".csv"))
  common$data_csv_sha256 <- data_hash
  common$data_rds <- file.path("data", paste0(plan$design, ".rds"))
  common$data_rds_sha256 <- rds_hash
  common$repo_stan_sha256 <- paste(names(hybrid_repo_hashes), hybrid_repo_hashes,
    sep = "=", collapse = ";")
  common$installed_stan_sha256 <- paste(names(hybrid_inst_hashes), hybrid_inst_hashes,
    sep = "=", collapse = ";")
  utils::write.csv(common, file.path(hybrid_out, "hybrid_manifest.csv"),
    row.names = FALSE, na = "")
  utils::write.csv(row_common, file.path(hybrid_out, "benchmark_run_manifest.csv"),
    row.names = FALSE, na = "")
  writeLines(c(
    "HYBRID BENCHMARK MANIFEST FROZEN",
    paste("frozen_utc:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste("baseline_commit:", hybrid_commit),
    paste("package:", as.character(packageVersion("antitrustBayes")), "from", hybrid_lib),
    "designs: seven mode-B primary panels, seeds/priors/settings retained",
    "change: fitted margin_likelihood=student_t with fixed df=5",
    "share_likelihood=gaussian; quality_parameterization=collapsed; use_rho=FALSE",
    "sampler: 4 chains x 1000 iterations, 500 warmup, adapt_delta=.95, depth=13, diag_e, cores=1",
    "gate: zero divergences/depth hits, BFMI>=.30, rank Rhat<=1.01, bulk/tail ESS>=100",
    "policy counterfactuals: intentionally not run in this bounded benchmark",
    "The reference B rows are copied to reference/primary_B_metrics.csv.",
    "All panel data are stored as both CSV inputs and RDS truth fixtures."),
    file.path(hybrid_out, "MANIFEST_FROZEN.txt"))
  invisible(list(plan = plan, panels = panels, reference = ref, manifest = common))
}

hybrid_assert_frozen <- function(plan) {
  path <- file.path(hybrid_out, "hybrid_manifest.csv")
  if (!file.exists(path)) stop("Missing frozen manifest; run --freeze first")
  man <- utils::read.csv(path, stringsAsFactors = FALSE, na.strings = c("", "NA"))
  if (nrow(man) != 7L || !all(man$baseline_commit == hybrid_commit) ||
      !all(man$package_version == as.character(packageVersion("antitrustBayes"))) ||
      !all(man$margin_likelihood == "student_t") || !all(man$share_likelihood == "gaussian") ||
      !all(man$quality_parameterization == "collapsed") ||
      !identical(as.integer(man$sampler_seed), as.integer(plan$sampler_seed)))
    stop("Frozen manifest does not match current requested baseline/API")
  for (i in seq_len(nrow(man))) {
    csv_path <- file.path(hybrid_out, man$data_csv[[i]])
    rds_path <- file.path(hybrid_out, man$data_rds[[i]])
    if (!identical(hybrid_digest_file(csv_path), man$data_csv_sha256[[i]]) ||
        !identical(hybrid_digest_file(rds_path), man$data_rds_sha256[[i]]))
      stop("Frozen data hash mismatch for ", man$design[[i]])
  }
  hybrid_validate_plan(plan)
  invisible(man)
}

hybrid_add_labels <- function(diag, case) {
  diag$source_mode <- "B"
  diag$mode <- "D"
  diag$share_likelihood <- "gaussian"
  diag$margin_likelihood <- "student_t"
  diag$quality_parameterization <- "collapsed"
  diag$comparison_label <- "D: collapsed Gaussian shares + fixed-df5 Student-t margins"
  diag$predictive_quality_conditioning <- "conditional Gaussian quality draws given shares"
  diag$gate_zero_divergences <- isTRUE(diag$divergences == 0)
  diag$gate_zero_depth_hits <- isTRUE(diag$treedepth_hits == 0)
  diag$gate_bfmi <- isTRUE(diag$min_bfmi >= .30)
  diag$gate_rhat <- isTRUE(diag$max_rhat <= 1.01)
  diag$gate_bulk_ess <- isTRUE(diag$min_bulk_ess >= 100)
  diag$gate_tail_ess <- isTRUE(diag$min_tail_ess >= 100)
  diag$gate_pass <- isTRUE(diag$usable)
  diag
}

## cga_fit_args derives the Stan seed from case_id.  D has a distinct result
## label, so this wrapper restores the exact seed used by the corresponding B
## case while changing no other estimator argument.
hybrid_estimator_for <- function(case) {
  sampler_seed <- as.integer(case$sampler_seed[[1L]])
  function(...) {
    args <- list(...)
    args$stan_options$seed <- sampler_seed
    do.call(getExportedValue("antitrustBayes", "bayes_est"), args)
  }
}

hybrid_failure_result <- function(panel, case, err) {
  d <- cga_fit_diagnostics(NULL, panel, case, NA_real_, err, "")
  list(diagnostics = hybrid_add_labels(d, case),
    bridge = data.frame(case_id = case$case_id, bridge_status = "not_attempted",
      bridge_log_marginal = NA_real_, bridge_error = err, stringsAsFactors = FALSE),
    policy = cga_pending_policy(case$case_id, "not run in bounded hybrid benchmark"),
    log_path = file.path(hybrid_out, "logs", paste0(case$case_id, ".log")), fit = NULL)
}

hybrid_bridge_with_id <- function(value, case_id) {
  if (is.data.frame(value)) {
    if (!"case_id" %in% names(value)) value$case_id <- case_id
    value <- value[, c("case_id", setdiff(names(value), "case_id")), drop = FALSE]
  } else {
    coalesce <- function(x, y) if (is.null(x) || !length(x) || is.na(x[[1L]])) y else x[[1L]]
    value <- data.frame(case_id = case_id,
      bridge_status = coalesce(value$status, "not_attempted"),
      bridge_log_marginal = coalesce(value$log_marginal, NA_real_),
      bridge_error = coalesce(value$error, ""), stringsAsFactors = FALSE)
  }
  value
}

hybrid_write_checkpoint <- function(state, diagnostics, bridges, policies, ledger) {
  utils::write.csv(state, file.path(hybrid_out, "benchmark_progress.csv"),
    row.names = FALSE, na = "")
  if (length(diagnostics)) utils::write.csv(do.call(rbind, diagnostics),
    file.path(hybrid_out, "benchmark.csv"), row.names = FALSE, na = "")
  if (length(bridges)) utils::write.csv(do.call(rbind, bridges),
    file.path(hybrid_out, "benchmark_bridges.csv"), row.names = FALSE, na = "")
  if (length(policies)) utils::write.csv(do.call(rbind, policies),
    file.path(hybrid_out, "benchmark_policy_hooks.csv"), row.names = FALSE, na = "")
  if (length(ledger)) utils::write.csv(do.call(rbind, ledger),
    file.path(hybrid_out, "benchmark_run_ledger.csv"), row.names = FALSE, na = "")
  invisible(TRUE)
}

hybrid_compare <- function(diag, ref) {
  metrics <- c("elapsed_seconds", "stan_seconds", "gradient_evaluations",
    "divergences", "treedepth_hits", "min_bfmi", "max_rhat", "min_bulk_ess",
    "min_tail_ess", "alpha_ess_per_sec", "gamma_ess_per_sec",
    "alpha_or_gamma_ess_per_sec", "s0_ess_per_sec", "tau_ess_per_sec",
    "sigma_share_ess_per_sec", "alpha_ess_per_stan_second",
    "gamma_ess_per_stan_second", "s0_ess_per_stan_second",
    "tau_ess_per_stan_second", "sigma_share_ess_per_stan_second",
    "alpha_q05", "alpha_q50", "alpha_q95", "gamma_q05", "gamma_q50", "gamma_q95",
    "s0_q05", "s0_q50", "s0_q95", "tau_q05", "tau_q50", "tau_q95",
    "sigma_share_q05", "sigma_share_q50", "sigma_share_q95", "sigma_margin_q05",
    "sigma_margin_q50", "sigma_margin_q95", "tau_sigma_correlation",
    "share_predictive_median_mae", "share_predictive_median_max_abs",
    "margin_predictive_median_mae", "margin_predictive_median_max_abs")
  keep_ref <- c("design", "case_id", "status", "usable", metrics[metrics %in% names(ref)])
  keep_h <- c("design", "case_id", "status", "usable", metrics[metrics %in% names(diag)])
  b <- ref[, keep_ref, drop = FALSE]
  h <- diag[, keep_h, drop = FALSE]
  names(b)[names(b) != "design"] <- paste0("B_", names(b)[names(b) != "design"])
  names(h)[names(h) != "design"] <- paste0("D_", names(h)[names(h) != "design"])
  out <- merge(b, h, by = "design", all = TRUE, sort = FALSE)
  out$comparison_note <- "B Gaussian margins vs D fixed-df5 Student-t margins; same data/prior/economics; no equality claim"
  utils::write.csv(out, file.path(hybrid_out, "benchmark_vs_primary_B.csv"),
    row.names = FALSE, na = "")
  out
}

hybrid_report <- function(diag, bridges, ref, compare) {
  usable <- sum(as.logical(diag$usable), na.rm = TRUE)
  status_lines <- paste(sprintf("- %s: status=%s, usable=%s, divergences=%s, depth_hits=%s, BFMI=%s, Rhat=%s, bulk_ESS=%s, tail_ESS=%s, elapsed_s=%s, native_Stan_s=%s",
    diag$design, diag$status, diag$usable, diag$divergences, diag$treedepth_hits,
    format(diag$min_bfmi, digits = 6), format(diag$max_rhat, digits = 6),
    format(diag$min_bulk_ess, digits = 6), format(diag$min_tail_ess, digits = 6),
    format(diag$elapsed_seconds, digits = 6), format(diag$stan_seconds, digits = 6)),
    collapse = "\n")
  bridge_lines <- paste(sprintf("- %s: %s; log marginal=%s; %s", bridges$case_id,
    bridges$bridge_status, format(bridges$bridge_log_marginal, digits = 8),
    ifelse(is.na(bridges$bridge_error) | !nzchar(bridges$bridge_error), "", bridges$bridge_error)),
    collapse = "\n")
  text <- c(
    "# D: collapsed Gaussian share + Student-t margin benchmark",
    "",
    "This bounded run uses the seven existing primary mode-B panels, seeds, priors, economic equations, and sampler settings from the baseline benchmark at commit `12b26c0`. The sole fitted likelihood change is `margin_likelihood=\"student_t\"` with fixed `df=5`; shares remain Gaussian and centered quality remains collapsed. D retains each original B-derived Stan sampler seed exactly.",
    "",
    "The usability gate is the prespecified gate: zero divergences, zero maximum-depth hits, minimum BFMI at least 0.30, maximum rank-normalized R-hat at most 1.01, and minimum monitored bulk and tail ESS at least 100. The diagnostics preserve wall elapsed time, native Stan time, gradient evaluations, diagnostics, ESS/sec, posterior summaries, predictive summaries, warnings, and all failed runs.",
    "",
    paste0("Usable fits: ", usable, "/", nrow(diag), "."),
    "",
    "## Per-fit gate results",
    "",
    status_lines,
    "",
    "## Bridge result",
    "",
    "Bridge sampling follows the baseline request and is attempted only for the noisy seed-271828 simultaneous-Bertrand cell. A bridge result is evidence for this D specification only.",
    "",
    bridge_lines,
    "",
    "## Comparison with primary B",
    "",
    "`benchmark_vs_primary_B.csv` aligns each D row with the corresponding prior B row. B uses Gaussian margins; D uses fixed-df5 Student-t margins. The comparison holds the realized panel data, structural priors, equations, and sampler settings fixed, and makes no claim that the two likelihoods are numerically identical.",
    "",
    "## Scope",
    "",
    "Policy counterfactuals were not run in this bounded benchmark. No defaults, priors, sampler settings, production source, or repository benchmark CSVs were changed.",
    "",
    paste0("Artifacts: `", hybrid_out, "` (durable archive: `/mnt/d/Projects/antitrustBayes/doc/statistical_toolkit_runs/20260916/hybrid`). The frozen input manifests are `frozen/hybrid_manifest.csv` and `frozen/benchmark_run_manifest.csv` in that archive; panel data are under `frozen/data/`; fit logs and RDS objects are under `fits/`.")
  )
  writeLines(text, file.path(hybrid_out, "hybrid_benchmark_report.md"))
  invisible(text)
}

plan <- hybrid_plan()
hybrid_validate_plan(plan)

if (hybrid_phase == "--freeze") {
  if (file.exists(file.path(hybrid_out, "hybrid_manifest.csv")))
    stop("Manifest already exists; preserve this frozen run and use a new output directory")
  hybrid_freeze(plan)
  cat("HYBRID MANIFEST FROZEN: ", hybrid_out, "\n", sep = "")
  cat("Seven panels/data inputs and baseline B reference metrics saved; no HMC run attempted.\n")
  quit(save = "no", status = 0L)
}

hybrid_assert_frozen(plan)
manifest <- utils::read.csv(file.path(hybrid_out, "hybrid_manifest.csv"),
  stringsAsFactors = FALSE, na.strings = c("", "NA"))
ref <- utils::read.csv(file.path(hybrid_out, "reference", "primary_B_metrics.csv"),
  stringsAsFactors = FALSE, na.strings = c("", "NA"))

state_path <- file.path(hybrid_out, "benchmark_progress.csv")
if (hybrid_resume) {
  if (!file.exists(state_path)) stop("--resume requires benchmark_progress.csv")
  state <- utils::read.csv(state_path, stringsAsFactors = FALSE)
  if (!identical(as.character(state$case_id), as.character(plan$case_id)))
    stop("Saved progress case order differs from frozen plan")
} else {
  if (file.exists(state_path) || file.exists(file.path(hybrid_out, "benchmark.csv")))
    stop("Existing run artifacts must be preserved; use --resume or a new output directory")
  state <- data.frame(case_id = plan$case_id, status = "unstarted", artifact = "",
    stringsAsFactors = FALSE)
}

diagnostics <- list(); bridges <- list(); policies <- list(); ledger <- list()
if (hybrid_resume && file.exists(file.path(hybrid_out, "benchmark.csv")))
  diagnostics <- list(utils::read.csv(file.path(hybrid_out, "benchmark.csv"),
    stringsAsFactors = FALSE, na.strings = c("", "NA")))
if (hybrid_resume && file.exists(file.path(hybrid_out, "benchmark_bridges.csv")))
  bridges <- list(utils::read.csv(file.path(hybrid_out, "benchmark_bridges.csv"),
    stringsAsFactors = FALSE, na.strings = c("", "NA")))
if (hybrid_resume && file.exists(file.path(hybrid_out, "benchmark_policy_hooks.csv")))
  policies <- list(utils::read.csv(file.path(hybrid_out, "benchmark_policy_hooks.csv"),
    stringsAsFactors = FALSE, na.strings = c("", "NA")))
if (hybrid_resume && file.exists(file.path(hybrid_out, "benchmark_run_ledger.csv")))
  ledger <- list(utils::read.csv(file.path(hybrid_out, "benchmark_run_ledger.csv"),
    stringsAsFactors = FALSE, na.strings = c("", "NA")))

if (length(bridges) && !"case_id" %in% names(bridges[[1L]])) {
  bridges[[1L]]$case_id <- plan$case_id
  bridges[[1L]] <- bridges[[1L]][, c("case_id", setdiff(names(bridges[[1L]]), "case_id")),
    drop = FALSE]
}

for (i in seq_len(nrow(plan))) {
  case <- plan[i, , drop = FALSE]
  if (hybrid_resume && state$status[[i]] != "unstarted") {
    cat("Reusing completed state without HMC rerun: ", case$case_id, "\n", sep = "")
    next
  }
  state$status[[i]] <- "running"
  hybrid_write_checkpoint(state, diagnostics, bridges, policies, ledger)
  cat("Starting ", case$case_id, "\n", sep = "")
  flush.console()
  panel <- readRDS(file.path(hybrid_out, "data", paste0(case$design, ".rds")))
  result <- tryCatch(cga_run_one(panel, case, hybrid_out,
      estimator = hybrid_estimator_for(case)),
    error = function(e) hybrid_failure_result(panel, case,
      paste("runner/postprocessing error:", conditionMessage(e))))
  d <- hybrid_add_labels(result$diagnostics, case)
  diagnostics <- c(diagnostics, list(d))
  bridges <- c(bridges, list(hybrid_bridge_with_id(result$bridge, case$case_id)))
  policies <- c(policies, list(result$policy))
  state$status[[i]] <- d$status[[1L]]
  state$artifact[[i]] <- result$log_path
  ledger <- c(ledger, list(data.frame(case_id = case$case_id,
    date = as.character(Sys.Date()), stage = "fit", status = d$status[[1L]],
    evidence = ifelse(is.na(d$error[[1L]]), "", d$error[[1L]]),
    artifact = result$log_path, stringsAsFactors = FALSE)))
  hybrid_write_checkpoint(state, diagnostics, bridges, policies, ledger)
  cat("Finished ", case$case_id, ": ", d$status[[1L]], "\n", sep = "")
  flush.console()
  rm(panel, result, d)
}

diag <- do.call(rbind, diagnostics)
bridge <- do.call(rbind, bridges)
policy <- do.call(rbind, policies)
ledger_df <- do.call(rbind, ledger)
utils::write.csv(diag, file.path(hybrid_out, "benchmark.csv"), row.names = FALSE, na = "")
utils::write.csv(bridge, file.path(hybrid_out, "benchmark_bridges.csv"), row.names = FALSE, na = "")
utils::write.csv(policy, file.path(hybrid_out, "benchmark_policy_hooks.csv"), row.names = FALSE, na = "")
utils::write.csv(ledger_df, file.path(hybrid_out, "benchmark_run_ledger.csv"), row.names = FALSE, na = "")
comparison <- hybrid_compare(diag, ref)
hybrid_report(diag, bridge, ref, comparison)
cat("HYBRID RUN COMPLETE: ", hybrid_out, "\n", sep = "")
cat("Usable fits: ", sum(diag$usable), "/", nrow(diag), "; policy hooks remain pending by scope.\n", sep = "")
