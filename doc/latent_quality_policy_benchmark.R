#!/usr/bin/env Rscript

## Frozen latent-quality policy comparison runner.
##
## This file intentionally has a preparation phase and an execution gate:
##
##   Rscript run_latent_quality_policy_benchmark.R --prepare
##   Rscript run_latent_quality_policy_benchmark.R --run --approved
##
## --prepare reads only the fourteen archived fits and writes the paired draw
## plans before any equilibrium is evaluated. --run refuses to start without
## --approved, the frozen manifest, and all fourteen plans. It never calls
## bayes_est(), Stan, or an installer.

options(warn = 1)

policy_args <- commandArgs(trailingOnly = TRUE)
if (any(policy_args %in% c("--help", "-h"))) {
  cat("Usage: Rscript run_latent_quality_policy_benchmark.R --prepare\n",
      "       Rscript run_latent_quality_policy_benchmark.R --run --approved [--resume]\n",
      sep = "")
  quit(save = "no", status = 0L)
}
policy_phases <- policy_args[policy_args %in% c("--prepare", "--run")]
if (length(policy_phases) != 1L) stop("Choose exactly one phase: --prepare or --run")
policy_phase <- policy_phases[[1L]]
policy_resume <- "--resume" %in% policy_args
if (policy_resume && policy_phase != "--run") stop("--resume requires --run")
if (policy_phase == "--run" && !("--approved" %in% policy_args)) {
  stop("The policy run is gated; pass --approved after the independent review.")
}

policy_file <- commandArgs(trailingOnly = FALSE)
policy_file <- policy_file[startsWith(policy_file, "--file=")]
policy_script_dir <- if (length(policy_file)) {
  dirname(normalizePath(sub("^--file=", "", policy_file[[1L]]), mustWork = FALSE))
} else getwd()
policy_repo <- Sys.getenv("ANTITRUSTBAYES_POLICY_REPO",
                          "/mnt/d/projects/antitrustbayes")
policy_archive <- Sys.getenv(
  "ANTITRUSTBAYES_POLICY_ARCHIVE",
  file.path(policy_repo, "doc/statistical_toolkit_runs/20260916/latent-quality-policy"))
policy_archive <- normalizePath(policy_archive, mustWork = FALSE)
policy_protocol <- file.path(policy_repo, "doc/latent_quality_policy_benchmark_protocol.md")
policy_b_dir <- Sys.getenv(
  "ANTITRUSTBAYES_POLICY_B_FITS",
  file.path(policy_repo, "doc/collapsed_gaussian_runs/20260916",
            "cga-primary-20260916-configfix/stanfit"))
policy_d_dir <- Sys.getenv(
  "ANTITRUSTBAYES_POLICY_D_FITS",
  file.path(policy_repo, "doc/statistical_toolkit_runs/20260916/hybrid/fits"))
policy_b_diag <- Sys.getenv(
  "ANTITRUSTBAYES_POLICY_B_DIAGNOSTICS",
  file.path(policy_repo, "doc/collapsed_gaussian_runs/20260916",
            "cga-primary-20260916-configfix/benchmark.csv"))
policy_d_plan <- Sys.getenv(
  "ANTITRUSTBAYES_POLICY_D_PLAN",
  file.path(policy_repo, "doc/statistical_toolkit_runs/20260916/hybrid",
            "frozen/benchmark_plan.csv"))
policy_b_plan <- Sys.getenv(
  "ANTITRUSTBAYES_POLICY_B_PLAN",
  file.path(policy_repo, "doc/collapsed_gaussian_runs/20260916",
            "cga-primary-20260916-configfix/benchmark_plan.csv"))
policy_run_out <- Sys.getenv("ANTITRUSTBAYES_POLICY_RUN_OUT",
                             "/tmp/latent-quality-policy-20260916")
policy_lib <- Sys.getenv("ANTITRUSTBAYES_POLICY_LIB", "")
if (nzchar(policy_lib)) .libPaths(c(policy_lib, .libPaths()))

if (!file.exists(policy_protocol)) stop("Missing frozen protocol: ", policy_protocol)
if (!file.exists(policy_b_diag)) stop("Missing B diagnostics: ", policy_b_diag)
if (!file.exists(policy_b_plan)) stop("Missing B plan: ", policy_b_plan)
if (!file.exists(policy_d_plan)) stop("Missing D plan: ", policy_d_plan)
if (!dir.exists(policy_b_dir)) stop("Missing B fit directory: ", policy_b_dir)
if (!dir.exists(policy_d_dir)) stop("Missing D fit directory: ", policy_d_dir)

suppressPackageStartupMessages({
  library(antitrustBayes)
  library(antitrust)
})
if (!exists("posterior_policy_draws", mode = "function")) {
  stop("The installed antitrustBayes does not contain posterior_policy_draws().",
       " Install the reviewed implementation into a dedicated policy library.")
}

policy_digest <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  digest::digest(file = path, algo = "sha256")
}

## This is the exact deterministic seed helper used by the frozen benchmark.
policy_derived_seed <- function(seed, ...) {
  bits <- paste(..., collapse = "|")
  bytes <- utf8ToInt(bits)
  z <- as.double(seed) * 104729 +
    sum(as.double(seq_along(bytes)) * as.double(bytes))
  as.integer(1 + z %% 2147483646)
}

policy_write <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE, na = "")
}

policy_append <- function(x, path) {
  if (!nrow(x)) return(invisible(FALSE))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, path, sep = ",", row.names = FALSE,
                     col.names = !file.exists(path), append = file.exists(path),
                     na = "", qmethod = "double")
  invisible(TRUE)
}

policy_read_plan <- function() {
  b <- utils::read.csv(policy_b_plan, stringsAsFactors = FALSE,
                       na.strings = c("", "NA"))
  d <- utils::read.csv(policy_d_plan, stringsAsFactors = FALSE,
                       na.strings = c("", "NA"))
  bd <- utils::read.csv(policy_b_diag, stringsAsFactors = FALSE,
                        na.strings = c("", "NA"))
  b <- b[b$mode == "B", , drop = FALSE]
  bd <- bd[bd$mode == "B" & as.logical(bd$usable), , drop = FALSE]
  d <- d[d$mode == "D", , drop = FALSE]
  designs <- c("logit_bertrand_exact_271828",
               "logit_bertrand_noisy_271828",
               "logit_bertrand_noisy_314159",
               "logit_cournot_noisy_271828",
               "logit_BF_noisy_271828",
               "logit_BS_noisy_271828",
               "ces_bertrand_noisy_271828")
  if (!identical(sort(unique(b$design)), sort(designs)) || nrow(b) != 7L ||
      nrow(bd) != 7L || !all(bd$design %in% designs) ||
      !identical(sort(unique(d$design)), sort(designs)) || nrow(d) != 7L) {
    stop("Frozen B/D inputs do not contain exactly the seven requested designs.")
  }
  b <- b[match(designs, b$design), , drop = FALSE]
  d <- d[match(designs, d$design), , drop = FALSE]
  bd <- bd[match(designs, bd$design), , drop = FALSE]
  if (anyNA(b$design) || anyNA(d$design) || anyNA(bd$design))
    stop("Frozen policy design mapping is incomplete.")
  out <- data.frame(
    design = designs,
    seed = as.integer(b$seed),
    demand_model = as.character(b$demand_model),
    supply_model = as.character(b$supply_model),
    structure = as.character(b$structure),
    partition_prior = as.character(b$partition_prior),
    b_case_id = paste0(designs, "__B"),
    d_case_id = paste0(designs, "__D"),
    b_fit = file.path(policy_b_dir, paste0(designs, "__B.rds")),
    d_fit = file.path(policy_d_dir, paste0(designs, "__D.rds")),
    stringsAsFactors = FALSE
  )
  if (any(!file.exists(out$b_fit)) || any(!file.exists(out$d_fit)))
    stop("At least one requested B or D fit RDS is missing.")
  out$policy_seed <- vapply(seq_len(nrow(out)), function(i)
    policy_derived_seed(out$seed[[i]], out$design[[i]],
                        "latent_quality_policy_plan"), integer(1))
  out
}

policy_fit_check <- function(fit, mode, case) {
  if (!methods::is(fit, "BayesFit")) stop("Not a BayesFit: ", case)
  if (!identical(as.character(fit@market_side), "output"))
    stop("Policy protocol requires output fits: ", case)
  if (!identical(as.character(fit@demand_model), case$demand_model))
    stop("Demand model mismatch for ", case)
  if (mode == "B" && !identical(as.character(fit@metadata$margin_likelihood), "gaussian"))
    stop("B fit does not have Gaussian margins: ", case)
  if (mode == "D" && !identical(as.character(fit@metadata$margin_likelihood), "student_t"))
    stop("D fit does not have Student-t margins: ", case)
  if (!identical(as.character(fit@metadata$quality_parameterization), "collapsed"))
    stop("Policy protocol requires collapsed quality: ", case)
  invisible(TRUE)
}

policy_prepare <- function(cases) {
  plan_dir <- file.path(policy_archive, "plans")
  manifest_path <- file.path(policy_archive, "policy_manifest.csv")
  if (file.exists(manifest_path) || dir.exists(plan_dir) && length(list.files(plan_dir)))
    stop("Existing policy plans must be preserved; choose a new archive.")
  dir.create(plan_dir, recursive = TRUE, showWarnings = FALSE)
  rows <- vector("list", 2L * nrow(cases))
  z <- 0L
  for (i in seq_len(nrow(cases))) {
    case <- cases[i, , drop = FALSE]
    for (mode in c("B", "D")) {
      fit_path <- case[[paste0(tolower(mode), "_fit")]][[1L]]
      fit <- readRDS(fit_path)
      case_id <- case[[paste0(tolower(mode), "_case_id")]][[1L]]
      policy_fit_check(fit, mode, case)
      plan <- posterior_policy_draws(
        fit, n_draws = 100L, seed = case$policy_seed[[1L]])
      if (length(plan$draw_indices) != 100L)
        stop("A policy plan did not contain exactly 100 draws: ", case_id)
      path <- file.path(plan_dir, paste0(case_id, ".rds"))
      saveRDS(plan, path)
      z <- z + 1L
      rows[[z]] <- data.frame(
        case_id = case_id, design = case$design, mode = mode,
        policy_seed = case$policy_seed, selected_draws = length(plan$draw_indices),
        fit_path = normalizePath(fit_path), fit_sha256 = policy_digest(fit_path),
        plan_path = normalizePath(path), plan_sha256 = policy_digest(path),
        quality_source = plan$quality_source,
        quality_conditioning = plan$quality_conditioning,
        exact_k = isTRUE(plan$metadata$exact_k), stringsAsFactors = FALSE)
    }
  }
  plans <- do.call(rbind, rows)
  policy_write(plans, file.path(policy_archive, "plans_manifest.csv"))
  manifest <- data.frame(
    prepared_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    protocol_path = normalizePath(policy_protocol),
    protocol_sha256 = policy_digest(policy_protocol),
    b_diagnostics = normalizePath(policy_b_diag),
    b_diagnostics_sha256 = policy_digest(policy_b_diag),
    b_plan = normalizePath(policy_b_plan), b_plan_sha256 = policy_digest(policy_b_plan),
    d_plan = normalizePath(policy_d_plan), d_plan_sha256 = policy_digest(policy_d_plan),
    package_version = as.character(packageVersion("antitrustBayes")),
    policy_library = paste(.libPaths(), collapse = ";"),
    policy_status = "plans_saved_no_equilibrium_run", stringsAsFactors = FALSE)
  policy_write(manifest, manifest_path)
  invisible(plans)
}

policy_assert_frozen <- function(cases) {
  manifest_path <- file.path(policy_archive, "policy_manifest.csv")
  plans_manifest_path <- file.path(policy_archive, "plans_manifest.csv")
  if (!file.exists(manifest_path) || !file.exists(plans_manifest_path))
    stop("Run --prepare first; no policy manifest or plan manifest found.")
  man <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
  if (nrow(man) != 1L || !identical(man$protocol_sha256[[1L]], policy_digest(policy_protocol)) ||
      !identical(man$b_diagnostics_sha256[[1L]], policy_digest(policy_b_diag)) ||
      !identical(man$b_plan_sha256[[1L]], policy_digest(policy_b_plan)) ||
      !identical(man$d_plan_sha256[[1L]], policy_digest(policy_d_plan)))
    stop("A frozen policy input changed after plan preparation.")
  pm <- utils::read.csv(plans_manifest_path, stringsAsFactors = FALSE)
  if (nrow(pm) != 14L) stop("Expected fourteen frozen B/D plans.")
  for (i in seq_len(nrow(pm))) {
    if (!file.exists(pm$fit_path[[i]]) || !file.exists(pm$plan_path[[i]]) ||
        !identical(pm$fit_sha256[[i]], policy_digest(pm$fit_path[[i]])) ||
        !identical(pm$plan_sha256[[i]], policy_digest(pm$plan_path[[i]])))
      stop("A frozen fit or policy plan changed: ", pm$case_id[[i]])
  }
  invisible(pm)
}

policy_scenarios <- function(fit) {
  d <- fit@data
  labels <- paste0("p", seq_len(nrow(d)))
  owner <- as.character(d$group)
  weight_size <- if (identical(as.character(fit@demand_model), "ces")) {
    as.numeric(d$quantities) * as.numeric(d$prices)
  } else as.numeric(d$quantities)
  size <- tapply(weight_size, owner, sum)
  largest <- names(size)[order(-as.numeric(size), names(size))][seq_len(2L)]
  merged <- owner
  merged[merged == largest[[2L]]] <- largest[[1L]]
  costs <- setNames(rep(.10, nrow(d)), labels)
  list(
    merger = antitrust::counterfactual(ownership = merged),
    cost10 = antitrust::counterfactual(costs = costs),
    merger_cost10 = antitrust::counterfactual(ownership = merged, costs = costs)
  )
}

policy_welfare <- function(model) {
  if (is.null(model)) return(c(consumer_welfare_change = NA_real_,
                               producer_surplus_change = NA_real_,
                               total_welfare_change = NA_real_))
  tryCatch({
    consumer <- -as.numeric(antitrust::CV(model))
    producer <- sum(antitrust::calcProducerSurplus(model, preMerger = FALSE)) -
      sum(antitrust::calcProducerSurplus(model, preMerger = TRUE))
    c(consumer_welfare_change = consumer,
      producer_surplus_change = producer,
      total_welfare_change = consumer + producer)
  }, error = function(e) {
    c(consumer_welfare_change = NA_real_, producer_surplus_change = NA_real_,
      total_welfare_change = NA_real_)
  })
}

policy_draw_metrics <- function(value, fit, case, mode, anchor, policy_name) {
  rows <- list(); failures <- list(); z <- 0L; f <- 0L
  weights <- as.numeric(fit@data$quantities)
  weights <- weights / sum(weights)
  for (u in seq_along(value@market_ids)) {
    unit_id <- value@market_ids[[u]]
    pre <- value@price_pre[[unit_id]]
    post <- value@price_post[[unit_id]]
    pchg <- value@price_change[[unit_id]]
    spre <- value@shares_pre[[unit_id]]
    spost <- value@shares_post[[unit_id]]
    if (ncol(pre) != length(weights))
      stop("The frozen policy runner currently expects one market-year unit per fit.")
    sims <- value@sims[[unit_id]]
    latent_diag <- value@metadata$latent_diagnostics[[unit_id]]
    for (i in seq_len(nrow(pre))) {
      draw_id <- value@metadata$draw_indices[[i]]
      err <- value@failures[[unit_id]][[i]]
      complete <- all(is.finite(pre[i, ])) && all(is.finite(post[i, ]))
      if (!complete) {
        f <- f + 1L
        failures[[f]] <- data.frame(
          case_id = paste0(case$design, "__", mode), design = case$design,
          mode = mode, anchor = anchor, policy = policy_name, unit_id = unit_id,
          draw_position = i, posterior_draw_index = draw_id,
          error = ifelse(is.na(err), "solver returned no finite result", err),
          stringsAsFactors = FALSE)
      }
      values <- c(
        setNames(as.numeric(pchg[i, ]), paste0("price_change:p", seq_len(ncol(pchg)))),
        setNames(as.numeric(spost[i, ] - spre[i, ]),
                 paste0("share_change:p", seq_len(ncol(spost)))),
        fixed_weight_price_index_change = if (complete)
          as.numeric(post[i, ] %*% weights) / as.numeric(pre[i, ] %*% weights) - 1 else NA_real_,
        inside_share_pre = sum(spre[i, ], na.rm = FALSE),
        inside_share_post = sum(spost[i, ], na.rm = FALSE),
        inside_share_change = sum(spost[i, ] - spre[i, ], na.rm = FALSE)
      )
      welfare <- policy_welfare(sims[[i]])
      values <- c(values, welfare)
      ld <- if (length(latent_diag) >= i) latent_diag[[i]] else NULL
      values <- c(values,
        latent_pre_foc_max = if (is.null(ld)) NA_real_ else ld$max_abs_foc_at_latent_pre,
        latent_pre_price_discrepancy_max = if (is.null(ld)) NA_real_ else
          ld$max_abs_pre_price_discrepancy)
      status <- if (complete) "complete" else "solver_failure"
      for (metric in names(values)) {
        z <- z + 1L
        rows[[z]] <- data.frame(
          case_id = paste0(case$design, "__", mode), design = case$design,
          mode = mode, anchor = anchor, policy = policy_name, unit_id = unit_id,
          draw_position = i, posterior_draw_index = draw_id, metric = metric,
          value = as.numeric(values[[metric]]), status = status,
          error = ifelse(is.na(err), "", err), stringsAsFactors = FALSE)
      }
    }
  }
  list(metrics = if (length(rows)) do.call(rbind, rows) else data.frame(),
       failures = if (length(failures)) do.call(rbind, failures) else data.frame())
}

policy_summary <- function(draw_metrics) {
  if (!nrow(draw_metrics)) return(data.frame())
  keys <- c("case_id", "design", "mode", "anchor", "policy", "unit_id", "metric")
  grouping_keys <- lapply(draw_metrics[, keys, drop = FALSE], function(x)
    addNA(as.factor(x), ifany = TRUE))
  groups <- interaction(grouping_keys, drop = TRUE, lex.order = TRUE, sep = "\r")
  out <- lapply(split(draw_metrics, groups), function(x) {
    good <- x$value[is.finite(x$value)]
    requested <- nrow(x)
    q <- if (length(good)) stats::quantile(good, c(.05, .50, .95), names = FALSE) else rep(NA_real_, 3L)
    data.frame(x[1L, keys, drop = FALSE], requested_draws = requested,
               valid_draws = length(good), solver_failures = sum(x$status != "complete"),
               solver_failure_rate = mean(x$status != "complete"),
               status = if (all(x$status == "simulation_error")) "simulation_error" else
                 if (!length(good)) "metric_unavailable" else
                 if (length(good) < requested) "partial_failure" else "complete",
               mean = if (length(good)) mean(good) else NA_real_,
               sd = if (length(good) > 1L) stats::sd(good) else NA_real_,
               q05 = q[[1L]], median = q[[2L]], q95 = q[[3L]], stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

policy_simulation_error_records <- function(pm, anchor, policy_name, plan,
                                            error, elapsed) {
  requested <- seq_along(plan$draw_indices)
  data.frame(case_id = pm$case_id, design = pm$design, mode = pm$mode,
             anchor = anchor, policy = policy_name, unit_id = NA_character_,
             draw_position = requested,
             posterior_draw_index = as.integer(plan$draw_indices),
             metric = "simulation_error", value = NA_real_,
             status = "simulation_error", error = error,
             elapsed_seconds = elapsed, stringsAsFactors = FALSE)
}

policy_paired <- function(draw_metrics) {
  if (!nrow(draw_metrics)) return(data.frame())
  ## B and D are different posterior distributions.  The design, mode, and
  ## case identity must remain in the join key even when their draw positions
  ## and policy labels coincide.
  key <- c("case_id", "design", "mode", "policy", "unit_id", "draw_position",
           "posterior_draw_index", "metric")
  x <- draw_metrics[draw_metrics$anchor %in% c("observed", "latent_quality"), , drop = FALSE]
  wide <- reshape(x[, c(key, "anchor", "value", "status"), drop = FALSE],
                  idvar = key, timevar = "anchor", direction = "wide")
  names(wide) <- sub("^value\\.", "value_", names(wide))
  if (!all(c("value_observed", "value_latent_quality") %in% names(wide)))
    return(data.frame())
  status_observed <- if ("status.observed" %in% names(wide)) {
    wide[["status.observed"]]
  } else if ("status_observed" %in% names(wide)) wide[["status_observed"]] else character()
  status_latent <- if ("status.latent_quality" %in% names(wide)) {
    wide[["status.latent_quality"]]
  } else if ("status_latent_quality" %in% names(wide)) wide[["status_latent_quality"]] else character()
  wide$joint_success <- rep(FALSE, nrow(wide))
  if (length(status_observed) == nrow(wide) && length(status_latent) == nrow(wide)) {
    wide$joint_success <- status_observed == "complete" &
      status_latent == "complete" & is.finite(wide$value_observed) &
      is.finite(wide$value_latent_quality)
  }
  wide$difference_latent_minus_observed <- ifelse(
    wide$joint_success,
    wide$value_latent_quality - wide$value_observed, NA_real_)
  wide[, c(key, "value_observed", "value_latent_quality",
           "difference_latent_minus_observed", "joint_success"), drop = FALSE]
}

## Bookkeeping self-checks run during --prepare before any equilibrium call.
## They protect the policy denominator and the B/D pairing key independently of
## the actual archived fit results.
policy_runner_self_check <- function() {
  pm <- data.frame(case_id = "case__B", design = "case", mode = "B",
                   stringsAsFactors = FALSE)
  fake_plan <- list(draw_indices = seq_len(100L))
  err <- policy_simulation_error_records(pm, "observed", "merger", fake_plan,
                                         "synthetic top-level failure", 0)
  if (nrow(err) != 100L || !identical(as.integer(err$draw_position), seq_len(100L)) ||
      !identical(as.integer(err$posterior_draw_index), seq_len(100L)))
    stop("Runner bookkeeping self-check failed: top-level errors must retain 100 draws.")
  err_summary <- policy_summary(err)
  if (nrow(err_summary) != 1L || err_summary$requested_draws != 100L ||
      err_summary$solver_failures != 100L || err_summary$solver_failure_rate != 1 ||
      err_summary$status != "simulation_error")
    stop("Runner bookkeeping self-check failed: summary must retain top-level errors.")
  rows <- list(); z <- 0L
  for (case in c("case__B", "case__D")) for (mode in c("B", "D")) {
    if ((case == "case__B") != (mode == "B")) next
    for (anchor in c("observed", "latent_quality")) for (dp in 1:2) {
      z <- z + 1L
      rows[[z]] <- data.frame(case_id=case, design="case", mode=mode,
        policy="merger", unit_id="u", draw_position=dp,
        posterior_draw_index=dp, metric="m", anchor=anchor,
        value=if (anchor == "observed") 1 else 3, status="complete")
    }
  }
  paired <- policy_paired(do.call(rbind, rows))
  if (nrow(paired) != 4L || !all(paired$joint_success) ||
      !all(paired$difference_latent_minus_observed == 2) ||
      length(unique(paste(paired$case_id, paired$mode))) != 2L)
    stop("Runner bookkeeping self-check failed: B/D pairing keys are not isolated.")
  rows <- do.call(rbind, rows)
  rows$status[rows$anchor == "latent_quality"] <- "solver_failure"
  paired <- policy_paired(rows)
  if (any(paired$joint_success) || any(is.finite(paired$difference_latent_minus_observed)))
    stop("Runner bookkeeping self-check failed: incomplete pairs must be excluded.")
  invisible(TRUE)
}

policy_plan <- policy_read_plan()
policy_runner_self_check()
if (policy_phase == "--prepare") {
  policy_prepare(policy_plan)
  cat("POLICY PLANS SAVED; no equilibrium calculations executed: ",
      policy_archive, "\n", sep = "")
  quit(save = "no", status = 0L)
}

policy_pm <- policy_assert_frozen(policy_plan)
dir.create(policy_run_out, recursive = TRUE, showWarnings = FALSE)
progress_path <- file.path(policy_archive, "policy_progress.csv")
summary_path <- file.path(policy_archive, "policy_summary.csv")
draw_path <- file.path(policy_archive, "policy_draw_metrics.csv")
failure_path <- file.path(policy_archive, "policy_failures.csv")
if (!policy_resume && any(file.exists(c(progress_path, summary_path, draw_path, failure_path))))
  stop("Existing policy results must be preserved; use --resume or a new archive.")
if (policy_resume && !file.exists(progress_path)) stop("--resume requires policy_progress.csv")
state <- if (policy_resume) utils::read.csv(progress_path, stringsAsFactors = FALSE) else
  data.frame(case_id = character(), anchor = character(), policy = character(),
             status = character(), result_path = character(), stringsAsFactors = FALSE)

for (i in seq_len(nrow(policy_pm))) {
  pm <- policy_pm[i, , drop = FALSE]
  case <- policy_plan[policy_plan$design == pm$design[[1L]], , drop = FALSE]
  if (nrow(case) != 1L) stop("Missing design in policy plan: ", pm$design[[1L]])
  fit <- readRDS(pm$fit_path[[1L]])
  policy_fit_check(fit, pm$mode[[1L]], case)
  plan <- readRDS(pm$plan_path[[1L]])
  scenarios <- policy_scenarios(fit)
  for (anchor in c("observed", "latent_quality")) {
    for (policy_name in names(scenarios)) {
      run_id <- paste(pm$case_id[[1L]], anchor, policy_name, sep = "__")
      done <- nrow(state) && any(state$case_id == pm$case_id[[1L]] &
                                 state$anchor == anchor & state$policy == policy_name &
                                 state$status == "complete")
      if (policy_resume && done) next
      result_path <- file.path(policy_run_out, paste0(run_id, ".rds"))
      started <- proc.time()[["elapsed"]]
      value <- tryCatch(
        antitrustBayes::bayes_sim(
          fit, ownerPost = scenarios[[policy_name]], draw_plan = plan,
          demand_anchoring = anchor, labels = paste0("p", seq_len(nrow(fit@data))),
          insideSize = 1, mc.cores = 1L),
        error = function(e) e)
      elapsed <- proc.time()[["elapsed"]] - started
      if (inherits(value, "condition")) {
        ## A top-level failure still represents all requested plan draws. Keep
        ## every row so summary denominators and paired failure accounting are
        ## 100 draws, rather than one synthetic error record.
        rec <- policy_simulation_error_records(
          pm, anchor, policy_name, plan, conditionMessage(value), elapsed)
        policy_append(rec, draw_path)
        policy_append(rec[, c("case_id", "design", "mode", "anchor", "policy",
                              "unit_id", "draw_position", "posterior_draw_index", "error")],
                      failure_path)
        state <- rbind(state, data.frame(case_id = pm$case_id, anchor = anchor,
                                         policy = policy_name, status = "simulation_error",
                                         result_path = "", stringsAsFactors = FALSE))
        policy_write(state, progress_path)
        next
      }
      ## The full BayesSim object is saved before metric extraction.
      saveRDS(value, result_path)
      extracted <- policy_draw_metrics(value, fit, case, pm$mode, anchor, policy_name)
      extracted$metrics$elapsed_seconds <- elapsed
      policy_append(extracted$metrics, draw_path)
      policy_append(extracted$failures, failure_path)
      state <- rbind(state, data.frame(case_id = pm$case_id, anchor = anchor,
                                       policy = policy_name, status = "complete",
                                       result_path = result_path, stringsAsFactors = FALSE))
      policy_write(state, progress_path)
      cat("Finished ", run_id, "\n", sep = "")
      flush.console()
    }
  }
}

all_draw <- utils::read.csv(draw_path, stringsAsFactors = FALSE,
                            na.strings = c("", "NA"))
policy_write(policy_summary(all_draw), summary_path)
policy_write(policy_paired(all_draw), file.path(policy_archive, "policy_paired.csv"))
cat("POLICY RUN COMPLETE: ", policy_archive, "\n", sep = "")
