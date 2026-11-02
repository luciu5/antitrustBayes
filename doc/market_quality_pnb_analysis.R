#!/usr/bin/env Rscript
# Bounded PNB postprocessing for the frozen market-quality recovery.
#
# This independently reviewed driver reads checkpointed serialized fits/results
# and calls only the installed antitrustBayes,
# antitrust, and posterior APIs.  It never fits a model, changes priors or
# solver settings, or edits the production checkout.
options(stringsAsFactors = FALSE)

mq_args <- commandArgs(trailingOnly = TRUE)
mq_repo <- normalizePath(Sys.getenv("MQ_REPO", "/mnt/d/projects/antitrustbayes"),
                         mustWork = FALSE)
mq_recovery <- file.path(mq_repo,
  "doc/statistical_toolkit_runs/20260916/performance-audit/recovery")
mq_default_root <- file.path(mq_recovery, "market-quality", "pnb-archived-retry-01")
mq_root <- normalizePath(Sys.getenv("MQ_PNB_ROOT", mq_default_root),
                         mustWork = FALSE)
mq_bundled_root <- normalizePath(Sys.getenv("MQ_PNB_BUNDLED_ROOT",
  file.path(mq_recovery, "market-quality", "pnb-bundled")), mustWork = FALSE)
mq_output <- normalizePath(Sys.getenv("MQ_PNB_ANALYSIS_ROOT",
  file.path(mq_root, "analysis")), mustWork = FALSE)
dir.create(mq_output, recursive = TRUE, showWarnings = FALSE)

mq_root_label <- basename(mq_root)
mq_plan_path <- file.path(mq_root, "frozen-plan.rds")
mq_data_path <- file.path(mq_root, "frozen-data.rds")
mq_limitations_path <- file.path(mq_output, "limitations.md")

mq_stop <- function(...) stop(paste0(...), call. = FALSE)
mq_now <- function() format(Sys.time(), tz = "UTC", usetz = TRUE)
mq_tmp_path <- function(path) paste0(path, ".tmp-", Sys.getpid())
mq_save <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- mq_tmp_path(path)
  saveRDS(x, tmp)
  if (!file.rename(tmp, path)) mq_stop("Unable to install output: ", path)
  invisible(path)
}
mq_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- mq_tmp_path(path)
  write.csv(x, tmp, row.names = FALSE, na = "")
  if (!file.rename(tmp, path)) mq_stop("Unable to install output: ", path)
  invisible(path)
}
mq_append_warning <- function(warnings, w) {
  c(warnings, conditionMessage(w))
}
mq_capture <- function(expr) {
  warnings <- character()
  value <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    warnings <<- mq_append_warning(warnings, w)
    invokeRestart("muffleWarning")
  }), error = function(e) e)
  list(value = value, warnings = unique(warnings),
       error = if (inherits(value, "error")) conditionMessage(value) else "")
}
mq_atomic_write_text <- function(lines, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- mq_tmp_path(path)
  writeLines(lines, tmp, useBytes = TRUE)
  if (!file.rename(tmp, path)) mq_stop("Unable to install output: ", path)
  invisible(path)
}
mq_bind_rows <- function(rows) {
  rows <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L, rows)
  if (!length(rows)) return(data.frame())
  columns <- unique(unlist(lapply(rows, names), use.names = FALSE))
  aligned <- lapply(rows, function(x) {
    for (nm in setdiff(columns, names(x))) x[[nm]] <- NA
    x[, columns, drop = FALSE]
  })
  do.call(rbind, aligned)
}
mq_complete_evidence <- function(table, expected, status) {
  is.data.frame(table) && identical(status, "complete") &&
    all(c("hypothesis", "logml_restricted", "posterior_probability") %in% names(table)) &&
    nrow(table) == length(expected) && !anyNA(table$hypothesis) &&
    !anyDuplicated(table$hypothesis) &&
    identical(sort(as.character(table$hypothesis)), sort(as.character(expected))) &&
    all(is.finite(table$logml_restricted)) &&
    all(is.finite(table$posterior_probability)) &&
    all(table$posterior_probability >= 0) &&
    abs(sum(table$posterior_probability) - 1) <= 1e-8
}
mq_validate_policy_pool <- function(arrays, case, reference_indices = NULL,
                                    expected_count = 200L) {
  indices <- arrays$draw_indices
  if (!is.numeric(indices) || length(indices) != expected_count ||
      any(!is.finite(indices)) || any(indices < 1) ||
      any(indices != as.integer(indices)) || anyDuplicated(indices))
    mq_stop("Malformed policy-pool posterior indices for ", case)
  indices <- as.integer(indices)
  if (!is.null(reference_indices) && !identical(indices, as.integer(reference_indices)))
    mq_stop("Policy-pool posterior draw grids differ across models: ", case)
  if (!identical(arrays$case, case) || !identical(arrays$anchor, "observed") ||
      !identical(arrays$product_labels, mq_groups))
    mq_stop("Policy pool has a different model, anchor or product order: ", case)
  x <- arrays$long
  fields <- c("case", "anchor", "draw_position", "product_position", "draw_index",
              "product_index", "product_label", "price_change", "failure")
  if (!is.data.frame(x) || !all(fields %in% names(x)) ||
      nrow(x) != length(indices) * length(mq_groups) ||
      anyDuplicated(x[c("draw_index", "product_index")]))
    mq_stop("Policy pool must have one row per draw/product: ", case)
  integer_positions <- function(z, maximum) is.numeric(z) &&
    all(is.finite(z)) && all(z == as.integer(z)) && all(z >= 1 & z <= maximum)
  if (!integer_positions(x$draw_position, length(indices)) ||
      !integer_positions(x$product_index, length(mq_groups)) ||
      !identical(x$product_position, x$product_index) ||
      !identical(as.integer(x$draw_index), indices[x$draw_position]) ||
      !identical(as.character(x$product_label), mq_groups[x$product_index]) ||
      anyNA(x$case) || anyNA(x$anchor) ||
      any(x$case != case) || any(x$anchor != "observed"))
    mq_stop("Policy pool violates the frozen draw/product indexing: ", case)
  if (!is.matrix(arrays$price_change) || !is.numeric(arrays$price_change) ||
      !identical(dim(arrays$price_change), c(length(indices), length(mq_groups))) ||
      !is.character(arrays$failures) || length(arrays$failures) != length(indices) ||
      any(!is.na(arrays$failures) & !nzchar(arrays$failures)) ||
      !identical(as.character(x$failure), arrays$failures[x$draw_position]) ||
      !identical(as.numeric(x$price_change),
        as.numeric(arrays$price_change[cbind(x$draw_position, x$product_index)])))
    mq_stop("Policy pool does not retain its saved values and failure labels: ", case)
  if (any(!is.finite(x$price_change[is.na(x$failure)])))
    mq_stop("A policy-pool cell marked successful is nonfinite: ", case)
  invisible(indices)
}

mq_usage <- function() {
  cat(paste(
    "Usage:",
    "  market_quality_pnb_analysis.R --summary",
    "  market_quality_pnb_analysis.R --bridge=canonical|rich",
    "  market_quality_pnb_analysis.R --policy=<case>",
    "  market_quality_pnb_analysis.R --bma",
    "",
    "MQ_PNB_ROOT defaults to recovery/market-quality/pnb-archived-retry-01.",
    "Set MQ_PNB_ROOT to pnb-bundled to analyze the separately labeled bundle.",
    sep = "\n"), "\n")
}

mq_action <- mq_args[grepl("^(--summary|--bridge=|--policy=|--bma$)", mq_args)]
if (length(mq_action) != 1L) {
  mq_usage()
  mq_stop("Exactly one explicit postprocessing action is required.")
}

suppressPackageStartupMessages({
  library(antitrustBayes)
  library(antitrust)
  library(posterior)
  library(rstan)
})

if (!file.exists(mq_plan_path) || !file.exists(mq_data_path)) {
  mq_stop("Frozen PNB inputs are missing under MQ_PNB_ROOT: ", mq_root)
}
mq_plan <- readRDS(mq_plan_path)
mq_data <- readRDS(mq_data_path)
if (!is.list(mq_plan) || !is.data.frame(mq_plan$models) ||
    !all(c("case", "supply", "structure") %in% names(mq_plan$models))) {
  mq_stop("frozen-plan.rds has no supported model manifest.")
}
if (!identical(as.integer(mq_plan$k_support), 1:5) ||
    !identical(as.numeric(mq_plan$k_prior), rep(.2, 5L)))
  mq_stop("Frozen exact ordered-k support/prior is not 1:5 with uniform mass.")
mq_models <- mq_plan$models
mq_all_cases <- as.character(mq_models$case)
if (!identical(mq_all_cases,
               c("B", "MC", "C", "A", "BF", "BS", "CF", "CS"))) {
  mq_stop("The frozen PNB model order/support changed unexpectedly.")
}
if (!is.data.frame(mq_data) || !all(c("group", "prices", "quantities", "margins",
                                      "margin_observed") %in% names(mq_data))) {
  mq_stop("frozen-data.rds has no supported PNB data contract.")
}
mq_groups <- as.character(mq_data$group)
if (anyNA(mq_groups) || any(!nzchar(mq_groups)) || anyDuplicated(mq_groups))
  mq_stop("PNB groups must be unique, non-empty character IDs.")
mq_acquirer <- as.character(attr(mq_data, "acquirer"))
mq_target <- as.character(attr(mq_data, "target"))
if (length(mq_acquirer) != 1L || length(mq_target) != 1L ||
    !mq_acquirer %in% mq_groups || !mq_target %in% mq_groups ||
    identical(mq_acquirer, mq_target)) {
  mq_stop("Frozen PNB data do not contain distinct validated acquirer/target groups.")
}
mq_bundled_data_path <- file.path(mq_bundled_root, "frozen-data.rds")
mq_bundled_groups <- NULL
if (file.exists(mq_bundled_data_path)) {
  mq_bundled_data <- readRDS(mq_bundled_data_path)
  mq_bundled_groups <- as.character(mq_bundled_data$group)
  if (!mq_acquirer %in% mq_bundled_groups || !mq_target %in% mq_bundled_groups)
    mq_stop("Acquirer/target are not both present in the bundled PNB group map.")
}
if (!isTRUE(antitrustBayes:::.antitrust_bayes_model_spec_version == 10L))
  mq_stop("Installed antitrustBayes is not specification 10.")

# The first archived BF attempt failed in driver validation before HMC began.
# Keep this operational configuration failure in its own ledger; it is not a
# sampler diagnostic and must never be relabeled as an unusable BF posterior.
mq_original_attempt <- data.frame(
  attempt_id = "pnb-archived-BF-01",
  source_root = file.path(mq_recovery, "market-quality", "pnb-archived"),
  case = "BF", stage = "pre_hmc_configuration", status = "failed",
  fit_present = FALSE, hmc_started = FALSE, usable_posterior = NA,
  configuration = "use_cutoff=FALSE; structure=fringe; ordered_k support=1:5",
  error = "use_cutoff = FALSE conflicts with structure = 'fringe'.",
  retained_artifact = "old_root/BF-started.rds and BF-result.rds",
  stringsAsFactors = FALSE)
mq_write_csv(mq_original_attempt, file.path(mq_output, "attempt-ledger.csv"))
cat("Retained pre-HMC configuration attempt:", mq_original_attempt$attempt_id,
    "status=failed; sampler not started; see attempt-ledger.csv\n")

mq_plan_prior <- mq_plan$prior
if (is.null(names(mq_plan_prior))) names(mq_plan_prior) <- mq_all_cases
mq_rich_prior <- c(MC = .25, B = 1 / 12, C = 1 / 12, A = .25,
                   BF = 1 / 12, CF = 1 / 12, BS = 1 / 12, CS = 1 / 12)
if (!identical(as.numeric(mq_plan_prior[mq_all_cases]),
               as.numeric(mq_rich_prior[mq_all_cases]))) {
  mq_stop("The frozen rich PNB prior differs from the required explicit prior.")
}
mq_canonical_prior <- c(MC = .25, B = .25, C = .25, A = .25)

mq_case_model_row <- function(case) {
  x <- mq_models[mq_models$case == case, , drop = FALSE]
  if (nrow(x) != 1L) mq_stop("Unknown PNB case: ", case)
  x
}
mq_fit_candidates <- function(case) {
  tag <- switch(case, MC = "monopolistic", case)
  c(file.path(mq_root, paste0(case, "-fit.rds")),
    file.path(mq_root, paste0("pnb_ab_fit_", tag, "_rich_main.rds")),
    file.path(mq_root, paste0("pnb_ab_fit_", tag, ".rds")))
}
mq_result_candidates <- function(case) {
  tag <- switch(case, MC = "monopolistic", case)
  c(file.path(mq_root, paste0(case, "-result.rds")),
    file.path(mq_root, paste0("pnb_ab_result_", tag, "_rich_main.rds")),
    file.path(mq_root, paste0("pnb_ab_result_", tag, ".rds")))
}
mq_first_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit)) hit[[1L]] else NA_character_
}
mq_result_error <- function(result) {
  if (!is.list(result) || is.null(result$error) || !length(result$error) ||
      is.na(result$error[[1L]])) "" else as.character(result$error[[1L]])
}
mq_result_num <- function(result, field) {
  z <- if (is.list(result)) result[[field]] else NULL
  if (!length(z) || !is.finite(as.numeric(z[[1L]]))) NA_real_ else as.numeric(z[[1L]])
}
mq_result_warnings <- function(result) {
  z <- if (is.list(result)) result$warnings else NULL
  if (!length(z)) character() else as.character(z)
}
mq_load_case <- function(case) {
  fp <- mq_first_existing(mq_fit_candidates(case))
  rp <- mq_first_existing(mq_result_candidates(case))
  fit <- NULL; result <- NULL; errors <- character()
  if (!is.na(fp)) fit <- tryCatch(readRDS(fp), error = function(e) {
    errors <<- c(errors, paste0("fit: ", conditionMessage(e))); NULL
  })
  if (!is.na(rp)) result <- tryCatch(readRDS(rp), error = function(e) {
    errors <<- c(errors, paste0("result: ", conditionMessage(e))); NULL
  })
  list(case = case, fit = fit, result = result, fit_path = fp,
       result_path = rp, errors = unique(errors))
}

mq_plan_row <- function(case) {
  p <- mq_case_model_row(case)
  list(case = case, supply = as.character(p$supply),
       structure = as.character(p$structure))
}
mq_sampler_diagnostics <- function(fit, case, result = NULL) {
  out <- list(case = case, fit = fit, parameter_table = data.frame(),
              posterior_k = NULL, posterior_k_error = "")
  sf <- fit@stanfit
  sampler <- rstan::get_sampler_params(sf, inc_warmup = FALSE)
  divs <- sum(vapply(sampler, function(x) sum(x[, "divergent__"]), numeric(1)))
  max_td <- max(vapply(sampler, function(x) max(x[, "treedepth__"]), numeric(1)))
  td_hits <- sum(vapply(sampler, function(x) sum(x[, "treedepth__"] >= 12), numeric(1)))
  bfmi <- vapply(sampler, function(x) {
    e <- x[, "energy__"]
    if (length(e) < 2L || !is.finite(var(e)) || var(e) <= 0) return(NA_real_)
    mean(diff(e)^2) / var(e)
  }, numeric(1))
  raw <- tryCatch(antitrustBayes:::.bayes_bridge_raw_parameter_names(sf),
                  error = function(e) character())
  model_pars <- as.character(sf@model_pars)
  raw <- intersect(raw, model_pars)
  positive_width <- function(nm) {
    dimensions <- sf@par_dims[[nm]]
    !is.null(dimensions) && (!length(dimensions) || prod(dimensions) > 0L)
  }
  raw <- raw[vapply(raw, positive_width, logical(1))]
  if (!length(raw)) mq_stop("Cannot identify sampled coordinates for diagnostics: ", case)
  economic <- intersect(c("alpha_structural_market", "gamma_structural_market",
    "s0_prob", "b_market", "sigma_b_strat", "sigma_share", "sigma_margin",
    "rho", "sigma_log_a", "sigma_s0_pool", "mu_s0_pool", "b_quality_all"),
    model_pars)
  economic <- economic[vapply(economic, positive_width, logical(1))]
  selected <- unique(c(raw, economic))
  arr <- rstan::extract(sf, pars = selected, permuted = FALSE,
                        inc_warmup = FALSE)
  darr <- posterior::as_draws_array(arr)
  tab <- as.data.frame(posterior::summarise_draws(darr))
  mat <- posterior::as_draws_matrix(darr)
  base <- sub("\\[.*$", "", tab$variable)
  raw_base <- unique(sub("\\[.*$", "", raw))
  finite_coord <- vapply(tab$variable, function(nm) {
    z <- as.numeric(mat[, nm]); all(is.finite(z))
  }, logical(1))
  constant <- vapply(tab$variable, function(nm) {
    z <- as.numeric(mat[, nm]); z <- z[is.finite(z)]
    length(z) > 0L && length(unique(z)) == 1L
  }, logical(1))
  tab$base_variable <- base
  tab$estimated_coordinate <- base %in% raw_base
  tab$economic_output <- base %in% unique(sub("\\[.*$", "", economic))
  tab$deterministic_output <- constant & !tab$estimated_coordinate
  tab$estimated_coordinate_failure <- tab$estimated_coordinate &
    (!finite_coord | constant)
  tab$finite_draws <- finite_coord
  out$parameter_table <- tab
  eligible <- !tab$deterministic_output
  finite_summary <- finite_coord & is.finite(tab$rhat) &
    is.finite(tab$ess_bulk) & is.finite(tab$ess_tail) & !constant
  all_eligible_valid <- any(eligible) && all(finite_summary[eligible])
  out$divergences <- as.integer(divs)
  out$treedepth_hits <- as.integer(td_hits)
  out$maximum_treedepth <- as.integer(max(vapply(sampler, function(x)
    max(x[, "treedepth__"]), numeric(1))))
  out$minimum_BFMI <- if (all(is.na(bfmi))) NA_real_ else min(bfmi, na.rm = TRUE)
  out$maximum_Rhat <- if (all_eligible_valid) max(tab$rhat[eligible]) else Inf
  out$minimum_bulk_ESS <- if (all_eligible_valid) min(tab$ess_bulk[eligible]) else 0
  out$minimum_tail_ESS <- if (all_eligible_valid) min(tab$ess_tail[eligible]) else 0
  out$estimated_coordinate_failures <- sum(tab$estimated_coordinate_failure)
  out$mean_step_size <- mean(vapply(sampler, function(x) mean(x[, "stepsize__"]), numeric(1)))
  out$sampler_failure <- !(out$divergences == 0L && out$treedepth_hits == 0L &&
    is.finite(out$minimum_BFMI) && out$minimum_BFMI > .3 &&
    out$maximum_Rhat <= 1.01 && out$minimum_bulk_ESS >= 100 &&
    out$minimum_tail_ESS >= 100 && out$estimated_coordinate_failures == 0L)
  out$wall_seconds <- mq_result_num(result, "elapsed")
  out$native_seconds <- mq_result_num(result, "native_seconds")
  out
}

mq_summary <- function() {
  rows <- vector("list", length(mq_all_cases)); parameter_rows <- list()
  k_rows <- list()
  for (j in seq_along(mq_all_cases)) {
    case <- mq_all_cases[[j]]; p <- mq_plan_row(case); z <- mq_load_case(case)
    result <- z$result; fit <- z$fit
    row <- data.frame(case = case, supply = p$supply, structure = p$structure,
      status = "not_run", fit_present = !is.null(fit), result_present = !is.null(result),
      wall_seconds = NA_real_, native_seconds = NA_real_, divergences = NA_integer_,
      treedepth_hits = NA_integer_, maximum_treedepth = NA_integer_, minimum_BFMI = NA_real_, maximum_Rhat = NA_real_,
      minimum_bulk_ESS = NA_real_, minimum_tail_ESS = NA_real_, mean_step_size = NA_real_,
      estimated_coordinate_failures = NA_integer_, sampler_failure = NA,
      posterior_k_status = "not_applicable", error = paste(z$errors, collapse = " | "),
      warnings = "", stringsAsFactors = FALSE)
    if (!is.null(result)) {
      row$wall_seconds <- mq_result_num(result, "elapsed")
      row$native_seconds <- mq_result_num(result, "native_seconds")
      row$error <- paste(c(row$error, mq_result_error(result)), collapse = " | ")
      row$warnings <- paste(unique(mq_result_warnings(result)), collapse = " | ")
      row$status <- if (nzchar(mq_result_error(result))) "failed" else "completed"
    } else if (file.exists(file.path(mq_root, paste0(case, "-started.rds")))) {
      row$status <- "unfinished"
    }
    if (!is.null(fit) && inherits(fit, "BayesFit") &&
        !nzchar(paste(z$errors, collapse = ""))) {
      d <- tryCatch(mq_sampler_diagnostics(fit, case, result), error = function(e) e)
      if (!inherits(d, "error")) {
        row$divergences <- d$divergences; row$treedepth_hits <- d$treedepth_hits
        row$maximum_treedepth <- d$maximum_treedepth
        row$minimum_BFMI <- d$minimum_BFMI; row$maximum_Rhat <- d$maximum_Rhat
        row$minimum_bulk_ESS <- d$minimum_bulk_ESS; row$minimum_tail_ESS <- d$minimum_tail_ESS
        row$mean_step_size <- d$mean_step_size
        row$estimated_coordinate_failures <- d$estimated_coordinate_failures
        row$sampler_failure <- d$sampler_failure
        if (row$status == "completed") row$status <- if (d$sampler_failure)
          "diagnostic_gate_failed" else "usable"
        parameter_rows[[length(parameter_rows) + 1L]] <- cbind(
          data.frame(case = case, stringsAsFactors = FALSE), d$parameter_table)
        mq_write_csv(d$parameter_table,
          file.path(mq_output, paste0("parameters-", case, ".csv")))
        structure_name <- p$structure
        if (structure_name %in% c("fringe", "stackelberg")) {
          pk <- tryCatch(posterior_k(fit), error = function(e) e)
          if (inherits(pk, "error")) {
            row$posterior_k_status <- "error"
            k_rows[[length(k_rows) + 1L]] <- data.frame(
              case = case, status = "error", error = conditionMessage(pk),
              stringsAsFactors = FALSE)
          } else {
            row$posterior_k_status <- "available"
            k_rows[[length(k_rows) + 1L]] <- cbind(
              data.frame(case = case, status = "available", stringsAsFactors = FALSE), pk)
            mq_write_csv(pk, file.path(mq_output, paste0("posterior-k-", case, ".csv")))
          }
        }
      } else {
        row$status <- if (row$status == "completed") "diagnostic_error" else row$status
        row$error <- paste(c(row$error, conditionMessage(d)), collapse = " | ")
      }
    } else if (!is.null(fit) && !inherits(fit, "BayesFit")) {
      row$status <- "failed"
      row$error <- paste(c(row$error, "serialized fit is not a BayesFit"), collapse = " | ")
    }
    rows[[j]] <- row
  }
  summary <- do.call(rbind, rows)
  parameters <- if (length(parameter_rows)) do.call(rbind, parameter_rows) else
    data.frame(case = character(), variable = character())
  posterior_k_table <- if (length(k_rows)) mq_bind_rows(k_rows) else
    data.frame(case = character(), status = character(), error = character())
  mq_write_csv(summary, file.path(mq_output, "model-summary.csv"))
  mq_write_csv(parameters, file.path(mq_output, "parameters-all.csv"))
  mq_write_csv(posterior_k_table, file.path(mq_output, "posterior-k-all.csv"))
  limitations <- c(
    "# PNB postprocessing limitations", "",
    paste0("Root label: `", mq_root_label, "` (`", mq_root, "`)."),
    "The driver consumes serialized fits/results only; it never runs HMC.",
    "Gaussian residual degrees of freedom are inactive in the frozen plan and are not estimated.",
    "Constant estimated coordinates remain in each parameter CSV and count as coordinate failures.",
    "The archived root is a tracked-input reconstruction whose historical rich_main summaries match; the original serialized process input is unavailable.",
    paste0("The bundled root is a separately labeled application at `", mq_bundled_root,
      "` with different margin inputs; it must not be pooled with the archived root."),
    "Bridge evidence is provisional whenever sampler diagnostics fail and is incomplete whenever any positive-prior model evidence is unavailable.",
    "Policy outputs retain all draw/product cells and failure labels; valid-cell summaries state their conditioning explicitly.",
    "Observed-anchor and latent-quality policy results use one shared posterior_policy_draws plan per fit. Auction latent-quality policy is explicitly unsupported.",
    "BMA is emitted only for a complete eight-model observed-anchor family and never normalizes over a successful subset."
  )
  mq_atomic_write_text(limitations, mq_limitations_path)
  invisible(summary)
}

mq_universe <- function(mode) {
  if (identical(mode, "canonical")) {
    fits <- list(MC = "MC", B = "B", C = "C", A = "A")
    return(bayes_model_universe(fits, mq_canonical_prior[names(fits)]))
  }
  fits <- list(MC = "MC", B = "B", C = "C", A = "A",
               BF = list(fit = "BF", k = "interior"),
               CF = list(fit = "CF", k = "interior"),
               BS = list(fit = "BS", k = "interior"),
               CS = list(fit = "CS", k = "interior"))
  bayes_model_universe(fits, mq_rich_prior[names(fits)])
}
mq_bridge <- function(mode) {
  expected <- if (identical(mode, "canonical")) c("MC", "B", "C", "A") else mq_all_cases
  loaded <- setNames(lapply(expected, mq_load_case), expected)
  fit_present <- vapply(loaded, function(z) !is.null(z$fit) && inherits(z$fit, "BayesFit"), logical(1))
  result_error <- vapply(loaded, function(z) {
    is.null(z$result) || nzchar(mq_result_error(z$result))
  }, logical(1))
  failure_rows <- do.call(rbind, lapply(expected, function(case) {
    z <- loaded[[case]]
    sampler_failure <- NA
    if (fit_present[[case]]) sampler_failure <- tryCatch(
      mq_sampler_diagnostics(z$fit, case, z$result)$sampler_failure,
      error = function(e) NA)
    data.frame(model = case, fit_present = fit_present[[case]],
      result_complete = !result_error[[case]], sampler_failure = sampler_failure,
      error = paste(z$errors, collapse = " | "), stringsAsFactors = FALSE)
  }))
  mq_write_csv(failure_rows, file.path(mq_output, paste0("bridge-", mode, "-failures.csv")))
  settings <- list(mode = mode, repetitions = 3L, maxiter = 5000L,
                   prior_weights = if (mode == "canonical") mq_canonical_prior else mq_rich_prior,
                   root = mq_root, started = mq_now())
  if (any(!fit_present) || any(result_error)) {
    artifact <- list(settings = settings, status = "incomplete_evidence",
                     error = "One or more positive-prior model fits/results are missing or failed.",
                     failures = failure_rows, comparison = NULL, completed = mq_now())
    mq_save(artifact, file.path(mq_output, paste0("bridge-", mode, "-result.rds")))
    return(invisible(artifact))
  }
  fits <- lapply(loaded, `[[`, "fit")
  names(fits) <- expected
  captured <- mq_capture(bayes_compare(fits, universe = mq_universe(mode), bridge = TRUE,
                                        bridge_repetitions = 3L, bridge_maxiter = 5000L))
  if (nzchar(captured$error)) {
    artifact <- list(settings = settings, status = "bridge_failed", error = captured$error,
                     warnings = captured$warnings, failures = failure_rows,
                     comparison = NULL, completed = mq_now())
    mq_save(artifact, file.path(mq_output, paste0("bridge-", mode, "-result.rds")))
    return(invisible(artifact))
  }
  comp <- captured$value
  provisional <- any(is.na(failure_rows$sampler_failure) |
                       failure_rows$sampler_failure)
  tab <- as.data.frame(comp)
  universe_status <- attr(comp, "status", exact = TRUE)
  evidence_complete <- mq_complete_evidence(tab, expected, universe_status)
  tab$evidence_provisional_sampler_failure <- provisional
  tab$evidence_complete <- evidence_complete
  mq_write_csv(tab, file.path(mq_output, paste0("bridge-", mode, "-comparison.csv")))
  source <- attr(comp, "source_comparison")
  if (!is.null(source)) {
    mq_save(source, file.path(mq_output, paste0("bridge-", mode, "-source-evidence.rds")))
    mq_write_csv(as.data.frame(source), file.path(mq_output,
      paste0("bridge-", mode, "-source-evidence.csv")))
    reps <- attr(source, "bridge_list")
    if (is.list(reps)) for (model in names(reps)) for (r in seq_along(reps[[model]]))
      mq_save(reps[[model]][[r]], file.path(mq_output,
        paste0("bridge-", mode, "-", model, "-rep", r, ".rds")))
  }
  artifact <- list(settings = settings,
    status = if (!evidence_complete) "incomplete_evidence" else
      if (provisional) "complete_provisional_sampler_failure" else "complete",
    universe_status = universe_status,
    error = "", warnings = captured$warnings, failures = failure_rows,
    comparison = comp, completed = mq_now())
  mq_save(artifact, file.path(mq_output, paste0("bridge-", mode, "-result.rds")))
  invisible(artifact)
}

mq_validate_merger <- function(groups = mq_groups) {
  groups <- as.character(groups)
  if (!identical(groups, mq_groups))
    mq_stop("Fit product groups do not match the frozen PNB data order.")
  if (!mq_acquirer %in% groups || !mq_target %in% groups)
    mq_stop("Acquirer/target validation failed against frozen groups.")
  if (!is.null(mq_bundled_groups) &&
      (!mq_acquirer %in% mq_bundled_groups || !mq_target %in% mq_bundled_groups))
    mq_stop("Acquirer/target validation failed against bundled PNB groups.")
  owner_post <- groups
  owner_post[owner_post == mq_target] <- mq_acquirer
  antitrust::counterfactual(ownership = owner_post)
}
mq_sim_arrays <- function(sim, case, anchor, draw_plan) {
  if (length(sim@price_change) != 1L || is.null(names(sim@price_change)) ||
      anyNA(names(sim@price_change)) || !nzchar(names(sim@price_change)[[1L]]))
    mq_stop("PNB simulation must contain exactly one identified market-period unit.")
  unit <- names(sim@price_change)[[1L]]
  pc <- sim@price_change[[unit]]
  labels <- sim@product_labels[[unit]]
  if (!is.matrix(pc) || !is.numeric(pc) || nrow(pc) < 1L ||
      ncol(pc) < 1L || nrow(pc) != sim@n_draws)
    mq_stop("Simulation price-change matrix violates its draw/product contract.")
  draw_indices <- sim@metadata$draw_indices
  if (!is.numeric(draw_indices) || length(draw_indices) != nrow(pc) ||
      any(!is.finite(draw_indices)) || any(draw_indices < 1) ||
      any(draw_indices != as.integer(draw_indices)) || anyDuplicated(draw_indices) ||
      !identical(as.integer(draw_indices), as.integer(draw_plan$draw_indices)))
    mq_stop("Simulation posterior indices do not match the shared draw plan.")
  draw_indices <- as.integer(draw_indices)
  if (!is.character(labels) || length(labels) != ncol(pc) || anyNA(labels) ||
      any(!nzchar(labels)) || anyDuplicated(labels) || !identical(labels, mq_groups))
    mq_stop("Simulation product labels do not match the frozen ordered products.")
  failures <- sim@failures[[unit]]
  if (!is.character(failures) || length(failures) != nrow(pc) ||
      any(!is.na(failures) & !nzchar(failures)))
    mq_stop("Simulation failure vector is malformed; no replacement was made.")
  if (any(!is.finite(pc[is.na(failures), , drop = FALSE])))
    mq_stop("A simulation row marked successful contains nonfinite price changes.")
  k <- rep(NA_integer_, nrow(pc)); cache <- sim@metadata$k_cache
  ordered_model <- case %in% c("BF", "BS", "CF", "CS")
  if (ordered_model && is.null(cache))
    mq_stop("Ordered-k simulation is missing its realized role cache: ", case)
  if (!is.null(cache) && (!is.list(cache) || is.null(names(cache)) ||
      anyNA(names(cache)) || anyDuplicated(names(cache))))
    mq_stop("Simulation role cache must be a posterior-draw-keyed list.")
  if (is.list(cache)) for (i in seq_along(draw_indices)) {
    one <- cache[[as.character(draw_indices[[i]])]]
    if (is.null(one) || is.null(names(one)) || !unit %in% names(one))
      mq_stop("Role cache is missing posterior draw/unit: ", draw_indices[[i]], "/", unit)
    realized <- one[[unit]]
    if (length(realized) != 1L || !is.numeric(realized) ||
        !is.finite(realized) || realized < 0 || realized != as.integer(realized))
      mq_stop("Malformed realized role state for posterior draw/unit: ", draw_indices[[i]], "/", unit)
    if (ordered_model && !realized %in% 1:5)
      mq_stop("Realized role state is outside the frozen PNB support: ", case, "/", realized)
    k[[i]] <- as.integer(realized)
  }
  failure_matrix <- matrix(NA_character_, nrow(pc), ncol(pc),
                           dimnames = dimnames(pc))
  failure_matrix[!is.na(failures), ] <- failures[!is.na(failures)]
  long <- expand.grid(draw_position = seq_len(nrow(pc)),
                      product_position = seq_len(ncol(pc)))
  long$case <- case; long$anchor <- anchor
  long$draw_index <- draw_indices[long$draw_position]
  long$product_index <- long$product_position
  long$product_label <- labels[long$product_position]
  long$k <- k[long$draw_position]
  long$price_change <- as.numeric(pc[cbind(long$draw_position, long$product_position)])
  long$failure <- as.character(failure_matrix[cbind(long$draw_position,
                                                    long$product_position)])
  valid <- is.finite(long$price_change) & is.na(long$failure)
  product_summary <- do.call(rbind, lapply(seq_len(ncol(pc)), function(j) {
    z <- pc[, j]; ok <- is.finite(z) & is.na(failures)
    data.frame(case = case, anchor = anchor, product_index = j,
      product_label = labels[[j]], n_requested = nrow(pc), n_valid = sum(ok),
      n_failed = sum(!ok), mean = if (any(ok)) mean(z[ok]) else NA_real_,
      sd = if (sum(ok) > 1L) sd(z[ok]) else NA_real_,
      q025 = if (any(ok)) unname(stats::quantile(z[ok], .025)) else NA_real_,
      q500 = if (any(ok)) unname(stats::quantile(z[ok], .5)) else NA_real_,
      q975 = if (any(ok)) unname(stats::quantile(z[ok], .975)) else NA_real_,
      stringsAsFactors = FALSE)
  }))
  flattened_summary <- data.frame(case = case, anchor = anchor,
    n_draws_requested = nrow(pc), n_draws_failed = sum(!vapply(seq_len(nrow(pc)),
      function(i) all(valid[long$draw_position == i]), logical(1))),
    n_valid_product_cells = sum(valid), mean = if (any(valid)) mean(long$price_change[valid]) else NA_real_,
    sd = if (sum(valid) > 1L) sd(long$price_change[valid]) else NA_real_,
    q025 = if (any(valid)) unname(stats::quantile(long$price_change[valid], .025)) else NA_real_,
    q500 = if (any(valid)) unname(stats::quantile(long$price_change[valid], .5)) else NA_real_,
    q975 = if (any(valid)) unname(stats::quantile(long$price_change[valid], .975)) else NA_real_,
    conditioning = "finite product cells only; failed draw/product cells retained in long table",
    stringsAsFactors = FALSE)
  list(case = case, anchor = anchor, unit = unit, draw_indices = draw_indices,
       product_labels = labels, price_change = pc, failures = failures, k = k,
       long = long, product_summary = product_summary,
       flattened_summary = flattened_summary,
       draw_plan_metadata = draw_plan$metadata)
}
mq_policy <- function(case) {
  if (!case %in% mq_all_cases) mq_stop("--policy case is not in the frozen 8-model plan: ", case)
  z <- mq_load_case(case)
  policy_dir <- file.path(mq_output, "policy"); dir.create(policy_dir, recursive = TRUE, showWarnings = FALSE)
  if (is.null(z$fit) || !inherits(z$fit, "BayesFit")) {
    out <- list(case = case, status = "failed", error = "serialized BayesFit missing",
                warnings = character(), completed = mq_now())
    mq_save(out, file.path(policy_dir, paste0(case, "-policy-status.rds")))
    return(invisible(out))
  }
  captured_plan <- mq_capture(posterior_policy_draws(z$fit, n_draws = 200L, seed = 1960L))
  if (nzchar(captured_plan$error)) {
    out <- list(case = case, status = "failed", error = captured_plan$error,
                warnings = captured_plan$warnings, completed = mq_now())
    mq_save(out, file.path(policy_dir, paste0(case, "-policy-status.rds")))
    return(invisible(out))
  }
  draw_plan <- captured_plan$value
  mq_save(draw_plan, file.path(policy_dir, paste0(case, "-draw-plan.rds")))
  fit_groups <- as.character(z$fit@data$group)
  if (anyNA(fit_groups) || anyDuplicated(fit_groups))
    mq_stop("Serialized fit has malformed product groups.")
  scenario <- mq_validate_merger(fit_groups)
  labels <- fit_groups
  inside <- as.numeric(attr(mq_data, "insideSize"))
  if (length(inside) != 1L || !is.finite(inside) || inside <= 0) inside <- NULL
  run_anchor <- function(anchor) {
    captured <- mq_capture(simulate(z$fit, scenario, n_draws = 200L,
      draw_plan = draw_plan, role_mapping = "mixture", demand_anchoring = anchor,
      mc.cores = 1L, insideSize = inside, labels = labels))
    if (nzchar(captured$error)) {
      failed <- list(case = case, anchor = anchor, status = "failed",
        error = captured$error, warnings = captured$warnings,
        simulation = NULL, completed = mq_now())
      mq_save(failed, file.path(policy_dir, paste0(case, "-", anchor, "-status.rds")))
      return(failed)
    }
    sim <- captured$value
    captured_arrays <- mq_capture(mq_sim_arrays(sim, case, anchor, draw_plan))
    if (nzchar(captured_arrays$error)) {
      failed <- list(case = case, anchor = anchor, status = "failed_postprocessing",
        error = captured_arrays$error,
        warnings = unique(c(captured$warnings, captured_arrays$warnings)),
        completed = mq_now())
      mq_save(sim, file.path(policy_dir, paste0(case, "-", anchor, "-invalid-sim.rds")))
      mq_save(failed, file.path(policy_dir, paste0(case, "-", anchor, "-status.rds")))
      return(failed)
    }
    arrays <- captured_arrays$value
    mq_save(sim, file.path(policy_dir, paste0(case, "-", anchor, "-sim.rds")))
    mq_save(arrays, file.path(policy_dir, paste0(case, "-", anchor, "-draws.rds")))
    mq_write_csv(arrays$long, file.path(policy_dir, paste0(case, "-", anchor, "-draws.csv")))
    mq_write_csv(arrays$product_summary, file.path(policy_dir,
      paste0(case, "-", anchor, "-product-summary.csv")))
    mq_write_csv(arrays$flattened_summary, file.path(policy_dir,
      paste0(case, "-", anchor, "-legacy-flattened-summary.csv")))
    valid <- is.finite(arrays$long$price_change) & is.na(arrays$long$failure)
    mq_write_csv(arrays$long[valid, c("case", "anchor", "draw_index", "product_index",
      "product_label", "k", "price_change")], file.path(policy_dir,
      paste0(case, "-", anchor, "-legacy-flattened.csv")))
    complete <- list(status = "complete", error = "", warnings = captured$warnings,
      simulation = sim, arrays = arrays, completed = mq_now())
    mq_save(complete[c("status", "error", "warnings", "completed")],
      file.path(policy_dir, paste0(case, "-", anchor, "-status.rds")))
    complete
  }
  observed <- run_anchor("observed")
  latent <- if (identical(case, "A")) {
    marked <- list(case = case, anchor = "latent_quality", status = "unsupported",
      reason = "Auction latent-quality anchoring is explicitly unsupported by the public API.",
      draw_plan_shared_with_observed = TRUE, completed = mq_now())
    mq_save(marked, file.path(policy_dir, paste0(case, "-latent_quality-unsupported.rds")))
    marked
  } else run_anchor("latent_quality")
  policy_complete <- identical(observed$status, "complete") &&
    (identical(latent$status, "complete") ||
      (identical(case, "A") && identical(latent$status, "unsupported")))
  status <- list(case = case, status = if (policy_complete) "complete" else "failed",
    observed = observed[c("status", "error", "warnings")],
    latent = latent[c("status", "error", "warnings", "reason")],
    draw_plan_metadata = draw_plan$metadata,
    same_draw_plan = TRUE, acquirer = mq_acquirer, target = mq_target,
    completed = mq_now())
  mq_save(status, file.path(policy_dir, paste0(case, "-policy-status.rds")))
  invisible(status)
}

mq_k_decomposition <- function(model, arrays) {
  x <- arrays$long; x <- x[is.finite(x$price_change) & is.na(x$failure) & is.finite(x$k), , drop = FALSE]
  if (!nrow(x)) return(data.frame(model = model, status = "no_realized_k_cells", stringsAsFactors = FALSE))
  products <- sort(unique(x$product_index)); out <- list()
  for (j in products) {
    y <- x[x$product_index == j, , drop = FALSE]
    byk <- split(y$price_change, y$k)
    kk <- as.integer(names(byk)); n <- vapply(byk, length, integer(1)); pk <- n / sum(n)
    means <- vapply(byk, mean, numeric(1)); vars <- vapply(byk, function(v) mean((v - mean(v))^2), numeric(1))
    mu <- sum(pk * means); within <- sum(pk * vars); between <- sum(pk * (means - mu)^2)
    for (i in seq_along(kk)) out[[length(out) + 1L]] <- data.frame(
      model = model, product_index = j, product_label = unique(y$product_label[y$product_index == j]),
      k = kk[[i]], n_valid = n[[i]], k_probability_conditioned = pk[[i]],
      mean = means[[i]], population_variance = vars[[i]],
      within_k_variance = within, between_k_mean_variance = between,
      total_variance = within + between, decomposition_residual = 0,
      conditioning = "valid finite price-change cells by realized posterior k",
      stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}
mq_bma <- function() {
  bma_dir <- file.path(mq_output, "bma"); dir.create(bma_dir, recursive = TRUE, showWarnings = FALSE)
  bridge_path <- file.path(mq_output, "bridge-rich-result.rds")
  if (!file.exists(bridge_path)) {
    out <- list(status = "incomplete_evidence", error = "Run explicit --bridge=rich first.", completed = mq_now())
    mq_save(out, file.path(bma_dir, "bma-status.rds")); return(invisible(out))
  }
  bridge <- readRDS(bridge_path); comp <- bridge$comparison
  tab <- if (!is.null(comp)) as.data.frame(comp) else NULL
  if (identical(as.character(bridge$status)[1L], "complete_provisional_sampler_failure")) {
    out <- list(status = "withheld_unusable_sampler",
      error = "Complete provisional evidence is retained; failed or unknown sampler diagnostics prevent a validated BMA posterior.",
      bridge_status = bridge$status, completed = mq_now())
    mq_save(out, file.path(bma_dir, "bma-status.rds")); return(invisible(out))
  }
  if (!identical(as.character(bridge$status)[1L], "complete") ||
      !mq_complete_evidence(tab, mq_all_cases, attr(comp, "status", exact = TRUE))) {
    out <- list(status = "incomplete_evidence", error = "Eight-model bridge evidence is missing, provisional, or incomplete.",
                bridge_status = bridge$status, completed = mq_now())
    mq_save(out, file.path(bma_dir, "bma-status.rds")); return(invisible(out))
  }
  weights <- tab$posterior_probability[match(mq_all_cases, tab$hypothesis)]
  if (any(!is.finite(weights)) || abs(sum(weights) - 1) > 1e-8)
    mq_stop("BMA posterior weights are not a complete probability vector.")
  arrays <- lapply(mq_all_cases, function(case) {
    p <- file.path(mq_output, "policy", paste0(case, "-observed-draws.rds"))
    if (!file.exists(p)) return(NULL)
    readRDS(p)
  }); names(arrays) <- mq_all_cases
  if (any(vapply(arrays, is.null, logical(1)))) {
    out <- list(status = "incomplete_draw_pool", error = "All eight observed-anchor policy draw pools are required; no subset was used.",
                missing = mq_all_cases[vapply(arrays, is.null, logical(1))], completed = mq_now())
    mq_save(out, file.path(bma_dir, "bma-status.rds")); return(invisible(out))
  }
  validation <- mq_capture({
    reference_indices <- NULL
    for (case in mq_all_cases) {
      indices <- mq_validate_policy_pool(arrays[[case]], case, reference_indices)
      if (is.null(reference_indices)) reference_indices <- indices
    }
    TRUE
  })
  if (nzchar(validation$error)) {
    out <- list(status = "incomplete_draw_pool", error = validation$error,
      warnings = validation$warnings, no_hidden_deletion = TRUE, completed = mq_now())
    mq_save(out, file.path(bma_dir, "bma-status.rds")); return(invisible(out))
  }
  pool <- do.call(rbind, lapply(mq_all_cases, function(case) {
    a <- arrays[[case]]; x <- a$long; n_draw <- length(unique(x$draw_index))
    x$model <- case; x$model_weight <- weights[match(case, mq_all_cases)]
    x$within_model_draw_weight <- 1 / n_draw
    x$cell_weight <- x$model_weight / n_draw
    x$included_in_complete_moments <- is.finite(x$price_change) & is.na(x$failure)
    x
  }))
  mq_save(list(status = "complete_draw_pools", pool = pool, weights = weights,
               conditioning = "common observed merger anchor; failed cells retained and excluded from finite moments",
               completed = mq_now()), file.path(bma_dir, "bma-8-observed-draw-pool.rds"))
  mq_write_csv(pool, file.path(bma_dir, "bma-8-observed-draw-pool.csv"))
  failure_counts <- do.call(rbind, lapply(mq_all_cases, function(case) {
    x <- pool[pool$model == case, , drop = FALSE]
    data.frame(model = case, posterior_weight = weights[match(case, mq_all_cases)],
      n_draws = length(unique(x$draw_index)), n_failed_cells = sum(!x$included_in_complete_moments),
      n_failed_draws = length(unique(x$draw_index[!x$included_in_complete_moments])),
      stringsAsFactors = FALSE)
  }))
  mq_write_csv(failure_counts, file.path(bma_dir, "bma-failure-counts.csv"))
  complete <- all(pool$included_in_complete_moments)
  if (!complete) {
    out <- list(status = "incomplete_draw_pool", error = "Failed draw/product cells prevent an exact complete-pool moment decomposition.",
                failure_counts = failure_counts, no_hidden_deletion = TRUE, completed = mq_now())
    mq_save(out, file.path(bma_dir, "bma-status.rds")); return(invisible(out))
  }
  products <- sort(unique(pool$product_index)); rows <- list()
  for (j in products) {
    per <- lapply(mq_all_cases, function(case) {
      v <- pool$price_change[pool$model == case & pool$product_index == j]
      c(mean = mean(v), variance = mean((v - mean(v))^2), n = length(v))
    }); names(per) <- mq_all_cases
    means <- vapply(per, `[[`, numeric(1), "mean"); vars <- vapply(per, `[[`, numeric(1), "variance")
    mu <- sum(weights * means); within <- sum(weights * vars); between <- sum(weights * (means - mu)^2)
    empirical <- sum(pool$cell_weight[pool$product_index == j] *
      (pool$price_change[pool$product_index == j] - mu)^2)
    rows[[length(rows) + 1L]] <- data.frame(product_index = j,
      product_label = unique(pool$product_label[pool$product_index == j]),
      bma_mean = mu, population_within_family_variance = within,
      between_family_mean_variance = between, total_variance = empirical,
      decomposition_residual = empirical - (within + between),
      conditioning = "complete observed-anchor draw pools; model weights unchanged",
      stringsAsFactors = FALSE)
  }
  mq_write_csv(do.call(rbind, rows), file.path(bma_dir, "bma-8-observed-variance-decomposition.csv"))
  kd <- mq_bind_rows(lapply(mq_all_cases, function(case) mq_k_decomposition(case, arrays[[case]])))
  mq_write_csv(kd, file.path(bma_dir, "bma-within-family-k-decomposition.csv"))
  out <- list(status = "complete", weights = weights, failure_counts = failure_counts,
              variance_file = "bma-8-observed-variance-decomposition.csv",
              k_file = "bma-within-family-k-decomposition.csv", latent = "separate policy outputs only",
              completed = mq_now())
  mq_save(out, file.path(bma_dir, "bma-status.rds"))
  invisible(out)
}

if (identical(mq_action, "--summary")) mq_summary()
if (startsWith(mq_action, "--bridge=")) {
  mode <- sub("^--bridge=", "", mq_action)
  if (!mode %in% c("canonical", "rich")) mq_stop("--bridge must be canonical or rich.")
  mq_bridge(mode)
}
if (startsWith(mq_action, "--policy=")) mq_policy(sub("^--policy=", "", mq_action))
if (identical(mq_action, "--bma")) mq_bma()
