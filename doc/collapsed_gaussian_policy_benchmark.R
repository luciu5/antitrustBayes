#!/usr/bin/env Rscript
## Reuse usable frozen-benchmark fits; never refit for a policy change.
cga_policy_benchmark <- function(run_dir, output = "doc/collapsed_gaussian_policy_comparison.csv") {
  e <- new.env()
  sys.source("doc/collapsed_gaussian_benchmark.R", e)
  diagnostics <- read.csv(file.path(run_dir, "benchmark.csv"), stringsAsFactors = FALSE)
  plan <- e$cga_make_plan()
  rows <- list()
  append <- function(case, policy, metric, values = numeric(), elapsed = NA_real_,
                     failure = NA_real_, status = "complete", error = "") {
    good <- values[is.finite(values)]
    if (identical(status, "complete") && length(good) < 100L)
      status <- if (length(good)) "partial_failure" else "metric_unavailable"
    q <- if (length(good)) quantile(good, c(.05, .5, .95), names = FALSE) else rep(NA_real_, 3)
    rows[[length(rows) + 1L]] <<- data.frame(
      case_id = case$case_id, design = case$design, mode = case$mode,
      policy = policy, metric = metric, status = status, error = error,
      requested_draws = 100L, valid_draws = length(good), q05 = q[1],
      median = q[2], q95 = q[3], mean = if (length(good)) mean(good) else NA_real_,
      sd = if (length(good) > 1) sd(good) else NA_real_,
      elapsed_seconds = elapsed, solver_failure_rate = failure,
      policy_ess_per_second = NA_real_,
      policy_ess_note = "not inferred from the prespecified 100-draw random subsample",
      reconstruction = "existing observed-share deterministic reconstruction",
      units = "insideSize=1; consumer change=-CV; total=producer change-CV",
      stringsAsFactors = FALSE)
    write.csv(do.call(rbind, rows), output, row.names = FALSE, na = "")
  }
  for (i in seq_len(nrow(plan))) {
    case <- plan[i, , drop = FALSE]
    diag <- diagnostics[diagnostics$case_id == case$case_id, , drop = FALSE]
    if (nrow(diag) != 1L || !isTRUE(diag$usable[[1]])) {
      append(case, "all", "not_run", status = "unusable_fit",
             error = if (nrow(diag)) diag$status[[1]] else "fit not completed")
      next
    }
    fit <- readRDS(file.path(run_dir, "stanfit", paste0(case$case_id, ".rds")))
    panel <- e$cga_make_dgp(seed = case$seed, demand_model = case$demand_model,
      supply_model = case$supply_model, structure = case$structure,
      share_scenario = case$share_scenario)
    d <- panel$data
    labels <- paste0("p", seq_len(nrow(d)))
    owner <- as.character(d$group)
    size <- tapply(d$quantities * if (case$demand_model == "ces") d$prices else 1,
                   owner, sum)
    largest <- names(size)[order(-size, names(size))][1:2]
    merged <- owner
    merged[merged == largest[2]] <- largest[1]
    costs <- setNames(rep(.10, nrow(d)), labels)
    policies <- list(
      merger = antitrust::counterfactual(ownership = merged),
      cost10 = antitrust::counterfactual(costs = costs),
      merger_cost10 = antitrust::counterfactual(ownership = merged, costs = costs))
    seed <- e$cga_derived_seed(case$seed, case$case_id, "paired_policies")
    for (policy in names(policies)) {
      cat(case$case_id, policy, "\n")
      started <- proc.time()[["elapsed"]]
      value <- tryCatch(e$cga_with_seed(seed, antitrustBayes::bayes_sim(
        fit, ownerPost = policies[[policy]], n_draws = 100L, labels = labels,
        insideSize = 1, mc.cores = 1L, role_mapping = "mixture")), error = function(x) x)
      elapsed <- proc.time()[["elapsed"]] - started
      if (inherits(value, "condition")) {
        append(case, policy, "not_run", elapsed = elapsed, failure = 1,
          status = "simulation_error", error = conditionMessage(value))
        next
      }
      failure <- 1 - sum(value@n_valid) / (value@n_draws * length(value@market_ids))
      unit <- value@market_ids[[1]]
      for (j in seq_along(labels)) {
        append(case, policy, paste0("price_change:", labels[j]),
               value@price_change[[unit]][, j], elapsed, failure)
        append(case, policy, paste0("share_change:", labels[j]),
               value@shares_post[[unit]][, j] - value@shares_pre[[unit]][, j], elapsed, failure)
      }
      weights <- d$quantities / sum(d$quantities)
      append(case, policy, "fixed_weight_price_index_change",
        as.numeric(value@price_post[[unit]] %*% weights) /
          as.numeric(value@price_pre[[unit]] %*% weights) - 1, elapsed, failure)
      append(case, policy, "inside_share_change",
        rowSums(value@shares_post[[unit]] - value@shares_pre[[unit]]), elapsed, failure)
      welfare_errors <- character()
      welfare <- lapply(value@sims[[unit]], function(model) {
        if (is.null(model)) return(rep(NA_real_, 3))
        tryCatch({
          consumer <- -as.numeric(antitrust::CV(model))
          producer <- sum(antitrust::calcProducerSurplus(model, preMerger = FALSE)) -
            sum(antitrust::calcProducerSurplus(model, preMerger = TRUE))
          c(consumer, producer, consumer + producer)
        }, error = function(x) {
          welfare_errors <<- c(welfare_errors, conditionMessage(x))
          rep(NA_real_, 3)
        })
      })
      welfare <- do.call(rbind, welfare)
      for (j in 1:3) append(case, policy,
        c("consumer_welfare_change", "producer_surplus_change", "total_welfare_change")[j],
        welfare[, j], elapsed, failure,
        error = paste(unique(welfare_errors), collapse = " | "))
    }
  }
  invisible(do.call(rbind, rows))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 1L || !startsWith(normalizePath(args[1], mustWork = FALSE), "/tmp/"))
    stop("Usage: Rscript doc/collapsed_gaussian_policy_benchmark.R /tmp/cga-run")
  cga_policy_benchmark(args[1])
}
