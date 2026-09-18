#!/usr/bin/env Rscript
# Summarize every frozen case; never omit failed/incomplete attempts.
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(antitrustBayes))
suppressPackageStartupMessages(library(posterior))
repo <- normalizePath(Sys.getenv("MQ_REPO", "/mnt/d/projects/antitrustbayes"))
root <- Sys.getenv("MQ_BENCH_ROOT", file.path(repo,
  "doc/statistical_toolkit_runs/20260916/performance-audit/recovery/market-quality/benchmark-v2"))
plan <- readRDS(file.path(root, "frozen-plan.rds"))
rows <- list()
for (j in seq_len(nrow(plan))) {
  p <- plan[j, ]; case <- p$case
  path <- file.path(root, paste0(case, "-result.rds"))
  row <- data.frame(case = case, fixture = p$fixture, mode = p$mode, seed = p$seed,
    status = "not run", wall_seconds = NA_real_, native_seconds = NA_real_,
    peak_RSS_kb = NA_real_, divergences = NA_integer_, treedepth_hits = NA_integer_,
    minimum_BFMI = NA_real_, maximum_Rhat = NA_real_, minimum_bulk_ESS = NA_real_,
    minimum_tail_ESS = NA_real_, minimum_bulk_ESS_second = NA_real_,
    minimum_tail_ESS_second = NA_real_, mean_step_size = NA_real_,
    quality_share_scale_correlation = NA_real_, estimated_coordinate_failures = NA_integer_,
    usable_gate = FALSE, error = "", warnings = "")
  if (!file.exists(path)) {
    if (file.exists(file.path(root, paste0(case, "-started.rds"))))
      row$status <- "unfinished (running or interrupted)"
    rows[[j]] <- row; next
  }
  result <- readRDS(path)
  row$wall_seconds <- result$elapsed; row$native_seconds <- result$native_seconds
  if (!is.null(result$peak_RSS_kb)) row$peak_RSS_kb <- result$peak_RSS_kb
  row$error <- result$error; row$warnings <- paste(unique(result$warnings), collapse = " | ")
  fit_path <- file.path(root, paste0(case, "-fit.rds"))
  if (nzchar(result$error) || !file.exists(fit_path)) {
    row$status <- "failed"; rows[[j]] <- row; next
  }
  fit <- readRDS(fit_path); stanfit <- fit@stanfit
  sampler <- rstan::get_sampler_params(stanfit, inc_warmup = FALSE)
  row$divergences <- sum(vapply(sampler, function(x) sum(x[, "divergent__"]), numeric(1)))
  row$treedepth_hits <- sum(vapply(sampler, function(x)
    sum(x[, "treedepth__"] >= p$max_treedepth), numeric(1)))
  bfmi <- vapply(sampler, function(x) mean(diff(x[, "energy__"])^2)/var(x[, "energy__"]), numeric(1))
  row$minimum_BFMI <- min(bfmi)
  row$mean_step_size <- mean(vapply(sampler, function(x) mean(x[, "stepsize__"]), numeric(1)))
  raw <- antitrustBayes:::.bayes_bridge_raw_parameter_names(stanfit)
  positive_width <- function(n) {
    dim <- stanfit@par_dims[[n]]
    if (is.null(dim)) return(FALSE)
    !length(dim) || prod(dim) > 0
  }
  raw <- raw[vapply(raw, positive_width, logical(1))]
  economic <- c("alpha_structural_market", "gamma_structural_market", "s0_prob",
    "b_market", "sigma_b_strat", "sigma_share", "sigma_margin", "rho",
    "sigma_log_a", "sigma_s0_pool", "mu_s0_pool", "b_quality_all")
  economic <- intersect(economic, stanfit@model_pars)
  economic <- economic[vapply(economic, positive_width, logical(1))]
  array <- rstan::extract(stanfit, pars = unique(c(raw, economic)), permuted = FALSE,
    inc_warmup = FALSE)
  summaries <- as.data.frame(posterior::summarise_draws(posterior::as_draws_array(array)))
  base <- sub("\\[.*$", "", summaries$variable)
  constant <- vapply(seq_len(dim(array)[[3L]]), function(k)
    length(unique(as.vector(array[, , k]))) == 1L, logical(1))
  # An estimated coordinate constant across all draws is a sampler failure.
  # Only truly deterministic transformed outputs may be excluded.
  deterministic <- constant & !(base %in% raw)
  summaries$estimated_coordinate <- base %in% raw
  summaries$deterministic_output <- deterministic
  eligible <- !deterministic
  finite <- is.finite(summaries$rhat) & is.finite(summaries$ess_bulk) & is.finite(summaries$ess_tail)
  row$estimated_coordinate_failures <- sum(eligible & !finite)
  row$maximum_Rhat <- if (all(finite[eligible])) max(summaries$rhat[eligible]) else Inf
  row$minimum_bulk_ESS <- if (all(finite[eligible])) min(summaries$ess_bulk[eligible]) else 0
  row$minimum_tail_ESS <- if (all(finite[eligible])) min(summaries$ess_tail[eligible]) else 0
  summaries$bulk_ESS_wall_second <- summaries$ess_bulk/result$elapsed
  summaries$tail_ESS_wall_second <- summaries$ess_tail/result$elapsed
  row$minimum_bulk_ESS_second <- row$minimum_bulk_ESS/result$elapsed
  row$minimum_tail_ESS_second <- row$minimum_tail_ESS/result$elapsed
  components <- rstan::extract(stanfit, pars = c("sigma_b_strat", "sigma_share"))
  row$quality_share_scale_correlation <- cor(components$sigma_b_strat, components$sigma_share)
  row$usable_gate <- row$divergences == 0 && row$treedepth_hits == 0 &&
    is.finite(row$minimum_BFMI) && row$minimum_BFMI > .3 && row$maximum_Rhat <= 1.01 &&
    row$minimum_bulk_ESS >= 100 && row$minimum_tail_ESS >= 100
  row$status <- if (row$usable_gate) "usable" else "diagnostic gate failed"
  write.csv(summaries, file.path(root, paste0(case, "-parameter-diagnostics.csv")), row.names = FALSE)
  write.csv(data.frame(chain = seq_along(bfmi), BFMI = bfmi),
    file.path(root, paste0(case, "-chain-diagnostics.csv")), row.names = FALSE)
  rows[[j]] <- row
}
write.csv(do.call(rbind, rows), file.path(root, "benchmark-summary.csv"), row.names = FALSE)
print(do.call(rbind, rows)[c("case", "status", "divergences", "minimum_BFMI", "maximum_Rhat",
  "minimum_bulk_ESS_second", "usable_gate")])
