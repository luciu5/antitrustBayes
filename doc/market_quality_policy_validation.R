#!/usr/bin/env Rscript
# Fixed structural states, not posterior recovery or a policy benchmark.
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(antitrustBayes))
suppressPackageStartupMessages(library(rstan))
repo <- normalizePath(Sys.getenv("MQ_REPO", "/mnt/d/projects/antitrustbayes"))
out <- Sys.getenv("MQ_POLICY_OUT", file.path(repo,
  "doc/statistical_toolkit_runs/20260916/performance-audit/recovery/market-quality/policy-validation-01"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
stopifnot(antitrustBayes:::.antitrust_bayes_model_spec_version == 10L)
env <- new.env(parent = globalenv())
env$baseline_mode <- FALSE
env$repo <- repo
env$utility <- new.env(parent = globalenv())
utilities <- c("cga_as_draw", "cga_init_one", "cga_set_all", "cga_common_init",
  "cga_model_for", "cga_fixed")
for (expr in parse(file.path(repo, "doc/collapsed_gaussian_production_validation.R"))) {
  if (is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
      is.symbol(expr[[2L]]) && as.character(expr[[2L]]) %in% utilities)
    eval(expr, env$utility)
}
for (expr in parse(file.path(repo, "doc/market_quality_native_validation.R"))) {
  if (is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
      is.symbol(expr[[2L]]) && as.character(expr[[2L]]) %in% c("fixture", "prepare", "bundle"))
    eval(expr, env)
}
rows <- list()
record <- function(case, check, error, tolerance = 2e-8, note = "") {
  rows[[length(rows)+1L]] <<- data.frame(case, check, error, tolerance,
    passed = is.finite(error) && error <= tolerance, note)
  write.csv(do.call(rbind, rows), file.path(out, "validation.csv"), row.names = FALSE)
}
shell <- function(fx, demand, conduct, structure, scope) {
  d <- env$prepare(fx, demand, conduct, structure, scope, "collapsed")
  b <- env$bundle(d, demand, conduct, structure)
  # Repeated fixed-parameter states expose public draw-selection paths without HMC.
  fit <- rstan::sampling(env$utility$cga_model_for(demand, structure), data = d,
    init = list(b$init), algorithm = "Fixed_param", chains = 1L,
    iter = 4L, warmup = 0L, seed = 41629L, refresh = 0L)
  data <- as.data.frame(fx)
  data$outside_price <- 1.3
  data$loo_cluster <- data$market
  methods::new("BayesLogit", stanfit = fit, stan_data = d, data = data,
    price_sd = if (demand == "logit") d$price_sd else d$log_price_sd,
    supply_model = as.integer(d$supply_model), demand_model = demand,
    market_side = "output", call = quote(fixed_policy_validation()),
    metadata = list(model_spec_version = 10L, quality_scope = scope,
      quality_parameterization = "collapsed", share_likelihood = "gaussian",
      margin_likelihood = "gaussian", posterior_draws = 4L,
      firm_market_ranking = attr(d, "firm_market_ranking", exact = TRUE),
      structure = structure, partition_prior = if (structure == "simultaneous")
        "cutoff" else "ordered_k", k_scope = "market"))
}
run <- function(fit, owner, plan, policy) {
  merged <- owner; merged[merged == "2"] <- "1"
  labels <- paste0("row", seq_len(nrow(fit@data)))
  units <- antitrustBayes:::.estimated_simulation_units(fit, market = NULL,
    year = c("1", "2"))
  cost_scenarios <- lapply(units$info, function(unit) {
    antitrust::counterfactual(
      ownership = if (policy == "merger_cost10") merged[unit$idx] else owner[unit$idx],
      costs = stats::setNames(rep(.1, length(unit$idx)), labels[unit$idx]))
  })
  names(cost_scenarios) <- vapply(units$info, function(unit) unit$unit_id, character(1))
  scenario <- switch(policy,
    unchanged = owner,
    merger = merged,
    cost10 = cost_scenarios,
    merger_cost10 = cost_scenarios)
  simulate(fit, ownerPost = scenario, draw_plan = plan, year = c("1", "2"),
    labels = labels, demand_anchoring = "latent_quality", insideSize = 1, mc.cores = 1L)
}
failures <- list()
replay <- commandArgs(trailingOnly = TRUE)
replay <- replay[startsWith(replay, "--replay=")]
if (length(replay)) {
  source <- sub("^--replay=", "", replay[[1L]])
  for (path in list.files(source, pattern = "_plan\\.rds$", full.names = TRUE)) {
    case <- sub("_plan\\.rds$", "", basename(path))
    plan <- readRDS(path)$plan
    for (policy in c("unchanged", "merger", "cost10", "merger_cost10")) {
      value <- readRDS(file.path(source, paste0(case, "_", policy, "_simulation.rds")))
      record(case, paste0(policy, "_solver_failures"), sum(value@n_draws-value@n_valid), 0)
      for (unit in value@market_ids) for (z in value@metadata$latent_diagnostics[[unit]]) {
        draw <- match(z$posterior_draw_index, plan$draw_indices)
        record(case, paste0(policy, "_", unit, "_quality_state_identity"),
          max(abs(z$quality-plan$quality_draws[draw, ])), 0,
          note = "replay of retained simulations; complete diagnostic quality vector")
      }
    }
  }
  # Independently recover the familiar welfare metrics from retained native models.
  for (path in list.files(source,
      pattern = "^(logit|ces)_(bertrand|cournot)_(unchanged|merger|cost10|merger_cost10)\\.rds$",
      full.names = TRUE)) {
    case <- sub("\\.rds$", "", basename(path))
    sims <- readRDS(path)$simulations
    for (unit in names(sims[[1L]]@sims)) {
      values <- lapply(sims, function(z) vapply(z@sims[[unit]], function(obj) {
        native <- if (methods::is(obj, "StructuralFit")) obj@model else obj
        consumer <- -as.numeric(antitrust::CV(native))
        producer <- sum(antitrust::calcProducerSurplus(native, preMerger = FALSE))-
          sum(antitrust::calcProducerSurplus(native, preMerger = TRUE))
        c(consumer = consumer, producer = producer, total = consumer+producer)
      }, numeric(3L)))
      for (metric in rownames(values[[1L]]))
        record(case, paste(unit, metric, "scope_identity", sep = "_"),
          max(abs(values[[1L]][metric, ]-values[[2L]][metric, ])))
    }
  }
  cat("Retained policy replay checks:", length(rows), "failed:",
    sum(!vapply(rows, function(x) x$passed, logical(1))), "\n")
  quit(status = as.integer(any(!vapply(rows, function(x) x$passed, logical(1)))))
}
if (!"--multi-only" %in% commandArgs(trailingOnly = TRUE)) {
for (demand in c("logit", "ces")) for (conduct in c("bertrand", "cournot")) {
  case <- paste(demand, conduct, sep = "_")
  fx <- env$fixture(M = 1L, T = 2L)
  fits <- lapply(c("global", "market"), function(scope)
    shell(fx, demand, conduct, "simultaneous", scope))
  plans <- lapply(fits, posterior_policy_draws, n_draws = 4L, seed = 812L)
  record(case, "single_market_quality_draw_identity",
    max(abs(plans[[1L]]$quality_draws-plans[[2L]]$quality_draws)), 0)
  owner <- as.character(fx$group)
  for (policy in c("unchanged", "merger", "cost10", "merger_cost10")) {
    sims <- lapply(seq_along(fits), function(i) tryCatch(
      run(fits[[i]], owner, plans[[i]], policy), error = function(e) e))
    for (i in seq_along(sims)) {
      if (inherits(sims[[i]], "condition")) {
        failures[[length(failures)+1L]] <- list(case = case, scope = c("global", "market")[i],
          policy = policy, error = conditionMessage(sims[[i]]))
        saveRDS(failures, file.path(out, "errors.rds"))
        record(case, paste0(policy, "_scope", i, "_execution"), Inf,
          note = conditionMessage(sims[[i]]))
      } else record(case, paste0(policy, "_scope", i, "_solver_failures"),
        sum(sims[[i]]@n_draws-sims[[i]]@n_valid), 0)
    }
    if (any(vapply(sims, inherits, logical(1), "condition"))) next
    for (field in c("price_pre", "price_post", "shares_pre", "shares_post"))
      record(case, paste(policy, field, "scope_identity", sep = "_"),
        max(abs(unlist(methods::slot(sims[[1L]], field))-
                  unlist(methods::slot(sims[[2L]], field)))))
    # Same realized persistent quality must be reused in every policy and period.
    saveRDS(list(plans = plans, simulations = sims),
      file.path(out, paste0(case, "_", policy, ".rds")))
  }
}
}
for (conduct in c("bertrand", "cournot")) for (structure in c("fringe", "stackelberg")) {
  case <- paste("multi_logit", conduct, structure, sep = "_")
  fx <- env$fixture(M = 2L, T = 2L)
  fit <- shell(fx, "logit", conduct, structure, "market")
  plan <- posterior_policy_draws(fit, n_draws = 4L, seed = 619L)
  record(case, "two_independent_quality_blocks",
    abs(ncol(plan$quality_draws)-sum(fit@stan_data$quality_block_n)), 0)
  for (m in 1:2) {
    d <- fit@stan_data
    jj <- d$quality_block_start[m]+seq_len(d$quality_block_n[m])-1L
    record(case, paste0("quality_center_market", m), max(abs(rowSums(plan$quality_draws[, jj]))), 2e-13)
    for (z in plan$k_cache) {
      states <- as.numeric(z[c(paste0(m, "_1"), paste0(m, "_2"))])
      record(case, paste0("persistent_k_market", m), abs(diff(states)), 0)
    }
  }
  for (policy in c("unchanged", "merger", "cost10", "merger_cost10")) {
    value <- tryCatch(run(fit, as.character(fx$group), plan, policy),
      error = function(e) e)
    if (inherits(value, "condition")) {
      record(case, paste0(policy, "_execution"), Inf, note = conditionMessage(value))
      next
    }
    record(case, paste0(policy, "_solver_failures"),
      sum(value@n_draws-value@n_valid), 0)
    for (unit in value@market_ids) for (z in value@metadata$latent_diagnostics[[unit]]) {
      # Diagnostics intentionally retain the complete scoped state, with the
      # observation map selecting the local block inside demand reconstruction.
      draw <- match(z$posterior_draw_index, plan$draw_indices)
      record(case, paste0(policy, "_", unit, "_quality_state_identity"),
        max(abs(z$quality-plan$quality_draws[draw, ])), 0)
    }
    saveRDS(value, file.path(out, paste0(case, "_", policy, "_simulation.rds")))
  }
  saveRDS(list(plan = plan, fit = fit), file.path(out, paste0(case, "_plan.rds")))
}
cat("Fixed policy checks:", length(rows), "failed:", sum(!vapply(rows, function(x) x$passed, logical(1))), "\n")
if (any(!vapply(rows, function(x) x$passed, logical(1)))) quit(status = 1L)
