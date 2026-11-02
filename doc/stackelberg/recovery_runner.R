#!/usr/bin/env Rscript

## Compact, prespecified recovery runner.
##
## This file is deliberately separate from the package and synthetic oracle.
## The default command prints a primary-run plan; HMC requires an explicit
## --run and the caller must provide an output directory under /tmp. There is
## no smoke sampling mode and no output written to the repository.

or_null <- function(x, y) if (is.null(x)) y else x

recovery_runner_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- a[startsWith(a, "--file=")]
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[[1L]]), mustWork = FALSE)) else
    getOption("stackelberg.recovery.runner.dir", getwd())
})

recovery_usage <- function() {
  cat(paste(
    "Usage:",
    "  Rscript recovery_runner.R --out=/tmp/native10-recovery [--run]",
    "",
    "Selectors accept comma-separated values:",
    "  --demand=logit,ces",
    "  --conduct=bertrand,cournot",
    "  --structure=simultaneous,fringe,stackelberg",
    "  --truth-structure=stackelberg,fringe",
    "  --seed=1103,2207,3301",
    "  --layout=single,multi_leader,multi_follower",
    "",
    "The default seed/layout mapping is 1103/single, 2207/multi_leader,",
    "and 3301/multi_follower. A single case can be screened with, e.g.,",
    "  --single-case=logit/bertrand/stackelberg/1103/single",
    "  --single-case=logit/bertrand/stackelberg/1103/single --truth-structure=fringe",
    "",
    "Primary sampler settings are fixed at two chains, 1000 iterations,",
    "500 warmup iterations, adapt_delta=.95, and max_treedepth=13.",
    "Without --run the script writes no fits or logs and makes no recovery claim.",
    sep = "\n"
  ), "\n")
}

recovery_option <- function(args, name, default = NULL) {
  key <- paste0("--", name, "=")
  hit <- args[startsWith(args, key)]
  if (!length(hit)) return(default)
  sub(key, "", hit[[1L]], fixed = TRUE)
}

recovery_values <- function(value, allowed, name) {
  if (is.null(value) || !nzchar(value)) return(allowed)
  ans <- strsplit(value, ",", fixed = TRUE)[[1L]]
  if (!length(ans) || any(!nzchar(ans)) || any(!ans %in% allowed))
    stop("--", name, " must contain only: ", paste(allowed, collapse = ","))
  unique(ans)
}

recovery_integer_values <- function(value, name) {
  if (is.null(value) || !nzchar(value)) return(NULL)
  ans <- suppressWarnings(as.integer(strsplit(value, ",", fixed = TRUE)[[1L]]))
  if (!length(ans) || anyNA(ans) || any(ans <= 0L))
    stop("--", name, " must contain positive integers")
  unique(ans)
}

recovery_parse_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (any(args %in% c("--help", "-h"))) return(list(help = TRUE))
  single <- recovery_option(args, "single-case")
  if (!is.null(single)) {
    bits <- strsplit(single, "/", fixed = TRUE)[[1L]]
    if (length(bits) != 5L)
      stop("--single-case must be demand/conduct/structure/seed/layout")
    if (any(vapply(bits, function(x) !nzchar(x), logical(1))))
      stop("--single-case has an empty field")
    if (length(recovery_values(bits[[1L]], c("logit", "ces"), "demand")) != 1L ||
        length(recovery_values(bits[[2L]], c("bertrand", "cournot"), "conduct")) != 1L ||
        length(recovery_values(bits[[3L]], c("simultaneous", "fringe", "stackelberg"), "structure")) != 1L)
      stop("--single-case fields are invalid")
    ss <- recovery_integer_values(bits[[4L]], "single-case seed")
    ll <- recovery_values(bits[[5L]], c("single", "multi_leader", "multi_follower"), "layout")
    tt <- recovery_values(recovery_option(args, "truth-structure", "stackelberg"),
                          c("stackelberg", "fringe"), "truth-structure")
    if (length(ss) != 1L || length(ll) != 1L) stop("--single-case seed/layout must be scalar")
    if (length(tt) != 1L) stop("--single-case truth-structure must be scalar")
    return(list(help = FALSE, out = recovery_option(args, "out"), run = any(args == "--run"),
                demand = bits[[1L]], conduct = bits[[2L]], structure = bits[[3L]],
                truth_structure = tt, seed = ss, layout = ll, single_case = TRUE))
  }
  list(
    help = FALSE, out = recovery_option(args, "out"), run = any(args == "--run"),
    demand = recovery_values(recovery_option(args, "demand"), c("logit", "ces"), "demand"),
    conduct = recovery_values(recovery_option(args, "conduct"), c("bertrand", "cournot"), "conduct"),
    structure = recovery_values(recovery_option(args, "structure"),
                                c("simultaneous", "fringe", "stackelberg"), "structure"),
    truth_structure = recovery_values(recovery_option(args, "truth-structure"),
                                      c("stackelberg", "fringe"), "truth-structure"),
    seed = recovery_integer_values(recovery_option(args, "seed"), "seed"),
    layout = {
      value <- recovery_option(args, "layout")
      if (is.null(value)) NULL else recovery_values(
        value, c("single", "multi_leader", "multi_follower"), "layout"
      )
    },
    single_case = FALSE
  )
}

recovery_make_plan <- function(demand = c("logit", "ces"),
                               conduct = c("bertrand", "cournot"),
                               structure = c("simultaneous", "fringe", "stackelberg"),
                               seed = c(1103L, 2207L, 3301L), layout = NULL,
                               truth_structure = c("stackelberg", "fringe"),
                               single_case = FALSE) {
  demand <- recovery_values(paste(demand, collapse = ","), c("logit", "ces"), "demand")
  conduct <- recovery_values(paste(conduct, collapse = ","), c("bertrand", "cournot"), "conduct")
  structure <- recovery_values(paste(structure, collapse = ","),
                               c("simultaneous", "fringe", "stackelberg"), "structure")
  truth_structure <- recovery_values(paste(truth_structure, collapse = ","),
                                     c("stackelberg", "fringe"), "truth_structure")
  seed <- as.integer(seed)
  if (!length(seed) || anyNA(seed) || any(seed <= 0L)) stop("seed must be positive")
  layout_supplied <- !is.null(layout)
  if (!layout_supplied) {
    mapped <- setNames(c("single", "multi_leader", "multi_follower"),
                       c("1103", "2207", "3301"))
    layout <- unname(mapped[as.character(seed)])
    layout[is.na(layout)] <- "single"
  } else {
    layout <- recovery_values(paste(layout, collapse = ","),
                              c("single", "multi_leader", "multi_follower"), "layout")
  }
  if (single_case && (length(demand) != 1L || length(conduct) != 1L ||
                      length(structure) != 1L || length(truth_structure) != 1L ||
                      length(seed) != 1L || length(layout) != 1L))
    stop("single-case selection must have one value for every selector")
  if (length(layout) == length(seed) && !layout_supplied) {
    rows <- expand.grid(demand = demand, conduct = conduct, structure = structure,
                        truth_structure = truth_structure,
                        case_index = seq_along(seed), stringsAsFactors = FALSE)
    rows$seed <- seed[rows$case_index]
    rows$layout <- layout[rows$case_index]
    rows$case_index <- NULL
  } else {
    rows <- expand.grid(demand = demand, conduct = conduct, structure = structure,
                        truth_structure = truth_structure, seed = seed,
                        layout = layout, stringsAsFactors = FALSE)
  }
  rows <- rows[order(match(rows$demand, c("logit", "ces")),
                     match(rows$conduct, c("bertrand", "cournot")),
                     match(rows$structure, c("simultaneous", "fringe", "stackelberg")),
                     match(rows$truth_structure, c("stackelberg", "fringe")),
                     rows$seed, match(rows$layout,
                                      c("single", "multi_leader", "multi_follower"))), , drop = FALSE]
  rownames(rows) <- NULL
  rows$case_id <- paste(rows$truth_structure, rows$demand, rows$conduct,
                        rows$structure, rows$seed, rows$layout, sep = "__")
  rows
}

recovery_primary_settings <- function() {
  list(chains = 2L, iter = 1000L, warmup = 500L, cores = 1L,
       refresh = 0L, control = list(adapt_delta = 0.95, max_treedepth = 13L),
       profile = "primary", smoke = FALSE)
}

recovery_require_tmp <- function(path) {
  if (is.null(path) || !nzchar(path)) stop("--out=/tmp/... is required")
  path <- normalizePath(path, mustWork = FALSE)
  if (!identical(path, "/tmp") && !startsWith(path, "/tmp/"))
    stop("--out must be a caller-supplied directory under /tmp")
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(path)) stop("Unable to create output directory: ", path)
  path
}

recovery_load_dgp <- function(path = file.path(recovery_runner_dir, "synthetic_dgp.R")) {
  if (!file.exists(path)) stop("Synthetic DGP not found: ", path)
  env <- new.env(parent = globalenv())
  sys.source(path, envir = env)
  if (!is.function(env$synthetic_stackelberg_panel))
    stop("Synthetic DGP did not define synthetic_stackelberg_panel()")
  env$synthetic_stackelberg_panel
}

recovery_make_dgp <- function(case, dgp_fun) {
  dgp_fun(demand = case$demand, conduct = case$conduct,
          seed = as.integer(case$seed), layout = case$layout,
          normalize_bayes = TRUE, partition = "ordered_k", k = 2L,
          structure = case$truth_structure)
}

recovery_validate_dgp <- function(panel, truth_structure = panel$structure) {
  if (!is.list(panel) || !is.data.frame(panel$data) || !length(panel$truth))
    stop("Synthetic DGP returned an invalid panel")
  if (!isTRUE(panel$parameters$normalize_bayes))
    stop("Recovery requires normalize_bayes=TRUE")
  if (!identical(panel$parameters$partition, "ordered_k") ||
      !identical(as.integer(panel$parameters$k), 2L))
    stop("Recovery DGP must be ordered-k with true k=2")
  if (!identical(panel$structure, truth_structure))
    stop("Recovery DGP truth structure does not match the selected case")
  if (length(panel$truth) != 8L)
    stop("The finalized recovery DGP must contain eight periods")
  k <- vapply(panel$truth, function(x) as.integer(x$k_realized), integer(1))
  if (!identical(k, rep(2L, 8L))) stop("DGP true k is not 2 in every period")
  fixed_rank <- panel$truth[[1L]]$fixed_rank
  if (!all(vapply(panel$truth, function(x) identical(x$fixed_rank, fixed_rank), logical(1))))
    stop("DGP firm ranking is not fixed across periods")
  if (any(!is.finite(panel$data$prices)) || any(!is.finite(panel$data$quantities)) ||
      any(!is.finite(panel$data$margins))) stop("DGP contains non-finite observations")
  invisible(TRUE)
}

recovery_fit_args <- function(panel, case, settings = recovery_primary_settings()) {
  out <- list(
    prices = panel$data$prices, quantities = panel$data$quantities,
    margins = panel$data$margins, market = panel$data$market,
    year = panel$data$year, group = panel$data$group,
    demand_model = case$demand, market_side = "output", supply_model = case$conduct,
    structure = case$structure, use_rho = FALSE, df = 5,
    use_market_resid_scale = FALSE, use_alpha_market_pool = FALSE,
    use_s0_market_year_pool = FALSE, use_s0_market_pool = FALSE,
    use_hmt = FALSE, prior_manifest = "baseline", likelihood = "joint",
    loo_cluster = panel$data$market, model_prior = 1 / 3,
    outside_price = panel$parameters$outside_price,
    stan_options = settings[c("chains", "iter", "warmup", "cores", "refresh", "control")]
  )
  bytes <- utf8ToInt(paste(case$truth_structure, case$demand, case$conduct,
                           case$structure, case$seed, case$layout, sep = "|"))
  out$stan_options$seed <- as.integer(1L + (as.numeric(case$seed) * 104729 +
    sum(seq_along(bytes) * bytes)) %% 2147483646)
  if (identical(case$structure, "fringe")) {
    out$partition_prior <- "ordered_k"
    out$k_support <- 0:5
    out$k_prior <- rep(1 / 6, 6)
    out$k_scope <- "market"
  } else if (identical(case$structure, "stackelberg")) {
    out$partition_prior <- "ordered_k"
    out$k_support <- 1:5
    out$k_prior <- rep(1 / 5, 5)
    out$k_scope <- "market"
  }
  out
}

recovery_write <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, na = "")
  invisible(path)
}

recovery_summary <- function(x, target) {
  z <- suppressWarnings(as.numeric(x)); z <- z[is.finite(z)]
  if (!length(z)) return(c(mean = NA_real_, median = NA_real_, q025 = NA_real_,
                           q975 = NA_real_, target = target, error = NA_real_, covered = NA_real_))
  q <- stats::quantile(z, c(.025, .5, .975), names = FALSE)
  c(mean = mean(z), median = q[[2L]], q025 = q[[1L]], q975 = q[[3L]],
    target = target, error = q[[2L]] - target,
    covered = as.numeric(target >= q[[1L]] && target <= q[[3L]]))
}

recovery_extract <- function(stanfit, name) {
  if (!requireNamespace("rstan", quietly = TRUE)) return(NULL)
  ans <- tryCatch(rstan::extract(stanfit, pars = name, permuted = TRUE),
                  error = function(e) NULL)
  if (is.null(ans) || !length(ans)) return(NULL)
  ans[[1L]]
}

recovery_matrix <- function(x) {
  if (is.null(x)) return(NULL)
  d <- dim(x)
  if (is.null(d)) return(matrix(as.numeric(x), nrow = length(x), ncol = 1L))
  if (!length(d) || d[[1L]] < 1L) return(NULL)
  matrix(as.numeric(x), nrow = d[[1L]], ncol = prod(d[-1L]))
}

recovery_margin_errors <- function(st, panel) {
  pm <- recovery_extract(st, "pred_margin")
  if (is.null(pm)) {
    inv <- recovery_extract(st, "pred_marginNormInv")
    if (!is.null(inv)) pm <- 1 / inv
  }
  mat <- recovery_matrix(pm)
  n <- nrow(panel$data)
  if (is.null(mat) || ncol(mat) != n) {
    return(c(margin_mae_overall = NA_real_, margin_mae_core = NA_real_,
             margin_mae_noncore = NA_real_))
  }
  med <- apply(mat, 2L, stats::median, na.rm = TRUE)
  err <- abs(med - as.numeric(panel$data$margins))
  core <- panel$data$group %in% as.character(panel$truth[[1L]]$core)
  c(margin_mae_overall = mean(err, na.rm = TRUE),
    margin_mae_core = if (any(core)) mean(err[core], na.rm = TRUE) else NA_real_,
    margin_mae_noncore = if (any(!core)) mean(err[!core], na.rm = TRUE) else NA_real_)
}

recovery_fit_diagnostics <- function(fit, case, panel, elapsed, error = "", warnings = "") {
  base <- data.frame(case_id = case$case_id, truth_structure = case$truth_structure,
                     demand = case$demand, conduct = case$conduct, structure = case$structure,
                     seed = case$seed, layout = case$layout,
                     run_profile = "primary", status = if (nzchar(error)) "fit_error" else "fit_complete",
                     error = error, warnings = warnings, elapsed_seconds = elapsed,
                     stringsAsFactors = FALSE)
  if (!is.null(fit) && !nzchar(error)) {
    st <- fit@stanfit
    demand_name <- if (case$demand == "logit") "alpha_structural_market" else "gamma_structural_market"
    demand_target <- if (case$demand == "logit") panel$parameters$alpha else panel$parameters$gamma
    ds <- recovery_summary(recovery_extract(st, demand_name), demand_target)
    os <- recovery_summary(recovery_extract(st, "s0_prob"), panel$parameters$s0)
    names(ds) <- paste0("demand_", names(ds))
    names(os) <- paste0("outside_", names(os))
    sm <- tryCatch(rstan::summary(st)$summary, error = function(e) NULL)
    monitored <- if (!is.null(sm)) grepl("^(alpha_structural_market|gamma_structural_market|s0_prob|sigma_|log_sigma)", rownames(sm)) else logical()
    if (!length(monitored) || !any(monitored)) monitored <- rep(TRUE, if (is.null(sm)) 0L else nrow(sm))
    rhat <- if (!is.null(sm) && length(monitored)) sm[monitored, "Rhat"] else numeric()
    ess <- if (!is.null(sm) && length(monitored)) sm[monitored, "n_eff"] else numeric()
    sp <- tryCatch(rstan::get_sampler_params(st, inc_warmup = FALSE), error = function(e) list())
    div <- sum(vapply(sp, function(x) if ("divergent__" %in% colnames(x)) sum(x[, "divergent__"] > 0) else 0, numeric(1)))
    depth <- recovery_primary_settings()$control$max_treedepth
    hits <- sum(vapply(sp, function(x) if ("treedepth__" %in% colnames(x)) sum(x[, "treedepth__"] >= depth) else 0, numeric(1)))
    bfmi <- vapply(sp, function(x) {
      if (!"energy__" %in% colnames(x)) return(NA_real_)
      e <- as.numeric(x[, "energy__"])
      if (length(e) < 2L || !is.finite(stats::var(e)) || stats::var(e) <= 0) NA_real_ else
        mean(diff(e)^2) / stats::var(e)
    }, numeric(1))
    gate <- isTRUE(length(rhat) && all(is.finite(rhat)) && max(rhat) <= 1.01 &&
      length(ess) && all(is.finite(ess)) && min(ess) >= 100 &&
      length(bfmi) && all(is.finite(bfmi)) && min(bfmi) >= .30 && div == 0 && hits == 0)
    vals_num <- c(ds, os, recovery_margin_errors(st, panel),
                  max_rhat = if (length(rhat)) max(rhat, na.rm = TRUE) else NA_real_,
                  min_ess = if (length(ess)) min(ess, na.rm = TRUE) else NA_real_,
                  min_bfmi = if (length(bfmi)) min(bfmi, na.rm = TRUE) else NA_real_,
                  divergences = div, treedepth_hits = hits)
    vals_chr <- c(gate_status = if (gate) "pass" else "fail",
                  recovery_claim = "not_claimed_by_runner")
  } else {
    vals_num <- c(setNames(rep(NA_real_, 17L),
                           c("demand_mean", "demand_median", "demand_q025", "demand_q975",
                             "demand_target", "demand_error", "demand_covered",
                             "outside_mean", "outside_median", "outside_q025", "outside_q975",
                             "outside_target", "outside_error", "outside_covered",
                             "margin_mae_overall", "margin_mae_core", "margin_mae_noncore")),
                  max_rhat = NA_real_, min_ess = NA_real_, min_bfmi = NA_real_,
                  divergences = NA_real_, treedepth_hits = NA_real_)
    vals_chr <- c(gate_status = "not_run", recovery_claim = "not_claimed_by_runner")
  }
  ## Keep numeric diagnostics numeric in the compact CSV; gate labels remain
  ## explicit character fields and never imply a recovery claim.
  for (nm in names(vals_num)) base[[nm]] <- unname(vals_num[[nm]])
  for (nm in names(vals_chr)) base[[nm]] <- unname(vals_chr[[nm]])
  base
}

recovery_posterior_k <- function(fit, case, panel) {
  if (identical(case$structure, "simultaneous")) {
    return(data.frame(case_id = case$case_id, truth_structure = case$truth_structure,
                      demand = case$demand, conduct = case$conduct,
                      structure = case$structure, seed = case$seed, layout = case$layout,
                      k = NA_integer_, posterior_probability = NA_real_, posterior_k_mean = NA_real_,
                      posterior_k_median = NA_real_, posterior_k_mode = NA_real_,
                      posterior_core_share = NA_real_, true_k_probability = NA_real_,
                      status = "not_applicable",
                      stringsAsFactors = FALSE))
  }
  support <- if (case$structure == "fringe") 0:5 else 1:5
  if (is.null(fit)) {
    return(data.frame(case_id = case$case_id, truth_structure = case$truth_structure,
                      demand = case$demand, conduct = case$conduct,
                      structure = case$structure, seed = case$seed, layout = case$layout,
                      k = support, posterior_probability = NA_real_, posterior_k_mean = NA_real_,
                      posterior_k_median = NA_real_, posterior_k_mode = NA_real_,
                      posterior_core_share = NA_real_, true_k_probability = NA_real_,
                      status = "fit_unavailable",
                      stringsAsFactors = FALSE))
  }
  st <- fit@stanfit
  p <- recovery_matrix(recovery_extract(st, "posterior_component_prob"))
  if (is.null(p) || ncol(p) != length(support)) {
    return(data.frame(case_id = case$case_id, truth_structure = case$truth_structure,
                      demand = case$demand, conduct = case$conduct,
                      structure = case$structure, seed = case$seed, layout = case$layout,
                      k = support, posterior_probability = NA_real_, posterior_k_mean = NA_real_,
                      posterior_k_median = NA_real_, posterior_k_mode = NA_real_,
                      posterior_core_share = NA_real_, true_k_probability = NA_real_,
                      status = "TODO_missing_generated_quantity",
                      stringsAsFactors = FALSE))
  }
  prob <- colMeans(p, na.rm = TRUE)
  kmean <- sum(support * prob)
  ord <- order(support)
  kmedian <- support[ord][which(cumsum(prob[ord]) >= .5)[1L]]
  kmode <- support[which.max(prob)]
  core <- recovery_extract(st, "posterior_core_share_conditional")
  core_value <- if (!is.null(core) && any(is.finite(core))) mean(as.numeric(core), na.rm = TRUE) else
    NA_real_
  true_prob <- prob[match(2L, support)]
  data.frame(case_id = case$case_id, truth_structure = case$truth_structure,
             demand = case$demand, conduct = case$conduct,
             structure = case$structure, seed = case$seed, layout = case$layout,
             k = support, posterior_probability = prob,
             posterior_k_mean = kmean, posterior_k_median = kmedian,
             posterior_k_mode = kmode, posterior_core_share = core_value,
             true_k_probability = true_prob,
             status = "complete", stringsAsFactors = FALSE)
}

recovery_run_one <- function(panel, case, out_dir, settings) {
  id <- case$case_id
  log_dir <- file.path(out_dir, "logs"); fit_dir <- file.path(out_dir, "stanfit")
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)
  log_path <- file.path(log_dir, paste0(id, ".log"))
  fit_args <- recovery_fit_args(panel, case, settings)
  fit <- NULL; err <- ""; warns <- character(); output <- character()
  started <- proc.time()[["elapsed"]]
  if (!requireNamespace("antitrustBayes", quietly = TRUE)) {
    err <- "antitrustBayes package is unavailable in the active runtime"
  } else {
    output <- capture.output({
      captured <- tryCatch({
        value <- withCallingHandlers(do.call(antitrustBayes::bayes_est, fit_args),
          warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
        list(value = value, error = NULL)
      }, error = function(e) list(value = NULL, error = conditionMessage(e)))
      fit <- captured$value
      if (!is.null(captured$error)) err <- captured$error
    }, type = "output")
  }
  writeLines(c(
    output,
    if (length(warns)) paste("WARNING:", warns) else character(),
    if (nzchar(err)) paste("ERROR:", err) else character()
  ), log_path)
  elapsed <- proc.time()[["elapsed"]] - started
  if (!is.null(fit) && !nzchar(err)) saveRDS(fit, file.path(fit_dir, paste0(id, ".rds")))
  diag <- recovery_fit_diagnostics(fit, case, panel, elapsed, err, paste(warns, collapse = " | "))
  list(fit = fit, diagnostics = diag,
       posterior_k = recovery_posterior_k(fit, case, panel), log_path = log_path)
}

recovery_run <- function(plan, out_dir, dgp_path = file.path(recovery_runner_dir, "synthetic_dgp.R")) {
  out_dir <- recovery_require_tmp(out_dir)
  settings <- recovery_primary_settings()
  recovery_write(plan, file.path(out_dir, "recovery_plan.csv"))
  manifest <- data.frame(run_profile = settings$profile, smoke = settings$smoke,
                         chains = settings$chains, iter = settings$iter, warmup = settings$warmup,
                         adapt_delta = settings$control$adapt_delta,
                         max_treedepth = settings$control$max_treedepth,
                         claim = "never claim recovery from smoke settings", stringsAsFactors = FALSE)
  recovery_write(manifest, file.path(out_dir, "recovery_run_manifest.csv"))
  dgp_fun <- recovery_load_dgp(dgp_path)
  diag_rows <- list(); k_rows <- list()
  for (i in seq_len(nrow(plan))) {
    case <- plan[i, , drop = FALSE]
    panel <- tryCatch(recovery_make_dgp(case, dgp_fun), error = function(e) e)
    if (inherits(panel, "condition")) {
      diag_rows[[length(diag_rows) + 1L]] <- recovery_fit_diagnostics(NULL, case, list(), 0,
        conditionMessage(panel), "")
      k_rows[[length(k_rows) + 1L]] <- recovery_posterior_k(NULL, case, list())
      next
    }
    valid <- tryCatch(recovery_validate_dgp(panel, case$truth_structure),
                      error = function(e) e)
    if (inherits(valid, "condition")) {
      diag_rows[[length(diag_rows) + 1L]] <- recovery_fit_diagnostics(NULL, case, panel, 0,
        paste0("DGP validation failed: ", conditionMessage(valid)), "")
      k_rows[[length(k_rows) + 1L]] <- recovery_posterior_k(NULL, case, panel)
      next
    }
    result <- recovery_run_one(panel, case, out_dir, settings)
    diag_rows[[length(diag_rows) + 1L]] <- result$diagnostics
    k_rows[[length(k_rows) + 1L]] <- result$posterior_k
  }
  diagnostics <- do.call(rbind, diag_rows)
  posterior_k <- do.call(rbind, k_rows)
  recovery_write(diagnostics, file.path(out_dir, "recovery_diagnostics.csv"))
  recovery_write(posterior_k, file.path(out_dir, "posterior_k.csv"))
  has_error <- !is.na(diagnostics$error) & nzchar(diagnostics$error)
  recovery_write(diagnostics[has_error, , drop = FALSE],
                 file.path(out_dir, "recovery_errors.csv"))
  invisible(list(plan = plan, diagnostics = diagnostics, posterior_k = posterior_k,
                 output = out_dir))
}

if (sys.nframe() == 0L) {
  cli <- recovery_parse_cli()
  if (isTRUE(cli$help)) {
    recovery_usage()
    quit(save = "no", status = 0L)
  }
  if (is.null(cli$out) && isTRUE(cli$run)) stop("--run requires --out=/tmp/...")
  seed <- or_null(cli$seed, c(1103L, 2207L, 3301L))
  layout <- cli$layout
  if (is.null(cli$layout) && !cli$single_case) layout <- NULL
  plan <- recovery_make_plan(cli$demand, cli$conduct, cli$structure, seed, layout,
                             truth_structure = cli$truth_structure,
                             single_case = cli$single_case)
  if (!isTRUE(cli$run)) {
    print(plan, row.names = FALSE)
    cat("PRIMARY PLAN ONLY: no HMC, fit object, log, or recovery claim. Use --run with --out=/tmp/... to sample.\n")
    quit(save = "no", status = 0L)
  }
  recovery_run(plan, cli$out)
  cat("PRIMARY RUN COMPLETE: compact diagnostics written under ", normalizePath(cli$out),
      "; recovery remains subject to all selected cases and gates.\n", sep = "")
}
