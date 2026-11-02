#!/usr/bin/env Rscript
# Frozen checkpointed market-quality comparison. Preparation does not fit.
options(stringsAsFactors = FALSE)
mq_repo <- normalizePath(Sys.getenv("MQ_REPO", "/mnt/d/projects/antitrustbayes"))
mq_root <- Sys.getenv("MQ_BENCH_ROOT", file.path(mq_repo,
  "doc/statistical_toolkit_runs/20260916/performance-audit/recovery/market-quality/benchmark-v2"))
mq_args <- commandArgs(trailingOnly = TRUE)
mq_ref <- new.env(parent = globalenv())
sys.source(file.path(mq_repo, "doc/market_quality_reference.R"), mq_ref)
mq_foc <- new.env(parent = globalenv())
mq_functions <- c("cga_supply_name", "cga_role_conduct_logit", "cga_role_inverse_lerner_ces")
for (expr in parse(file.path(mq_repo, "doc/collapsed_gaussian_benchmark.R"))) {
  if (is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
      is.symbol(expr[[2L]]) && as.character(expr[[2L]]) %in% mq_functions)
    eval(expr, mq_foc)
}
mq_hash <- function(x) digest::digest(x, algo = "sha256")
mq_save_once <- function(x, path) {
  if (file.exists(path)) {
    if (!identical(mq_hash(readRDS(path)), mq_hash(x)))
      stop("Existing frozen artifact differs: ", path)
    return(invisible(path))
  }
  tmp <- paste0(path, ".tmp")
  saveRDS(x, tmp); stopifnot(file.rename(tmp, path)); invisible(path)
}
mq_fixture <- function(M, T = 2L, F = 6L, demand = "logit", conduct = "bertrand",
                       structure = "simultaneous", seed = 271828L) {
  set.seed(seed)
  C <- mq_ref$helmert_basis_local(F)
  quality <- matrix(0, M, F)
  for (m in seq_len(M)) quality[m, ] <- as.numeric(C %*% rnorm(F-1L, 0, .15))
  market <- rep(seq_len(M), each = T*F)
  year <- rep(rep(seq_len(T), each = F), M)
  firm <- rep(seq_len(F), M*T)
  # Bounded variation across markets: no price trend growing with M.
  prices <- 1.06 + .085*firm + .025*year + .04*sin(market/3) +
    .035*cos(firm*year + market/5)
  alpha <- 6; gamma <- 3.5; p0 <- 1.3
  # Keep the true outside share constant within each market, so the routine
  # market-level outside-share pool is correctly specified as well. Preserve
  # product-specific price variation; offset only the common period level.
  for (m in seq_len(M)) {
    first <- (m-1L)*T*F+seq_len(F)
    reference <- mq_ref$logsumexp(quality[m, ] - if (demand == "logit")
      alpha*prices[first] else gamma*log(prices[first]))
    if (T > 1L) for (tt in 2:T) {
      ii <- (m-1L)*T*F+(tt-1L)*F+seq_len(F)
      current <- mq_ref$logsumexp(quality[m, ] - if (demand == "logit")
        alpha*prices[ii] else gamma*log(prices[ii]))
      if (demand == "logit") prices[ii] <- prices[ii]+(current-reference)/alpha
      else prices[ii] <- prices[ii]*exp((current-reference)/gamma)
    }
  }
  std <- if (demand == "logit") (prices-mean(prices))/sd(prices) else
    (log(prices)-mean(log(prices)))/sd(log(prices))
  a <- if (demand == "logit") alpha*sd(prices) else (gamma-1)*sd(log(prices))
  latent <- numeric(length(prices)); observed <- latent; s0 <- numeric(M*T)
  share_noise <- latent
  for (m in seq_len(M)) for (tt in seq_len(T)) {
    j <- (m-1L)*T+tt; ii <- (j-1L)*F+seq_len(F)
    score <- if (demand == "logit") quality[m, ]-a*std[ii] else
      quality[m, ]-(a+sd(log(prices)))*std[ii]+gamma*log(p0)
    s0[j] <- plogis(-mq_ref$logsumexp(score))
    latent[ii] <- score-mq_ref$logsumexp(score)
    noise <- rnorm(F, 0, .04)
    observed[ii] <- score+noise-mq_ref$logsumexp(score+noise)
    share_noise[ii] <- observed[ii]-latent[ii]
  }
  inv_margin <- numeric(length(prices)); fixed_rank <- matrix(0L, M, F)
  for (m in seq_len(M)) {
    rows <- (m-1L)*T*F+seq_len(T*F)
    metric <- exp(observed[rows])
    if (demand == "ces") {
      metric <- prices[rows]*metric
      metric <- as.numeric(apply(matrix(metric, T, F, byrow = TRUE), 1L,
                                  function(v) v/sum(v)))
    }
    mean_metric <- colMeans(matrix(metric, T, F, byrow = TRUE))
    ordering <- order(-mean_metric, seq_len(F))
    fixed_rank[m, ordering] <- seq_len(F)
    for (tt in seq_len(T)) {
      j <- (m-1L)*T+tt; ii <- (j-1L)*F+seq_len(F)
      r <- exp(observed[ii]); p <- prices[ii]
      if (demand == "logit") inv_margin[ii] <- p*alpha*
        mq_foc$cga_role_conduct_logit(r, s0[j], conduct, structure,
                                    fixed_rank[m, ], 2L)
      else {
        rev <- p*r/sum(p*r)
        R0 <- plogis(qlogis(s0[j])+log(p0)-log(sum(p*r)))
        inv_margin[ii] <- mq_foc$cga_role_inverse_lerner_ces(rev, R0, gamma,
          conduct, structure, fixed_rank[m, ], 2L)
      }
    }
  }
  # Known compositional perturbation; recorded log shares are not independent
  # Gaussian outcomes. Margin innovations have prespecified correlation with
  # the realized normalized share disturbance, not unnormalized noise.
  margin_noise <- .30*(.06/.04)*share_noise + rnorm(length(prices), 0, .06*sqrt(1-.3^2))
  margins <- exp(-(log(inv_margin)+margin_noise))
  outside_time_error <- max(vapply(seq_len(M), function(m) {
    values <- s0[(m-1L)*T+seq_len(T)]
    max(abs(values-values[[1L]]))
  }, numeric(1)))
  stopifnot(all(is.finite(margins)), all(margins > 0), all(margins < 1),
            all(abs(rowSums(quality)) < 1e-13), outside_time_error < 1e-13)
  dat <- data.frame(prices = prices, quantities = 1000*exp(observed), margins = margins,
    market = paste0("m", market), year = as.character(year),
    group = paste0("firm", firm), margin_observed = seq_along(prices) %% 11L != 0L)
  list(data = dat, truth = list(quality = quality, alpha = alpha, gamma = gamma,
    s0 = s0, tau = .15, innovation_scale = .04, sigma_margin = .06, rho = .30,
    latent_log_share = latent, observed_log_share = observed, fixed_rank = fixed_rank,
    true_k = if (structure == "simultaneous") NULL else 2L,
    outside_time_error = outside_time_error),
    M = M, T = T, F = F, demand = demand, conduct = conduct,
    structure = structure, seed = seed, outside_price = p0)
}
mq_plan <- function() {
  cases <- data.frame(fixture = c("B", "C", "BF", "BS", "CES_B", "routine_auto_B"),
    demand = c(rep("logit", 4L), "ces", "logit"),
    conduct = c("bertrand", "cournot", rep("bertrand", 4L)),
    structure = c("simultaneous", "simultaneous", "fringe", "stackelberg",
                  "simultaneous", "simultaneous"),
    automatic_pooling = c(rep(FALSE, 5L), TRUE))
  plan <- merge(merge(cases, data.frame(seed = c(271828L, 314159L))),
    data.frame(mode = c("legacy", "expanded_market", "collapsed_market")))
  plan$case <- paste(plan$fixture, plan$mode, plan$seed, sep = "_")
  plan$M <- 3L; plan$T <- 2L; plan$F <- 6L
  plan$chains <- 4L; plan$iter <- 1000L; plan$warmup <- 500L; plan$cores <- 1L
  plan$adapt_delta <- .95; plan$max_treedepth <- 13L
  plan$market_s0_parameterization <- "legacy"
  plan$use_rho <- TRUE; plan$use_market_resid_scale <- FALSE
  plan$prior_manifest <- "baseline"; plan$df <- 5
  plan[order(plan$fixture, plan$seed, plan$mode), ]
}
mq_prepare <- function() {
  dir.create(mq_root, recursive = TRUE, showWarnings = FALSE)
  plan <- mq_plan()
  for (j in seq_len(nrow(plan))) {
    p <- plan[j, ]
    fixture <- mq_fixture(p$M, p$T, p$F, p$demand, p$conduct, p$structure, p$seed)
    path <- file.path(mq_root, paste0(p$fixture, "_", p$seed, "-fixture.rds"))
    mq_save_once(fixture, path)
    plan$data_sha256[j] <- mq_hash(fixture$data)
  }
  mq_save_once(plan, file.path(mq_root, "frozen-plan.rds"))
  write.csv(plan, file.path(mq_root, "frozen-plan.csv"), row.names = FALSE)
  cat("Prepared", nrow(plan), "frozen fits; no sampling performed.\n")
}
mq_fit_one <- function(case) {
  suppressPackageStartupMessages(library(antitrustBayes))
  stopifnot(antitrustBayes:::.antitrust_bayes_model_spec_version == 10L)
  gate <- file.path(dirname(mq_root), "native-validation", "approved-gate.rds")
  if (!file.exists(gate) || !isTRUE(readRDS(gate)$approved))
    stop("Native density/gradient and independent review gate not approved.")
  plan <- readRDS(file.path(mq_root, "frozen-plan.rds"))
  p <- plan[plan$case == case, ]
  if (nrow(p) != 1L) stop("Unknown or ambiguous case: ", case)
  path <- file.path(mq_root, paste0(case, "-fit.rds"))
  if (file.exists(path)) stop("Fit already saved; inspect and preserve it: ", path)
  started_path <- file.path(mq_root, paste0(case, "-started.rds"))
  result_path <- file.path(mq_root, paste0(case, "-result.rds"))
  if (file.exists(started_path) || file.exists(result_path))
    stop("An earlier attempt exists. Preserve and assess it; run any operational retry in a separately labeled artifact directory.")
  fx <- readRDS(file.path(mq_root, paste0(p$fixture, "_", p$seed, "-fixture.rds")))
  stopifnot(identical(mq_hash(fx$data), p$data_sha256))
  legacy <- p$mode == "legacy"
  args <- c(as.list(fx$data), list(demand_model = p$demand, supply_model = p$conduct,
    market_side = "output", use_rho = TRUE,
    share_likelihood = if (legacy) "student_t" else "gaussian",
    margin_likelihood = if (legacy) "student_t" else "gaussian",
    quality_parameterization = if (p$mode == "collapsed_market") "collapsed" else "sampled",
    quality_scope = if (legacy) "global" else "market", df = 5,
    use_market_resid_scale = FALSE, use_alpha_market_pool = p$automatic_pooling,
    s0_pooling = if (p$automatic_pooling) "auto" else "none",
    market_s0_parameterization = p$market_s0_parameterization,
    prior_manifest = "baseline", outside_price = fx$outside_price,
    structure = p$structure, likelihood = "joint",
    stan_options = list(chains = p$chains, iter = p$iter, warmup = p$warmup,
      cores = p$cores, seed = p$seed, init = 0, refresh = 0,
      control = list(adapt_delta = p$adapt_delta, max_treedepth = p$max_treedepth,
                     metric = "diag_e"))))
  if (p$structure != "simultaneous") {
    args$partition_prior <- "ordered_k"
    args$k_support <- if (p$structure == "fringe") 0:4 else 1:4
    args$k_prior <- rep(1/length(args$k_support), length(args$k_support))
  }
  # Retain partially written native chains if WSL stops during sampling.
  args$stan_options$sample_file <- file.path(mq_root, paste0(case, "-stan.csv"))
  driver_snapshot <- file.path(mq_root, paste0(case, "-driver.R"))
  stopifnot(file.copy(file.path(mq_repo, "doc/market_quality_benchmark.R"), driver_snapshot,
    overwrite = FALSE))
  kernels <- get("stanmodels", envir = asNamespace("antitrustBayes"))
  model_hashes <- vapply(kernels, function(m) digest::digest(m@model_code,
    algo = "sha256", serialize = FALSE), character(1))
  mq_save_once(list(plan = p, args = args, started = Sys.time(), PID = Sys.getpid(),
    driver_sha256 = digest::digest(driver_snapshot, algo = "sha256", file = TRUE),
    model_code_sha256 = model_hashes, R_version = R.version.string), started_path)
  warnings <- character(); error <- ""
  start <- proc.time()[["elapsed"]]
  fit <- tryCatch(withCallingHandlers(do.call(bayes_est, args), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  }), error = function(e) {error <<- conditionMessage(e); NULL})
  elapsed <- proc.time()[["elapsed"]]-start
  if (!is.null(fit)) saveRDS(fit, path)
  peak <- readLines("/proc/self/status", warn = FALSE)
  peak_kb <- as.numeric(sub("^VmHWM:\\s+([0-9]+).*$", "\\1", peak[grepl("^VmHWM:", peak)]))
  saveRDS(list(case = case, elapsed = elapsed, warnings = warnings, error = error,
    native_seconds = if (!is.null(fit)) sum(rstan::get_elapsed_time(fit@stanfit)) else NA_real_,
    peak_RSS_kb = peak_kb, completed = Sys.time()), result_path)
  if (nzchar(error)) stop(error)
  cat(case, "saved; wall seconds:", elapsed, "\n")
}
if ("--prepare" %in% mq_args) mq_prepare()
case_arg <- mq_args[startsWith(mq_args, "--case=")]
if (length(case_arg)) mq_fit_one(sub("^--case=", "", case_arg[[1L]]))
