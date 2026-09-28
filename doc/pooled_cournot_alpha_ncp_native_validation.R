#!/usr/bin/env Rscript

## Validation gate for non-centering the pooled Cournot alpha hierarchy in
## inst/stan/logit_supply.stan and inst/stan/logit_stackelberg.stan.
##
## Part A (pure R, no Stan): the OLD and NEW Stan formulas are transcribed
## exactly and checked against the hand-derived exact identity
##   log_prob_new(r) - log_prob_old(kappa) == N_market * log(sigma_log_a)
## at a matched economic state (same a_market), for a grid of sigma_log_a
## values approaching zero -- the regime where the centered OLD coordinate
## funnels and the noncentered NEW coordinate does not.
##
## Part B (real CmdStanR): compiles the NEW working-tree kernels and runs a
## fixed_param evaluation with explicit constrained inits, to (1) confirm
## they compile and execute without error, and (2) cross-check Stan's own
## computed a_market/log_a_market_derived against the Part-A R formulas, so
## the analytic identity is tied to the actual compiled arithmetic and not
## just to a hand transcription. Also checks use_centered_alpha == 1 gives
## the same a_market as use_centered_alpha == 0 at a matched state -- an
## internal self-consistency check within the NEW kernel alone.
##
## No HMC; the posterior-equivalence smoke test against the OLD kernel is a
## separate script, doc/pooled_cournot_alpha_ncp_hmc_smoke.R.

options(stringsAsFactors = FALSE)
lib <- Sys.getenv("PCA_LIB", unset = "")
if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(antitrustBayes))
suppressPackageStartupMessages(library(cmdstanr))

script_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
outdir <- if (length(script_arg) && nzchar(script_arg)) dirname(normalizePath(script_arg)) else "doc"

## ---------------------------------------------------------------------
## Part A: pure-R analytic identity check.
## ---------------------------------------------------------------------

## OLD Stan formulas (inst/stan/logit_supply.stan before this change):
##   log_a_market_derived = log_kappa_market - mean_log_s0_market
##   a_market = exp(log_a_market_derived)
##   target += normal_lpdf(log_a_market_derived | mu_log_a, sigma_log_a)
## kappa (log_kappa_market) is the free sampled parameter; log_a_derived is a
## deterministic (unit-Jacobian) shift of it by the data-dependent constant
## c = mean_log_s0_market, so the density formula above is, term for term,
## the density of kappa under kappa ~ N(mu_log_a + c, sigma_log_a).
old_log_prob <- function(kappa, mu, sigma, c) {
  log_a <- kappa - c
  list(a_market = exp(log_a),
       log_a_market_derived = log_a,
       target = sum(dnorm(log_a, mean = mu, sd = sigma, log = TRUE)))
}

## NEW Stan formulas (this change), noncentered branch:
##   log_a_market_derived = mu_log_a + sigma_log_a * r_market_a_raw
##   a_market = exp(log_a_market_derived)
##   target += std_normal_lpdf(r_market_a_raw)
new_log_prob_noncentered <- function(r, mu, sigma) {
  log_a <- mu + sigma * r
  list(a_market = exp(log_a),
       log_a_market_derived = log_a,
       target = sum(dnorm(r, mean = 0, sd = 1, log = TRUE)))
}

## NEW Stan formulas, centered branch (use_centered_alpha == 1):
##   log_a_market_derived = log_a_market_pool
##   target += normal_lpdf(log_a_market_pool | mu_log_a, sigma_log_a)
new_log_prob_centered <- function(log_a_pool, mu, sigma) {
  list(a_market = exp(log_a_pool),
       log_a_market_derived = log_a_pool,
       target = sum(dnorm(log_a_pool, mean = mu, sd = sigma, log = TRUE)))
}

rows <- list()
add <- function(scenario, quantity, error, tolerance, detail) {
  rows[[length(rows) + 1L]] <<- data.frame(
    scenario = scenario, quantity = quantity, max_abs_error = error,
    tolerance = tolerance, passed = is.finite(error) && error <= tolerance,
    method = "analytic_r_transcription", detail = detail
  )
}

grid <- expand.grid(
  sigma = c(0.6, 0.15, 0.02, 0.005, 0.0005),
  mu = c(log(0.84), -0.2, 0.35)
)
r_vec <- c(-1.5, 0.3, -0.7, 1.1, -0.4, 0.9)
c_vec <- c(0.38, -0.07, 0.21, 0.15, -0.31, 0.44)

for (i in seq_len(nrow(grid))) {
  sigma <- grid$sigma[i]; mu <- grid$mu[i]
  ## Matched economic state: choose kappa so that log_a matches exactly.
  log_a <- mu + sigma * r_vec
  kappa <- log_a + c_vec
  old <- old_log_prob(kappa, mu, sigma, c_vec)
  new <- new_log_prob_noncentered(r_vec, mu, sigma)
  scenario <- sprintf("sigma=%.4g_mu=%.3f", sigma, mu)

  add(scenario, "a_market_matched", max(abs(old$a_market - new$a_market)), 1e-12,
      "matched-state a_market must agree between the two parameterizations")
  expected_gap <- length(r_vec) * log(sigma)
  add(scenario, "log_prob_gap", abs((new$target - old$target) - expected_gap),
      max(1e-9, 1e-9 * abs(expected_gap)),
      "log_prob_new(r) - log_prob_old(kappa) == N_market * log(sigma_log_a)")

  ## Centered branch: log_a_market_pool set to the same matched log_a must
  ## reproduce OLD's density exactly (zero gap, no Jacobian -- both are
  ## centered coordinates related by only a fixed unit-Jacobian shift).
  centered <- new_log_prob_centered(log_a, mu, sigma)
  add(scenario, "centered_log_prob_gap", abs(centered$target - old$target), 1e-9,
      "use_centered_alpha == 1 at the same log_a must equal OLD's density exactly")
}

target_rows <- do.call(rbind, rows)

## ---------------------------------------------------------------------
## Part B: compile and run the real NEW kernels.
## ---------------------------------------------------------------------

kernels <- c("logit_supply", "logit_stackelberg")
repo <- normalizePath(file.path(outdir, ".."), mustWork = TRUE)

make_fixture <- function() {
  M <- 3L; T <- 2L; J <- 3L; L <- M * T
  list(
    prices = rep(c(1.05, 1.28, 1.52), L) + rep(c(0, .01, -.01), L),
    quantities = rep(c(32, 21, 14), L) + rep(c(0, 1, -1), L),
    margins = rep(c(.30, .27, .33), L),
    market = rep(rep(1:M, each = T), each = J),
    year = rep(rep(1:T, times = M), each = J),
    group = rep(c(1L, 1L, 2L), L),
    margin_observed = rep(c(TRUE, FALSE, TRUE), L),
    loo_cluster = rep(rep(1:M, each = T), each = J)
  )
}
fx <- make_fixture()

prepare_data <- function(structure) {
  common <- list(
    prices = fx$prices, quantities = fx$quantities, margins = fx$margins,
    market = fx$market, year = fx$year, group = fx$group,
    supply_model = 3L, use_cutoff = FALSE, market_side = "output",
    use_rho = TRUE, use_market_resid_scale = FALSE,
    use_alpha_market_pool = TRUE, use_s0_market_year_pool = FALSE,
    use_s0_market_pool = TRUE, use_centered_alpha = FALSE,
    use_centered_s0_pool = FALSE, market_s0_parameterization = "legacy",
    use_hmt = FALSE, priors = list(), X_s0 = NULL,
    margin_observed = fx$margin_observed, outside_price = 1,
    likelihood_mode = 1L, loo_cluster = fx$loo_cluster,
    estimate_df = FALSE, resid_df_fixed = 5,
    structure = structure,
    partition_prior = if (identical(structure, "stackelberg")) "ordered_k" else "cutoff",
    k_support = if (identical(structure, "stackelberg")) 1:2 else NULL,
    k_prior = if (identical(structure, "stackelberg")) c(.4, .6) else NULL,
    k_scope = "market", share_likelihood = "gaussian",
    margin_likelihood = "gaussian",
    quality_parameterization = "collapsed", quality_scope = "market"
  )
  do.call(antitrustBayes:::prepare_stan_data, common)
}

structures <- list(logit_supply = "simultaneous", logit_stackelberg = "stackelberg")
mu <- log(0.84); sigma <- 0.043; r <- c(-1.2, 0.6, -0.3)
log_a_pool <- mu + sigma * r  ## matched centered init

compiled_rows <- list()
addc <- function(kernel, quantity, error, tolerance, detail) {
  compiled_rows[[length(compiled_rows) + 1L]] <<- data.frame(
    scenario = paste0(kernel, ":compiled"), quantity = quantity,
    max_abs_error = error, tolerance = tolerance,
    passed = is.finite(error) && error <= tolerance,
    method = "native_cmdstanr_fixed_param", detail = detail
  )
}

for (kernel in kernels) {
  data <- prepare_data(structures[[kernel]])
  model <- cmdstanr::cmdstan_model(
    file.path(repo, "inst", "stan", paste0(kernel, ".stan")),
    force_recompile = TRUE,
    cpp_options = list(PRECOMPILED_HEADERS = "false", stan_threads = TRUE)
  )

  init_noncentered <- list(mu_log_a_free = mu, sigma_log_a_vec = sigma,
                           r_market_a_raw = r)
  fit_nc <- model$sample(data = data, chains = 1L, threads_per_chain = 1L,
                         iter_warmup = 0L,
                         iter_sampling = 1L, fixed_param = TRUE,
                         init = list(init_noncentered), seed = 190921L,
                         sig_figs = 18, refresh = 0)
  a_market_nc <- as.numeric(fit_nc$draws(variables = "a_market"))
  log_a_nc <- as.numeric(fit_nc$draws(variables = "log_a_market_derived"))
  expected_log_a <- mu + sigma * r
  addc(kernel, "noncentered_log_a_market_derived",
       max(abs(log_a_nc - expected_log_a)), 1e-9,
       "compiled kernel's log_a_market_derived must equal mu + sigma*r")
  addc(kernel, "noncentered_a_market", max(abs(a_market_nc - exp(expected_log_a))),
       1e-9, "compiled kernel's a_market must equal exp(mu + sigma*r)")

  data$use_centered_alpha <- 1L
  init_centered <- list(mu_log_a_free = mu, sigma_log_a_vec = sigma,
                        log_a_market_pool = log_a_pool)
  fit_c <- model$sample(data = data, chains = 1L, threads_per_chain = 1L,
                        iter_warmup = 0L,
                        iter_sampling = 1L, fixed_param = TRUE,
                        init = list(init_centered), seed = 190922L,
                        sig_figs = 18, refresh = 0)
  a_market_c <- as.numeric(fit_c$draws(variables = "a_market"))
  addc(kernel, "centered_matches_noncentered_a_market",
       max(abs(a_market_c - a_market_nc)), 1e-9,
       "use_centered_alpha == 1 at a matched log_a_market_pool must give the same a_market")

  ## Non-Cournot pooled state must be unaffected (regression on the
  ## untouched B/A/MC branch). Use this kernel's own structure so the data
  ## block matches what was compiled (logit_stackelberg always needs
  ## exact-k data regardless of conduct model).
  data_bertrand <- prepare_data(structures[[kernel]])
  data_bertrand$supply_model <- 1L
  fit_b <- model$sample(data = data_bertrand, chains = 1L, threads_per_chain = 1L,
                        iter_warmup = 0L,
                        iter_sampling = 1L, fixed_param = TRUE,
                        init = list(init_noncentered), seed = 190923L,
                        sig_figs = 18, refresh = 0)
  a_market_b <- as.numeric(fit_b$draws(variables = "a_market"))
  addc(kernel, "bertrand_regression_a_market",
       max(abs(a_market_b - exp(expected_log_a))), 1e-9,
       "the untouched is_cournot == 0 branch must be unaffected by this change")
}

compiled_target_rows <- do.call(rbind, compiled_rows)
all_rows <- rbind(target_rows, compiled_target_rows)

write.csv(all_rows, file.path(outdir, "pooled_cournot_alpha_ncp_target_checks.csv"),
          row.names = FALSE)

if (any(!all_rows$passed)) {
  print(subset(all_rows, !passed))
  stop("pooled-Cournot noncentering validation failed")
}
cat(sprintf("passed %d analytic and %d compiled-kernel checks\n",
            nrow(target_rows), nrow(compiled_target_rows)))
