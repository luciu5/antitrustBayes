#!/usr/bin/env Rscript

## Posterior-equivalence HMC smoke test for the pooled Cournot alpha
## noncentering (see doc/pooled_cournot_alpha_ncp.md). Fits the OLD (centered
## log_kappa_market) and NEW (noncentered r_market_a_raw) kernels on the same
## pooled-Cournot fixture, same seed, same short HMC budget, and compares
## economic posterior means within combined MCSE (protocol of
## doc/market_intercept_enterprise/local_pilot_results.md). Divergences in
## the OLD run are expected -- that is the motivation for this change -- and
## are reported, not gated on.

options(stringsAsFactors = FALSE)
lib <- Sys.getenv("PCA_LIB", unset = "")
if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(antitrustBayes))
suppressPackageStartupMessages(library(cmdstanr))
suppressPackageStartupMessages(library(posterior))

script_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
outdir <- if (length(script_arg) && nzchar(script_arg)) dirname(normalizePath(script_arg)) else "doc"
repo <- normalizePath(file.path(outdir, ".."), mustWork = TRUE)
old_commit <- Sys.getenv("PCA_OLD_COMMIT", "84f970d230053b7b3e7a0d6b98a8c49c2c761b7b")
build_dir <- Sys.getenv("PCA_BUILD_DIR", file.path(tempdir(), "pooled-cournot-ncp-hmc"))
dir.create(build_dir, recursive = TRUE, showWarnings = FALSE)

checkout_old <- function(kernel) {
  dest <- file.path(build_dir, "old", paste0(kernel, ".stan"))
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  src <- system2("git", c("-C", shQuote(repo), "show",
                          paste0(old_commit, ":inst/stan/", kernel, ".stan")),
                stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(src, "status"))) {
    stop("Cannot check out old ", kernel, ".stan at ", old_commit)
  }
  writeLines(src, dest)
  dest
}
compile_kernel <- function(path) {
  cmdstanr::cmdstan_model(path, force_recompile = TRUE,
                          cpp_options = list(PRECOMPILED_HEADERS = "false",
                                            stan_threads = TRUE))
}

kernel <- "logit_supply"
model_old <- compile_kernel(checkout_old(kernel))
model_new <- compile_kernel(file.path(repo, "inst", "stan", paste0(kernel, ".stan")))

## Slightly larger fixture than the native-validation one, so the alpha
## hierarchy is identified but still weakly enough to be diagnostic (small
## sigma_log_a should be visitable).
make_fixture <- function() {
  M <- 6L; T <- 1L; J <- 4L; L <- M * T
  set.seed(190924L)
  list(
    prices = rep(c(1.02, 1.15, 1.30, 1.48), L) + rnorm(L * J, 0, 0.01),
    quantities = rep(c(30, 22, 18, 11), L) + rnorm(L * J, 0, 0.5),
    margins = rep(c(.28, .24, .31, .22), L) + rnorm(L * J, 0, 0.005),
    market = rep(seq_len(M), each = J),
    year = rep(1L, L * J),
    group = rep(1:J, L),
    margin_observed = rep(TRUE, L * J),
    loo_cluster = rep(seq_len(M), each = J)
  )
}
fx <- make_fixture()

prepare_data <- function() {
  antitrustBayes:::prepare_stan_data(
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
    structure = "simultaneous", partition_prior = "cutoff",
    k_scope = "market", share_likelihood = "gaussian",
    margin_likelihood = "gaussian",
    quality_parameterization = "collapsed", quality_scope = "market"
  )
}
data <- prepare_data()

fit_one <- function(model, seed) {
  model$sample(data = data, chains = 4L, threads_per_chain = 1L,
              parallel_chains = 4L, iter_warmup = 500L, iter_sampling = 500L,
              adapt_delta = 0.95, max_treedepth = 12L, init = 0,
              seed = seed, refresh = 0)
}

cat("Fitting OLD kernel (centered log_kappa_market)...\n")
fit_old <- fit_one(model_old, 20260928L)
cat("Fitting NEW kernel (noncentered r_market_a_raw)...\n")
fit_new <- fit_one(model_new, 20260928L)

diag <- function(fit, label) {
  sd <- fit$sampler_diagnostics(format = "draws_df")
  cat(sprintf("%s: divergences=%d treedepth_hits=%d\n", label,
              sum(sd$divergent__), sum(sd$treedepth__ >= 12L)))
}
diag(fit_old, "OLD")
diag(fit_new, "NEW")

vars <- c("a_market", "s0_prob", "sigma_log_a_vec", "mu_log_a")
old_summary <- summarise_draws(fit_old$draws(vars), "mean", "sd", "ess_bulk")
new_summary <- summarise_draws(fit_new$draws(vars), "mean", "sd", "ess_bulk")

old_mcse <- old_summary$sd / sqrt(old_summary$ess_bulk)
new_mcse <- new_summary$sd / sqrt(new_summary$ess_bulk)
combined_mcse <- sqrt(old_mcse^2 + new_mcse^2)
delta <- abs(old_summary$mean - new_summary$mean)
z <- delta / combined_mcse

report <- data.frame(
  variable = old_summary$variable,
  old_mean = old_summary$mean, new_mean = new_summary$mean,
  combined_mcse = combined_mcse, combined_mcse_z = z
)
print(report)

write.csv(report, file.path(outdir, "pooled_cournot_alpha_ncp_hmc_smoke.csv"),
          row.names = FALSE)

if (any(z > 3, na.rm = TRUE)) {
  cat("WARNING: at least one economic quantity differs by more than 3 combined MCSE.\n")
}
cat("HMC smoke test complete.\n")
