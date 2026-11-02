#!/usr/bin/env Rscript
# Prespecified prior diagnostics; no fitting and no changes to priors.
suppressPackageStartupMessages(library(antitrustBayes))
root <- Sys.getenv("MQ_REPO", "/mnt/d/projects/antitrustbayes")
out <- Sys.getenv("MQ_PRIOR_OUT", file.path(root,
  "doc/statistical_toolkit_runs/20260916/performance-audit/recovery/market-quality/prior-predictive"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
prior <- bayes_prior_manifest("baseline")
set.seed(1729L)
n <- 20000L
# Exact inverse-CDF draws from the existing positive-truncated Normal(1, sd).
lower <- pnorm(0, 1, prior$prior_sigma_meanval_strat)
tau <- qnorm(lower+(1-lower)*runif(n), 1, prior$prior_sigma_meanval_strat)
sigma_share <- exp(rnorm(n, log(prior$prior_sigma_share), prior$prior_log_sigma_share_sd))
summarize <- function(x, groups, quantity) {
  q <- quantile(x, c(.025, .5, .975))
  odds_mean_infinite <- identical(quantity, "relative_share_odds_pair") &&
    prior$prior_sigma_meanval_strat >= 1/sqrt(2)
  data.frame(groups = groups, quantity = quantity,
    mean = if (odds_mean_infinite) NA_real_ else mean(x),
    mean_status = if (odds_mean_infinite) "infinite under baseline prior" else "Monte Carlo estimate",
    lower_95 = q[[1L]], median = q[[2L]], upper_95 = q[[3L]])
}
rows <- list(summarize(tau, NA_integer_, "common_quality_scale"),
  summarize(sigma_share, NA_integer_, "share_residual_scale"))
for (G in c(1L, 2L, 6L, 39L)) {
  independent <- matrix(rnorm(n*G), n, G)*tau
  quality <- independent-rowMeans(independent)
  stopifnot(max(abs(rowSums(quality))) < 1e-12)
  rows[[length(rows)+1L]] <- summarize(quality[, 1L], G, "quality_of_one_group")
  if (G > 1L) {
    odds <- exp(quality[, 1L]-quality[, 2L])
    rows[[length(rows)+1L]] <- summarize(odds, G, "relative_share_odds_pair")
  }
  # Equal price/other utility terms isolate quality and observation dispersion.
  shares <- exp(quality-apply(quality, 1L, max))
  shares <- shares/rowSums(shares)
  observed_score <- quality+matrix(rnorm(n*G), n, G)*sigma_share
  observed <- exp(observed_score-apply(observed_score, 1L, max))
  observed <- observed/rowSums(observed)
  rows[[length(rows)+1L]] <- summarize(shares[, 1L], G, "latent_relative_share")
  rows[[length(rows)+1L]] <- summarize(observed[, 1L]-shares[, 1L], G,
    "compositional_perturbation_share_difference")
}
write.csv(do.call(rbind, rows), file.path(out, "prior-predictive.csv"), row.names = FALSE)
saveRDS(list(seed = 1729L, n = n, prior = prior,
  interpretation = "Equal prices and other demand terms. Compositional perturbation is an interpretable diagnostic, not the production log-share likelihood."),
  file.path(out, "provenance.rds"))
cat("Saved prior-predictive diagnostics. No priors were changed.\n")
