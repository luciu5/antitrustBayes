#!/usr/bin/env Rscript

# Extended statistical validation against an independently generated market.
# Run with an installed current antitrustBayes and CmdStan:
#   Rscript --vanilla doc/validation/known_truth_recovery.R
#
# Optional paths:
#   AB_RECOVERY_CMDSTAN=/writable/cmdstan
#   AB_RECOVERY_CACHE=/writable/stan-cache
#   AB_RECOVERY_OUTPUT=/path/to/results
# Set AB_RECOVERY_SCOPE=recovery for the shorter parameter-recovery gate.
# The default "full" scope also checks independent conduct model evidence
# and requires a writable CmdStan installation.
#
# This is deliberately outside ordinary testthat discovery. The fixed design
# and HMC budget make this a release/nightly validation, not a commit-time test.

if (!requireNamespace("antitrustBayes", quietly = TRUE) ||
    !requireNamespace("cmdstanr", quietly = TRUE)) {
  stop("Install antitrustBayes and CmdStanR before recovery validation.")
}

cmdstan_path <- Sys.getenv("AB_RECOVERY_CMDSTAN", "")
if (nzchar(cmdstan_path)) cmdstanr::set_cmdstan_path(cmdstan_path)
if (!nzchar(cmdstanr::cmdstan_path()))
  stop("Install CmdStan or set AB_RECOVERY_CMDSTAN.")
cache <- Sys.getenv("AB_RECOVERY_CACHE",
                    file.path(tempdir(), "antitrustBayes-stan-cache"))
dir.create(cache, recursive = TRUE, showWarnings = FALSE)
options(antitrustBayes.cmdstan_cache = cache)
out_dir <- Sys.getenv("AB_RECOVERY_OUTPUT",
                      file.path(tempdir(), "known-truth-recovery"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
scope <- Sys.getenv("AB_RECOVERY_SCOPE", "full")
if (!scope %in% c("full", "recovery"))
  stop("AB_RECOVERY_SCOPE must be 'full' or 'recovery'.")

log_sum_exp <- function(x) {
  peak <- max(x)
  peak + log(sum(exp(x - peak)))
}

# No antitrustBayes demand, conduct, or simulation helper enters this DGP.
# Firms own one product each. For Logit Bertrand output pricing the
# independently specified FOC is
#   1 / Lerner_f = price_f * alpha * (1 - (1 - s0) * inside_share_f).
make_fixture <- function(seed = 817L, n_market = 6L, n_year = 3L,
                         n_firm = 5L, alpha = 2.5) {
  set.seed(seed)
  quality <- matrix(0, n_market, n_firm)
  for (m in seq_len(n_market)) {
    effects <- c(-.52, -.27, 0, .25, .54) +
      rnorm(n_firm, 0, .035)
    quality[m, ] <- effects - mean(effects)
  }
  price <- array(0, c(n_market, n_year, n_firm))
  for (m in seq_len(n_market))
    for (tt in seq_len(n_year))
      for (f in seq_len(n_firm))
        price[m, tt, f] <- 1.14 + .16 * f + .025 * m + .055 * tt +
          .035 * sin(f * tt + m) + runif(1, -.025, .025)

  # A common within-market/year price shift holds the true outside share
  # fixed across years while product-relative price variation remains.
  for (m in seq_len(n_market)) for (tt in 2:n_year) {
    reference <- log_sum_exp(quality[m, ] - alpha * price[m, 1L, ])
    current <- log_sum_exp(quality[m, ] - alpha * price[m, tt, ])
    price[m, tt, ] <- price[m, tt, ] + (current - reference) / alpha
  }
  centered_price <- mean(price)
  rows <- vector("list", n_market * n_year)
  outside_by_year <- matrix(NA_real_, n_market, n_year)
  for (m in seq_len(n_market)) for (tt in seq_len(n_year)) {
    score <- quality[m, ] - alpha * (price[m, tt, ] - centered_price)
    outside_by_year[m, tt] <- stats::plogis(-log_sum_exp(score))
    share_noise <- rnorm(n_firm, 0, .025)
    inside_share <- exp(score + share_noise -
                          log_sum_exp(score + share_noise))
    inverse_margin <- price[m, tt, ] * alpha *
      (1 - (1 - outside_by_year[m, tt]) * inside_share)
    margin <- exp(-log(inverse_margin) + rnorm(n_firm, 0, .035))
    rows[[(m - 1L) * n_year + tt]] <- data.frame(
      prices = price[m, tt, ], quantities = 1000 * inside_share,
      margins = margin, market = paste0("m", m),
      year = paste0("t", tt), group = paste0("f", seq_len(n_firm))
    )
  }
  data <- do.call(rbind, rows)
  if (max(abs(outside_by_year - outside_by_year[, 1L])) > 1e-12)
    stop("Recovery market outside share changed across years.")
  true_s0 <- outside_by_year[, 1L]
  if (any(!is.finite(as.matrix(data[c("prices", "quantities", "margins")]))) ||
      any(data$margins <= 0 | data$margins >= 1) ||
      any(true_s0 <= .1 | true_s0 >= .3)) {
    stop("The prespecified synthetic market is invalid.")
  }
  cell <- interaction(data$market, data$year, drop = TRUE)
  if (max(abs(as.numeric(tapply(data$quantities, cell, sum)) - 1000)) > 1e-8)
    stop("Synthetic conditional inside shares do not sum to one.")
  list(data = data, alpha = alpha, s0 = true_s0,
       share_noise_sd = .025, margin_noise_sd = .035, seed = seed)
}

# A lower-dimensional evidence fixture isolates model selection. With one
# market the sampled market-intercept parameterization is a no-op, so Cournot
# cannot mimic Bertrand by raising a separate market outside share. The
# within-year price pattern still gives distinct observed firm shares.
make_evidence_fixture <- function(seed = 183L, n_year = 8L,
                                  n_firm = 5L, alpha = 2.5) {
  set.seed(seed)
  quality <- c(-.55, -.30, -.05, .22, .68)
  quality <- quality - mean(quality)
  price <- matrix(0, n_year, n_firm)
  for (tt in seq_len(n_year))
    for (f in seq_len(n_firm))
      price[tt, f] <- 1.12 + .16 * f + .05 * tt +
        .04 * sin(f * tt) + runif(1, -.02, .02)
  for (tt in 2:n_year)
    price[tt, ] <- price[tt, ] +
      (log_sum_exp(quality - alpha * price[tt, ]) -
         log_sum_exp(quality - alpha * price[1L, ])) / alpha

  center <- mean(price)
  rows <- vector("list", n_year)
  outside <- numeric(n_year)
  for (tt in seq_len(n_year)) {
    score <- quality - alpha * (price[tt, ] - center)
    outside[[tt]] <- stats::plogis(-log_sum_exp(score))
    share_noise <- rnorm(n_firm, 0, .02)
    inside_share <- exp(score + share_noise -
                          log_sum_exp(score + share_noise))
    inverse_margin <- price[tt, ] * alpha *
      (1 - (1 - outside[[tt]]) * inside_share)
    rows[[tt]] <- data.frame(
      prices = price[tt, ], quantities = 1000 * inside_share,
      margins = exp(-log(inverse_margin) + rnorm(n_firm, 0, .025)),
      market = "m1", year = paste0("t", tt),
      group = paste0("f", seq_len(n_firm))
    )
  }
  if (max(abs(outside - outside[[1L]])) > 1e-12)
    stop("Evidence market outside share changed across years.")
  list(data = do.call(rbind, rows), alpha = alpha, s0 = outside[[1L]],
       seed = seed)
}

fit_conduct <- function(fixture, supply_model, seed) {
  message("Fitting ", supply_model, " to the same known-truth market")
  do.call(antitrustBayes::bayes_est, c(as.list(fixture$data), list(
    demand_model = "logit", market_side = "output",
    supply_model = supply_model,
    share_likelihood = "gaussian", margin_likelihood = "gaussian",
    quality_parameterization = "collapsed", quality_scope = "market",
    market_intercept_parameterization = "sampled", use_rho = FALSE,
    use_alpha_market_pool = FALSE, s0_pooling = "market",
    stan_options = list(chains = 2L, iter = 500L, warmup = 250L,
                        cores = 2L, seed = seed, refresh = 0L)
  )))
}

fit_diagnostics <- function(fit, model) {
  summary <- fit@stanfit$summary(variables = c(
    "alpha_structural_market", "s0_prob", "sigma_share", "sigma_margin"
  ))
  sampler <- fit@stanfit$sampler_diagnostics(format = "matrix")
  data.frame(model = model, max_rhat = max(summary$rhat),
             min_bulk_ess = min(summary$ess_bulk),
             divergences = sum(sampler[, "divergent__"]),
             transitions = nrow(sampler))
}

fixture <- make_fixture()
bertrand <- fit_conduct(fixture, "bertrand", 867L)
diagnostics <- fit_diagnostics(bertrand, "Bertrand")

draw_matrix <- function(name, columns) {
  draws <- as.matrix(bertrand@stanfit$draws(variables = name,
                                          format = "matrix"))
  expected <- if (columns == 1L) name else
    paste0(name, "[", seq_len(columns), "]")
  if (!identical(colnames(draws), expected))
    stop("CmdStan output changed the expected ", name, " coordinates.")
  matrix(as.numeric(draws), nrow = nrow(draws), ncol = ncol(draws))
}
alpha_draws <- draw_matrix("alpha_structural_market", length(fixture$s0))
s0_draws <- draw_matrix("s0_prob", length(bertrand@stan_data$mky_to_market))
margin_scale_draws <- as.numeric(draw_matrix("sigma_margin", 1L))
market_map <- bertrand@stan_data$mky_to_market
if (!is.matrix(alpha_draws) || ncol(alpha_draws) != length(fixture$s0) ||
    !is.matrix(s0_draws) || ncol(s0_draws) != length(market_map) ||
    !identical(sort(unique(as.integer(market_map))),
               seq_along(fixture$s0))) {
  stop("Posterior dimensions no longer match the known-truth market.")
}
if (max(abs(alpha_draws - alpha_draws[, 1L])) > 1e-10)
  stop("A common structural alpha changed across markets.")
for (m in seq_along(fixture$s0)) {
  columns <- which(market_map == m)
  if (max(abs(s0_draws[, columns, drop = FALSE] -
              s0_draws[, columns[[1L]]])) > 1e-10)
    stop("Market-level outside share changed across years.")
}

summarize_truth <- function(draws, parameter, truth) {
  q <- as.numeric(stats::quantile(draws, c(.025, .5, .975)))
  data.frame(parameter = parameter, truth = truth, q025 = q[[1L]],
             median = q[[2L]], q975 = q[[3L]],
             abs_error = abs(q[[2L]] - truth),
             covered_95 = q[[1L]] <= truth && truth <= q[[3L]])
}
recovery <- summarize_truth(alpha_draws[, 1L], "alpha", fixture$alpha)
for (m in seq_along(fixture$s0)) {
  column <- which(market_map == m)[[1L]]
  recovery <- rbind(recovery, summarize_truth(
    s0_draws[, column], paste0("s0_m", m), fixture$s0[[m]]
  ))
}
recovery <- rbind(recovery, summarize_truth(
  margin_scale_draws, "sigma_margin", fixture$margin_noise_sd
))

utils::write.csv(recovery, file.path(out_dir, "recovery.csv"), row.names = FALSE)
utils::write.csv(diagnostics, file.path(out_dir, "diagnostics.csv"),
                 row.names = FALSE)
print(recovery, row.names = FALSE)
print(diagnostics, row.names = FALSE)
cat("Results:", normalizePath(out_dir), "\n")

# These thresholds are deliberately broad enough to allow posterior and
# bridge Monte Carlo error, but narrow enough to reject economically material
# drift in the structural price effect, outside shares, or model preference.
recovery_ok <- recovery$covered_95[[1L]] &&
  recovery$abs_error[[1L]] <= .08 * fixture$alpha &&
  sum(recovery$covered_95[2:7]) >= 5L &&
  mean(recovery$abs_error[2:7]) <= .06 &&
  recovery$median[[8L]] >= .02 && recovery$median[[8L]] <= .06
diagnostics_ok <- diagnostics$max_rhat[[1L]] <= 1.05 &&
  diagnostics$min_bulk_ess[[1L]] >= 100 &&
  diagnostics$divergences[[1L]] == 0
if (!isTRUE(recovery_ok) || !isTRUE(diagnostics_ok)) {
  stop("Known-truth recovery gate failed: recovery=", recovery_ok,
       ", diagnostics=", diagnostics_ok,
       ". Inspect the written results.", call. = FALSE)
}
cat("Known-truth parameter recovery gate PASSED.\n")
if (scope == "recovery") quit(save = "no", status = 0L)

evidence_fixture <- make_evidence_fixture()
evidence_bertrand <- fit_conduct(evidence_fixture, "bertrand", 367L)
evidence_cournot <- fit_conduct(evidence_fixture, "cournot", 368L)
evidence_diagnostics <- rbind(
  fit_diagnostics(evidence_bertrand, "Bertrand evidence"),
  fit_diagnostics(evidence_cournot, "Cournot evidence")
)
utils::write.csv(evidence_diagnostics,
                 file.path(out_dir, "evidence_diagnostics.csv"),
                 row.names = FALSE)

evidence_draws <- function(name) {
  draws <- as.matrix(evidence_bertrand@stanfit$draws(
    variables = name, format = "matrix"))
  matrix(as.numeric(draws), nrow = nrow(draws), ncol = ncol(draws))
}
evidence_s0_draws <- evidence_draws("s0_prob")
if (ncol(evidence_s0_draws) != length(unique(evidence_fixture$data$year)) ||
    max(abs(evidence_s0_draws - evidence_s0_draws[, 1L])) > 1e-10)
  stop("Evidence-fit market outside share changed across years.")
evidence_recovery <- rbind(
  summarize_truth(evidence_draws("alpha_structural_market")[, 1L],
                  "alpha", evidence_fixture$alpha),
  summarize_truth(evidence_s0_draws[, 1L],
                  "s0", evidence_fixture$s0)
)
utils::write.csv(evidence_recovery,
                 file.path(out_dir, "evidence_recovery.csv"),
                 row.names = FALSE)

message("Comparing normalized marginal likelihoods by bridge sampling")
set.seed(18815L)
comparison <- antitrustBayes::bayes_compare(
  Bertrand = evidence_bertrand, Cournot = evidence_cournot, bridge = TRUE,
  bridge_repetitions = 2L, bridge_maxiter = 2000L,
  prior_weights = c(.5, .5)
)
bridge_context <- attr(comparison, "bridge_context_info")
evidence <- comparison[, c("model", "logml", "logml_sd", "post_prob")]
utils::write.csv(evidence, file.path(out_dir, "evidence.csv"), row.names = FALSE)
print(evidence_recovery, row.names = FALSE)
print(evidence_diagnostics, row.names = FALSE)
print(evidence, row.names = FALSE)

# An independent, identifiable conduct design must show both a meaningful
# Bayes factor and stable repeated numerical estimates of that factor.
correct <- match("Bertrand", evidence$model)
wrong <- match("Cournot", evidence$model)
evidence_ok <- all(evidence_recovery$covered_95) &&
  evidence_recovery$abs_error[[1L]] <= .08 * evidence_fixture$alpha &&
  evidence_recovery$abs_error[[2L]] <= .05 &&
  all(evidence_diagnostics$max_rhat <= 1.05) &&
  all(evidence_diagnostics$min_bulk_ess >= 100) &&
  all(evidence_diagnostics$divergences == 0) &&
  length(bridge_context) == nrow(evidence) &&
  all(vapply(bridge_context, function(x)
    identical(x$status, "validated") &&
      identical(x$saved_hash, x$current_hash), logical(1))) &&
  all(is.finite(evidence$logml)) &&
  all(is.finite(evidence$logml_sd) & evidence$logml_sd <= .5) &&
  is.finite(correct) && is.finite(wrong) &&
  evidence$logml[[correct]] - evidence$logml[[wrong]] >= log(10) &&
  evidence$post_prob[[correct]] >= .9
if (!isTRUE(evidence_ok))
  stop("Known-truth model-evidence gate failed. Inspect evidence outputs.",
       call. = FALSE)
cat("Known-truth model-evidence gate PASSED.\n")
