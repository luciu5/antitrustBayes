#!/usr/bin/env Rscript

## Conditional identification audit for the prespecified Logit-Bertrand
## Stackelberg seed-1103/single-product DGP. This fixes the continuous demand
## parameters at truth and varies only the exact ordered-k state.

here <- local({
  args <- commandArgs(trailingOnly = FALSE)
  f <- args[startsWith(args, "--file=")]
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[[1L]]))) else getwd()
})
out_dir <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[[1L]] else here
source(file.path(here, "synthetic_dgp.R"))

log_sum_exp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

predict_log_inverse_margin <- function(panel, k) {
  alpha <- panel$parameters$alpha
  s0 <- panel$parameters$s0
  unlist(lapply(seq_along(panel$truth), function(period) {
    tr <- panel$truth[[period]]
    unconditional <- tr$size * (1 - s0)
    follower <- seq_along(unconditional) > k
    phi <- 1 / (1 - sum(
      unconditional[follower]^2 /
        (1 - unconditional[follower] + unconditional[follower]^2)
    ))
    h <- ifelse(
      follower,
      1 / (1 - unconditional),
      1 / (1 - phi * unconditional)
    )
    rows <- which(panel$data$year == as.character(period))
    firm <- match(panel$data$group[rows], tr$firm_ids)
    log(alpha * tr$prices / h[firm])
  }))
}

panel <- synthetic_stackelberg_panel(
  demand = "logit", conduct = "bertrand", seed = 1103L,
  layout = "single", normalize_bayes = TRUE,
  partition = "ordered_k", k = 2L, structure = "stackelberg"
)
support <- 1:5
mu <- vapply(support, function(k) predict_log_inverse_margin(panel, k),
             numeric(nrow(panel$data)))
y <- log(1 / panel$data$margins)
df <- panel$parameters$margin_noise_df
margin_sd <- panel$parameters$margin_noise_sd
student_scale <- margin_sd * sqrt((df - 2) / df)
score <- colSums(apply(mu, 2L, function(pred) {
  stats::dt((y - pred) / student_scale, df = df, log = TRUE) -
    log(student_scale)
}))
oracle_probability <- exp(score - log_sum_exp(score))

score_table <- data.frame(
  k = support,
  noiseless_rmse_log_inverse_margin = sqrt(colMeans((mu - mu[, 2L])^2)),
  noiseless_rmse_over_margin_sd =
    sqrt(colMeans((mu - mu[, 2L])^2)) / margin_sd,
  realized_log_likelihood = score,
  realized_log_likelihood_relative_to_k2 = score - score[[2L]],
  oracle_posterior_probability = oracle_probability
)

mc_seed <- 20260915L
n_replication <- 5000L
set.seed(mc_seed)
probability <- matrix(NA_real_, n_replication, length(support))
for (b in seq_len(n_replication)) {
  simulated_y <- mu[, 2L] + stats::rt(nrow(mu), df = df) * student_scale
  simulated_score <- colSums(apply(mu, 2L, function(pred) {
    stats::dt((simulated_y - pred) / student_scale, df = df, log = TRUE) -
      log(student_scale)
  }))
  probability[b, ] <- exp(simulated_score - log_sum_exp(simulated_score))
}

power_table <- data.frame(
  replication_seed = mc_seed,
  replications = n_replication,
  true_k = 2L,
  margin_noise_sd = margin_sd,
  probability_true_k_is_modal = mean(max.col(probability) == 2L),
  mean_p_true_k = mean(probability[, 2L]),
  median_p_true_k = stats::median(probability[, 2L]),
  mean_posterior_k = mean(probability %*% support)
)

zero_noise <- do.call(rbind, lapply(c(.04, .02), function(sd) {
  scale <- sd * sqrt((df - 2) / df)
  zero_score <- colSums(apply(mu, 2L, function(pred) {
    stats::dt((mu[, 2L] - pred) / scale, df = df, log = TRUE) - log(scale)
  }))
  data.frame(margin_noise_sd = sd, k = support,
             posterior_probability = exp(zero_score - log_sum_exp(zero_score)))
}))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(score_table, file.path(out_dir, "bertrand_k_oracle_scores.csv"),
                 row.names = FALSE)
utils::write.csv(power_table, file.path(out_dir, "bertrand_k_oracle_power.csv"),
                 row.names = FALSE)
utils::write.csv(zero_noise, file.path(out_dir, "bertrand_k_zero_noise.csv"),
                 row.names = FALSE)
cat("PASS: exact-k conditional identification audit; outputs in ", out_dir, "\n", sep = "")
