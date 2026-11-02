#!/usr/bin/env Rscript

## Independent target and chain-rule oracle for the exact contrast-total-scale
## / share-variance-fraction coordinates.  This file does not load package
## code and keeps every prior normalizer.

options(stringsAsFactors = FALSE)

log_half_normal_shifted <- function(tau, mean, sd) {
  dnorm(tau, mean, sd, log = TRUE) - pnorm(0, mean, sd,
    lower.tail = FALSE, log.p = TRUE)
}

legacy_economic_target <- function(ell, tau, state) {
  ## A small normalized Gaussian observation block supplies nontrivial scale
  ## curvature without duplicating the production structured likelihood.
  covariance <- tau^2 * state$K + exp(2 * ell) * diag(nrow(state$K))
  R <- chol(covariance)
  z <- forwardsolve(t(R), state$residual)
  likelihood <- -0.5 * (
    length(z) * log(2 * pi) + 2 * sum(log(diag(R))) + sum(z^2)
  )
  prior <- dnorm(ell, log(state$prior_sigma_share),
                 state$prior_log_sigma_share_sd, log = TRUE) +
    log_half_normal_shifted(tau, 1, state$prior_tau_sd)
  likelihood + prior
}

legacy_unconstrained_target <- function(aell, state) {
  ell <- aell[[1L]]
  tau <- exp(aell[[2L]])
  legacy_economic_target(ell, tau, state) + log(tau)
}

total_fraction_to_economic <- function(xy) {
  x <- xy[[1L]]
  y <- xy[[2L]]
  lambda <- plogis(y)
  list(
    ell = x + 0.5 * plogis(y, log.p = TRUE),
    tau = exp(x + 0.5 * plogis(-y, log.p = TRUE)),
    lambda = lambda
  )
}

total_fraction_target <- function(xy, state) {
  economic <- total_fraction_to_economic(xy)
  legacy_economic_target(economic$ell, economic$tau, state) +
    log(economic$tau) - log(2)
}

finite_gradient <- function(fun, x, step = 1e-6) {
  vapply(seq_along(x), function(j) {
    h <- step * max(1, abs(x[[j]]))
    xp <- xm <- x
    xp[[j]] <- xp[[j]] + h
    xm[[j]] <- xm[[j]] - h
    (fun(xp) - fun(xm)) / (2 * h)
  }, numeric(1))
}

## Avoid a package dependency for a two-term log-sum-exp.
log_add_exp <- function(a, b) {
  m <- max(a, b)
  m + log(exp(a - m) + exp(b - m))
}
economic_to_total_fraction <- function(ell, tau) {
  log_tau <- log(tau)
  list(x = 0.5 * log_add_exp(2 * ell, 2 * log_tau),
       y = 2 * (ell - log_tau))
}

G <- 5L
P <- diag(G) - matrix(1 / G, G, G)
state <- list(
  K = P,
  residual = c(-0.31, 0.17, 0.42, -0.19, -0.09),
  prior_sigma_share = 0.5,
  prior_log_sigma_share_sd = 0.5,
  prior_tau_sd = 0.75
)

scenarios <- list(
  ordinary = c(sigma_share = 0.45, tau = 0.8),
  share_small = c(sigma_share = 0.01, tau = 1.3),
  quality_small = c(sigma_share = 0.6, tau = 0.01),
  both_small = c(sigma_share = 0.012, tau = 0.008),
  share_large = c(sigma_share = 2.2, tau = 0.3)
)

target_rows <- list()
gradient_rows <- list()
for (name in names(scenarios)) {
  sigma <- scenarios[[name]][["sigma_share"]]
  tau <- scenarios[[name]][["tau"]]
  ell <- log(sigma)
  mapped <- economic_to_total_fraction(ell, tau)
  old_u <- c(ell, log(tau))
  new_u <- c(mapped$x, mapped$y)
  old <- legacy_unconstrained_target(old_u, state)
  new <- total_fraction_target(new_u, state)
  reconstructed <- total_fraction_to_economic(new_u)
  target_rows[[length(target_rows) + 1L]] <- data.frame(
    scenario = name,
    sigma_share = sigma,
    sigma_quality = tau,
    reconstruction_error = max(abs(c(
      reconstructed$ell - ell, reconstructed$tau - tau
    ))),
    target_error_after_coordinate_constant = abs(new - old + log(2)),
    tolerance = 2e-10,
    passed = max(abs(c(reconstructed$ell - ell,
                       reconstructed$tau - tau)),
                 abs(new - old + log(2))) <= 2e-10
  )

  old_grad <- finite_gradient(
    function(x) legacy_unconstrained_target(x, state), old_u
  )
  new_grad <- finite_gradient(
    function(x) total_fraction_target(x, state), new_u
  )
  lambda <- reconstructed$lambda
  ## Chain rule from legacy unconstrained coordinates (ell, log tau) to
  ## (x, y).  The -log(2) coordinate constant has zero gradient.
  expected <- c(
    old_grad[[1L]] + old_grad[[2L]],
    0.5 * (1 - lambda) * old_grad[[1L]] -
      0.5 * lambda * old_grad[[2L]]
  )
  for (j in 1:2) {
    gradient_rows[[length(gradient_rows) + 1L]] <- data.frame(
      scenario = name,
      coordinate = c("log_contrast_total_scale",
                     "logit_share_variance_fraction")[[j]],
      finite_gradient = new_grad[[j]],
      chain_rule_gradient = expected[[j]],
      relative_error = abs(new_grad[[j]] - expected[[j]]) /
        (1 + abs(expected[[j]])),
      tolerance = 2e-6,
      passed = abs(new_grad[[j]] - expected[[j]]) /
        (1 + abs(expected[[j]])) <= 2e-6
    )
  }
}

targets <- do.call(rbind, target_rows)
gradients <- do.call(rbind, gradient_rows)
args <- commandArgs(trailingOnly = FALSE)
script <- args[startsWith(args, "--file=")]
doc_dir <- if (length(script)) {
  dirname(normalizePath(sub("^--file=", "", script[[1L]])))
} else "doc"
write.csv(targets, file.path(doc_dir, "quality_scale_rotation_target_checks.csv"),
          row.names = FALSE)
write.csv(gradients,
          file.path(doc_dir, "quality_scale_rotation_gradient_checks.csv"),
          row.names = FALSE)
cat(sprintf("targets %d/%d; gradients %d/%d\n",
            sum(targets$passed), nrow(targets),
            sum(gradients$passed), nrow(gradients)))
if (any(!targets$passed) || any(!gradients$passed)) {
  stop("quality-scale rotation reference checks failed")
}
