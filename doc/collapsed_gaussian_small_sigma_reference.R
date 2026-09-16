#!/usr/bin/env Rscript

## Preserve the failed subtractive algorithm and test its stable replacement.
## A centered, one-product-per-firm residual has an exact spectral density:
## the centered subspace has variance sigma^2 + tau^2 and its orthogonal
## constant direction has variance sigma^2. No ill-conditioned dense oracle
## is required for this particular stress fixture.
args <- commandArgs(trailingOnly = FALSE)
script <- args[startsWith(args, "--file=")]
doc_dir <- if (length(script)) dirname(normalizePath(sub("^--file=", "", script[1]))) else "doc"
reference <- new.env(parent = globalenv())
sys.source(file.path(doc_dir, "collapsed_gaussian_structured_reference.R"), reference)

initial_subtractive_density <- function(z, v, group, x, tau, G) {
  inv_v <- 1 / v
  a <- reference$group_sum(x * x * inv_v, group, G)
  c_stat <- reference$group_sum(x * z * inv_v, group, G)
  e <- sum(z * z * inv_v)
  w <- 1 / (1 + tau * tau * a)
  quadratic <- e - tau * tau * sum(c_stat * c_stat * w)
  quadratic <- quadratic + tau * tau * sum(c_stat * w)^2 / sum(w)
  logdet <- sum(log(v)) + sum(log1p(tau * tau * a)) + log(mean(w))
  -0.5 * (length(z) * log(2 * pi) + logdet + quadratic)
}

N <- 24L
tau <- 0.5
z <- seq(-1, 1, length.out = N)
z <- z - mean(z)
rows <- list()
for (sigma in c(1e-2, 1e-3, 1e-4, 1e-6, 1e-8)) {
  v <- rep(sigma^2, N)
  eigen_variance <- sigma^2 + tau^2
  oracle <- -0.5 * (N * log(2 * pi) + log(sigma^2) +
    (N - 1) * log(eigen_variance) + sum(z^2) / eigen_variance)
  values <- c(initial_subtraction = initial_subtractive_density(z, v, seq_len(N), rep(1, N), tau, N),
    stable_within_residual = reference$structured_gaussian_log_density(z, v, seq_len(N), rep(1, N), tau, N)$log_density)
  for (method in names(values)) {
    error <- abs(values[[method]] - oracle)
    rows[[length(rows) + 1L]] <- data.frame(method, N, sigma_share = sigma,
      error, tolerance = 1e-8, status = if (error <= 1e-8) "passed" else "failed")
  }
}
result <- do.call(rbind, rows)
write.csv(result, file.path(doc_dir, "collapsed_gaussian_small_sigma_validation.csv"), row.names = FALSE)
print(result, row.names = FALSE)
stopifnot(all(result$status[result$method == "stable_within_residual"] == "passed"))
