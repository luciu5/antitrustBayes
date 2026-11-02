#!/usr/bin/env Rscript
# Fixed-parameter covariance/information identities, not HMC performance claims.
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[[1L]] else file.path("doc", "statistical_toolkit_runs",
  "20260916", "performance-audit", "recovery", "market-quality", "information-reference")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
checks <- list()
record <- function(label, actual, expected) {
  error <- max(abs(actual - expected))
  checks[[length(checks) + 1L]] <<- data.frame(check = label,
    maximum_absolute_error = error, tolerance = 2e-10, passed = error <= 2e-10)
}
for (G in c(1L, 2L, 5L)) for (T in c(1L, 2L, 7L)) {
  group <- rep(seq_len(G), each = T)
  Z <- diag(G)[group, , drop = FALSE]
  P <- diag(G) - matrix(1 / G, G, G)
  sigma <- .37; tau <- .83
  V <- sigma^2 * diag(G * T) + tau^2 * Z %*% P %*% t(Z)
  expected <- sort(c(rep(sigma^2 + T * tau^2, G - 1L),
    rep(sigma^2, G * T - G + 1L)))
  record(paste("balanced eigenvalues", G, T),
    sort(eigen(V, symmetric = TRUE, only.values = TRUE)$values), expected)
}
# Independent direct residual-covariance inverse determines quality curvature.
# Include a single-product group, repeated observations, missing margins and
# positive/negative correlations; dense matrices are small reference-only.
for (rho in c(-.7, 0, .31)) {
  group <- c(1L, 1L, 2L, 3L, 3L, 3L)
  matched <- c(1L, 3L, 4L, 6L)
  N <- length(group); J <- length(matched); G <- max(group)
  Z <- diag(G)[group, , drop = FALSE]
  R <- diag(N)[, matched, drop = FALSE]
  sigma_s <- .37; sigma_m <- .61
  V <- rbind(cbind(sigma_s^2 * diag(N), rho * sigma_s * sigma_m * R),
    cbind(rho * sigma_s * sigma_m * t(R), sigma_m^2 * diag(J)))
  X <- rbind(Z, matrix(0, J, G))
  curvature <- crossprod(X, solve(V, X))
  observed <- tabulate(group[matched], G)
  missing <- tabulate(group, G) - observed
  a <- missing / sigma_s^2 + observed / (sigma_s^2 * (1 - rho^2))
  record(paste("joint data curvature", rho), curvature, diag(a))
}
result <- do.call(rbind, checks)
write.csv(result, file.path(out, "validation.csv"), row.names = FALSE)
print(result)
if (!all(result$passed)) stop("Information reference identity failed.")
