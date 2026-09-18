#!/usr/bin/env Rscript
# Small-N diagnostic of the existing Logit normalization, not a model change.
# Dense covariance is oracle-only here (39 PNB observations).
options(stringsAsFactors = FALSE)
repo <- normalizePath(Sys.getenv("MQ_REPO", "/mnt/d/projects/antitrustbayes"))
base <- file.path(repo, "doc/statistical_toolkit_runs/20260916/performance-audit/recovery/market-quality")
data <- readRDS(file.path(base, "pnb-archived/frozen-data.rds"))
out <- file.path(base, "mean-channel-reference")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
G <- nrow(data)
stopifnot(length(unique(data$group)) == G)
y <- log(data$quantities/sum(data$quantities))
x <- (data$prices-mean(data$prices))/sd(data$prices)
C <- contr.helmert(G)
C <- sweep(C, 2L, sqrt(colSums(C^2)), "/")
quality <- as.numeric(C %*% seq(-.13, .17, length.out = G-1L))
a <- rep(1/G, G)
eta <- qlogis(.21)
sigma_s <- .35; sigma_m <- .22
margin_idx <- which(seq_len(G) %% 7L != 0L)
R <- diag(G)[, margin_idx, drop = FALSE]
rows <- list()
check <- function(name, error, tolerance = 2e-12) {
  rows[[length(rows)+1L]] <<- data.frame(check = name, error = error,
    tolerance = tolerance, passed = is.finite(error) && error <= tolerance)
}
check("centered_price_mean", abs(mean(x)))
check("centered_quality_mean", abs(mean(quality)))
for (sign in c(-1, 1)) for (alpha in c(.4, 1.7, 6))
  check(paste("share_mean_projection", sign, alpha, sep = "_"),
    abs(mean(eta + sign*alpha*x + quality)-eta))
for (tau in c(0, .43, 3)) for (rho in c(0, .31)) {
  S <- tau^2 * tcrossprod(C) + diag(sigma_s^2, G)
  B <- rho*sigma_s*sigma_m*R
  V <- rbind(cbind(S, B), cbind(t(B), diag(sigma_m^2, length(margin_idx))))
  stopifnot(all(is.finite(chol(V))))
  projection <- c(a, rep(0, length(margin_idx)))
  variance <- as.numeric(t(projection) %*% V %*% projection)
  check(paste("share_average_noise_variance", tau, rho, sep = "_"),
    abs(variance - sigma_s^2/G))
}
check("Jensen_bound", max(0, mean(y)+log(G)))
validation <- do.call(rbind, rows)
write.csv(validation, file.path(out, "validation.csv"), row.names = FALSE)
summary <- data.frame(groups = G, mean_log_inside_share = mean(y),
  geometric_mean_inside_share = exp(mean(y)),
  small_share_noise_outside_point = plogis(mean(y)),
  uniform_share_upper_bound = 1/(G+1),
  interpretation = "Normalization diagnostic; not an estimated outside share or an explanation of every Cournot anomaly")
write.csv(summary, file.path(out, "summary.csv"), row.names = FALSE)
print(summary)
print(validation)
stopifnot(all(validation$passed))
