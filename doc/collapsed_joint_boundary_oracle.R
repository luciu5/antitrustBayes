#!/usr/bin/env Rscript

## Independent boundary audit for the production collapsed joint-Gaussian
## likelihood.  The dense block covariance is the small-N oracle.  The
## structured evaluator is sourced from the package-independent O(N + G)
## reference used by the original collapsed-Gaussian validation.

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = FALSE)
script <- args[startsWith(args, "--file=")]
doc_dir <- if (length(script)) {
  dirname(normalizePath(sub("^--file=", "", script[[1L]])))
} else {
  "doc"
}
reference <- new.env(parent = globalenv())
sys.source(file.path(doc_dir, "collapsed_gaussian_structured_reference.R"),
           reference)

helmert_basis <- function(G) {
  if (G == 1L) return(matrix(numeric(), nrow = 1L, ncol = 0L))
  C <- matrix(0, G, G - 1L)
  for (j in seq_len(G - 1L)) {
    C[seq_len(j), j] <- 1 / sqrt(j * (j + 1))
    C[j + 1L, j] <- -j / sqrt(j * (j + 1))
  }
  C
}

make_fixture <- function(kind = c("sparse", "panel")) {
  kind <- match.arg(kind)
  if (kind == "sparse") {
    group <- seq_len(8L)
    matched <- c(1L, 3L, 4L, 7L, 8L)
  } else {
    ## Four persistent groups, two products per group, repeated for three
    ## periods.  This exercises repeated quality loadings and missing margins.
    group <- rep(rep(seq_len(4L), each = 2L), times = 3L)
    matched <- setdiff(seq_along(group), c(2L, 5L, 9L, 14L, 20L))
  }
  list(kind = kind, group = group, matched = matched, G = max(group))
}

make_residuals <- function(fixture, tau, sigma_share, sigma_margin, rho,
                           seed) {
  set.seed(seed)
  C <- helmert_basis(fixture$G)
  beta <- if (ncol(C)) rnorm(ncol(C), sd = tau) else numeric()
  quality <- if (ncol(C)) as.numeric(C %*% beta) else 0
  e_share <- rnorm(length(fixture$group))
  share <- quality[fixture$group] + sigma_share * e_share
  e_margin <- rnorm(length(fixture$matched))
  margin <- sigma_margin * (
    rho * e_share[fixture$matched] + sqrt(1 - rho^2) * e_margin
  )
  list(share = share, margin = margin)
}

evaluate_pair <- function(fixture, tau, sigma_share, sigma_margin, rho,
                          seed) {
  residual <- make_residuals(
    fixture, tau, sigma_share, sigma_margin, rho, seed
  )
  structured <- reference$rho_structured_log_density(
    residual$share, fixture$group, residual$margin, fixture$matched,
    tau, sigma_share, sigma_margin, rho, G = fixture$G
  )
  dense <- reference$rho_dense_original_log_density(
    residual$share, fixture$group, residual$margin, fixture$matched,
    tau, sigma_share, sigma_margin, rho, G = fixture$G
  )
  c(
    structured = structured$log_density,
    dense = dense$log_density,
    abs_error = abs(structured$log_density - dense$log_density),
    logdet_error = abs(structured$logdet - dense$logdet),
    quadratic_error = abs(structured$quadratic - dense$quadratic),
    slope = structured$augmented$slope
  )
}

finite_gradient <- function(fun, x, step = 2e-6) {
  vapply(seq_along(x), function(j) {
    h <- step * max(1, abs(x[[j]]))
    xp <- xm <- x
    xp[[j]] <- xp[[j]] + h
    xm[[j]] <- xm[[j]] - h
    (fun(xp) - fun(xm)) / (2 * h)
  }, numeric(1))
}

rows <- list()
gradient_rows <- list()
fixtures <- lapply(c("sparse", "panel"), make_fixture)
sigma_values <- c(0.5, 0.1, 0.03, 0.01, 0.003, 0.001)
rho_values <- c(0, 0.8, -0.8, 0.94, -0.94)
tau_values <- c(0.2, 1.5)
case <- 0L

for (fixture in fixtures) {
  for (tau in tau_values) for (sigma_share in sigma_values) {
    for (rho in rho_values) {
      case <- case + 1L
      seed <- 91000L + case
      ans <- evaluate_pair(
        fixture, tau, sigma_share, 0.35, rho, seed
      )
      ## The dense Cholesky oracle itself loses precision as its condition
      ## number approaches machine precision.  The primary posterior-relevant
      ## range uses sigma_share >= .01 and a 1e-7 gate; smaller stress states
      ## retain a transparent, scale-aware tolerance.
      tolerance <- if (sigma_share >= 0.01) 1e-7 else 2e-5
      rows[[length(rows) + 1L]] <- data.frame(
        fixture = fixture$kind,
        tau = tau,
        sigma_share = sigma_share,
        sigma_margin = 0.35,
        rho = rho,
        slope = unname(ans[["slope"]]),
        structured_log_density = unname(ans[["structured"]]),
        dense_log_density = unname(ans[["dense"]]),
        abs_error = unname(ans[["abs_error"]]),
        logdet_error = unname(ans[["logdet_error"]]),
        quadratic_error = unname(ans[["quadratic_error"]]),
        tolerance = tolerance,
        passed = is.finite(ans[["abs_error"]]) &&
          ans[["abs_error"]] <= tolerance
      )

      if (sigma_share %in% c(0.1, 0.01, 0.001) &&
          rho %in% c(-0.94, 0.94)) {
        fixture_fixed <- fixture
        seed_fixed <- seed
        ## Freeze the realized residual vector while differentiating the
        ## likelihood, as Stan does for observed data.  Coordinates are
        ## log(tau), log(sigma_share), and atanh(rho / .95).
        residual <- make_residuals(
          fixture_fixed, tau, sigma_share, 0.35, rho, seed_fixed
        )
        eval_method <- function(theta, method) {
          tau_i <- exp(theta[[1L]])
          sigma_i <- exp(theta[[2L]])
          rho_i <- 0.95 * tanh(theta[[3L]])
          if (method == "structured") {
            reference$rho_structured_log_density(
              residual$share, fixture_fixed$group, residual$margin,
              fixture_fixed$matched, tau_i, sigma_i, 0.35, rho_i,
              G = fixture_fixed$G
            )$log_density
          } else {
            reference$rho_dense_original_log_density(
              residual$share, fixture_fixed$group, residual$margin,
              fixture_fixed$matched, tau_i, sigma_i, 0.35, rho_i,
              G = fixture_fixed$G
            )$log_density
          }
        }
        theta <- c(log(tau), log(sigma_share), atanh(rho / 0.95))
        gs <- finite_gradient(function(x) eval_method(x, "structured"), theta)
        gd <- finite_gradient(function(x) eval_method(x, "dense"), theta)
        for (j in seq_along(theta)) {
          tolerance_g <- if (sigma_share >= 0.01) 2e-5 else 2e-3
          gradient_rows[[length(gradient_rows) + 1L]] <- data.frame(
            fixture = fixture$kind,
            tau = tau,
            sigma_share = sigma_share,
            rho = rho,
            coordinate = c("log_tau", "log_sigma_share", "rho_raw")[[j]],
            structured_gradient = gs[[j]],
            dense_gradient = gd[[j]],
            relative_error = abs(gs[[j]] - gd[[j]]) /
              (1 + abs(gd[[j]])),
            tolerance = tolerance_g,
            passed = is.finite(gs[[j]]) && is.finite(gd[[j]]) &&
              abs(gs[[j]] - gd[[j]]) / (1 + abs(gd[[j]])) <= tolerance_g
          )
        }
      }
    }
  }
}

results <- do.call(rbind, rows)
gradients <- do.call(rbind, gradient_rows)
write.csv(results, file.path(doc_dir, "collapsed_joint_boundary_oracle.csv"),
          row.names = FALSE)
write.csv(gradients,
          file.path(doc_dir, "collapsed_joint_boundary_gradient_oracle.csv"),
          row.names = FALSE)

cat(sprintf(
  "density: %d/%d passed, max error %.6g; gradients: %d/%d passed, max relative error %.6g\n",
  sum(results$passed), nrow(results), max(results$abs_error),
  sum(gradients$passed), nrow(gradients), max(gradients$relative_error)
))
if (any(!results$passed) || any(!gradients$passed)) {
  print(subset(results, !passed), row.names = FALSE)
  print(subset(gradients, !passed), row.names = FALSE)
  stop("collapsed joint boundary oracle failed")
}
