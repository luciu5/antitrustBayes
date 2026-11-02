## Independent base-R reference for the collapsed Gaussian quality block.
##
## This file deliberately does not source package code. It implements:
##   (i) the current orthonormal Helmert quality coordinates;
##   (ii) fully normalized expanded and collapsed Gaussian densities;
##   (iii) direct low-dimensional integration;
##   (iv) Gaussian quality posterior moments;
##   (v) the collapsed jointly Gaussian rho block with missing margins.
##
## It is a mathematical/reference check, not a production likelihood.

helmert_basis <- function(n_group) {
  out <- matrix(0, nrow = n_group, ncol = max(0L, n_group - 1L))
  if (n_group > 1L) {
    for (k in seq_len(n_group - 1L)) {
      out[seq_len(k), k] <- 1 / sqrt(k * (k + 1))
      out[k + 1L, k] <- -k / sqrt(k * (k + 1))
    }
  }
  out
}

one_hot <- function(group) {
  lev <- sort(unique(group))
  z <- matrix(0, nrow = length(group), ncol = length(lev))
  z[cbind(seq_along(group), match(group, lev))] <- 1
  z
}

solve_spd <- function(a, b) {
  if (!length(a)) return(matrix(numeric(), nrow = 0L, ncol = ncol(as.matrix(b))))
  r <- chol(a)
  b <- as.matrix(b)
  z <- forwardsolve(t(r), b)
  backsolve(r, z)
}

logdet_spd <- function(a) {
  if (!length(a)) return(0)
  2 * sum(log(diag(chol(a))))
}

log_mvnorm <- function(x, mean, covariance) {
  x <- as.numeric(x)
  mean <- as.numeric(mean)
  if (!length(x)) return(0)
  covariance <- as.matrix(covariance)
  r <- x - mean
  chol_a <- chol(covariance)
  z <- forwardsolve(t(chol_a), r)
  -0.5 * (length(x) * log(2 * pi) + 2 * sum(log(diag(chol_a))) +
            sum(z * z))
}

make_quality_system <- function(group) {
  z <- one_hot(group)
  c_basis <- helmert_basis(ncol(z))
  list(Z = z, C = c_basis, H = z %*% c_basis)
}

collapsed_share_covariance <- function(Z, C, tau, sigma) {
  h <- Z %*% C
  sigma^2 * diag(nrow(Z)) + tau^2 * tcrossprod(h)
}

expanded_log_density <- function(y, mu, beta, H, tau, sigma) {
  log_mvnorm(y, mu + as.numeric(H %*% beta),
             sigma^2 * diag(length(y))) +
    log_mvnorm(beta, rep(0, ncol(H)), tau^2 * diag(ncol(H)))
}

collapsed_share_log_density <- function(y, mu, Z, C, tau, sigma) {
  covariance <- collapsed_share_covariance(Z, C, tau, sigma)
  log_mvnorm(y, mu, covariance)
}

quality_posterior <- function(y, mu, H, tau, sigma,
                              residual_covariance = NULL) {
  d <- ncol(H)
  if (is.null(residual_covariance)) {
    residual_covariance <- sigma^2 * diag(length(y))
  }
  if (!d) return(list(mean = numeric(), covariance = matrix(numeric(), 0L, 0L)))
  prior_precision <- diag(1 / tau^2, d)
  k_inv_h <- solve_spd(residual_covariance, H)
  k_inv_r <- solve_spd(residual_covariance, y - mu)
  posterior_precision <- prior_precision + crossprod(H, k_inv_h)
  v <- solve_spd(posterior_precision, diag(d))
  m <- as.numeric(v %*% crossprod(H, k_inv_r))
  list(mean = m, covariance = v)
}

collapsed_joint_covariance <- function(Z, C, tau, sigma_share, sigma_margin,
                                       rho, R) {
  n <- nrow(Z)
  q <- ncol(R)
  share_cov <- collapsed_share_covariance(Z, C, tau, sigma_share)
  if (!q) return(share_cov)
  b <- rho * sigma_share * sigma_margin * R
  rbind(cbind(share_cov, b),
        cbind(t(b), sigma_margin^2 * diag(q)))
}

collapsed_joint_log_density <- function(y_share, mu, y_margin, nu,
                                         Z, C, tau, sigma_share,
                                         sigma_margin, rho, R) {
  y <- c(y_share, y_margin)
  mean <- c(mu, nu)
  covariance <- collapsed_joint_covariance(
    Z, C, tau, sigma_share, sigma_margin, rho, R
  )
  log_mvnorm(y, mean, covariance)
}

collapsed_joint_conditional_log_density <- function(y_share, mu, y_margin,
                                                     nu, Z, C, tau,
                                                     sigma_share,
                                                     sigma_margin, rho, R) {
  share_cov <- collapsed_share_covariance(Z, C, tau, sigma_share)
  share_resid <- y_share - mu
  if (!length(y_margin)) return(log_mvnorm(y_share, mu, share_cov))
  b <- rho * sigma_share * sigma_margin * R
  cond_mean <- nu + as.numeric(t(b) %*% solve_spd(share_cov, share_resid))
  cond_cov <- sigma_margin^2 * diag(length(y_margin)) -
    t(b) %*% solve_spd(share_cov, b)
  log_mvnorm(y_share, mu, share_cov) +
    log_mvnorm(y_margin, cond_mean, cond_cov)
}

joint_quality_posterior <- function(y_share, mu, y_margin, nu, H,
                                    tau, sigma_share, sigma_margin,
                                    rho, R) {
  n <- length(y_share)
  q <- length(y_margin)
  d <- ncol(H)
  if (!d) return(list(mean = numeric(), covariance = matrix(numeric(), 0L, 0L)))
  k <- if (!q) {
    sigma_share^2 * diag(n)
  } else {
    b <- rho * sigma_share * sigma_margin * R
    rbind(cbind(sigma_share^2 * diag(n), b),
          cbind(t(b), sigma_margin^2 * diag(q)))
  }
  x <- rbind(H, matrix(0, nrow = q, ncol = d))
  residual <- c(y_share - mu, y_margin - nu)
  prior_precision <- diag(1 / tau^2, d)
  k_inv_x <- solve_spd(k, x)
  k_inv_residual <- solve_spd(k, residual)
  v <- solve_spd(prior_precision + crossprod(x, k_inv_x), diag(d))
  m <- as.numeric(v %*% crossprod(x, k_inv_residual))
  list(mean = m, covariance = v)
}

direct_integral_1d <- function(log_integrand) {
  anchor <- log_integrand(0)
  scaled <- function(b) vapply(
    b, function(one_b) exp(log_integrand(one_b) - anchor), numeric(1)
  )
  z <- integrate(scaled, -Inf, Inf, rel.tol = 1e-10,
                 subdivisions = 1000L)$value
  first <- integrate(function(b) b * scaled(b), -Inf, Inf,
                     rel.tol = 1e-10, subdivisions = 1000L)$value
  second <- integrate(function(b) b * b * scaled(b), -Inf, Inf,
                      rel.tol = 1e-10, subdivisions = 1000L)$value
  list(log_density = anchor + log(z), mean = first / z,
       variance = second / z - (first / z)^2)
}

direct_integral_2d <- function(log_integrand) {
  anchor <- log_integrand(c(0, 0))
  scaled <- function(b1, b2) vapply(
    b1, function(one_b1) exp(log_integrand(c(one_b1, b2)) - anchor),
    numeric(1)
  )
  inner_one <- function(b2) {
    integrate(function(b1) scaled(b1, b2), -Inf, Inf,
              rel.tol = 2e-8, subdivisions = 300L)$value
  }
  inner <- function(b2) vapply(b2, inner_one, numeric(1))
  z <- integrate(inner, -Inf, Inf, rel.tol = 2e-8,
                 subdivisions = 300L)$value
  anchor + log(z)
}

run_checks <- function(verbose = TRUE) {
  checks <- list()
  add_check <- function(name, value) {
    checks[[name]] <<- as.numeric(value)
  }

  ## Helmert identities, including the zero-dimensional G = 1 boundary.
  for (g in c(1L, 2L, 3L, 5L)) {
    c_basis <- helmert_basis(g)
    p <- diag(g) - matrix(1 / g, nrow = g, ncol = g)
    if (g > 1L) {
      add_check(paste0("helmert_orthonormal_G", g),
                max(abs(crossprod(c_basis) - diag(g - 1L))))
      add_check(paste0("helmert_projection_G", g),
                max(abs(tcrossprod(c_basis) - p)))
      add_check(paste0("helmert_center_G", g),
                max(abs(colSums(c_basis))))
    } else {
      add_check("helmert_zero_dimension_G1", length(c_basis))
    }
  }

  ## A multiproduct firm: rows one and two share a group-quality effect.
  group <- c(1L, 1L, 2L, 2L)
  system <- make_quality_system(group)
  z <- system$Z
  c_basis <- system$C
  h <- system$H
  stopifnot(nrow(z) == 4L, ncol(z) == 2L, ncol(h) == 1L,
            identical(z[1, ], z[2, ]), max(abs(z - diag(4)[, seq_len(2)])) > 0)
  tau <- 0.37
  sigma <- 0.23
  mu <- c(-0.25, 0.18, -0.08, 0.35)
  y <- c(-0.47, 0.34, -0.12, 0.27)
  expanded_integrand <- function(beta) {
    expanded_log_density(y, mu, beta, h, tau, sigma)
  }
  collapsed <- collapsed_share_log_density(y, mu, z, c_basis, tau, sigma)
  integrated <- direct_integral_1d(expanded_integrand)
  add_check("expanded_collapsed_G2_multiproduct_log_density",
            abs(collapsed - integrated$log_density))
  add_check("expanded_collapsed_G2_multiproduct_density",
            abs(exp(collapsed) - exp(integrated$log_density)))

  beta_test <- 0.41
  group_quality <- as.numeric(c_basis %*% beta_test)
  add_check("group_quality_zero_sum", abs(sum(group_quality)))
  add_check("multiproduct_quality_identity",
            max(abs(as.numeric(h %*% beta_test) -
                      group_quality[match(group, seq_along(group_quality))])))

  posterior <- quality_posterior(y, mu, h, tau, sigma)
  add_check("posterior_mean_G2",
            abs(posterior$mean - integrated$mean))
  add_check("posterior_variance_G2",
            abs(posterior$covariance - integrated$variance))
  add_check("posterior_group_mean_G2",
            max(abs(as.numeric(c_basis %*% posterior$mean) -
                      group_quality * (posterior$mean / beta_test))))
  add_check("posterior_group_covariance_centered_G2",
            max(abs(c_basis %*% posterior$covariance %*% t(c_basis) -
                      posterior$covariance[1, 1] * tcrossprod(c_basis))))

  ## A three-group two-dimensional integral gives a separate direct check.
  group3 <- c(1L, 1L, 2L, 3L, 3L)
  system3 <- make_quality_system(group3)
  mu3 <- c(-0.3, 0.1, -0.2, 0.4, 0.05)
  y3 <- c(-0.5, 0.32, -0.15, 0.31, 0.1)
  tau3 <- 0.29
  sigma3 <- 0.19
  expanded3 <- function(beta) {
    expanded_log_density(y3, mu3, beta, system3$H, tau3, sigma3)
  }
  collapsed3 <- collapsed_share_log_density(
    y3, mu3, system3$Z, system3$C, tau3, sigma3
  )
  integrated3 <- direct_integral_2d(expanded3)
  add_check("expanded_collapsed_G3_direct_2d_log_density",
            abs(collapsed3 - integrated3))

  ## Low-rank covariance form equals the transparent covariance form.
  u <- tau3 * system3$H
  low_rank <- sigma3^2 * diag(length(y3)) + tcrossprod(u)
  direct_cov <- collapsed_share_covariance(
    system3$Z, system3$C, tau3, sigma3
  )
  add_check("low_rank_covariance_identity", max(abs(low_rank - direct_cov)))

  ## Optional rho block with q = 3 observed margins out of n = 5 shares.
  n <- 5L
  q <- 3L
  matched_rows <- c(2L, 4L, 5L)
  r_match <- matrix(0, nrow = n, ncol = q)
  r_match[cbind(matched_rows, seq_len(q))] <- 1
  sigma_margin <- 0.31
  rho <- 0.43
  y_margin <- c(0.55, 0.62, 0.48)
  nu <- c(0.5, 0.6, 0.45)
  joint_cov <- collapsed_joint_covariance(
    system3$Z, system3$C, tau3, sigma3, sigma_margin, rho, r_match
  )
  eigen_min <- min(eigen(joint_cov, symmetric = TRUE,
                         only.values = TRUE)$values)
  add_check("rho_missing_margin_min_eigenvalue", eigen_min)
  joint_direct <- collapsed_joint_log_density(
    y3, mu3, y_margin, nu, system3$Z, system3$C, tau3, sigma3,
    sigma_margin, rho, r_match
  )
  joint_conditional <- collapsed_joint_conditional_log_density(
    y3, mu3, y_margin, nu, system3$Z, system3$C, tau3, sigma3,
    sigma_margin, rho, r_match
  )
  add_check("rho_block_conditional_factorization",
            abs(joint_direct - joint_conditional))

  ## The rho posterior formula is independently evaluated from the expanded
  ## residual covariance, and reduces to the share-only formula at rho = 0.
  joint_post <- joint_quality_posterior(
    y3, mu3, y_margin, nu, system3$H, tau3, sigma3,
    sigma_margin, rho, r_match
  )
  zero_post <- joint_quality_posterior(
    y3, mu3, y_margin, nu, system3$H, tau3, sigma3,
    sigma_margin, 0, r_match
  )
  share_post <- quality_posterior(y3, mu3, system3$H, tau3, sigma3)
  add_check("rho_zero_posterior_mean_reduction",
            max(abs(zero_post$mean - share_post$mean)))
  add_check("rho_zero_posterior_covariance_reduction",
            max(abs(zero_post$covariance - share_post$covariance)))
  add_check("rho_joint_posterior_finite",
            max(abs(c(joint_post$mean, joint_post$covariance))))

  ## q = 0 is the missing-margin boundary, and the collapsed density remains
  ## exactly the share density. The share term is therefore common to all k.
  r_empty <- matrix(numeric(), nrow = nrow(system3$Z), ncol = 0L)
  no_margin_joint <- collapsed_joint_log_density(
    y3, mu3, numeric(), numeric(), system3$Z, system3$C, tau3,
    sigma3, sigma_margin, rho, r_empty
  )
  add_check("zero_margin_block_reduction", abs(no_margin_joint - collapsed3))
  share_by_k <- vapply(0:5, function(k) {
    collapsed_share_log_density(y3, mu3, system3$Z, system3$C, tau3, sigma3)
  }, numeric(1))
  add_check("share_density_common_across_k", max(abs(share_by_k - share_by_k[1])))

  ## Exact-k likelihood algebra: for any continuous parameter draw, one
  ## market factorizes as common share density plus log-sum-exp of complete
  ## margin-state densities. Independent market states multiply after their
  ## own marginalization, rather than enumerating K^M states.
  margin_scores_m1 <- c(-2.1, -1.4, -3.0)
  margin_scores_m2 <- c(-0.8, -2.5)
  prior_m1 <- c(0.2, 0.3, 0.5)
  prior_m2 <- c(0.6, 0.4)
  lse <- function(x) {
    max_x <- max(x)
    max_x + log(sum(exp(x - max_x)))
  }
  state_m1 <- lse(log(prior_m1) + margin_scores_m1)
  state_m2 <- lse(log(prior_m2) + margin_scores_m2)
  ## Each market has its own share term. Enumerating the Cartesian product of
  ## state indices gives the same value because the state prior factorizes.
  separate_with_second_share <- 2 * collapsed3 + state_m1 + state_m2
  enumerated_with_second_share <- 2 * collapsed3 + lse(
    rep(log(prior_m1), each = length(prior_m2)) +
      rep(margin_scores_m1, each = length(prior_m2)) +
      rep(log(prior_m2), times = length(prior_m1)) +
      rep(margin_scores_m2, times = length(prior_m1))
  )
  add_check("exact_k_market_factorization",
            abs(separate_with_second_share - enumerated_with_second_share))

  ## Bounds for assertions. The eigenvalue is a positive quantity rather than
  ## an error, so assert it separately.
  stopifnot(
    checks$helmert_orthonormal_G2 < 1e-14,
    checks$helmert_projection_G2 < 1e-14,
    checks$helmert_center_G2 < 1e-14,
    checks$expanded_collapsed_G2_multiproduct_log_density < 1e-8,
    checks$expanded_collapsed_G3_direct_2d_log_density < 1e-6,
    checks$low_rank_covariance_identity < 1e-14,
    checks$group_quality_zero_sum < 1e-14,
    checks$multiproduct_quality_identity < 1e-14,
    checks$posterior_mean_G2 < 1e-8,
    checks$posterior_variance_G2 < 1e-8,
    checks$rho_missing_margin_min_eigenvalue > 0,
    checks$rho_block_conditional_factorization < 1e-8,
    checks$rho_zero_posterior_mean_reduction < 1e-10,
    checks$rho_zero_posterior_covariance_reduction < 1e-10,
    is.finite(checks$rho_joint_posterior_finite),
    checks$zero_margin_block_reduction < 1e-10,
    checks$share_density_common_across_k < 1e-14,
    checks$exact_k_market_factorization < 1e-10
  )

  out <- data.frame(
    check = names(checks),
    value = unname(unlist(checks)),
    row.names = NULL
  )
  if (verbose) {
    print(out, row.names = FALSE)
    cat("\nAll independent collapsed-Gaussian reference checks passed.\n")
    cat("The expanded and collapsed densities agree including Normalizing constants.\n")
  }
  invisible(out)
}

if (identical(environment(), globalenv())) run_checks()
