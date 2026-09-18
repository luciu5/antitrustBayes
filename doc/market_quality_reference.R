## Independent reference for market-block centered quality.
##
## This file intentionally does not source package code or reproduce the
## production Stan kernel.  It implements the proposed statistical block in
## base R.  A market m has local group qualities
##
##   b_m = C_m beta_m,       beta_m ~ N(0, tau^2 I_(G_m - 1)),
##
## independently across markets, with one common tau.  C_m is an orthonormal
## Helmert basis.  Thus a group label reused in two markets denotes two
## independent quality units; repeated rows within one market share their
## market-local quality.  The observed residual contract is the package
## contract: z_s = log_shareIn - share_mean and z_m =
## log_marginNormInv - pred_log_margin.
##
## The main structured likelihood conditions shares on the observed margins:
##   z_s[i] - rho*sigma_s/sigma_m*z_m[i] | beta_m
##       ~ N(b_m[group[i]], sigma_s^2*(1-rho^2)) for paired rows,
##   z_s[i] | beta_m ~ N(b_m[group[i]], sigma_s^2) otherwise,
## and adds the normalized N(0, sigma_m^2) density for observed margins.
## The dense oracle uses the original joint covariance and is only intended for
## small fixtures.  It verifies that this conditional factorization is exact.

options(stringsAsFactors = FALSE)

`%||%` <- function(x, y) if (is.null(x)) y else x

tol <- list(
  algebra = 2e-10,
  moments = 2e-9,
  quadrature = 2e-7,
  tiny_sigma = 2e-5,
  state = 2e-10
)

finite_scalar <- function(x, name, lower = -Inf, strict_lower = FALSE) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) ||
      (strict_lower && x <= lower) || (!strict_lower && x < lower)) {
    stop(name, " must be a finite scalar in the requested range")
  }
  as.numeric(x)
}

validate_integer <- function(x, name, n = NULL, allow_empty = FALSE) {
  if (!is.numeric(x) || (!allow_empty && !length(x)) || any(!is.finite(x)) ||
      any(x != floor(x)) || any(x < 1L) || (!is.null(n) && any(x > n))) {
    stop(name, " must contain positive integer labels")
  }
  as.integer(x)
}

helmert_basis_local <- function(g) {
  g <- as.integer(finite_scalar(g, "G", lower = 1))
  out <- matrix(0, nrow = g, ncol = max(0L, g - 1L))
  if (g > 1L) {
    for (k in seq_len(g - 1L)) {
      out[seq_len(k), k] <- 1 / sqrt(k * (k + 1))
      out[k + 1L, k] <- -k / sqrt(k * (k + 1))
    }
  }
  out
}

one_hot_local <- function(group, labels = sort(unique(group))) {
  group <- as.integer(group)
  z <- matrix(0, nrow = length(group), ncol = length(labels))
  if (length(group)) z[cbind(seq_along(group), match(group, labels))] <- 1
  z
}

## Convert market/group labels into independent market-local blocks.  Groups
## are sorted within each market only to give deterministic contrast columns.
make_market_layout <- function(market, group, n_market = NULL) {
  if (length(market) != length(group) || !length(market))
    stop("market and group must be nonempty vectors of equal length")
  market <- validate_integer(market, "market", n_market)
  n_market <- n_market %||% max(market)
  n_market <- as.integer(finite_scalar(n_market, "n_market", lower = 1))
  if (any(market > n_market)) stop("market index exceeds n_market")
  if (!is.numeric(group) || any(!is.finite(group)) || any(group != floor(group)) ||
      any(group < 1L)) stop("group must contain positive integer labels")
  group <- as.integer(group)
  blocks <- vector("list", n_market)
  group_local <- integer(length(group))
  group_labels <- vector("list", n_market)
  for (m in seq_len(n_market)) {
    rows <- which(market == m)
    if (!length(rows)) {
      blocks[[m]] <- integer()
      group_labels[[m]] <- integer()
      next
    }
    labels <- sort(unique(group[rows]))
    group_local[rows] <- match(group[rows], labels)
    blocks[[m]] <- rows
    group_labels[[m]] <- labels
  }
  list(market = market, group = group, group_local = group_local,
       n_market = n_market, rows = blocks,
       group_labels = group_labels,
       n_group = vapply(group_labels, length, integer(1)))
}

validate_quality_inputs <- function(z, v, group, x, tau, G = max(group)) {
  if (!length(z) || length(v) != length(z) || length(group) != length(z) ||
      length(x) != length(z)) stop("z, v, group, x must have equal nonzero lengths")
  if (any(!is.finite(z)) || any(!is.finite(v)) || any(v <= 0) ||
      any(!is.finite(x))) stop("z, x finite and v positive required")
  group <- validate_integer(group, "group", as.integer(G))
  G <- as.integer(finite_scalar(G, "G", lower = 1))
  tau <- finite_scalar(tau, "tau", lower = 0)
  invisible(list(z = as.numeric(z), v = as.numeric(v), group = group,
                x = as.numeric(x), tau = tau, G = G))
}

## Stable enough for small reference fixtures and entirely independent of the
## production one-pass update.  The group residual is calculated after the
## weighted mean, which avoids subtracting large global quadratic forms.
group_stats <- function(z, v, group, x = rep(1, length(z)), G = max(group)) {
  validate_quality_inputs(z, v, group, x, 0, G)
  G <- as.integer(G)
  w <- x * x / v
  cstat <- numeric(G)
  a <- numeric(G)
  for (i in seq_along(z)) {
    a[group[i]] <- a[group[i]] + w[i]
    cstat[group[i]] <- cstat[group[i]] + x[i] * z[i] / v[i]
  }
  means <- numeric(G)
  positive <- a > 0
  means[positive] <- cstat[positive] / a[positive]
  within <- 0
  for (i in seq_along(z)) {
    within <- within + (z[i] - x[i] * means[group[i]])^2 / v[i]
  }
  list(a = a, c = cstat, within = within,
       e = sum(z * z / v), logdet_v = sum(log(v)))
}

## Exact normalized N(X b, diag(v)) integrated under b ~ N(0,tau^2 P_G),
## where each row has a single group loading x[i].
structured_centered_log_density <- function(z, v, group, x = rep(1, length(z)),
                                             tau, G = max(group)) {
  s <- group_stats(z, v, group, x, G)
  G <- as.integer(G)
  L <- length(z)
  tau <- finite_scalar(tau, "tau", lower = 0)
  tau2 <- tau * tau
  if (G == 1L || tau == 0) {
    ## P_1 = 0 and tau=0 both leave only the observation noise.
    return(list(log_density = -0.5 * (L * log(2 * pi) + s$logdet_v + s$e),
                logdet = s$logdet_v, quadratic = s$e, stats = s,
                w = rep(1, G), sum_w = G))
  }
  w <- 1 / (1 + tau2 * s$a)
  sum_w <- sum(w)
  q0 <- s$within
  positive <- s$a > 0
  q0 <- q0 + sum((s$c[positive]^2 / s$a[positive]) * w[positive])
  quadratic <- q0 + tau2 * sum(s$c * w)^2 / sum_w
  logdet <- s$logdet_v + sum(log1p(tau2 * s$a)) + log(sum_w / G)
  list(log_density = -0.5 * (L * log(2 * pi) + logdet + quadratic),
       logdet = logdet, quadratic = quadratic, stats = s,
       w = w, sum_w = sum_w)
}

## Build the requested conditional-share rows.  `margin_idx` indexes the
## share rows, as in the current Stan data contract; each entry corresponds
## to one observed margin residual.
conditional_share_rows <- function(share_resid, market, group_local,
                                   margin_resid, margin_idx, sigma_share,
                                   sigma_margin, rho) {
  n <- length(share_resid)
  if (length(market) != n || length(group_local) != n)
    stop("share residual and layout lengths disagree")
  market <- as.integer(market); group_local <- as.integer(group_local)
  margin_idx <- as.integer(margin_idx)
  if (length(margin_resid) != length(margin_idx) ||
      anyDuplicated(margin_idx) || any(margin_idx < 1L | margin_idx > n))
    stop("margin residuals and margin_idx must be valid one-to-one rows")
  sigma_share <- finite_scalar(sigma_share, "sigma_share", lower = 0,
                               strict_lower = TRUE)
  sigma_margin <- finite_scalar(sigma_margin, "sigma_margin", lower = 0,
                                strict_lower = TRUE)
  rho <- finite_scalar(rho, "rho")
  if (abs(rho) >= 1) stop("rho must satisfy |rho| < 1")
  if (any(!is.finite(share_resid)) || any(!is.finite(margin_resid)))
    stop("residuals must be finite")
  slope <- rho * sigma_share / sigma_margin
  variance_paired <- sigma_share^2 * (1 - rho^2)
  z <- as.numeric(share_resid)
  v <- rep(sigma_share^2, n)
  paired <- logical(n)
  paired[margin_idx] <- TRUE
  z[margin_idx] <- z[margin_idx] - slope * margin_resid
  v[margin_idx] <- variance_paired
  list(z = z, v = v, paired = paired, slope = slope,
       margin_idx = margin_idx, margin_resid = as.numeric(margin_resid),
       market = market, group_local = group_local)
}

## O(N+G) market-block likelihood.  It uses local group sufficient statistics
## and never constructs a dense observation covariance matrix.
market_structured_log_likelihood <- function(share_resid, margin_resid,
                                             margin_idx, layout, tau,
                                             sigma_share, sigma_margin, rho) {
  if (!is.list(layout) || is.null(layout$market) ||
      is.null(layout$group_local) || is.null(layout$rows))
    stop("layout must be returned by make_market_layout")
  n <- length(share_resid)
  if (length(layout$market) != n) stop("layout and residual lengths disagree")
  rows <- conditional_share_rows(share_resid, layout$market,
                                 layout$group_local, margin_resid, margin_idx,
                                 sigma_share, sigma_margin, rho)
  block <- vector("list", layout$n_market)
  total <- 0
  ## The margin factor is fully normalized and has no quality loading.
  margin_log <- if (length(margin_resid)) {
    sum(-0.5 * (log(2 * pi * sigma_margin^2) +
                (margin_resid / sigma_margin)^2))
  } else 0
  for (m in seq_len(layout$n_market)) {
    ii <- layout$rows[[m]]
    if (!length(ii)) {
      block[[m]] <- list(log_density = 0, rows = ii)
      next
    }
    ans <- structured_centered_log_density(
      rows$z[ii], rows$v[ii], layout$group_local[ii], rep(1, length(ii)),
      tau = tau, G = layout$n_group[m]
    )
    ans$rows <- ii
    block[[m]] <- ans
    total <- total + ans$log_density
  }
  list(log_density = total + margin_log, margin_log_density = margin_log,
       conditional_share_log_density = total, blocks = block,
       rows = rows, layout = layout)
}

## Independent cross-check in the current Stan orientation: condition the
## margin on the share and apply the unit-Jacobian triangular transform.  This
## is mathematically the same density as market_structured_log_likelihood(),
## but has different sufficient statistics (margin rows load -lambda on the
## quality).  Keeping it here makes the reverse-conditioning identity
## auditable without calling the production implementation.
stan_style_augmented_rows <- function(share_resid, margin_resid, margin_idx,
                                      sigma_share, sigma_margin, rho) {
  n <- length(share_resid); margin_idx <- as.integer(margin_idx)
  if (length(margin_resid) != length(margin_idx) || anyDuplicated(margin_idx) ||
      any(margin_idx < 1L | margin_idx > n)) stop("invalid margin_idx")
  sigma_share <- finite_scalar(sigma_share, "sigma_share", lower = 0, strict_lower = TRUE)
  sigma_margin <- finite_scalar(sigma_margin, "sigma_margin", lower = 0, strict_lower = TRUE)
  rho <- finite_scalar(rho, "rho")
  if (abs(rho) >= 1) stop("rho must satisfy |rho| < 1")
  if (any(!is.finite(share_resid)) || any(!is.finite(margin_resid)))
    stop("residuals must be finite")
  lambda <- rho * sigma_margin / sigma_share
  list(z = c(share_resid, margin_resid - lambda * share_resid[margin_idx]),
       x = c(rep(1, n), rep(-lambda, length(margin_idx))),
       v = c(rep(sigma_share^2, n),
             rep(sigma_margin^2 * (1 - rho^2), length(margin_idx))),
       margin_idx = margin_idx, lambda = lambda)
}

stan_style_augmented_structured_log_likelihood <- function(share_resid,
                                                            margin_resid,
                                                            margin_idx,
                                                            layout, tau,
                                                            sigma_share,
                                                            sigma_margin, rho) {
  a <- stan_style_augmented_rows(share_resid, margin_resid, margin_idx,
                                 sigma_share, sigma_margin, rho)
  total <- 0
  blocks <- vector("list", layout$n_market)
  for (m in seq_len(layout$n_market)) {
    ii <- layout$rows[[m]]; jj <- which(layout$market[margin_idx] == m)
    if (length(jj)) {
      z <- c(a$z[ii], a$z[length(share_resid) + jj])
      x <- c(a$x[ii], a$x[length(share_resid) + jj])
      v <- c(a$v[ii], a$v[length(share_resid) + jj])
      groups <- c(layout$group_local[ii], layout$group_local[margin_idx[jj]])
    } else {
      z <- a$z[ii]; x <- a$x[ii]; v <- a$v[ii]
      groups <- layout$group_local[ii]
    }
    if (!length(ii)) {
      blocks[[m]] <- list(log_density = 0, rows = ii)
      next
    }
    ans <- structured_centered_log_density(z, v, groups, x, tau,
                                           G = layout$n_group[m])
    ans$rows <- ii; blocks[[m]] <- ans; total <- total + ans$log_density
  }
  list(log_density = total, blocks = blocks, rows = a, layout = layout)
}

posterior_quality_fixed_state_augmented_structured <- function(share_resid,
                                                               margin_resid,
                                                               margin_idx,
                                                               layout, tau,
                                                               sigma_share,
                                                               sigma_margin,
                                                               rho) {
  a <- stan_style_augmented_rows(share_resid, margin_resid, margin_idx,
                                 sigma_share, sigma_margin, rho)
  n <- length(share_resid); tau <- finite_scalar(tau, "tau", lower = 0)
  out <- vector("list", layout$n_market)
  for (m in seq_len(layout$n_market)) {
    ii <- layout$rows[[m]]; jj <- which(layout$market[margin_idx] == m)
    G <- layout$n_group[m]
    if (length(jj)) {
      z <- c(a$z[ii], a$z[n + jj]); x <- c(a$x[ii], a$x[n + jj]);
      v <- c(a$v[ii], a$v[n + jj]); groups <- c(layout$group_local[ii],
                                                  layout$group_local[margin_idx[jj]])
    } else {
      z <- a$z[ii]; x <- a$x[ii]; v <- a$v[ii]; groups <- layout$group_local[ii]
    }
    if (!length(ii) || tau == 0) {
      out[[m]] <- list(group_mean = numeric(G), group_covariance = matrix(0, G, G), rows = ii)
      next
    }
    s <- group_stats(z, v, groups, x, G)
    d <- tau^2 / (1 + tau^2 * s$a); den <- sum(d); raw_mean <- d * s$c
    means <- raw_mean - d * sum(raw_mean) / den
    covariance <- diag(d) - tcrossprod(d) / den
    out[[m]] <- list(group_mean = means, group_covariance = covariance,
                     rows = ii, stats = s, d = d)
  }
  out
}

make_global_quality_system <- function(group) {
  labels <- sort(unique(as.integer(group)))
  Z <- one_hot_local(as.integer(group), labels)
  C <- helmert_basis_local(length(labels))
  list(labels = labels, Z = Z, C = C, H = Z %*% C)
}

log_mvnorm_dense <- function(x, mean = rep(0, length(x)), covariance) {
  x <- as.numeric(x); mean <- as.numeric(mean)
  if (!length(x)) return(0)
  R <- chol(covariance)
  d <- x - mean
  sol <- backsolve(R, forwardsolve(t(R), d))
  -0.5 * (length(x) * log(2 * pi) +
           2 * sum(log(diag(R))) + sum(d * sol))
}

## Explicit completion-of-the-square evaluator.  Keeping the normalized
## likelihood, normalized N(0,tau^2 I) prior, and (2*pi)^(d/2)|A|^-1/2
## Gaussian integral visible is useful for auditing all constants in the
## sampled/collapsed equality.
exact_completion_log_likelihood <- function(share_resid, margin_resid,
                                            margin_idx, layout, tau,
                                            sigma_share, sigma_margin, rho) {
  rows <- conditional_share_rows(share_resid, layout$market,
                                 layout$group_local, margin_resid, margin_idx,
                                 sigma_share, sigma_margin, rho)
  tau <- finite_scalar(tau, "tau", lower = 0)
  margin_log <- if (length(margin_resid)) {
    sum(-0.5 * (log(2 * pi * sigma_margin^2) +
                (margin_resid / sigma_margin)^2))
  } else 0
  total <- margin_log
  blocks <- vector("list", layout$n_market)
  for (m in seq_len(layout$n_market)) {
    ii <- layout$rows[[m]]; G <- layout$n_group[m]; d <- max(0L, G - 1L)
    if (!length(ii)) {
      blocks[[m]] <- list(log_density = 0, rows = ii)
      next
    }
    z <- rows$z[ii]; v <- rows$v[ii]
    V_inv <- diag(1 / v, nrow = length(ii))
    if (!d || tau == 0) {
      q <- sum(z * z / v)
      ld <- sum(log(v))
      value <- -0.5 * (length(ii) * log(2 * pi) + ld + q)
      blocks[[m]] <- list(log_density = value, quadratic = q, rows = ii)
      total <- total + value
      next
    }
    Z <- one_hot_local(layout$group_local[ii], seq_len(G))
    H <- Z %*% helmert_basis_local(G)
    A <- diag(1 / tau^2, d) + crossprod(H, V_inv %*% H)
    h <- as.numeric(crossprod(H, V_inv %*% z))
    A_inv_h <- solve(A, h)
    q <- sum(z * z / v) - sum(h * A_inv_h)
    logdet_A <- 2 * sum(log(diag(chol(A))))
    logdet_V <- sum(log(v))
    ## Prior normalization and Gaussian integral are retained explicitly.
    value <- -0.5 * (length(ii) * log(2 * pi) + logdet_V +
                      d * log(tau^2) + logdet_A + q)
    blocks[[m]] <- list(log_density = value, quadratic = q,
                        logdet_V = logdet_V, logdet_A = logdet_A, rows = ii)
    total <- total + value
  }
  list(log_density = total, margin_log_density = margin_log, blocks = blocks,
       rows = rows, layout = layout)
}

## Dense market-block conditional oracle.  The market block covariance is
## diag(v_m) + tau^2 H_m H_m', and the observed margin factor is added once.
dense_market_conditional_log_likelihood <- function(share_resid, margin_resid,
                                                    margin_idx, layout, tau,
                                                    sigma_share, sigma_margin,
                                                    rho) {
  rows <- conditional_share_rows(share_resid, layout$market,
                                 layout$group_local, margin_resid, margin_idx,
                                 sigma_share, sigma_margin, rho)
  margin_log <- if (length(margin_resid)) {
    sum(-0.5 * (log(2 * pi * sigma_margin^2) +
                (margin_resid / sigma_margin)^2))
  } else 0
  total <- margin_log
  block <- vector("list", layout$n_market)
  for (m in seq_len(layout$n_market)) {
    ii <- layout$rows[[m]]
    if (!length(ii)) {
      block[[m]] <- list(log_density = 0, covariance = matrix(numeric(), 0, 0), rows = ii)
      next
    }
    Z <- one_hot_local(layout$group_local[ii], seq_len(layout$n_group[m]))
    C <- helmert_basis_local(layout$n_group[m])
    H <- Z %*% C
    covariance <- diag(rows$v[ii], nrow = length(ii)) + tau^2 * tcrossprod(H)
    value <- log_mvnorm_dense(rows$z[ii], covariance = covariance)
    block[[m]] <- list(log_density = value, covariance = covariance, H = H, rows = ii)
    total <- total + value
  }
  list(log_density = total, margin_log_density = margin_log,
       blocks = block, rows = rows, layout = layout)
}

## Original joint covariance oracle: share error/margin error covariance is
## rho*sigma_s*sigma_m on matched rows.  Quality loads only into shares.
dense_original_joint_log_likelihood <- function(share_resid, margin_resid,
                                                margin_idx, layout, tau,
                                                sigma_share, sigma_margin, rho) {
  n <- length(share_resid); q <- length(margin_resid)
  if (length(layout$market) != n) stop("layout and residual lengths disagree")
  margin_idx <- as.integer(margin_idx)
  if (length(margin_idx) != q || anyDuplicated(margin_idx) ||
      any(margin_idx < 1L | margin_idx > n)) stop("invalid margin_idx")
  sigma_share <- finite_scalar(sigma_share, "sigma_share", lower = 0, strict_lower = TRUE)
  sigma_margin <- finite_scalar(sigma_margin, "sigma_margin", lower = 0, strict_lower = TRUE)
  rho <- finite_scalar(rho, "rho")
  if (abs(rho) >= 1) stop("rho must satisfy |rho| < 1")
  q_error <- matrix(0, nrow = n + q, ncol = n + q)
  q_error[seq_len(n), seq_len(n)] <- diag(sigma_share^2, n)
  if (q) {
    q_error[n + seq_len(q), n + seq_len(q)] <- diag(sigma_margin^2, q)
    q_error[cbind(margin_idx, n + seq_len(q))] <- rho * sigma_share * sigma_margin
    q_error[cbind(n + seq_len(q), margin_idx)] <- rho * sigma_share * sigma_margin
  }
  H <- matrix(0, nrow = n + q, ncol = sum(pmax(layout$n_group - 1L, 0L)))
  offset <- 0L
  for (m in seq_len(layout$n_market)) {
    ii <- layout$rows[[m]]
    d <- max(0L, layout$n_group[m] - 1L)
    if (length(ii) && d) {
      Z <- one_hot_local(layout$group_local[ii], seq_len(layout$n_group[m]))
      H[ii, offset + seq_len(d)] <- Z %*% helmert_basis_local(layout$n_group[m])
    }
    offset <- offset + d
  }
  covariance <- q_error + tau^2 * tcrossprod(H)
  log_mvnorm_dense(c(share_resid, margin_resid), covariance = covariance)
}

## Direct contrast posterior moments for fixed k/state.  `margin_resid` is
## state-specific, but the resulting share posterior is the same conditional
## Gaussian algebra used by the likelihood.  The returned group moments are
## local to each market and satisfy sum_g E[b_mg]=0.
posterior_quality_fixed_state_dense <- function(share_resid, margin_resid,
                                                margin_idx, layout, tau,
                                                sigma_share, sigma_margin, rho) {
  rows <- conditional_share_rows(share_resid, layout$market,
                                 layout$group_local, margin_resid, margin_idx,
                                 sigma_share, sigma_margin, rho)
  tau <- finite_scalar(tau, "tau", lower = 0)
  out <- vector("list", layout$n_market)
  for (m in seq_len(layout$n_market)) {
    ii <- layout$rows[[m]]; G <- layout$n_group[m]; d <- max(0L, G - 1L)
    C <- helmert_basis_local(G)
    if (!length(ii) || !d || tau == 0) {
      out[[m]] <- list(contrast_mean = numeric(d),
                       contrast_covariance = matrix(0, d, d),
                       group_mean = numeric(G), group_covariance = matrix(0, G, G),
                       rows = ii)
      next
    }
    Z <- one_hot_local(layout$group_local[ii], seq_len(G))
    H <- Z %*% C
    W <- diag(1 / rows$v[ii], nrow = length(ii))
    precision <- diag(1 / tau^2, d) + crossprod(H, W %*% H)
    V <- solve(precision)
    mm <- as.numeric(V %*% crossprod(H, W %*% rows$z[ii]))
    out[[m]] <- list(contrast_mean = mm, contrast_covariance = V,
                     group_mean = as.numeric(C %*% mm),
                     group_covariance = C %*% V %*% t(C), rows = ii)
  }
  out
}

## Structured posterior reconstruction conditions independent group normals on
## the one linear zero-sum constraint.  This is algebraically equal to the
## direct Helmert contrast calculation above, but uses only O(N+G) statistics.
posterior_quality_fixed_state_structured <- function(share_resid, margin_resid,
                                                     margin_idx, layout, tau,
                                                     sigma_share, sigma_margin,
                                                     rho) {
  rows <- conditional_share_rows(share_resid, layout$market,
                                 layout$group_local, margin_resid, margin_idx,
                                 sigma_share, sigma_margin, rho)
  tau <- finite_scalar(tau, "tau", lower = 0)
  out <- vector("list", layout$n_market)
  for (m in seq_len(layout$n_market)) {
    ii <- layout$rows[[m]]; G <- layout$n_group[m]
    if (!length(ii) || tau == 0) {
      out[[m]] <- list(group_mean = numeric(G), group_covariance = matrix(0, G, G), rows = ii)
      next
    }
    s <- group_stats(rows$z[ii], rows$v[ii], layout$group_local[ii])
    tau2 <- tau^2
    d <- tau2 / (1 + tau2 * s$a)
    raw_mean <- d * s$c
    den <- sum(d)
    if (den == 0) {
      means <- numeric(G); covariance <- matrix(0, G, G)
    } else {
      means <- raw_mean - d * sum(raw_mean) / den
      covariance <- diag(d) - tcrossprod(d) / den
    }
    out[[m]] <- list(group_mean = means, group_covariance = covariance, rows = ii,
                     stats = s, d = d)
  }
  out
}

sum_market_moments <- function(x, y) {
  if (length(x) != length(y)) stop("moment block lengths disagree")
  max(c(0, vapply(seq_along(x), function(m) {
    max(abs(x[[m]]$group_mean - y[[m]]$group_mean),
        abs(x[[m]]$group_covariance - y[[m]]$group_covariance))
  }, numeric(1))))
}

## A self-contained expanded density with normalized quality priors.  This is
## the object sampled quality would use for a fixed state.
expanded_fixed_state_log_density <- function(share_resid, margin_resid,
                                             margin_idx, layout, beta_blocks,
                                             tau, sigma_share, sigma_margin, rho) {
  rows <- conditional_share_rows(share_resid, layout$market,
                                 layout$group_local, margin_resid, margin_idx,
                                 sigma_share, sigma_margin, rho)
  if (length(beta_blocks) != layout$n_market)
    stop("one beta block is required per market")
  total <- if (length(margin_resid)) {
    sum(-0.5 * (log(2 * pi * sigma_margin^2) +
                (margin_resid / sigma_margin)^2))
  } else 0
  for (m in seq_len(layout$n_market)) {
    ii <- layout$rows[[m]]; G <- layout$n_group[m]; d <- max(0L, G - 1L)
    beta <- as.numeric(beta_blocks[[m]])
    if (length(beta) != d) stop("beta block has incorrect dimension")
    C <- helmert_basis_local(G)
    b <- if (d) as.numeric(C %*% beta) else numeric(G)
    total <- total + sum(-0.5 * (log(2 * pi * rows$v[ii]) +
                                (rows$z[ii] - b[layout$group_local[ii]])^2 /
                                  rows$v[ii]))
    if (d && tau > 0) {
      total <- total + sum(-0.5 * (log(2 * pi * tau^2) + beta^2 / tau^2))
    }
  }
  total
}

gauss_hermite_rule <- function(n = 20L) {
  n <- as.integer(n)
  if (n < 3L) stop("Gauss-Hermite rule requires at least three nodes")
  ## Golub-Welsch eigenvalue construction for integral exp(-x^2) f(x) dx.
  off <- sqrt(seq_len(n - 1L) / 2)
  J <- matrix(0, n, n)
  J[cbind(seq_len(n - 1L), seq_len(n - 1L) + 1L)] <- off
  J[cbind(seq_len(n - 1L) + 1L, seq_len(n - 1L))] <- off
  e <- eigen(J, symmetric = TRUE)
  ord <- order(e$values)
  list(nodes = e$values[ord], weights = sqrt(pi) * e$vectors[1L, ord]^2)
}

## Numerical integration of a small expanded model, retained as an
## independent all-constants check.  Gauss-Hermite quadrature integrates over
## the normalized quality priors; the expanded density still supplies every
## Normal likelihood/prior constant before the prior term is removed to change
## the integration measure.  It is deliberately restricted to small contrast
## dimensions; the production/reference likelihood is analytic.
integrate_expanded_fixed_state <- function(share_resid, margin_resid, margin_idx,
                                            layout, tau, sigma_share, sigma_margin,
                                            rho, rel.tol = 2e-7) {
  dims <- pmax(layout$n_group - 1L, 0L)
  total_dim <- sum(dims)
  if (total_dim > 3L) stop("numerical integration fixture exceeds 3 dimensions")
  if (tau == 0 || total_dim == 0L) {
    return(expanded_fixed_state_log_density(
      share_resid, margin_resid, margin_idx, layout,
      lapply(dims, function(d) numeric(d)), tau, sigma_share, sigma_margin, rho
    ))
  }
  ## Thirty-five nodes per coordinate is still a small (<=42,875 point)
  ## check here, while reducing the quadrature remainder below the frozen
  ## 2e-7 comparison tolerance for the deliberately sharp share likelihood.
  rule <- gauss_hermite_rule(35L)
  grid <- expand.grid(rep(list(seq_along(rule$nodes)), total_dim))
  log_terms <- numeric(nrow(grid))
  prior_constant <- -0.5 * total_dim * log(2 * pi * tau^2)
  for (row in seq_len(nrow(grid))) {
    std <- sqrt(2) * rule$nodes[as.integer(grid[row, ])]
    beta_all <- tau * std
    blocks <- vector("list", length(dims)); off <- 0L
    log_weight <- -0.5 * total_dim * log(pi)
    for (j in seq_len(total_dim)) log_weight <- log_weight +
      log(rule$weights[as.integer(grid[row, j])])
    for (m in seq_along(dims)) {
      d <- dims[m]
      blocks[[m]] <- if (d) beta_all[off + seq_len(d)] else numeric()
      off <- off + d
    }
    expanded <- expanded_fixed_state_log_density(
      share_resid, margin_resid, margin_idx, layout, blocks,
      tau, sigma_share, sigma_margin, rho
    )
    ## Expanded = likelihood * normalized prior.  Quadrature weights here are
    ## for the normalized N(0,tau^2 I) measure, so divide that prior density.
    prior_log <- prior_constant - 0.5 * sum(beta_all^2 / tau^2)
    log_terms[row] <- log_weight + expanded - prior_log
  }
  max(log_terms) + log(sum(exp(log_terms - max(log_terms))))
}

logsumexp <- function(x) {
  if (!length(x)) return(-Inf)
  a <- max(x)
  a + log(sum(exp(x - a)))
}

## Exact K^M state oracle versus the factorized O(MK) expression.  Each
## market/state supplies its own margin prediction, so each state changes only
## that market's margin residuals.
enumerated_state_log_likelihood <- function(share_resid, margin_observed,
                                            margin_idx, layout, margin_mean,
                                            state_log_prior, tau, sigma_share,
                                            sigma_margin, rho) {
  M <- layout$n_market; K <- ncol(margin_mean)
  if (nrow(margin_mean) != length(margin_idx) || length(state_log_prior) != K)
    stop("state matrices have incompatible dimensions")
  if (nrow(margin_observed) != length(margin_idx) || ncol(margin_observed) != K)
    stop("margin_observed must be length(margin_idx) by K")
  choices <- expand.grid(rep(list(seq_len(K)), M))
  scores <- numeric(nrow(choices))
  for (r in seq_len(nrow(choices))) {
    states <- as.integer(choices[r, ])
    margin_resid <- numeric(length(margin_idx))
    for (j in seq_along(margin_idx)) {
      margin_resid[j] <- margin_observed[j, states[layout$market[margin_idx[j]]]] -
        margin_mean[j, states[layout$market[margin_idx[j]]]]
    }
    score <- market_structured_log_likelihood(
      share_resid, margin_resid, margin_idx, layout, tau,
      sigma_share, sigma_margin, rho
    )$log_density
    score <- score + sum(state_log_prior[states])
    scores[r] <- score
  }
  list(log_density = logsumexp(scores), scores = scores, choices = choices)
}

factorized_state_log_likelihood <- function(share_resid, margin_observed,
                                            margin_idx, layout, margin_mean,
                                            state_log_prior, tau, sigma_share,
                                            sigma_margin, rho) {
  M <- layout$n_market; K <- ncol(margin_mean)
  per_market <- numeric(M)
  for (m in seq_len(M)) {
    ii <- which(layout$market == m)
    jj <- which(layout$market[margin_idx] == m)
    local_scores <- numeric(K)
    for (k in seq_len(K)) {
      local_margin <- margin_observed[jj, k] - margin_mean[jj, k]
      local <- market_structured_log_likelihood(
        share_resid[ii], local_margin,
        match(margin_idx[jj], ii),
        make_market_layout(rep(1L, length(ii)), layout$group_local[ii], 1L),
        tau, sigma_share, sigma_margin, rho
      )$log_density
      local_scores[k] <- state_log_prior[k] + local
    }
    per_market[m] <- logsumexp(local_scores)
  }
  list(log_density = sum(per_market), per_market = per_market)
}

assert_close <- function(name, value, tolerance, table, notes = "") {
  value <- as.numeric(value)
  status <- if (is.finite(value) && value <= tolerance) "passed" else "failed"
  table[[length(table) + 1L]] <- data.frame(
    check = name, error = value, tolerance = tolerance, status = status,
    notes = notes, stringsAsFactors = FALSE
  )
  table
}

fixture_panel <- function() {
  ## Unequal local group counts, multiproduct rows, entry/exit over years, and
  ## repeated local groups across markets.  The time index documents the
  ## persistence interpretation but is not a likelihood input.
  market <- c(rep(1L, 8L), rep(2L, 6L), rep(3L, 9L))
  year <- c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L,
            1L, 1L, 2L, 2L, 2L, 2L,
            1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L)
  group <- c(1L, 1L, 2L, 3L, 1L, 2L, 3L, 4L,
             1L, 2L, 2L, 3L, 1L, 3L,
             1L, 1L, 2L, 1L, 2L, 3L, 1L, 4L, 4L)
  list(market = market, year = year, group = group,
       layout = make_market_layout(market, group))
}

run_validation <- function(output_dir = getwd(),
                           write_artifacts = TRUE, verbose = TRUE) {
  set.seed(20260916L)
  f <- fixture_panel(); layout <- f$layout; n <- length(f$market)
  share_resid <- rnorm(n, sd = 0.7)
  sigma_share <- 0.31; sigma_margin <- 0.47
  margin_idx <- c(1L, 4L, 6L, 10L, 14L, 18L, 22L)
  margin_resid <- rnorm(length(margin_idx), sd = 0.6)
  errors <- list()
  add <- function(name, value, t, notes = "") {
    errors[[length(errors) + 1L]] <<- data.frame(
      check = name, error = as.numeric(value), tolerance = t,
      status = if (is.finite(value) && value <= t) "passed" else "failed",
      notes = notes, stringsAsFactors = FALSE
    )
  }

  ## Helmert identities and market-local dimension accounting.
  for (m in seq_len(layout$n_market)) {
    G <- layout$n_group[m]; C <- helmert_basis_local(G)
    P <- diag(G) - matrix(1 / G, G, G)
    add(paste0("helmert_orthonormal_m", m), if (G > 1) max(abs(crossprod(C) - diag(G - 1))) else 0,
        tol$algebra)
    add(paste0("helmert_projection_m", m), max(abs(tcrossprod(C) - P)), tol$algebra)
    add(paste0("helmert_center_m", m), if (G > 1) max(abs(colSums(C))) else 0, tol$algebra)
  }

  ## Structured versus dense conditional market blocks at rho=0 and rho != 0,
  ## with tau=0 and ordinary positive tau.
  for (rho in c(0, 0.37, -0.42)) {
    for (tau in c(0, 0.43, 1.1)) {
      s <- market_structured_log_likelihood(
        share_resid, margin_resid, margin_idx, layout, tau,
        sigma_share, sigma_margin, rho
      )
      d <- dense_market_conditional_log_likelihood(
        share_resid, margin_resid, margin_idx, layout, tau,
        sigma_share, sigma_margin, rho
      )
      add(sprintf("structured_dense_rho_%s_tau_%s", rho, tau),
          abs(s$log_density - d$log_density), tol$algebra,
          "market-block conditional density")
      j <- dense_original_joint_log_likelihood(
        share_resid, margin_resid, margin_idx, layout, tau,
        sigma_share, sigma_margin, rho
      )
      add(sprintf("joint_conditional_factorization_rho_%s_tau_%s", rho, tau),
          abs(j - s$log_density), tol$algebra,
          "original joint covariance equals margin x conditional share")
      a <- stan_style_augmented_structured_log_likelihood(
        share_resid, margin_resid, margin_idx, layout, tau,
        sigma_share, sigma_margin, rho
      )
      add(sprintf("reverse_vs_stan_style_rho_%s_tau_%s", rho, tau),
          abs(a$log_density - s$log_density), tol$algebra,
          "reverse conditioning equals current margin-given-share transform")
    }
  }

  ## Missing-margin compatibility and all-unpaired rho invariance.
  s0 <- market_structured_log_likelihood(
    share_resid, numeric(), integer(), layout, 0.61,
    sigma_share, sigma_margin, 0
  )
  s1 <- market_structured_log_likelihood(
    share_resid, numeric(), integer(), layout, 0.61,
    sigma_share, sigma_margin, 0.81
  )
  add("missing_margins_rho_invariance", abs(s0$log_density - s1$log_density), tol$algebra)

  ## Single-market global quality and market-block quality are identical.
  one <- which(f$market == 1L)
  one_layout <- make_market_layout(rep(1L, length(one)), f$group[one], 1L)
  one_margin <- which(match(margin_idx, one) %in% seq_along(one))
  ## The selected indices above intentionally preserve a sparse missing mask.
  one_idx <- match(margin_idx[one_margin], one)
  one_mres <- margin_resid[one_margin]
  global <- make_global_quality_system(f$group[one])
  ## Compare the same conditional covariance using the global ZC map directly.
  rr <- conditional_share_rows(share_resid[one], rep(1L, length(one)),
                               match(f$group[one], global$labels), one_mres,
                               one_idx, sigma_share, sigma_margin, 0.37)
  cov_global <- diag(rr$v, nrow = length(one)) +
    0.43^2 * tcrossprod(global$H)
  global_value <- log_mvnorm_dense(rr$z, covariance = cov_global) +
    sum(-0.5 * (log(2 * pi * sigma_margin^2) + (one_mres / sigma_margin)^2))
  market_value <- market_structured_log_likelihood(
    share_resid[one], one_mres, one_idx, one_layout, 0.43,
    sigma_share, sigma_margin, 0.37
  )$log_density
  add("single_market_global_equals_market", abs(global_value - market_value), tol$algebra)

  ## Dense direct-contrast posterior moments equal O(N+G) reconstruction for
  ## a fixed margin state (nonzero rho and missing margins included).
  pd <- posterior_quality_fixed_state_dense(
    share_resid, margin_resid, margin_idx, layout, 0.53,
    sigma_share, sigma_margin, 0.37
  )
  ps <- posterior_quality_fixed_state_structured(
    share_resid, margin_resid, margin_idx, layout, 0.53,
    sigma_share, sigma_margin, 0.37
  )
  add("fixed_state_direct_contrast_mean_cov", sum_market_moments(pd, ps), tol$moments,
      "group-coordinate posterior moments")
  completion <- exact_completion_log_likelihood(
    share_resid, margin_resid, margin_idx, layout, 0.53,
    sigma_share, sigma_margin, 0.37
  )
  collapsed_for_completion <- market_structured_log_likelihood(
    share_resid, margin_resid, margin_idx, layout, 0.53,
    sigma_share, sigma_margin, 0.37
  )
  add("sampled_collapsed_exact_completion", abs(
    completion$log_density - collapsed_for_completion$log_density
  ), tol$algebra, "explicit prior and Gaussian integral constants")
  pa <- posterior_quality_fixed_state_augmented_structured(
    share_resid, margin_resid, margin_idx, layout, 0.53,
    sigma_share, sigma_margin, 0.37
  )
  add("fixed_state_reverse_vs_stan_style_moments", sum_market_moments(ps, pa), tol$moments,
      "posterior moments agree under both triangular orientations")
  add("fixed_state_zero_sum_means", max(vapply(ps, function(x) abs(sum(x$group_mean)), numeric(1))),
      tol$moments)
  add("fixed_state_zero_sum_covariance", max(vapply(ps, function(x) max(abs(rowSums(x$group_covariance))), numeric(1))),
      tol$moments)

  ## tau=0 and tiny sigma stress cases are finite and match the dense diagonal
  ## oracle within a relaxed, pre-frozen scale-aware tolerance.
  tiny <- market_structured_log_likelihood(
    share_resid, margin_resid, margin_idx, layout, 0,
    1e-4, 0.47, 0.37
  )
  tiny_dense <- dense_market_conditional_log_likelihood(
    share_resid, margin_resid, margin_idx, layout, 0,
    1e-4, 0.47, 0.37
  )
  add("tau_zero_tiny_sigma", abs(tiny$log_density - tiny_dense$log_density), tol$tiny_sigma)
  add("tau_zero_quality_moments", max(vapply(
    posterior_quality_fixed_state_structured(share_resid, margin_resid, margin_idx,
      layout, 0, sigma_share, sigma_margin, 0.37),
    function(x) max(abs(x$group_mean), abs(x$group_covariance)), numeric(1)
  )), tol$moments)

  ## Expanded/sampled versus collapsed integration with every Normal constant.
  ## Use one market with three contrast coordinates so base-R quadrature stays
  ## a small, transparent independent check.
  small_rows <- which(f$market == 1L)
  small_layout <- make_market_layout(rep(1L, length(small_rows)), f$group[small_rows], 1L)
  small_midx <- which(f$market[margin_idx] == 1L)
  small_idx <- match(margin_idx[small_midx], small_rows)
  small_mres <- margin_resid[small_midx]
  collapsed_small <- market_structured_log_likelihood(
    share_resid[small_rows], small_mres, small_idx, small_layout, 0.29,
    sigma_share, sigma_margin, -0.31
  )$log_density
  integrated_small <- integrate_expanded_fixed_state(
    share_resid[small_rows], small_mres, small_idx, small_layout, 0.29,
    sigma_share, sigma_margin, -0.31, rel.tol = 2e-7
  )
  add("expanded_collapsed_exact_integration", abs(collapsed_small - integrated_small),
      tol$quadrature, "all likelihood and prior constants retained")

  ## K^M enumeration versus the market-factorized logsumexp.  The margin
  ## means are state/market-specific; shares and common parameters are fixed.
  K <- 3L; M <- layout$n_market
  margin_mean <- matrix(0, nrow = length(margin_idx), ncol = K)
  margin_observed <- matrix(0, nrow = length(margin_idx), ncol = K)
  for (j in seq_along(margin_idx)) {
    margin_mean[j, ] <- c(-0.11, 0.08, 0.24) + 0.03 * layout$market[margin_idx[j]]
    margin_observed[j, ] <- margin_resid[j] + margin_mean[j, 2L]
  }
  state_log_prior <- log(c(0.2, 0.5, 0.3))
  enumerated <- enumerated_state_log_likelihood(
    share_resid, margin_observed, margin_idx, layout, margin_mean,
    state_log_prior, 0.43, sigma_share, sigma_margin, 0.37
  )
  factored <- factorized_state_log_likelihood(
    share_resid, margin_observed, margin_idx, layout, margin_mean,
    state_log_prior, 0.43, sigma_share, sigma_margin, 0.37
  )
  add("K_to_M_enumeration_equals_market_logsumexp",
      abs(enumerated$log_density - factored$log_density), tol$state,
      sprintf("K=%d, M=%d; enumerated=%d states", K, M, nrow(enumerated$choices)))

  ## A direct conditional/original joint comparison also covers nonzero rho
  ## when one market has no observed margin rows.
  no_margin_layout <- make_market_layout(rep(1L, 4L), c(1L, 1L, 2L, 3L), 1L)
  no_margin_share <- rnorm(4L)
  no_margin_s <- market_structured_log_likelihood(
    no_margin_share, numeric(), integer(), no_margin_layout, 0.7,
    0.25, 0.44, -0.5
  )$log_density
  no_margin_j <- dense_original_joint_log_likelihood(
    no_margin_share, numeric(), integer(), no_margin_layout, 0.7,
    0.25, 0.44, -0.5
  )
  add("all_missing_margin_joint_factorization", abs(no_margin_s - no_margin_j), tol$algebra)

  result <- do.call(rbind, errors)
  if (write_artifacts) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(result, file.path(output_dir, "error_table.csv"), row.names = FALSE)
    failed <- result[result$status == "failed", , drop = FALSE]
    write.csv(failed, file.path(output_dir, "failed_attempts.csv"), row.names = FALSE)
    write.csv(result, file.path(output_dir, "attempt_history.csv"), row.names = FALSE)
  }
  if (verbose) {
    print(result, row.names = FALSE)
    cat("\nchecks:", nrow(result), " failed:", sum(result$status == "failed"), "\n")
  }
  if (any(result$status == "failed")) {
    stop("market-quality oracle validation failed; see failed_attempts.csv")
  }
  invisible(result)
}

if (identical(environment(), globalenv())) {
  args <- commandArgs(trailingOnly = TRUE)
  out_arg <- args[startsWith(args, "--out=")]
  out <- if (length(out_arg)) sub("^--out=", "", out_arg[[1L]]) else
    file.path("doc", "statistical_toolkit_runs", "20260916", "performance-audit",
              "recovery", "market-quality", "reference-repetition")
  run_validation(output_dir = out)
}
