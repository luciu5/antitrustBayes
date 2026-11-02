#!/usr/bin/env Rscript

## Independent fixed-target reference for the market/outside-share rotation.
##
## This file is deliberately self-contained.  It mirrors the live Stan
## contracts that matter for this coordinate experiment:
##   * bounded_inv_logit() followed by logit() in the demand mean;
##   * unpooled, market-year, and market-level s0 bases;
##   * centered and noncentered s0 hierarchies;
##   * market-intercept hierarchy and the single-market b=0 convention;
##   * Logit Cournot's kappa/s0 alpha construction;
##   * CES's quantity-s0 to revenue-R0 conversion; and
##   * sampled and collapsed Gaussian quality likelihoods.
##
## No production code is sourced or modified.  Running this script writes only
## /tmp/market_s0_rotation_reference_results.csv and prints its checks.

options(stringsAsFactors = FALSE)

EPS <- 1e-6
PI <- pi

bounded_inv_logit <- function(u) EPS + (1 - 2 * EPS) * plogis(u)
effective_eta <- function(u) qlogis(bounded_inv_logit(u))
effective_eta_derivative <- function(u) {
  p <- plogis(u)
  h <- bounded_inv_logit(u)
  (1 - 2 * EPS) * p * (1 - p) / (h * (1 - h))
}
log_norm <- function(x, mean = 0, sd = 1) {
  sum(dnorm(x, mean = mean, sd = sd, log = TRUE))
}
log_half_norm <- function(x, sd) {
  sum(dnorm(x, mean = 0, sd = sd, log = TRUE) - log(0.5))
}
log_trunc_norm <- function(x, sd, lower, upper) {
  log_norm(x, 0, sd) - length(x) *
    log(pnorm(upper, 0, sd) - pnorm(lower, 0, sd))
}

make_Q <- function(G) {
  if (G <= 1L) return(matrix(0, nrow = G, ncol = 0L))
  Q <- matrix(0, nrow = G, ncol = G - 1L)
  for (k in seq_len(G - 1L)) {
    Q[seq_len(k), k] <- 1 / sqrt(k * (k + 1))
    Q[k + 1L, k] <- -k / sqrt(k * (k + 1))
  }
  Q
}

make_fixture <- function(M, T, J = 4L, K = 2L, seed = 20260916L) {
  set.seed(seed + 100L * M + 10L * T + J)
  L <- M * T
  mky_market <- rep(seq_len(M), each = T)
  mky_year <- rep(seq_len(T), times = M)
  market <- rep(mky_market, each = J)
  mky <- rep(seq_len(L), each = J)
  group <- rep(seq_len(J), times = L)
  N <- length(market)
  X_s0 <- cbind(intercept = seq(-0.15, 0.15, length.out = L),
                trend = (seq_len(L) - mean(seq_len(L))) / max(1, L))
  if (K == 0L) X_s0 <- matrix(0, nrow = L, ncol = 0L)
  if (K != ncol(X_s0)) stop("reference fixture only supports K in {0,2}")
  price <- rep(seq(0.8, 1.35, length.out = J), times = L) + rnorm(N, 0, 0.01)
  price_std <- as.numeric(scale(price))
  log_price_std <- as.numeric(scale(log(price)))
  outside_price <- seq(0.75, 1.4, length.out = L)
  inside_price_index <- seq(0.9, 1.5, length.out = L)
  firm_share <- rep(seq(0.10, 0.40, length.out = J), times = L)
  firm_rev_share <- rep(seq(0.12, 0.42, length.out = J), times = L)
  eta_truth <- seq(-1.5, -0.45, length.out = L)
  b_truth <- if (M > 1L) seq(-0.18, 0.22, length.out = M) else 0
  beta_truth <- if (K) c(0.12, -0.08) else numeric(0)
  q_truth <- if (J > 1L) seq(-0.15, 0.18, length.out = J - 1L) else numeric(0)
  u_truth <- eta_truth + drop(X_s0 %*% beta_truth)
  s0_truth <- bounded_inv_logit(u_truth)
  y_mean <- effective_eta(u_truth) + b_truth[mky_market] +
    rep(drop(make_Q(J) %*% q_truth), times = L) + 0.08 * price_std
  y_share <- y_mean + rnorm(N, 0, 0.045)
  ## A fixed k state is enough to exercise the s0/alpha/FOC dependency.  The
  ## production exact-k code supplies the same k to each row in a market.
  k_by_market <- if (M > 1L) seq_len(M) - 1L else 1L
  y_margin <- 0.3 + 0.05 * sin(seq_len(N)) + rnorm(N, 0, 0.03)
  margin_obs <- rep(TRUE, N)
  if (N >= 2L) margin_obs[c(2L, N)] <- FALSE
  list(
    M = M, T = T, L = L, J = J, K = K, N = N,
    mky_market = mky_market, mky_year = mky_year,
    market = market, mky = mky, group = group,
    market_sign = -1L,
    X_s0 = X_s0, price = price, price_std = price_std,
    log_price_std = log_price_std, price_sd = sd(price),
    outside_price = outside_price, inside_price_index = inside_price_index,
    firm_share = firm_share, firm_rev_share = firm_rev_share,
    y_share = y_share, y_margin = y_margin, margin_obs = margin_obs,
    k_by_market = k_by_market, Q = make_Q(J),
    s0_prior_mean = -1.446786403764666, s0_prior_sd = 0.456234014550216,
    prior_sigma_s0_pool = 0.5, prior_sigma_beta_s0 = 1,
    prior_alpha_mean = 0, prior_alpha_sd = 1,
    prior_sigma_alpha = 0.5, share_sigma = 0.35,
    margin_sigma = 0.22, sigma_b_strat = 0.28,
    eta_truth = eta_truth, b_truth = b_truth,
    beta_truth = beta_truth, q_truth = q_truth
  )
}

make_parameters <- function(fx, s0_structure, centered_s0,
                            demand = c("logit", "ces"),
                            supply = c("bertrand", "cournot"),
                            quality = c("sampled", "collapsed"),
                            alpha_pool = FALSE) {
  demand <- match.arg(demand)
  supply <- match.arg(supply)
  quality <- match.arg(quality)
  if (!s0_structure %in% c("none", "market_year", "market"))
    stop("bad s0 structure")
  M <- fx$M
  ## The fixture stores economically meaningful eta levels.  A hierarchy's
  ## location and scale are held fixed across the paired evaluations.
  eta_lat <- switch(s0_structure,
    none = fx$eta_truth,
    market_year = fx$eta_truth,
    market = vapply(seq_len(M), function(m) mean(fx$eta_truth[fx$mky_market == m]), numeric(1))
  )
  mu_s0 <- -1.1
  sigma_s0 <- 0.37
  mu_b <- if (M > 1L) -0.04 else 0
  sigma_b <- if (M > 1L) 0.43 else 1
  b_raw <- if (M > 1L) (fx$b_truth - mu_b) / sigma_b else numeric(0)
    log_a <- log(0.9)
  log_kappa <- log(0.65)
  log_kappa_market <- if (M > 1L) rep(log_kappa, M) + seq(-0.08, 0.08, length.out = M) else numeric(0)
    mu_log_a_free <- log(0.85)
    sigma_log_a <- 0.26
    alpha_raw <- if (alpha_pool) seq(-0.25, 0.25, length.out = M) else numeric(0)
  if (s0_structure == "market" && M == 1L) {
    ## The live single-market market pool has only a direct mu_s0 parameter;
    ## there is no between-market sigma or market raw effect.
    eta_lat <- mu_s0
  }
  list(
    s0_structure = s0_structure, centered_s0 = centered_s0,
    demand = demand, supply = supply, quality = quality,
    alpha_pool = alpha_pool,
    eta_lat = eta_lat, mu_s0 = mu_s0, sigma_s0 = sigma_s0,
    mu_b = mu_b, sigma_b = sigma_b, b_raw = b_raw,
    beta_s0 = fx$beta_truth, q_centered = fx$q_truth,
    log_a = log_a, log_kappa = log_kappa,
    log_kappa_market = log_kappa_market,
    mu_log_a_free = mu_log_a_free, sigma_log_a = sigma_log_a,
    alpha_raw = alpha_raw,
    log_sigma_share = log(fx$share_sigma),
    log_sigma_margin = log(fx$margin_sigma), rho = 0.31,
    nu = 5
  )
}

expand_s0_latent <- function(lat, fx, s0_structure) {
  if (s0_structure %in% c("none", "market_year")) return(as.numeric(lat))
  as.numeric(lat[fx$mky_market])
}

expand_b <- function(par, fx) {
  if (fx$M == 1L) return(rep(0, fx$M))
  par$mu_b + par$sigma_b * par$b_raw
}

latent_b_for_s0 <- function(b_market, fx, s0_structure) {
  if (s0_structure %in% c("none", "market_year")) return(b_market[fx$mky_market])
  b_market
}

old_eta_and_raw <- function(par, fx) {
  if (par$s0_structure == "market" && fx$M == 1L) {
    return(list(eta_lat = par$eta_lat, raw = numeric(0)))
  }
  eta_lat <- par$eta_lat
  raw <- if (par$s0_structure == "none") numeric(0) else
    (eta_lat - par$mu_s0) / par$sigma_s0
  list(eta_lat = eta_lat, raw = raw)
}

## Coordinate map for the proposed exact raw affine rotation.
##
## centered/unpooled: c = eta + b, with unit Jacobian.
## noncentered actual-level option: c = mu + sigma*z + b; the inverse has
##   |dz/dc| = sigma^(-n), hence dnorm(z)-n*log(sigma).
## noncentered shifted-raw option: w = z + b/sigma.  Then c = mu+sigma*w,
##   eta = mu+sigma*w-b and |dz/dw|=1.  This is the exact rotation that
##   preserves the noncentered standard-normal coordinate style.
rotation_coordinates <- function(par, fx, option = c("legacy", "actual", "shifted")) {
  option <- match.arg(option)
  b_market <- expand_b(par, fx)
  old <- old_eta_and_raw(par, fx)
  b_lat <- latent_b_for_s0(b_market, fx, par$s0_structure)
  if (par$s0_structure == "market" && fx$M == 1L) {
    return(list(option = option, b_market = b_market, b_lat = b_lat,
                eta_lat = old$eta_lat, raw = numeric(0), c_lat = old$eta_lat,
                jacobian_log = 0))
  }
  if (option == "legacy") {
    if (par$centered_s0 || par$s0_structure == "none") {
      return(list(option = option, b_market = b_market, b_lat = b_lat,
                  eta_lat = old$eta_lat, raw = old$raw, c_lat = old$eta_lat,
                  jacobian_log = 0))
    }
    return(list(option = option, b_market = b_market, b_lat = b_lat,
                eta_lat = old$eta_lat, raw = old$raw, c_lat = old$eta_lat,
                jacobian_log = 0))
  }
  c_lat <- old$eta_lat + b_lat
  if (option == "actual") {
    raw <- if (par$centered_s0 || par$s0_structure == "none") numeric(0) else
      (c_lat - b_lat - par$mu_s0) / par$sigma_s0
    jac <- if (par$centered_s0 || par$s0_structure == "none") 0 else
      -length(c_lat) * log(par$sigma_s0)
  } else {
    if (par$centered_s0 || par$s0_structure == "none")
      stop("shifted option is for a noncentered hierarchy")
    raw <- (old$eta_lat - par$mu_s0) / par$sigma_s0
    c_lat <- par$mu_s0 + par$sigma_s0 * (raw + b_lat / par$sigma_s0)
    raw <- raw + b_lat / par$sigma_s0
    jac <- 0
  }
  list(option = option, b_market = b_market, b_lat = b_lat,
       eta_lat = c_lat - b_lat, raw = raw, c_lat = c_lat,
       jacobian_log = jac)
}

state_from_coordinates <- function(par, fx, coord) {
  ## coord$eta_lat is already reconstructed by rotation_coordinates().
  eta_mky <- expand_s0_latent(coord$eta_lat, fx, par$s0_structure)
  u <- eta_mky + drop(fx$X_s0 %*% par$beta_s0)
  s0 <- bounded_inv_logit(u)
  eta_eff <- effective_eta(u)
  b_market <- coord$b_market
  b_mky <- b_market[fx$market]
  if (par$demand == "logit") {
    a <- structural_a(par, fx, s0, coord)
    share_mean <- eta_eff + b_mky + fx$market_sign * a[fx$market] * fx$price_std
  } else {
    a <- structural_a(par, fx, s0, coord)
    ## Live CES uses one a_market (or one per market) and gamma=1+a/log_price_sd.
    log_price_sd <- sd(log(fx$price))
    gamma_market <- 1 + a / log_price_sd
    gamma_i <- gamma_market[fx$market]
    ## Use the actual product row's standardized log price.  The CES source
    ## uses market_sign=-1 in the output-market contract.
    share_mean <- eta_eff + b_mky +
      fx$market_sign *
      (a[fx$market] + log_price_sd) * fx$log_price_std +
      gamma_i * log(fx$outside_price[fx$mky])
  }
  ## The historical derived rotation, including the one-period convention.
  n_per_market <- tabulate(fx$mky_market, nbins = fx$M)
  gamma_market_derived <- ifelse(n_per_market > 1,
    b_market + vapply(seq_len(fx$M), function(m) mean(eta_eff[fx$mky_market == m]), numeric(1)),
    b_market)
  list(
    eta_mky = eta_mky, u = u, s0 = s0, eta_eff = eta_eff,
    effective_c_mky = eta_eff + b_mky, b_market = b_market,
    a_market = a, alpha_structural_market = a / fx$price_sd,
    share_mean = share_mean, gamma_market_derived = gamma_market_derived,
    R0 = plogis(qlogis(s0) + log(fx$outside_price) - log(fx$inside_price_index))
  )
}

structural_a <- function(par, fx, s0, coord) {
  M <- fx$M
  if (par$demand == "logit" && par$supply == "cournot") {
    mean_log_s0 <- vapply(seq_len(M), function(m) {
      mean(log(s0[fx$mky_market == m]))
    }, numeric(1))
    if (!par$alpha_pool) return(rep(exp(par$log_kappa - mean(mean_log_s0)), M))
    return(exp(par$log_kappa_market - mean_log_s0))
  }
  if (!par$alpha_pool) return(rep(exp(par$log_a), M))
  ## The non-Cournot market-alpha hierarchy is included to make the reference
  ## useful for an exact fixed-target comparison beyond the s0 block.
  exp(par$mu_log_a_free + par$sigma_log_a * par$alpha_raw)
}

margin_prediction <- function(par, fx, st) {
  m <- fx$market
  mky <- fx$mky
  k <- fx$k_by_market[fx$mky_market[mky]]
  if (par$demand == "logit") {
    if (par$supply == "cournot") conduct <- st$s0[mky] /
      (st$s0[mky] + fx$firm_share)
    else conduct <- 1 - fx$firm_share
    invm <- st$alpha_structural_market[m] * conduct * (1 + 0.04 * k)
    return(log(fx$price * invm + 1e-9))
  }
  gamma <- 1 + st$a_market / sd(log(fx$price))
  R0 <- st$R0[mky]
  r <- fx$firm_rev_share * (1 - R0)
  if (par$supply == "cournot") {
    lerner <- 1 / gamma[m] + (gamma[m] - 1) * r /
      (gamma[m] * (1 + (gamma[m] - 1) * R0))
    invm <- 1 / lerner
  } else invm <- gamma[m] - (gamma[m] - 1) * r
  log(invm * (1 + 0.04 * k) + 1e-9)
}

collapsed_normal <- function(y, mu, group, Q, sigma, tau) {
  V <- diag(sigma^2, length(y))
  if (ncol(Q)) {
    A <- Q[group, , drop = FALSE]
    V <- V + tau^2 * tcrossprod(A)
  }
  r <- y - mu
  ch <- chol(V)
  logdet <- 2 * sum(log(diag(ch)))
  -0.5 * (length(y) * log(2 * PI) + logdet +
    sum(backsolve(ch, forwardsolve(t(ch), r)) * r))
}

collapsed_joint_normal <- function(y_share, share_mean, y_margin, pred_margin,
                                   margin_idx, group, Q, sigma_share,
                                   sigma_margin, rho, tau) {
  z_s <- y_share - share_mean
  slope <- rho * sigma_margin / sigma_share
  z <- c(z_s, y_margin[margin_idx] - pred_margin[margin_idx] -
    slope * z_s[margin_idx])
  x <- c(rep(1, length(y_share)), rep(-slope, length(margin_idx)))
  g <- c(group, group[margin_idx])
  v <- c(rep(sigma_share^2, length(y_share)),
         rep(sigma_margin^2 * (1 - rho^2), length(margin_idx)))
  V <- diag(v, length(v))
  if (ncol(Q)) {
    A <- Q[g, , drop = FALSE] * x
    V <- V + tau^2 * tcrossprod(A)
  }
  ch <- chol(V)
  logdet <- 2 * sum(log(diag(ch)))
  -0.5 * (length(z) * log(2 * PI) + logdet +
    sum(backsolve(ch, forwardsolve(t(ch), z)) * z))
}

target <- function(par, fx, coord_option = c("legacy", "actual", "shifted"),
                   share_path = c("gaussian", "student_t"),
                   include_margin = TRUE) {
  share_path <- match.arg(share_path)
  coord_option <- match.arg(coord_option)
  coord <- rotation_coordinates(par, fx, coord_option)
  st <- state_from_coordinates(par, fx, coord)
  ans <- 0
  M <- fx$M

  ## b hierarchy, exactly normalized (the b=0 single-market branch is empty).
  if (M > 1L) {
    ans <- ans + log_norm(par$mu_b, 0, 1) +
      log_half_norm(par$sigma_b, 0.5) + log_norm(par$b_raw)
  }

  ## The prior is evaluated on reconstructed eta.  For the actual-level c
  ## coordinate, the change-of-variables Jacobian is -n log(sigma), which is
  ## added once at the end of this function.  Algebraically this is the same
  ## as the normalized Normal(c-b | mu,sigma) density.  The shifted raw
  ## coordinate has a unit Jacobian and retains std_normal(w-b/sigma).
  if (par$s0_structure == "market" && M == 1L) {
    ans <- ans + log_norm(coord$eta_lat, fx$s0_prior_mean, fx$s0_prior_sd)
  } else if (par$s0_structure == "none") {
    ans <- ans + log_norm(coord$eta_lat, fx$s0_prior_mean, fx$s0_prior_sd)
  } else {
    n_lat <- length(coord$eta_lat)
    ans <- ans + log_norm(par$mu_s0, fx$s0_prior_mean, fx$s0_prior_sd) +
      log_half_norm(par$sigma_s0, fx$prior_sigma_s0_pool)
    if (coord_option == "shifted" && !par$centered_s0) {
      raw <- coord$raw - coord$b_lat / par$sigma_s0
      ans <- ans + log_norm(raw)
    } else if (!par$centered_s0) {
      raw <- if (coord_option == "legacy") {
        coord$raw
      } else {
        (coord$c_lat - coord$b_lat - par$mu_s0) / par$sigma_s0
      }
      ## For actual c the -n log(sigma) term is coord$jacobian_log below;
      ## keeping it here too would apply the change-of-variables twice.
      ans <- ans + log_norm(raw)
    } else {
      eta_for_prior <- if (coord_option == "actual") {
        coord$c_lat - coord$b_lat
      } else {
        coord$eta_lat
      }
      ans <- ans + log_norm(eta_for_prior, par$mu_s0, par$sigma_s0)
    }
  }
  ans <- ans + log_norm(par$beta_s0, 0, fx$prior_sigma_beta_s0)

  ## alpha/Cournot terms.  These are functions of reconstructed s0 only and
  ## are therefore intentionally identical under every coordinate map.
  if (par$demand == "logit" && par$supply == "cournot" && !par$alpha_pool) {
    mu_log_a <- par$log_kappa - mean(log(st$s0))
    ans <- ans + log_norm(mu_log_a, fx$prior_alpha_mean, fx$prior_alpha_sd)
  } else if (par$alpha_pool && par$demand == "logit" && par$supply == "cournot") {
    mean_log_s0 <- vapply(seq_len(M), function(m) mean(log(st$s0[fx$mky_market == m])), numeric(1))
    log_a_derived <- par$log_kappa_market - mean_log_s0
    ans <- ans + log_norm(par$mu_log_a_free, fx$prior_alpha_mean, fx$prior_alpha_sd) +
      log_half_norm(par$sigma_log_a, fx$prior_sigma_alpha) +
      log_norm(log_a_derived, par$mu_log_a_free, par$sigma_log_a)
  } else if (!par$alpha_pool) {
    ans <- ans + log_norm(par$log_a, fx$prior_alpha_mean, fx$prior_alpha_sd)
  } else {
    ans <- ans + log_norm(par$mu_log_a_free, fx$prior_alpha_mean, fx$prior_alpha_sd) +
      log_half_norm(par$sigma_log_a, fx$prior_sigma_alpha) +
      log_norm(par$alpha_raw)
  }
  ans <- ans + log_norm(par$log_sigma_share, log(fx$share_sigma), 0.5) +
    log_norm(par$log_sigma_margin, log(fx$margin_sigma), 0.5)

  q_i <- if (ncol(fx$Q)) rep(drop(fx$Q %*% par$q_centered), times = fx$L) else rep(0, fx$N)
  sigma_share <- exp(par$log_sigma_share)
  sigma_margin <- exp(par$log_sigma_margin)
  if (par$quality == "sampled") {
    ans <- ans + log_norm(par$q_centered, 0, fx$sigma_b_strat)
    share_mu <- st$share_mean + q_i
    if (share_path == "gaussian") {
      ans <- ans + log_norm(fx$y_share, share_mu, sigma_share)
    } else {
      ans <- ans + sum(dt((fx$y_share - share_mu) /
        (sigma_share * sqrt((par$nu - 2) / par$nu)), par$nu, log = TRUE) -
        log(sigma_share * sqrt((par$nu - 2) / par$nu)))
    }
    if (include_margin && any(fx$margin_obs)) {
      idx <- which(fx$margin_obs)
      pred <- margin_prediction(par, fx, st)
      err_share <- fx$y_share - share_mu
      if (share_path == "gaussian" && abs(par$rho) > 0) {
        ans <- ans + log_norm(fx$y_margin[idx],
          pred[idx] + par$rho * sigma_margin / sigma_share * err_share[idx],
          sigma_margin * sqrt(1 - par$rho^2))
      } else if (share_path == "gaussian") {
        ans <- ans + log_norm(fx$y_margin[idx], pred[idx], sigma_margin)
      } else {
        ans <- ans + log_norm(fx$y_margin[idx], pred[idx], sigma_margin)
      }
    }
  } else {
    if (share_path != "gaussian") stop("collapsed quality requires Gaussian shares")
    idx <- if (include_margin) which(fx$margin_obs) else integer(0)
    pred <- margin_prediction(par, fx, st)
    if (length(idx) && abs(par$rho) > 0) {
      ans <- ans + collapsed_joint_normal(
        fx$y_share, st$share_mean, fx$y_margin, pred, idx,
        fx$group, fx$Q, sigma_share, sigma_margin, par$rho,
        fx$sigma_b_strat)
    } else {
      ans <- ans + collapsed_normal(fx$y_share, st$share_mean, fx$group,
        fx$Q, sigma_share, fx$sigma_b_strat)
      if (length(idx)) ans <- ans + log_norm(fx$y_margin[idx], pred[idx], sigma_margin)
    }
  }
  ans + coord$jacobian_log
}

finite_gradient <- function(fun, x, h = 1e-5) {
  out <- numeric(length(x))
  for (j in seq_along(x)) {
    xp <- xm <- x
    xp[j] <- xp[j] + h
    xm[j] <- xm[j] - h
    out[j] <- (fun(xp) - fun(xm)) / (2 * h)
  }
  out
}

## Evaluate a target after replacing the old or new s0 coordinate.  This is
## intentionally written as a coordinate map rather than an algebraic
## shortcut, so the finite-difference checks exercise the complete target,
## including Cournot's derived alpha and CES's R0.
target_at_s0_coordinate <- function(x, par, fx, coordinate,
                                    include_margin = TRUE) {
  p <- par
  b_market <- expand_b(p, fx)
  b_lat <- latent_b_for_s0(b_market, fx, p$s0_structure)
  if (coordinate == "old") {
    if (p$s0_structure == "market" && fx$M == 1L) p$eta_lat <- x
    else if (p$s0_structure == "none") p$eta_lat <- x
    else if (p$centered_s0) p$eta_lat <- x
    else p$eta_lat <- p$mu_s0 + p$sigma_s0 * x
    return(target(p, fx, "legacy", "gaussian", include_margin))
  }
  if (coordinate == "actual") {
    p$eta_lat <- x - b_lat
    return(target(p, fx, "actual", "gaussian", include_margin))
  }
  if (coordinate == "shifted") {
    p$eta_lat <- p$mu_s0 + p$sigma_s0 * (x - b_lat / p$sigma_s0)
    return(target(p, fx, "shifted", "gaussian", include_margin))
  }
  stop("unknown s0 coordinate")
}

gradient_checks <- function() {
  rows <- list()
  add <- function(name, value, detail) {
    rows[[length(rows) + 1L]] <<- data.frame(
      check = name, max_abs_error = as.numeric(value), tolerance = 3e-5,
      status = if (is.finite(value) && value < 3e-5) "passed" else "failed",
      detail = detail
    )
  }
  fx <- make_fixture(3L, 2L)
  for (s0_structure in c("none", "market_year", "market")) {
    for (centered in c(FALSE, TRUE)) {
      if (s0_structure == "none" && centered) next
      if (s0_structure == "market" && fx$M == 1L && !centered) next
      for (demand in c("logit", "ces")) {
        for (supply in c("bertrand", "cournot")) {
          p <- make_parameters(fx, s0_structure, centered, demand, supply, "sampled")
          p$rho <- 0.23
          old_coord <- if (s0_structure == "market" && fx$M == 1L) p$eta_lat else
            if (s0_structure == "none" || centered) p$eta_lat else
              (p$eta_lat - p$mu_s0) / p$sigma_s0
          b_lat <- latent_b_for_s0(expand_b(p, fx), fx, s0_structure)
          actual_coord <- old_coord
          if (!(s0_structure == "market" && fx$M == 1L)) actual_coord <-
            p$eta_lat + b_lat
          old_grad <- finite_gradient(function(z) target_at_s0_coordinate(
            z, p, fx, "old"), old_coord)
          actual_grad <- finite_gradient(function(z) target_at_s0_coordinate(
            z, p, fx, "actual"), actual_coord)
          expected_actual <- if (!centered && s0_structure != "none" &&
                                 !(s0_structure == "market" && fx$M == 1L))
            old_grad / p$sigma_s0 else old_grad
          nm <- paste("gradient_actual", s0_structure, centered, demand,
                      supply, sep = "_")
          add(nm, max(abs(actual_grad - expected_actual)),
              "finite-difference chain rule; actual c")
          if (!centered && s0_structure != "none" &&
              !(s0_structure == "market" && fx$M == 1L)) {
            shifted_coord <- old_coord + b_lat / p$sigma_s0
            shifted_grad <- finite_gradient(function(z) target_at_s0_coordinate(
              z, p, fx, "shifted"), shifted_coord)
            add(paste("gradient_shifted", s0_structure, demand, supply, sep = "_"),
                max(abs(shifted_grad - old_grad)),
                "finite-difference chain rule; shifted noncentered w")
          }
        }
      }
    }
  }
  do.call(rbind, rows)
}

run_checks <- function() {
  rows <- list()
  add <- function(name, value, tolerance, status, detail = "") {
    rows[[length(rows) + 1L]] <<- data.frame(
      check = name, max_abs_error = as.numeric(value), tolerance = tolerance,
      status = status, detail = detail
    )
  }
  ## Dimension sweep covers one-cell, one-market panels, several markets with
  ## one period, and a genuinely repeated multi-market panel.
  dims <- list(c(1L, 1L), c(1L, 3L), c(3L, 1L), c(3L, 3L))
  for (dm in dims) {
    fx <- make_fixture(dm[1], dm[2])
    for (s0_structure in c("none", "market_year", "market")) {
      if (s0_structure == "market" && fx$M == 1L) centered_grid <- TRUE else centered_grid <- c(FALSE, TRUE)
      for (centered in centered_grid) {
        for (demand in c("logit", "ces")) {
          for (supply in c("bertrand", "cournot")) {
            for (quality in c("sampled", "collapsed")) {
              par <- make_parameters(fx, s0_structure, centered, demand, supply, quality)
              share_paths <- if (quality == "collapsed") "gaussian" else c("gaussian", "student_t")
              for (share_path in share_paths) {
                ## rho=0 and rho!=0 are both tested where the corresponding
                ## Gaussian path exists.  Student-t is tested at rho=0.
                for (rho in if (share_path == "gaussian") c(0, 0.31) else 0) {
                  par$rho <- rho
                  nm <- paste(dm[1], dm[2], s0_structure, centered,
                              demand, supply, quality, share_path, rho, sep = "_")
                  old <- target(par, fx, "legacy", share_path, include_margin = TRUE)
                  actual <- target(par, fx, "actual", share_path, include_margin = TRUE)
                  ## Densities with respect to coordinates differ by the
                  ## change-of-variables Jacobian.  For a noncentered
                  ## actual-level c, that Jacobian is -n*log(sigma); compare
                  ## the transformed legacy density to the rotated density.
                  jac_actual <- rotation_coordinates(par, fx, "actual")$jacobian_log
                  err <- abs((old + jac_actual) - actual)
                  add(paste0("actual_", nm), err, 2e-10,
                      if (is.finite(err) && err < 2e-10) "passed" else "failed",
                      "legacy plus change-of-variables Jacobian versus actual-level c")
                  if (!centered && s0_structure != "none" && !(s0_structure == "market" && fx$M == 1L)) {
                    shifted <- target(par, fx, "shifted", share_path, include_margin = TRUE)
                    err2 <- abs(old - shifted)
                    add(paste0("shifted_", nm), err2, 2e-10,
                        if (is.finite(err2) && err2 < 2e-10) "passed" else "failed",
                        "legacy versus unit-Jacobian shifted noncentered coordinate")
                  }
                  ## Missing margins must leave the share target untouched and
                  ## the paired target equal after the same rotation.
                  old_nom <- target(par, fx, "legacy", share_path, include_margin = FALSE)
                  new_nom <- target(par, fx, "actual", share_path, include_margin = FALSE)
                  err3 <- abs((old_nom + jac_actual) - new_nom)
                  add(paste0("missing_margin_", nm), err3, 2e-10,
                      if (is.finite(err3) && err3 < 2e-10) "passed" else "failed",
                      "margin mask excludes two synthetic rows")
                }
              }
            }
          }
        }
      }
    }
  }
  ## Nonzero bounded-link correction: this should fail at an extreme raw eta,
  ## and the failed control is deliberately retained in the output.
  fx <- make_fixture(2L, 2L)
  par <- make_parameters(fx, "none", TRUE, "logit", "bertrand", "sampled")
  par$eta_lat <- c(-18, -12, 8, 15)
  legacy_st <- state_from_coordinates(par, fx, rotation_coordinates(par, fx, "legacy"))
  old_mean <- legacy_st$share_mean
  ## Calculate the actual and naive Gaussian share targets on the same y.
  y0 <- fx$y_share
  actual_ll <- log_norm(y0, old_mean, exp(par$log_sigma_share))
  coord <- rotation_coordinates(par, fx, "actual")
  st <- state_from_coordinates(par, fx, coord)
  c_mky <- expand_s0_latent(coord$c_lat, fx, par$s0_structure)
  ## This expression replaces effective_eta(u)+b with c+X beta.
  naive_mean <- c_mky[fx$mky] + drop(fx$X_s0[fx$mky, , drop = FALSE] %*% par$beta_s0) +
    (st$share_mean - st$eta_eff[fx$mky] - st$b_market[fx$market]) +
    st$b_market[fx$market]
  naive_ll <- log_norm(y0, naive_mean, exp(par$log_sigma_share))
  failed_err <- abs(actual_ll - naive_ll)
  add("expected_failed_control_raw_c_without_bounded_link_correction",
      failed_err, 1e-8,
      if (is.finite(failed_err) && failed_err > 1e-8) "expected_fail_preserved" else "unexpected_pass",
      "raw c shortcut is not exact at extreme u")
  grad <- gradient_checks()
  out <- rbind(do.call(rbind, rows), grad)
  write.csv(out, "/tmp/market_s0_rotation_reference_results.csv", row.names = FALSE)
  print(out)
  if (any(out$status == "failed")) stop("one or more exact target checks failed")
  invisible(out)
}

if (identical(environment(), globalenv())) run_checks()
