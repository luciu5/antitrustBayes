#!/usr/bin/env Rscript

## Independent dense reference for the proposed joint market/s0 rotation.
## This is deliberately separate from the package implementation.  It uses a
## fully normalized dense Gaussian likelihood with market-centred persistent
## quality, observation-matched share/margin covariance, unequal market sizes,
## and missing margins.  Production code must continue to use the structured
## O(N + G) likelihood.

options(stringsAsFactors = FALSE)

eps <- 1e-6
bounded_inv_logit <- function(x) eps + (1 - 2 * eps) * plogis(x)
effective_logit <- function(x) qlogis(bounded_inv_logit(x))
stable_hypot <- function(x, y) {
  scale <- max(abs(c(x, y)))
  if (scale == 0) return(0)
  scale * sqrt((x / scale)^2 + (y / scale)^2)
}

log_normal <- function(x, mean = 0, sd = 1) {
  sum(dnorm(x, mean, sd, log = TRUE))
}

log_half_normal <- function(x, sd) {
  sum(dnorm(x, 0, sd, log = TRUE) - log(0.5))
}

log_mvn <- function(y, mean, Sigma) {
  R <- chol(Sigma)
  residual <- y - mean
  whitened <- forwardsolve(t(R), residual)
  -0.5 * (length(y) * log(2 * pi) +
            2 * sum(log(diag(R))) + sum(whitened^2))
}

make_fixture <- function() {
  n_market <- c(3L, 5L, 4L)
  market <- rep(seq_along(n_market), n_market)
  local_group <- unlist(lapply(n_market, seq_len), use.names = FALSE)
  group_offset <- c(0L, cumsum(n_market)[-length(n_market)])
  group <- local_group + group_offset[market]
  N <- length(market)
  price <- seq(0.82, 1.37, length.out = N)
  firm_share <- unlist(lapply(n_market, function(n) {
    w <- seq(1, 2, length.out = n)
    0.72 * w / sum(w)
  }), use.names = FALSE)
  margin_observed <- rep(TRUE, N)
  margin_observed[c(2L, 7L, 12L)] <- FALSE
  list(
    M = length(n_market), N = N, n_market = n_market,
    market = market, group = group, local_group = local_group,
    price = price, firm_share = firm_share,
    margin_observed = margin_observed,
    y_share = log(firm_share) + 0.017 * sin(seq_len(N)),
    y_margin = log(4.1 + 0.15 * cos(seq_len(N))),
    alpha = 0.84, sigma_share = 0.19, sigma_margin = 0.14,
    rho = 0.37, tau = 0.21
  )
}

joint_to_legacy <- function(v, h, sigma_b, sigma_s0) {
  S <- stable_hypot(sigma_b, sigma_s0)
  list(
    u = (sigma_b * v - sigma_s0 * h) / S,
    z = (sigma_s0 * v + sigma_b * h) / S,
    S = S
  )
}

legacy_to_joint <- function(u, z, sigma_b, sigma_s0) {
  S <- stable_hypot(sigma_b, sigma_s0)
  list(
    v = (sigma_b * u + sigma_s0 * z) / S,
    h = (-sigma_s0 * u + sigma_b * z) / S,
    S = S
  )
}

economic_state <- function(u, z, mu_b, sigma_b, mu_s0, sigma_s0) {
  b <- mu_b + sigma_b * u
  eta <- mu_s0 + sigma_s0 * z
  list(b = b, eta = eta, s0 = bounded_inv_logit(eta))
}

quality_covariance <- function(fx) {
  out <- matrix(0, fx$N, fx$N)
  first <- 1L
  for (n in fx$n_market) {
    idx <- first:(first + n - 1L)
    P <- diag(n) - matrix(1 / n, n, n)
    out[idx, idx] <- fx$tau^2 * P
    first <- first + n
  }
  out
}

likelihood <- function(fx, state) {
  eta_effective <- effective_logit(state$eta)
  share_mean <- eta_effective[fx$market] + state$b[fx$market] -
    fx$alpha * fx$price
  ## Smooth supply mean that exercises both the outside share and the same
  ## market mapping as the production Cournot path.
  margin_mean <- log(fx$alpha * fx$price *
    (fx$firm_share + state$s0[fx$market]) /
    pmax(state$s0[fx$market], 1e-12))

  midx <- which(fx$margin_observed)
  y <- c(fx$y_share, fx$y_margin[midx])
  mean <- c(share_mean, margin_mean[midx])
  Sigma <- matrix(0, length(y), length(y))
  Sigma[seq_len(fx$N), seq_len(fx$N)] <-
    quality_covariance(fx) + diag(fx$sigma_share^2, fx$N)
  margin_pos <- fx$N + seq_along(midx)
  Sigma[margin_pos, margin_pos] <- diag(fx$sigma_margin^2, length(midx))
  for (j in seq_along(midx)) {
    Sigma[midx[j], margin_pos[j]] <- fx$rho * fx$sigma_share * fx$sigma_margin
    Sigma[margin_pos[j], midx[j]] <- Sigma[midx[j], margin_pos[j]]
  }
  list(
    value = log_mvn(y, mean, Sigma), b = state$b, eta = state$eta,
    s0 = state$s0, share_mean = share_mean, margin_mean = margin_mean,
    Sigma = Sigma
  )
}

common_prior <- function(mu_b, sigma_b, mu_s0, sigma_s0) {
  log_normal(mu_b, 0, 1) + log_half_normal(sigma_b, 0.5) +
    log_normal(mu_s0, -1.45, 0.46) + log_half_normal(sigma_s0, 0.5)
}

legacy_target <- function(x, fx) {
  M <- fx$M
  u <- x[seq_len(M)]
  z <- x[M + seq_len(M)]
  mu_b <- x[2 * M + 1L]
  sigma_b <- exp(x[2 * M + 2L])
  mu_s0 <- x[2 * M + 3L]
  sigma_s0 <- exp(x[2 * M + 4L])
  state <- economic_state(u, z, mu_b, sigma_b, mu_s0, sigma_s0)
  ll <- likelihood(fx, state)
  ## The last two coordinates are unconstrained log-scales.  Add their Stan
  ## transform Jacobians so this is the adjust_transform=TRUE target.
  target <- common_prior(mu_b, sigma_b, mu_s0, sigma_s0) +
    log_normal(u) + log_normal(z) + ll$value +
    log(sigma_b) + log(sigma_s0)
  c(list(target = target, u = u, z = z, mu_b = mu_b,
         sigma_b = sigma_b, mu_s0 = mu_s0, sigma_s0 = sigma_s0), ll)
}

shifted_target <- function(x, fx) {
  M <- fx$M
  u <- x[seq_len(M)]
  w <- x[M + seq_len(M)]
  mu_b <- x[2 * M + 1L]
  sigma_b <- exp(x[2 * M + 2L])
  mu_s0 <- x[2 * M + 3L]
  sigma_s0 <- exp(x[2 * M + 4L])
  b <- mu_b + sigma_b * u
  z <- w - b / sigma_s0
  state <- economic_state(u, z, mu_b, sigma_b, mu_s0, sigma_s0)
  ll <- likelihood(fx, state)
  target <- common_prior(mu_b, sigma_b, mu_s0, sigma_s0) +
    log_normal(u) + log_normal(z) + ll$value +
    log(sigma_b) + log(sigma_s0)
  c(list(target = target, u = u, z = z, w = w, mu_b = mu_b,
         sigma_b = sigma_b, mu_s0 = mu_s0, sigma_s0 = sigma_s0), ll)
}

joint_target <- function(x, fx) {
  M <- fx$M
  v <- x[seq_len(M)]
  h <- x[M + seq_len(M)]
  mu_b <- x[2 * M + 1L]
  sigma_b <- exp(x[2 * M + 2L])
  mu_s0 <- x[2 * M + 3L]
  sigma_s0 <- exp(x[2 * M + 4L])
  legacy <- joint_to_legacy(v, h, sigma_b, sigma_s0)
  state <- economic_state(legacy$u, legacy$z, mu_b, sigma_b,
                          mu_s0, sigma_s0)
  ll <- likelihood(fx, state)
  target <- common_prior(mu_b, sigma_b, mu_s0, sigma_s0) +
    log_normal(v) + log_normal(h) + ll$value +
    log(sigma_b) + log(sigma_s0)
  c(list(target = target, v = v, h = h, u = legacy$u, z = legacy$z,
         S = legacy$S, mu_b = mu_b, sigma_b = sigma_b,
         mu_s0 = mu_s0, sigma_s0 = sigma_s0), ll)
}

finite_gradient <- function(fun, x, rel_step = 2e-6) {
  vapply(seq_along(x), function(j) {
    h <- rel_step * max(1, abs(x[j]))
    xp <- xm <- x
    xp[j] <- xp[j] + h
    xm[j] <- xm[j] - h
    (fun(xp) - fun(xm)) / (2 * h)
  }, numeric(1))
}

finite_jacobian <- function(fun, x, rel_step = 2e-6) {
  base <- fun(x)
  out <- matrix(NA_real_, length(base), length(x))
  for (j in seq_along(x)) {
    h <- rel_step * max(1, abs(x[j]))
    xp <- xm <- x
    xp[j] <- xp[j] + h
    xm[j] <- xm[j] - h
    out[, j] <- (fun(xp) - fun(xm)) / (2 * h)
  }
  out
}

map_legacy_to_shifted <- function(x, fx) {
  st <- legacy_target(x, fx)
  c(st$u, st$z + st$b / st$sigma_s0,
    st$mu_b, log(st$sigma_b), st$mu_s0, log(st$sigma_s0))
}

map_legacy_to_joint <- function(x, fx) {
  st <- legacy_target(x, fx)
  q <- legacy_to_joint(st$u, st$z, st$sigma_b, st$sigma_s0)
  c(q$v, q$h, st$mu_b, log(st$sigma_b), st$mu_s0, log(st$sigma_s0))
}

map_joint_to_legacy <- function(x, fx) {
  st <- joint_target(x, fx)
  c(st$u, st$z, st$mu_b, log(st$sigma_b), st$mu_s0, log(st$sigma_s0))
}

run_checks <- function() {
  fx <- make_fixture()
  M <- fx$M
  scenarios <- list(
    ordinary = c(sigma_b = 0.43, sigma_s0 = 0.37),
    s0_near_zero = c(sigma_b = 0.43, sigma_s0 = 1e-8),
    b_near_zero = c(sigma_b = 1e-8, sigma_s0 = 0.37),
    both_small = c(sigma_b = 2e-8, sigma_s0 = 1e-8)
  )
  rows <- list()
  grows <- list()
  add <- function(scenario, quantity, error, tolerance, detail) {
    rows[[length(rows) + 1L]] <<- data.frame(
      scenario = scenario, quantity = quantity, max_abs_error = error,
      tolerance = tolerance, passed = is.finite(error) && error <= tolerance,
      method = "independent_dense_reference",
      detail = detail
    )
  }
  addg <- function(scenario, quantity, error, tolerance, detail) {
    grows[[length(grows) + 1L]] <<- data.frame(
      scenario = scenario, quantity = quantity, max_abs_error = error,
      tolerance = tolerance, passed = is.finite(error) && error <= tolerance,
      method = "finite_difference_chain_rule",
      detail = detail
    )
  }

  for (nm in names(scenarios)) {
    scales <- scenarios[[nm]]
    x_old <- c(
      u = c(-0.7, 0.25, 1.1), z = c(0.45, -0.8, 0.2),
      mu_b = 0.17, log_sigma_b = log(scales[["sigma_b"]]),
      mu_s0 = -1.12, log_sigma_s0 = log(scales[["sigma_s0"]])
    )
    x_shifted <- map_legacy_to_shifted(x_old, fx)
    x_joint <- map_legacy_to_joint(x_old, fx)
    old <- legacy_target(x_old, fx)
    shifted <- shifted_target(x_shifted, fx)
    joint <- joint_target(x_joint, fx)

    add(nm, "orthogonality",
        max(abs(crossprod(matrix(c(
          joint$sigma_b / joint$S, joint$sigma_s0 / joint$S,
          -joint$sigma_s0 / joint$S, joint$sigma_b / joint$S
        ), 2, 2)) - diag(2))), 2e-15,
        "Q'Q versus identity")
    add(nm, "determinant",
        abs((joint$sigma_b^2 + joint$sigma_s0^2) / joint$S^2 - 1),
        2e-15, "absolute determinant of the orthogonal map")
    for (field in c("b", "eta", "s0", "share_mean", "margin_mean", "Sigma")) {
      add(nm, paste0("joint_", field),
          max(abs(old[[field]] - joint[[field]])), 2e-10,
          "legacy economic state versus joint coordinate")
      add(nm, paste0("shifted_", field),
          max(abs(old[[field]] - shifted[[field]])), 2e-10,
          "legacy economic state versus current shifted coordinate")
    }
    add(nm, "joint_normalized_target", abs(old$target - joint$target),
        2e-9, "all normalized priors, transform Jacobians, and dense likelihood")
    add(nm, "shifted_normalized_target", abs(old$target - shifted$target),
        2e-9, "current shifted rotation target")
    add(nm, "bounded_coefficients",
        max(0, abs(c(joint$sigma_b / joint$S,
                     joint$sigma_s0 / joint$S)) - 1), 2e-15,
        "joint-rotation scale ratios remain bounded")

    ## Chain-rule gradient oracle.  Both targets are independently evaluated;
    ## the Jacobian maps the full joint coordinate (including log-scales) back
    ## to the legacy coordinate.  Its scale columns therefore test the
    ## scale-dependent orthogonal rotation and block-triangular argument.
    grad_old <- finite_gradient(function(x) legacy_target(x, fx)$target, x_old)
    grad_joint <- finite_gradient(function(x) joint_target(x, fx)$target, x_joint)
    J <- finite_jacobian(function(x) map_joint_to_legacy(x, fx), x_joint)
    expected <- drop(crossprod(J, grad_old))
    addg(nm, "full_target_chain_rule", max(abs(grad_joint - expected)),
         if (nm == "ordinary") 2e-4 else 2e-3,
         "finite-difference targets with full scale-dependent coordinate Jacobian")

    ## Analytic coordinate derivatives of u,z at fixed scales provide an
    ## independent check of the two raw-coordinate columns.
    Q <- matrix(c(joint$sigma_b / joint$S,
                  joint$sigma_s0 / joint$S,
                  -joint$sigma_s0 / joint$S,
                  joint$sigma_b / joint$S), 2, 2)
    ## Q acts separately for every market, so reorder the stacked result.
    expected_raw <- c(
      joint$sigma_b / joint$S * grad_old[seq_len(M)] +
        joint$sigma_s0 / joint$S * grad_old[M + seq_len(M)],
      -joint$sigma_s0 / joint$S * grad_old[seq_len(M)] +
        joint$sigma_b / joint$S * grad_old[M + seq_len(M)]
    )
    addg(nm, "raw_coordinate_chain_rule",
         max(abs(grad_joint[seq_len(2 * M)] - expected_raw)),
         if (nm == "ordinary") 8e-5 else 8e-4,
         "analytic orthogonal-map chain rule for v and h")
  }

  target_rows <- do.call(rbind, rows)
  gradient_rows <- do.call(rbind, grows)
  write.csv(target_rows, "doc/s0_joint_rotation_target_checks.csv",
            row.names = FALSE)
  write.csv(gradient_rows, "doc/s0_joint_rotation_gradient_checks.csv",
            row.names = FALSE)
  if (any(!target_rows$passed) || any(!gradient_rows$passed)) {
    print(target_rows)
    print(gradient_rows)
    stop("joint-rotation reference checks failed")
  }
  cat(sprintf("passed %d target and %d gradient checks\n",
              nrow(target_rows), nrow(gradient_rows)))
  invisible(list(target = target_rows, gradient = gradient_rows))
}

if (identical(environment(), globalenv())) run_checks()
