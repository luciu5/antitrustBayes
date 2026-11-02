#!/usr/bin/env Rscript
## Independent CES Cournot check.  This script is deliberately outside the
## repository's production tree and imports no package code.
suppressPackageStartupMessages({
  library(numDeriv)
  library(nleqslv)
})

make_market <- function(gamma) {
  p0 <- 1.7; E <- 1
  p <- c(2.13, 2.46, 2.79, 3.12, 3.45, 3.78)
  owner <- c("A", "A", "B", "B", "C", "D")
  leaders <- c("A", "C")
  r0 <- .17
  r <- c(.18, .16, .14, .13, .12, .10)
  stopifnot(abs(sum(r) + r0 - 1) < 1e-14)
  A <- (r / r0) * (p / p0)^(gamma - 1)
  q <- E * r / p
  firms <- unique(owner)
  R <- vapply(firms, function(f) sum(r[owner == f]), numeric(1))
  names(R) <- firms
  H <- 1 + (gamma - 1) * r0
  b <- (gamma - 1) / (gamma * H)
  mf <- 1 / gamma + b * R
  af <- vapply(firms, function(f) {
    m <- mf[f]; aa <- 1 - m
    gamma * b * (aa - b * R[f] + R[f] * (gamma - 1)^2 * r0 / H^2) /
      (aa + (gamma - 1) * b * R[f])
  }, numeric(1))
  names(af) <- firms
  C <- sum(R[setdiff(firms, leaders)] * af[setdiff(firms, leaders)])
  ml <- 1 / gamma + b * R / (1 + C)
  m <- ifelse(firms %in% leaders, ml, mf)
  names(m) <- firms
  cst <- p * (1 - m[owner])
  list(gamma = gamma, p0 = p0, E = E, p = p, q = q, r = r, r0 = r0,
       A = A, owner = owner, leaders = leaders, firms = firms, R = R,
       H = H, b = b, mf = mf, af = af, C = C, ml = ml, m = m, c = cst)
}

inv_p <- function(q, mk) {
  stopifnot(all(q > 0))
  g <- mk$gamma; A <- mk$A; p0 <- mk$p0; E <- mk$E
  bb <- A^(1 / g) * q^(-1 / g)
  B <- sum(bb * q); outside <- p0^(1 - g)
  fun <- function(t) outside * t^g + B * t - E
  hi <- 1
  while (fun(hi) < 0) hi <- hi * 2
  t <- uniroot(fun, c(0, hi), tol = 1e-13)$root
  t * bb
}

inv_jac <- function(q, p, mk) {
  g <- mk$gamma; E <- mk$E
  r0 <- 1 - sum(p * q) / E
  H <- 1 + (g - 1) * r0
  b <- (g - 1) / (g * H)
  -diag(p / q) / g - b * outer(p, p) / E
}

q_from_p <- function(p, mk) {
  g <- mk$gamma
  den <- mk$p0^(1 - g) + sum(mk$A * p^(1 - g))
  mk$E * mk$A * p^(-g) / den
}

raw_foc <- function(q, mk) {
  p <- inv_p(q, mk); K <- inv_jac(q, p, mk)
  vapply(seq_along(q), function(k) {
    own <- mk$owner == mk$owner[k]
    p[k] - mk$c[k] + sum(q[own] * K[own, k])
  }, numeric(1))
}

solve_f <- function(q_l, q_start, mk, li, fi) {
  base <- q_start; base[li] <- q_l
  ff <- function(x) {
    z <- base; z[fi] <- x
    raw_foc(z, mk)[fi]
  }
  sol <- nleqslv(base[fi], ff, control = list(ftol = 1e-12, xtol = 1e-12,
                                               maxit = 300))
  if (max(abs(sol$fvec)) > 2e-8) stop("follower solve failed")
  base[fi] <- sol$x
  base
}

check <- function(gamma) {
  mk <- make_market(gamma); owner <- mk$owner
  li <- which(owner %in% mk$leaders); fi <- which(!(owner %in% mk$leaders))
  q0 <- mk$q; qsol <- solve_f(q0[li], q0, mk, li, fi)
  follower_foc <- max(abs(raw_foc(qsol, mk)[fi]))

  ## Claimed formula dq_j^F/dq_k^L = -q_j a_f p_k/[E(1+C)].
  claimed <- outer(q0[fi], mk$p[li], function(qj, pk) {
    -qj * mk$af[owner[fi]][1] * pk / (mk$E * (1 + mk$C))
  })
  ## Fix the row-specific firm coefficient (outer above uses first row).
  claimed <- sweep(claimed, 1, mk$af[owner[fi]], "*") /
    rep(mk$af[owner[fi]][1], each = length(li))
  numerical <- matrix(0, length(fi), length(li))
  for (kk in seq_along(li)) {
    h <- 2e-5 * q0[li[kk]]
    up <- q0[li]; dn <- q0[li]
    up[kk] <- up[kk] + h; dn[kk] <- dn[kk] - h
    numerical[, kk] <- (solve_f(up, q0, mk, li, fi)[fi] -
                          solve_f(dn, q0, mk, li, fi)[fi]) / (2 * h)
  }

  ## Reduced leader finite-difference gradients of the actual profit.
  profit <- function(x, firm) {
    z <- solve_f(x, q0, mk, li, fi); pp <- inv_p(z, mk)
    ix <- owner == firm
    sum((pp[ix] - mk$c[ix]) * z[ix])
  }
  grad <- numeric()
  for (firm in mk$leaders) {
    own <- which(owner == firm); pos <- match(own, li)
    for (kk in pos) {
      h <- 2e-5 * q0[li[kk]]; up <- q0[li]; dn <- q0[li]
      up[kk] <- up[kk] + h; dn[kk] <- dn[kk] - h
      grad <- c(grad, (profit(up, firm) - profit(dn, firm)) / (2 * h))
    }
  }
  ## Check the proposed K against a numerical derivative of inverse demand,
  ## whose values are obtained from the CES primitive rather than a copied
  ## Jacobian formula.  Also verify inverse-demand round trip.
  K_num <- jacobian(function(x) inv_p(x, mk), q0, method = "Richardson",
                    method.args = list(eps = 1e-5, d = 1e-5))
  inverse_err <- max(abs(inv_jac(q0, mk$p, mk) - K_num))
  roundtrip_err <- max(abs(q_from_p(inv_p(q0, mk), mk) - q0))
  c(gamma = gamma, follower_foc = follower_foc,
    reaction_error = max(abs(claimed - numerical)),
    leader_gradient = max(abs(grad)), inverse_error = inverse_err,
    roundtrip_error = roundtrip_err,
    C = mk$C, min_af = min(mk$af), max_af = max(mk$af))
}

ans <- rbind(check(1.5), check(2), check(6))
print(ans, digits = 8)
stopifnot(max(ans[, "follower_foc"]) < 2e-8,
          max(ans[, "reaction_error"]) < 2e-5,
          max(ans[, "leader_gradient"]) < 2e-5,
          max(ans[, "inverse_error"]) < 2e-7,
          max(ans[, "roundtrip_error"]) < 2e-12)
cat("PASS: independent CES Cournot rank-one/reaction check.\n")
