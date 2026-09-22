#!/usr/bin/env Rscript
# Independent Gaussian oracle for the market-intercept collapse extension.
# All dense covariance calculations here are small-fixture references only.
options(stringsAsFactors = FALSE)
args <- commandArgs(FALSE)
script <- sub("^--file=", "", args[startsWith(args, "--file=")])
outdir <- if (length(script)) dirname(normalizePath(script[1])) else "doc"

log_gaussian <- function(x, V) {
  if (!length(x)) return(0)
  U <- chol(V)
  u <- forwardsolve(t(U), x)
  -.5 * (length(x) * log(2*pi) + 2 * sum(log(diag(U))) + sum(u*u))
}

log_sum_exp <- function(x) {
  xmax <- max(x)
  xmax + log(sum(exp(x - xmax)))
}

contrast_basis <- function(G) {
  if (G <= 1L) return(matrix(numeric(), nrow = G, ncol = 0L))
  P <- diag(G) - matrix(1/G, G, G)
  ee <- eigen(P, symmetric = TRUE)
  ee$vectors[, ee$values > .5, drop = FALSE]
}

incidence <- function(group) {
  G <- max(group)
  Z <- matrix(0, length(group), G)
  Z[cbind(seq_along(group), group)] <- 1
  Z
}

# Structured normalized density for one market, after margins are
# conditionally removed. z has the prior center mu_b removed; v can vary by
# observation, as it does when each market has its own residual scales.
collapse_market <- function(z, v, group, tau, market_sd) {
  G <- max(group)
  if (G == 1L) tau <- 0
  a <- cstat <- numeric(G)
  for (i in seq_along(z)) {
    a[group[i]] <- a[group[i]] + 1/v[i]
    cstat[group[i]] <- cstat[group[i]] + z[i]/v[i]
  }
  group_mean <- cstat/a
  within <- sum((z - group_mean[group])^2/v)
  h <- 1/(1 + tau^2*a)
  ah <- a*h
  A <- sum(ah)
  cbar <- sum(ah*group_mean)/A
  d <- mean(h)
  den <- d + market_sd^2*A
  logdet <- sum(log(v)) + sum(log1p(tau^2*a)) + log(den)
  quadratic <- within + sum(ah*(group_mean-cbar)^2) + A*cbar^2/den
  list(log_density = -.5*(length(z)*log(2*pi) + logdet + quadratic),
       intercept_mean_centered = market_sd^2*A*cbar/den,
       intercept_variance = market_sd^2*d/den,
       sufficient = list(a = a, c = cstat, h = h, A = A,
                         cbar = cbar, d = d))
}

condition_margins <- function(y_share, y_margin, eta_share, eta_margin,
                              paired, share_sd, margin_sd, rho, mu_b) {
  N <- length(y_share)
  share_sd <- rep_len(share_sd, N)
  margin_sd <- rep_len(margin_sd, length(paired))
  r_share <- y_share - eta_share - mu_b
  r_margin <- y_margin - eta_margin
  z <- r_share
  v <- share_sd^2
  if (length(paired)) {
    z[paired] <- z[paired] - rho*share_sd[paired]/margin_sd*r_margin
    v[paired] <- share_sd[paired]^2*(1-rho^2)
  }
  list(z = z, v = v, margin_log_density =
         sum(dnorm(r_margin, mean = 0, sd = margin_sd, log = TRUE)))
}

structured_target <- function(y_share, y_margin, eta_share, eta_margin,
                              group, paired, tau, share_sd, margin_sd,
                              rho, market_sd, mu_b) {
  adjusted <- condition_margins(y_share, y_margin, eta_share, eta_margin,
                                paired, share_sd, margin_sd, rho, mu_b)
  ans <- collapse_market(adjusted$z, adjusted$v, group, tau, market_sd)
  ans$log_density <- ans$log_density + adjusted$margin_log_density
  ans$adjusted <- adjusted
  ans
}

# Direct joint Gaussian density of observed shares and only observed margins.
# The full normalized constants are included by log_gaussian().
dense_observed_covariance <- function(group, paired, tau, share_sd,
                                      margin_sd, rho, market_sd) {
  N <- length(group)
  Q <- length(paired)
  G <- max(group)
  Z <- incidence(group)
  P <- diag(G) - matrix(1/G, G, G)
  share_sd <- rep_len(share_sd, N)
  margin_sd <- rep_len(margin_sd, Q)
  cross <- matrix(0, N, Q)
  if (Q) for (j in seq_len(Q))
    cross[paired[j], j] <- rho*share_sd[paired[j]]*margin_sd[j]
  share_cov <- tau^2*Z%*%P%*%t(Z) + diag(share_sd^2, N) +
    matrix(market_sd^2, N, N)
  rbind(cbind(share_cov, cross),
        cbind(t(cross), diag(margin_sd^2, Q)))
}

dense_target <- function(y_share, y_margin, eta_share, eta_margin,
                         group, paired, tau, share_sd, margin_sd, rho,
                         market_sd, mu_b) {
  residual <- c(y_share - eta_share - mu_b,
                y_margin - eta_margin)
  log_gaussian(residual, dense_observed_covariance(
    group, paired, tau, share_sd, margin_sd, rho, market_sd))
}

posterior_map <- function(G, U) {
  D <- 1L + ncol(U)
  M <- matrix(0, G + 1L, D)
  M[1L, 1L] <- 1
  if (ncol(U)) M[2:(G+1L), 2:D] <- U
  M
}

# Posterior from sufficient Gaussian linear-model precision. Output order is
# (b, q_1, ..., q_G), with b on its original nonzero prior-mean scale.
joint_posterior_structured <- function(z, v, group, tau, market_sd, mu_b) {
  G <- max(group)
  U <- contrast_basis(G)
  Z <- incidence(group)
  H <- cbind(1, Z %*% U)
  D <- 1L + ncol(U)
  prior_var <- c(market_sd^2, rep(tau^2, ncol(U)))
  prior_precision <- diag(1/prior_var, D)
  weighted_H <- H/v
  precision <- prior_precision + crossprod(H, weighted_H)
  covariance_x <- solve(precision)
  mean_x <- covariance_x %*% crossprod(H, z/v)
  M <- posterior_map(G, U)
  mean_out <- as.vector(M %*% mean_x)
  mean_out[1L] <- mean_out[1L] + mu_b
  list(mean = mean_out,
       covariance = M %*% covariance_x %*% t(M),
       U = U, H = H)
}

# Independent dense conditioning identity using Cov(x,z) and Cov(z).
joint_posterior_dense <- function(z, v, group, tau, market_sd, mu_b) {
  G <- max(group)
  U <- contrast_basis(G)
  H <- cbind(1, incidence(group) %*% U)
  D <- 1L + ncol(U)
  prior <- diag(c(market_sd^2, rep(tau^2, ncol(U))), D)
  cross <- prior %*% t(H)
  z_cov <- H %*% prior %*% t(H) + diag(v)
  mean_x <- cross %*% solve(z_cov, z)
  covariance_x <- prior - cross %*% solve(z_cov, t(cross))
  M <- posterior_map(G, U)
  mean_out <- as.vector(M %*% mean_x)
  mean_out[1L] <- mean_out[1L] + mu_b
  list(mean = mean_out,
       covariance = M %*% covariance_x %*% t(M),
       U = U, H = H)
}

# Direct conditional q | b, data moments. The draw algorithm first draws b
# from its marginal, then uses these moments for q, preserving their coupling.
quality_given_intercept <- function(z, v, group, tau, intercept_centered) {
  G <- max(group)
  U <- contrast_basis(G)
  if (!ncol(U)) return(list(mean = rep(0, G), covariance = matrix(0, G, G)))
  X <- incidence(group) %*% U
  W_X <- X/v
  covariance_theta <- solve(diag(1/tau^2, ncol(U)) + crossprod(X, W_X))
  mean_theta <- covariance_theta %*% crossprod(X, (z-intercept_centered)/v)
  list(mean = as.vector(U %*% mean_theta),
       covariance = U %*% covariance_theta %*% t(U))
}

sequential_draws <- function(z, v, group, tau, market_sd, mu_b,
                             b_mean_centered, b_variance, n_draw) {
  G <- max(group)
  U <- contrast_basis(G)
  b0 <- rnorm(n_draw, b_mean_centered, sqrt(b_variance))
  if (!ncol(U)) return(cbind(b = b0 + mu_b, matrix(0, n_draw, G)))
  X <- incidence(group) %*% U
  covariance_theta <- solve(diag(1/tau^2, ncol(U)) + crossprod(X, X/v))
  base <- as.vector(covariance_theta %*% crossprod(X, z/v))
  intercept_slope <- as.vector(covariance_theta %*% crossprod(X, rep(1, length(z))/v))
  theta_mean <- outer(base, rep(1, n_draw)) -
    outer(intercept_slope, b0)
  ee <- eigen(covariance_theta, symmetric = TRUE)
  root <- ee$vectors %*% diag(sqrt(pmax(ee$values, 0)), ncol(U))
  theta <- theta_mean + root %*% matrix(rnorm(ncol(U)*n_draw), ncol(U), n_draw)
  q <- t(U %*% theta)
  cbind(b = b0 + mu_b, q)
}

scenario <- function(name, group, year, paired, y_share, y_margin,
                     eta_share, eta_margin, share_sd, margin_sd, rho,
                     tau, market_sd, mu_b) {
  list(name = name, group = as.integer(group), year = as.integer(year),
       paired = as.integer(paired), y_share = y_share, y_margin = y_margin,
       eta_share = eta_share, eta_margin = eta_margin,
       share_sd = share_sd, margin_sd = margin_sd, rho = rho,
       tau = tau, market_sd = market_sd, mu_b = mu_b)
}

set.seed(20260920L)
scenarios <- list(
  scenario("unbalanced_three_year_partial_margins",
    group = c(1,1,2,1,2,1,2,3,3), year = c(1,1,1,2,2,3,3,3,3),
    paired = c(1,4,6,8), y_share = rnorm(9, .1, .3),
    y_margin = rnorm(4, -.05, .25), eta_share = seq(-.10,.12,length.out=9),
    eta_margin = c(-.04,.03,.07,-.02), share_sd = .23, margin_sd = .41,
    rho = .73, tau = .48, market_sd = .36, mu_b = .27),
  scenario("all_margins_missing",
    group = c(1,1,2,3,1,2), year = c(1,1,2,2,3,3), paired = integer(),
    y_share = rnorm(6, -.1, .35), y_margin = numeric(),
    eta_share = c(.02,-.08,.04,.11,-.03,.07), eta_margin = numeric(),
    share_sd = .31, margin_sd = numeric(), rho = -.81,
    tau = .62, market_sd = .29, mu_b = -.34),
  scenario("single_product_repeated_years",
    group = rep(1, 6), year = c(1,1,2,2,3,3), paired = c(2,5),
    y_share = rnorm(6, .05, .25), y_margin = rnorm(2, .03, .2),
    eta_share = seq(-.06,.09,length.out=6), eta_margin = c(.01,-.05),
    share_sd = .19, margin_sd = .33, rho = -.68,
    tau = .51, market_sd = .43, mu_b = .31)
)

density_rows <- list()
joint_rows <- list()
for (fx in scenarios) {
  st <- structured_target(fx$y_share, fx$y_margin, fx$eta_share,
    fx$eta_margin, fx$group, fx$paired, fx$tau, fx$share_sd,
    fx$margin_sd, fx$rho, fx$market_sd, fx$mu_b)
  dn <- dense_target(fx$y_share, fx$y_margin, fx$eta_share,
    fx$eta_margin, fx$group, fx$paired, fx$tau, fx$share_sd,
    fx$margin_sd, fx$rho, fx$market_sd, fx$mu_b)
  err <- abs(st$log_density - dn)
  density_rows[[length(density_rows)+1L]] <- data.frame(
    fixture = fx$name, years = length(unique(fx$year)),
    observations = length(fx$group), groups = max(fx$group),
    observed_margins = length(fx$paired), missing_margins =
      length(fx$group)-length(fx$paired), mu_b = fx$mu_b,
    share_sd = fx$share_sd,
    margin_sd = if (length(fx$margin_sd)) fx$margin_sd[1L] else NA_real_,
    rho = fx$rho,
    structured_log_density = st$log_density, dense_log_density = dn,
    absolute_error = err, tolerance = 1e-7, passed = err < 1e-7)

  # The all-missing fixture still has a valid share posterior; the one-group
  # fixture is handled by its zero-dimensional centered-quality subspace.
  jp <- joint_posterior_structured(st$adjusted$z, st$adjusted$v,
    fx$group, fx$tau, fx$market_sd, fx$mu_b)
  jd <- joint_posterior_dense(st$adjusted$z, st$adjusted$v,
    fx$group, fx$tau, fx$market_sd, fx$mu_b)
  collapsed_b_mean <- fx$mu_b + st$intercept_mean_centered
  collapsed_b_var <- st$intercept_variance
  q_direct <- quality_given_intercept(st$adjusted$z, st$adjusted$v,
    fx$group, fx$tau, st$intercept_mean_centered)
  bq_error <- if (max(fx$group) > 1L)
    max(abs(jp$covariance[1L, -1L] - jd$covariance[1L, -1L])) else 0
  qcond_mean <- quality_given_intercept(st$adjusted$z, st$adjusted$v,
    fx$group, fx$tau, st$intercept_mean_centered)
  # Conditional Gaussian identity evaluated at the posterior mean of b.
  qcond_from_joint <- if (max(fx$group) > 1L) {
    cov_b <- jd$covariance[1L, 1L]
    jd$mean[-1L] + jd$covariance[-1L, 1L]/cov_b *
      (collapsed_b_mean-jd$mean[1L])
  } else rep(0, max(fx$group))
  qcond_cov_from_joint <- if (max(fx$group) > 1L) {
    jd$covariance[-1L, -1L, drop=FALSE] -
      tcrossprod(jd$covariance[-1L, 1L])/jd$covariance[1L,1L]
  } else matrix(0, max(fx$group), max(fx$group))
  qcond_err <- max(abs(qcond_mean$mean-qcond_from_joint),
                   abs(qcond_mean$covariance-qcond_cov_from_joint))

  draws <- sequential_draws(st$adjusted$z, st$adjusted$v, fx$group,
    fx$tau, fx$market_sd, fx$mu_b, st$intercept_mean_centered,
    st$intercept_variance, n_draw = 30000L)
  draw_mean <- colMeans(draws)
  draw_cov <- cov(draws)
  mean_se <- sqrt(pmax(diag(jd$covariance), 1e-30)/nrow(draws))
  cov_se <- sqrt((outer(diag(jd$covariance), diag(jd$covariance)) +
    jd$covariance^2)/(nrow(draws)-1))
  max_mean_z <- max(abs(draw_mean-jd$mean)/mean_se)
  max_cov_z <- max(abs(draw_cov-jd$covariance)/pmax(cov_se,1e-30))
  joint_rows[[length(joint_rows)+1L]] <- data.frame(
    fixture = fx$name, max_structured_dense_mean_error = max(abs(jp$mean-jd$mean)),
    max_structured_dense_covariance_error = max(abs(jp$covariance-jd$covariance)),
    collapsed_vs_joint_b_mean_error = abs(collapsed_b_mean-jd$mean[1L]),
    collapsed_vs_joint_b_variance_error = abs(collapsed_b_var-jd$covariance[1L,1L]),
    max_bq_cross_covariance_error = bq_error,
    q_given_b_formula_error = qcond_err,
    sequential_draw_max_mean_mc_z = max_mean_z,
    sequential_draw_max_covariance_mc_z = max_cov_z,
    monte_carlo_z_tolerance = 6,
    passed = max(abs(jp$mean-jd$mean), abs(jp$covariance-jd$covariance),
      abs(collapsed_b_mean-jd$mean[1L]), abs(collapsed_b_var-jd$covariance[1L,1L]),
      bq_error, qcond_err) < 1e-8 && max_mean_z < 6 && max_cov_z < 6)
}

# Exact-k oracle with three markets, three states each, market-specific scales,
# repeated years, unequal group counts, and one market with every margin absent.
market_data <- list(
  list(group=c(1,1,2,1), year=c(1,1,2,2), paired=c(1,4),
       ys=c(.18,-.08,.24,.31), ym=c(-.13,.17), ss=.22, sm=.29),
  list(group=c(1,2,2,3,2), year=c(1,1,2,2,3), paired=integer(),
       ys=c(-.19,.11,.06,.28,-.03), ym=numeric(), ss=.37, sm=numeric()),
  list(group=c(1,2,1,3,3,2), year=c(1,1,2,2,3,3), paired=c(2,3,6),
       ys=c(.03,.21,-.14,.33,.12,-.05), ym=c(.04,-.22,.16), ss=.18, sm=.53)
)
K <- 3L
log_prior <- rbind(log(c(.52,.31,.17)), log(c(.22,.50,.28)),
                   log(c(.39,.26,.35)))
tau_mix <- .57
sb_mix <- .34
mu_mix <- .26
rho_mix <- -.64
state_scores <- matrix(NA_real_, length(market_data), K)
for (m in seq_along(market_data)) {
  md <- market_data[[m]]
  N <- length(md$group)
  for (k in seq_len(K)) {
    eta_s <- seq(-.07,.08,length.out=N) + (k-1L)*seq(.018,.043,length.out=N)
    eta_m <- if (length(md$paired)) seq(-.03,.05,length.out=length(md$paired)) + .06*(k-1L) else numeric()
    fit <- structured_target(md$ys, md$ym, eta_s, eta_m, md$group,
      md$paired, tau_mix, md$ss, md$sm, rho_mix, sb_mix, mu_mix)
    state_scores[m,k] <- log_prior[m,k] + fit$log_density
  }
}
market_log_marginal <- apply(state_scores, 1L, log_sum_exp)
factorized_log_evidence <- sum(market_log_marginal)
assignments <- expand.grid(rep(list(seq_len(K)), length(market_data)))
joint_scores <- apply(assignments, 1L, function(kvec)
  sum(vapply(seq_along(kvec), function(m) state_scores[m,kvec[m]], numeric(1))))
enumerated_log_evidence <- log_sum_exp(joint_scores)
enumerated_probs <- exp(joint_scores-enumerated_log_evidence)
factor_marginal_probs <- lapply(seq_along(market_data), function(m)
  exp(state_scores[m,]-market_log_marginal[m]))
enum_marginal_probs <- lapply(seq_along(market_data), function(m) {
  vapply(seq_len(K), function(k)
    sum(enumerated_probs[assignments[[m]] == k]), numeric(1))
})
prob_error <- max(unlist(Map(function(a,b) max(abs(a-b)),
                             factor_marginal_probs, enum_marginal_probs)))
factorization_row <- data.frame(
  markets = length(market_data), states_per_market = K,
  enumerated_joint_states = nrow(assignments),
  factorized_log_evidence = factorized_log_evidence,
  enumerated_log_evidence = enumerated_log_evidence,
  absolute_log_evidence_error = abs(factorized_log_evidence-enumerated_log_evidence),
  max_market_state_probability_error = prob_error,
  normalized_state_densities = TRUE, tolerance = 1e-10,
  passed = abs(factorized_log_evidence-enumerated_log_evidence) < 1e-10 &&
    prob_error < 1e-10)

tables <- list(
  density = do.call(rbind, density_rows),
  joint = do.call(rbind, joint_rows),
  exact_k = factorization_row
)
for (name in names(tables))
  write.csv(tables[[name]], file.path(outdir,
    paste0("market_intercept_collapse_extended_", name, "_checks.csv")),
    row.names = FALSE)

for (name in names(tables)) {
  x <- tables[[name]]
  cat(name, sum(x$passed), "/", nrow(x), "passed\n")
}
cat("maximum dense likelihood error:", max(tables$density$absolute_error), "\n")
cat("maximum joint posterior mean error:",
    max(tables$joint$max_structured_dense_mean_error), "\n")
cat("maximum joint posterior covariance error:",
    max(tables$joint$max_structured_dense_covariance_error), "\n")
cat("maximum q|b formula error:", max(tables$joint$q_given_b_formula_error), "\n")
cat("maximum sequential draw mean MC z:",
    max(tables$joint$sequential_draw_max_mean_mc_z), "\n")
cat("maximum sequential draw covariance MC z:",
    max(tables$joint$sequential_draw_max_covariance_mc_z), "\n")
cat("exact-k enumeration states:", nrow(assignments), "\n")
cat("exact-k log evidence error:", factorization_row$absolute_log_evidence_error, "\n")
cat("exact-k state probability error:", prob_error, "\n")
stopifnot(all(vapply(tables, function(x) all(x$passed), logical(1))))
