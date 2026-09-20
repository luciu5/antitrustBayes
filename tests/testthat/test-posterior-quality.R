posterior_quality_fixture <- function(group, log_share, labels = NULL,
                                      quality = "collapsed", share = "gaussian",
                                      likelihood_mode = 2L, n_draws = NULL) {
  group <- as.integer(group)
  n_group <- max(group)
  if (is.null(labels)) labels <- paste0("g", seq_len(n_group))
  stan_data <- list(
    N = length(group), N_group = n_group, group = group,
    log_shareIn = as.numeric(log_share), likelihood_mode = likelihood_mode
  )
  attr(stan_data, "group_levels") <- labels
  attr(stan_data, "quality_parameterization") <- quality
  attr(stan_data, "share_likelihood") <- share
  metadata <- list(quality_parameterization = quality,
                   share_likelihood = share)
  if (!is.null(n_draws)) metadata$posterior_draws <- n_draws
  methods::new(
    "BayesLogit", stanfit = NULL, data = data.frame(group = labels[group]),
    stan_data = stan_data, price_sd = 1, supply_model = 1L,
    demand_model = "logit", market_side = "input", call = quote(fake()),
    metadata = metadata
  )
}

posterior_quality_draws_fixture <- function(share_mean, tau, sigma) {
  list(share_mean = share_mean, sigma_b_strat = tau, sigma_share = sigma)
}

posterior_quality_dense_moments <- function(group, y, share_mean, tau, sigma) {
  g <- max(group)
  z <- matrix(0, nrow = length(group), ncol = g)
  z[cbind(seq_along(group), group)] <- 1
  q <- matrix(0, nrow = g, ncol = max(0L, g - 1L))
  if (g > 1L) {
    for (k in seq_len(g - 1L)) {
      q[seq_len(k), k] <- 1 / sqrt(k * (k + 1))
      q[k + 1L, k] <- -k / sqrt(k * (k + 1))
    }
  }
  h <- z %*% q
  if (g == 1L) return(list(mean = 0, covariance = matrix(0, 1, 1)))
  v <- solve(diag(1 / tau^2, g - 1L) + crossprod(h) / sigma^2)
  m <- as.numeric(v %*% crossprod(h, y - share_mean) / sigma^2)
  list(mean = as.numeric(q %*% m), covariance = q %*% v %*% t(q))
}

test_that("posterior_quality matches dense expanded Gaussian conditioning", {
  group <- c(1L, 2L, 1L, 3L, 2L, 3L)
  y <- c(.21, .45, .36, .12, .49, .2)
  mu <- c(.1, .25, .2, .15, .23, .18)
  tau <- .73
  sigma <- .41
  fit <- posterior_quality_fixture(group, y, c("firm_a", "firm_b", "firm_c"))
  got <- posterior_quality(
    fit, seed = 19,
    draws = posterior_quality_draws_fixture(matrix(mu, nrow = 1L), tau, sigma)
  )
  oracle <- posterior_quality_dense_moments(group, y, mu, tau, sigma)
  expect_equal(as.numeric(got$conditional_means[1L, ]), oracle$mean,
               tolerance = 1e-12)
  expect_equal(as.numeric(got$conditional_sds[1L, ]),
               sqrt(diag(oracle$covariance)), tolerance = 1e-12)
  expect_equal(sum(got$draws[1L, ]), 0, tolerance = 1e-14)
  expect_equal(colnames(got$draws), c("firm_a", "firm_b", "firm_c"))
})

test_that("posterior_quality draws have the dense conditional covariance", {
  group <- c(1L, 2L, 1L, 3L, 2L, 3L)
  y <- c(.21, .45, .36, .12, .49, .2)
  mu <- c(.1, .25, .2, .15, .23, .18)
  tau <- .73
  sigma <- .41
  n_mc <- 3000L
  fit <- posterior_quality_fixture(group, y, c("firm_a", "firm_b", "firm_c"))
  got <- posterior_quality(
    fit, seed = 29,
    draws = posterior_quality_draws_fixture(
      matrix(rep(mu, n_mc), nrow = n_mc, byrow = TRUE),
      rep(tau, n_mc), rep(sigma, n_mc)
    )
  )
  oracle <- posterior_quality_dense_moments(group, y, mu, tau, sigma)
  expect_true(max(abs(unname(cov(got$draws)) - unname(oracle$covariance))) < .01)
})

test_that("posterior_quality handles one global group as zero quality", {
  fit <- posterior_quality_fixture(rep(1L, 4L), c(.1, .2, .3, .4), "only")
  got <- posterior_quality(
    fit, draws = posterior_quality_draws_fixture(
      matrix(c(.0, .1, .2, .3), nrow = 1L), .8, .4
    ), seed = 2
  )
  expect_equal(got$draws, matrix(0, 1L, 1L,
                                  dimnames = list("1", "only")))
  expect_equal(got$conditional_means, got$draws)
  expect_equal(got$conditional_sds, got$draws)
  expect_equal(got$summary$posterior_mean, 0)
  expect_equal(got$summary$posterior_sd, 0)
})

test_that("posterior_quality uses repeated multiproduct panel groups", {
  group <- c(1L, 1L, 2L, 2L, 1L, 2L)
  y <- c(.4, .5, .1, .2, .45, .15)
  mu <- c(.2, .3, .12, .18, .25, .17)
  fit <- posterior_quality_fixture(group, y, c("A", "B"))
  dr <- posterior_quality_draws_fixture(rbind(mu, mu + .03), c(.6, .6), c(.5, .5))
  got <- posterior_quality(fit, draws = dr, seed = 5)
  expect_identical(got$group_labels, c("A", "B"))
  expect_identical(got$draw_indices, 1:2)
  expect_equal(nrow(got$draws_long), 4L)
  expect_equal(got$draws_long$draw_index, c(1L, 1L, 2L, 2L))
  expect_equal(got$draws_long$group, c("A", "B", "A", "B"))
  expect_equal(unname(rowSums(got$draws)), c(0, 0), tolerance = 1e-14)
})

test_that("posterior_quality seed is deterministic and preserves RNG state", {
  fit <- posterior_quality_fixture(c(1L, 2L, 1L, 3L), c(.2, .3, .1, .4))
  dr <- posterior_quality_draws_fixture(
    rbind(c(.1, .2, .3, .4), c(.2, .3, .4, .5), c(.3, .4, .5, .6)),
    c(.5, .6, .7), c(.4, .4, .4)
  )
  set.seed(71)
  before <- .Random.seed
  a <- posterior_quality(fit, n_draws = 2L, seed = 123, draws = dr)
  expect_identical(.Random.seed, before)
  b <- posterior_quality(fit, n_draws = 2L, seed = 123, draws = dr)
  expect_identical(a$draw_indices, b$draw_indices)
  expect_identical(a$draws, b$draws)
  expect_identical(a$conditional_means, b$conditional_means)
  expect_identical(a$conditional_sds, b$conditional_sds)
})

test_that("posterior_quality rejects unsupported modes and malformed draws", {
  dr <- posterior_quality_draws_fixture(matrix(c(.1, .2), nrow = 1L), .5, .4)
  sampled <- posterior_quality_fixture(c(1L, 2L), c(.2, .3), quality = "sampled")
  expect_error(posterior_quality(sampled, draws = dr), "sampled")
  student <- posterior_quality_fixture(c(1L, 2L), c(.2, .3), share = "student_t")
  expect_error(posterior_quality(student, draws = dr), "non-Gaussian")
  fit <- posterior_quality_fixture(c(1L, 2L), c(.2, .3))
  expect_error(posterior_quality(fit, draws = list()), "missing field")
  expect_error(
    posterior_quality(fit, draws = posterior_quality_draws_fixture(
      matrix(c(.1, .2, .3), nrow = 1L), .5, .4
    )),
    "draws by observations"
  )
  expect_error(
    posterior_quality(fit, draws = posterior_quality_draws_fixture(dr$share_mean,
                                                                    c(.5, .6), .4)),
    "expected 1"
  )
  dr$share_mean[1L, 1L] <- NA_real_
  expect_error(posterior_quality(fit, draws = dr), "non-finite")
})

test_that("Gaussian rho quality recovery matches original joint conditioning", {
  group <- c(1L, 1L, 2L, 3L, 2L, 3L)
  y <- c(.4, .2, -.1, .6, .3, .1)
  mu <- c(.2, .1, -.2, .4, .2, .2)
  idx <- c(1L, 3L, 6L)
  fit <- posterior_quality_fixture(group, y, likelihood_mode = 1L)
  fit@stan_data$use_rho <- 1L
  fit@stan_data$margin_idx <- idx
  fit@stan_data$log_marginNormInv <- c(.8, NA, 1.1, NA, NA, .9)
  tau <- .6; sigma_s <- .4; sigma_m <- .35
  z <- matrix(0, 6L, 3L); z[cbind(seq_along(group), group)] <- 1
  contrasts <- qr.Q(qr(cbind(rep(1, 3L), diag(3L))))[, -1L, drop = FALSE]
  h <- rbind(z %*% contrasts, matrix(0, length(idx), 2L))
  r <- matrix(0, 6L, length(idx)); r[cbind(idx, seq_along(idx))] <- 1
  prediction <- rep(2, 6L)
  residual <- c(y - mu, fit@stan_data$log_marginNormInv[idx] - log(prediction[idx] + 1e-9))
  for (rho in c(0, .3, -.4)) {
    covariance <- rbind(cbind(diag(sigma_s^2, 6L), rho * sigma_s * sigma_m * r),
      cbind(rho * sigma_s * sigma_m * t(r), diag(sigma_m^2, length(idx))))
    precision <- diag(1 / tau^2, 2L) + crossprod(h, solve(covariance, h))
    beta_covariance <- solve(precision)
    expected_mean <- contrasts %*% beta_covariance %*% crossprod(h, solve(covariance, residual))
    expected_covariance <- contrasts %*% beta_covariance %*% t(contrasts)
    dr <- posterior_quality_draws_fixture(matrix(mu, 1L), tau, sigma_s)
    dr$rho <- rho; dr$sigma_margin <- sigma_m
    dr$pred_marginNormInv <- matrix(prediction, 1L)
    got <- posterior_quality(fit, draws = dr, seed = 9)
    expect_equal(as.numeric(got$conditional_means), as.numeric(expected_mean), tolerance = 1e-12)
    expect_equal(as.numeric(got$conditional_sds), sqrt(diag(expected_covariance)), tolerance = 1e-12)
    expect_equal(sum(got$draws), 0, tolerance = 1e-14)
    expect_match(got$metadata$conditioning, "matched Gaussian margins")
  }
  fit@metadata$partition_prior <- "ordered_k"
  expect_error(posterior_quality(fit, draws = dr), "component-specific margin means")
  share_only <- posterior_quality(fit, draws = dr, conditioning = "share", seed = 9)
  expected <- posterior_quality_dense_moments(group, y, mu, tau, sigma_s)
  expect_equal(as.numeric(share_only$conditional_means), expected$mean, tolerance = 1e-12)
  expect_identical(share_only$metadata$conditioning, "observed log shares only")
})

test_that("scalar rstan arrays are accepted without changing draw alignment", {
  fit <- posterior_quality_fixture(c(1L, 2L), c(.2, .3))
  plain <- posterior_quality_draws_fixture(rbind(c(.1, .2), c(.2, .3)), c(.5, .6), c(.4, .4))
  arrays <- plain
  arrays$sigma_b_strat <- array(plain$sigma_b_strat, dim = 2L)
  arrays$sigma_share <- array(plain$sigma_share, dim = 2L)
  expect_identical(posterior_quality(fit, draws = arrays, seed = 31)$draws,
    posterior_quality(fit, draws = plain, seed = 31)$draws)
})

test_that("posterior_quality uses economic quality scale under total/fraction coordinates", {
  group <- c(1L, 2L, 1L, 3L)
  y <- c(.2, .3, .1, .4)
  mu <- matrix(c(.1, .2, .15, .3), nrow = 1L)
  tau <- .6
  sigma <- .4
  total <- sqrt(tau^2 + sigma^2)
  fraction_logit <- 2 * (log(sigma) - log(tau))
  fit <- posterior_quality_fixture(group, y)
  fit@metadata$quality_scale_parameterization <- "total_fraction"
  fit@metadata$quality_scale_parameterization_effective <- "total_fraction"
  fit@stan_data$quality_scale_parameterization <- 2L

  contract <- antitrustBayes:::.posterior_contract(
    fit, purpose = "posterior_quality", conditioning = "share"
  )
  expect_identical(
    contract$actual[contract$canonical == "sigma_b_strat"],
    "sigma_b_strat_economic"
  )
  expect_error(antitrustBayes:::.posterior_normalize_draw_names(list(
    share_mean = mu, sigma_b_strat = total, sigma_share = sigma
  ), contract), "sigma_b_strat")

  legacy <- posterior_quality(fit, seed = 47, draws = list(
    share_mean = mu, sigma_b_strat_economic = tau,
    ## Deliberately wrong raw value confirms that the economic alias wins.
    sigma_b_strat = total * 3, log_sigma_share = fraction_logit,
    sigma_share = sigma
  ))
  reconstructed <- posterior_quality(fit, seed = 47, draws = list(
    share_mean = mu, sigma_b_strat = total,
    log_sigma_share = fraction_logit, sigma_share = sigma
  ))
  expect_identical(reconstructed$draws, legacy$draws)
  expect_identical(reconstructed$conditional_means, legacy$conditional_means)
  expect_identical(reconstructed$conditional_sds, legacy$conditional_sds)

  expect_error(posterior_quality(fit, draws = list(
    share_mean = mu, sigma_b_strat = total, sigma_share = sigma
  )), "raw pair")
})
