relative_math_fixture <- function(markets = 1L, periods = 1L,
                                  ragged = FALSE, global = FALSE) {
  d <- expand.grid(product = 1:4, period = seq_len(periods),
                   market = seq_len(markets))
  if (ragged) {
    keep <- !(d$product == 4L & d$period %% 2L == 0L)
    keep <- keep & !(d$market == 1L & d$period == 4L & d$product != 2L)
    d <- d[keep, , drop = FALSE]
  }
  d <- d[rev(seq_len(nrow(d))), , drop = FALSE]
  group <- if (global) d$product else (d$market - 1L) * 4L + d$product
  blocks <- if (global) rep(1L, 5L) else c(rep(seq_len(markets), each = 4L), markets)
  list(residual = sin(seq_len(nrow(d)) / 3) + 2 * d$market - d$period,
       variance = (.3 + .07 * d$market)^2,
       quality_index = group, quality_block = blocks,
       cell = paste(d$market, d$period, sep = ":"), tau = .8)
}

relative_math_oracle <- function(residual, variance, quality_index,
                                 quality_block, cell, tau, rotation = NULL) {
  n <- length(residual); g <- length(quality_block)
  ee <- outer(cell, unique(cell), "==") * 1
  nc <- ncol(ee)
  q <- if (n == nc) matrix(0, n, 0L) else
    qr.Q(qr(ee), complete = TRUE)[, seq.int(nc + 1L, n), drop = FALSE]
  if (!is.null(rotation)) q <- q %*% rotation
  z <- diag(g)[quality_index, , drop = FALSE]
  p <- diag(g)
  for (ii in split(seq_len(g), quality_block)) p[ii, ii] <- p[ii, ii] - 1 / length(ii)
  covariance <- tau^2 * p
  if (!ncol(q)) return(list(log_density = 0, mean = rep(0, g),
                            covariance = covariance, precision = matrix(0, n, n)))
  v <- crossprod(q, (diag(rep(variance, length.out = n)) + z %*% covariance %*% t(z)) %*% q)
  ch <- chol(v)
  observed <- as.numeric(crossprod(q, residual))
  whitened <- forwardsolve(t(ch), observed)
  log_density <- -.5 * (length(observed) * log(2 * pi) +
                         2 * sum(log(diag(ch))) + sum(whitened^2))
  gain <- covariance %*% t(z) %*% q
  list(log_density = log_density,
       mean = as.numeric(gain %*% solve(v, observed)),
       covariance = covariance - gain %*% solve(v, t(gain)),
       precision = q %*% solve(v, t(q)))
}

test_that("relative likelihood and quality moments match expanded panel Gaussians", {
  for (markets in c(1L, 4L)) for (periods in c(1L, 4L)) {
    for (global in c(FALSE, TRUE)) for (rho in c(0, .63, -.95)) {
      a <- relative_math_fixture(markets, periods, ragged = TRUE, global = global)
      paired <- seq_along(a$residual) %% 3L != 0L
      a$residual[paired] <- a$residual[paired] - rho * .4 / .7 * cos(which(paired))
      a$variance[paired] <- a$variance[paired] * (1 - rho^2)
      for (tau in c(0, .02, .8, 3)) {
        a$tau <- tau
        expected <- do.call(relative_math_oracle, a)
        got <- do.call(.relative_share_gaussian, c(a, list(return_moments = TRUE)))
        expect_equal(got$log_density, expected$log_density, tolerance = 2e-10)
        expect_equal(got$conditional_mean, expected$mean, tolerance = 2e-10)
        expect_equal(got$conditional_variance, diag(expected$covariance), tolerance = 2e-10)
        expect_equal(got$dimension, length(a$residual) - length(unique(a$cell)))
      }
    }
  }
})

test_that("cell levels carry no relative-share information", {
  a <- relative_math_fixture(4L, 4L, ragged = TRUE)
  got <- do.call(.relative_share_gaussian, c(a, list(return_moments = TRUE)))
  shifted <- a
  shifted$residual <- a$residual + 100 * match(a$cell, unique(a$cell))
  other <- do.call(.relative_share_gaussian, c(shifted, list(return_moments = TRUE)))
  expect_equal(other, got, tolerance = 2e-10)
  expected <- do.call(relative_math_oracle, a)
  x <- cbind(a$residual, model.matrix(~ factor(a$cell) - 1))
  args <- a[names(a) != "residual"]
  precision <- do.call(.relative_share_precision, c(list(x = x), args))
  expect_equal(precision, unname(expected$precision %*% x), tolerance = 1e-10)
  expect_equal(unname(precision[, -1L]), matrix(0, nrow(x), ncol(x) - 1L),
               tolerance = 1e-12)
  expect_equal(do.call(.relative_share_precision, c(list(x = a$residual), args)),
               as.numeric(precision[, 1L]), tolerance = 1e-12)
})

test_that("normalization and inference do not depend on contrast orientation", {
  a <- relative_math_fixture(1L, 4L, ragged = TRUE)
  q <- .relative_share_basis(a$cell)
  expect_equal(crossprod(q), diag(ncol(q)), tolerance = 1e-13)
  e <- model.matrix(~ factor(a$cell) - 1)
  expect_equal(unname(crossprod(q, e)), matrix(0, ncol(q), ncol(e)), tolerance = 1e-13)
  rotate <- qr.Q(qr(matrix(sin(seq_len(ncol(q)^2)), ncol(q))))
  expected <- do.call(relative_math_oracle, c(a, list(rotation = rotate)))
  got <- do.call(.relative_share_gaussian, c(a, list(return_moments = TRUE)))
  expect_equal(got$log_density, expected$log_density, tolerance = 1e-11)
  expect_equal(got$conditional_mean, expected$mean, tolerance = 1e-11)
  a$quality_index <- rep(1L, length(a$residual))
  a$quality_block <- 1L
  got <- do.call(.relative_share_gaussian, c(a, list(return_moments = TRUE)))
  expected <- do.call(relative_math_oracle, a)
  expect_equal(got$log_density, expected$log_density, tolerance = 1e-11)
  expect_identical(got$conditional_mean, 0)
  expect_identical(got$conditional_variance, 0)
})

test_that("singleton cells retain proper quality uncertainty but no share dimension", {
  a <- list(residual = c(40, -20, 5), variance = c(.1, .2, .7),
            quality_index = c(1L, 2L, 1L), quality_block = c(1L, 1L, 1L),
            cell = c("first", "second", "third"), tau = .8)
  got <- do.call(.relative_share_gaussian, c(a, list(return_moments = TRUE)))
  expect_identical(got$log_density, 0)
  expect_identical(got$dimension, 0L)
  expect_identical(got$conditional_mean, rep(0, 3L))
  expect_equal(got$conditional_variance, rep(.8^2 * 2 / 3, 3L))
  expect_equal(dim(.relative_share_basis(a$cell)), c(3L, 0L))
  a$tau <- 0
  expect_identical(do.call(.relative_share_quality_draw, c(a, list(n_draws = 3L, seed = 5))),
                   matrix(0, 3L, 3L))
})

test_that("positive square-root factorization retains the thin-noise ridge", {
  delta <- 2^-34
  a <- list(residual = c(.25, -.25, .25 + delta, -.25 - delta),
            variance = delta^2, quality_index = c(1L, 2L, 1L, 2L),
            quality_block = c(1L, 1L), cell = c(1L, 1L, 2L, 2L), tau = 1)
  expected <- dnorm(.5 + delta, 0, sqrt(2 + delta^2), log = TRUE) +
    dnorm(-delta, 0, delta, log = TRUE)
  got <- do.call(.relative_share_gaussian, c(a, list(return_moments = TRUE)))
  expect_equal(got$log_density, expected, tolerance = 1e-8)
  expect_equal(got$conditional_mean, c(.25 + delta / 2, -.25 - delta / 2),
               tolerance = 1e-12)
  expect_equal(got$conditional_variance, rep(delta^2 / (2 * (2 + delta^2)), 2L),
               tolerance = 1e-30)
})

test_that("conditional quality draws preserve centered posterior correlations", {
  a <- relative_math_fixture(2L, 2L, ragged = TRUE)
  expected <- do.call(relative_math_oracle, a)
  n <- 16000L
  set.seed(219)
  before <- .Random.seed
  draws <- do.call(.relative_share_quality_draw, c(a, list(n_draws = n, seed = 191)))
  expect_identical(.Random.seed, before)
  expect_identical(draws, do.call(.relative_share_quality_draw,
                                 c(a, list(n_draws = n, seed = 191))))
  expect_equal(dim(draws), c(n, length(a$quality_block)))
  for (ii in split(seq_along(a$quality_block), a$quality_block))
    expect_true(max(abs(rowSums(draws[, ii, drop = FALSE]))) < 1e-13)
  expect_true(all(abs(colMeans(draws) - expected$mean) <
                    6 * sqrt(diag(expected$covariance) / n)))
  covariance_mcse <- sqrt((outer(diag(expected$covariance), diag(expected$covariance)) +
                            expected$covariance^2) / (n - 1L))
  expect_true(all(abs(cov(draws) - expected$covariance) < 6 * covariance_mcse))
})

test_that("relative Gaussian rejects malformed incidence and scale contracts", {
  a <- relative_math_fixture()
  wrong <- a; wrong$variance[1L] <- 0
  expect_error(do.call(.relative_share_gaussian, wrong), "positive variances")
  wrong <- a; wrong$tau <- -1
  expect_error(do.call(.relative_share_gaussian, wrong), "nonnegative")
  wrong <- a; wrong$quality_index[1L] <- 1.5
  expect_error(do.call(.relative_share_gaussian, wrong), "quality_index")
  wrong <- a; wrong$quality_block <- seq_along(a$quality_block)
  expect_error(do.call(.relative_share_gaussian, wrong), "one quality-centering block")
  wrong <- a; wrong$residual[1L] <- Inf
  expect_error(do.call(.relative_share_gaussian, wrong), "finite length-N")
  wrong <- a; wrong$cell <- wrong$cell[-1L]
  expect_error(do.call(.relative_share_gaussian, wrong), "cell")
  expect_error(do.call(.relative_share_quality_draw, c(a, list(seed = NA_real_))), "seed")
  expect_error(do.call(.relative_share_quality_draw, c(a, list(n_draws = 0))), "positive integer")
})
