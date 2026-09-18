helmert_quality_basis <- function(n_group) {
  out <- matrix(0, nrow = n_group, ncol = max(0L, n_group - 1L))
  if (n_group > 1L) {
    for (k in seq_len(n_group - 1L)) {
      out[seq_len(k), k] <- 1 / sqrt(k * (k + 1))
      out[k + 1L, k] <- -k / sqrt(k * (k + 1))
    }
  }
  out
}

test_that("centered quality basis is orthonormal and zero-sum", {
  for (n_group in c(1L, 3L, 5L, 8L)) {
    Q <- helmert_quality_basis(n_group)
    P <- diag(n_group) - matrix(1 / n_group, nrow = n_group, ncol = n_group)
    expect_equal(dim(Q), c(n_group, max(0L, n_group - 1L)))
    if (n_group == 1L) {
      expect_equal(Q, matrix(numeric(), nrow = 1L, ncol = 0L))
      expect_equal(as.numeric(Q %*% numeric()), 0)
    } else {
      expect_equal(crossprod(Q), diag(n_group - 1L), tolerance = 1e-14)
      expect_equal(Q %*% t(Q), P, tolerance = 1e-14)
      expect_equal(as.numeric(colSums(Q)), rep(0, n_group - 1L), tolerance = 1e-14)
      beta <- seq_len(n_group - 1L) / 7
      expect_equal(sum(as.numeric(Q %*% beta)), 0, tolerance = 1e-14)
    }

    for (sigma in c(0.02, 0.1, 1, 5)) {
      beta <- if (n_group > 1L) seq_len(n_group - 1L) * sigma / 3 else numeric()
      # The old raw-Gaussian decomposition has an orthonormal mean coordinate
      # plus these G-1 coordinates.  Its sigma^(G-1) Jacobian gives the same
      # density as the centered parameter's Normal(0, sigma) prior.
      old_coordinate_log_density <- sum(dnorm(beta / sigma, log = TRUE)) -
        length(beta) * log(sigma)
      new_coordinate_log_density <- sum(dnorm(beta, sd = sigma, log = TRUE))
      expect_equal(old_coordinate_log_density, new_coordinate_log_density,
                   tolerance = 1e-12)
      expect_equal(sigma^2 * (Q %*% t(Q)), sigma^2 * P, tolerance = 1e-14)
    }
  }
})

test_that("all Stan graphs use the centered quality coordinate contract", {
  package_root <- normalizePath(file.path(testthat::test_path(), "..", ".."),
                                mustWork = TRUE)
  source_stan_root <- file.path(package_root, "inst", "stan")
  stan_root <- if (dir.exists(source_stan_root)) {
    source_stan_root
  } else {
    system.file("stan", package = "antitrustBayes", mustWork = TRUE)
  }
  files <- file.path(stan_root, c(
    "logit_supply.stan", "ces_supply.stan",
    "logit_stackelberg.stan", "ces_stackelberg.stan"
  ))
  for (f in files) {
    code <- paste(readLines(f, warn = FALSE), collapse = "\n")
    expect_false(grepl("b_group_raw", code, fixed = TRUE), info = f)
    expect_true(grepl("vector[quality_parameterization == 1 ? N_quality_coord : 0] b_group_centered", code,
                      fixed = TRUE), info = f)
    expect_false(grepl("matrix[N_group, quality_parameterization == 1 ? N_group - 1 : 0] Q", code,
                       fixed = TRUE), info = f)
    expect_true(grepl("for (block in 1:N_quality_block)", code, fixed = TRUE), info = f)
    expect_true(grepl("real suffix = 0.0;", code, fixed = TRUE), info = f)
    expect_true(grepl("b_quality_all[quality_block_start[block] + j - 1] = suffix - k * v_k;", code, fixed = TRUE), info = f)
    expect_true(grepl("normal_lpdf(b_group_centered | 0, sigma_b_strat)",
                      code, fixed = TRUE), info = f)
  }
})
