test_that("collapsed Gaussian Stan paths are sparse, normalized, and scalar-accounted", {
  stan_dir <- .ab_test_path("inst", "stan")
  files <- file.path(stan_dir, c(
    "logit_supply.stan", "ces_supply.stan",
    "logit_stackelberg.stan", "ces_stackelberg.stan"
  ))
  expect_true(all(file.exists(files)))
  code <- lapply(files, function(f) paste(readLines(f, warn = FALSE), collapse = "\n"))

  for (src in code) {
    expect_match(src, "centered_gaussian_from_stats", fixed = TRUE)
    expect_match(src, "log(den / G)", fixed = TRUE)
    expect_match(src, "L * log(2.0 * pi())", fixed = TRUE)
    expect_match(src, "if (G == 1)", fixed = TRUE)
    expect_match(src, "real q0 = within;", fixed = TRUE)
    expect_match(src, "square(c[j] / a[j]) * a[j] * w[j]", fixed = TRUE)
    expect_false(grepl("q0 -=", src, fixed = TRUE))
    expect_match(src, "log_lik_share_collapsed", fixed = TRUE)
    expect_match(src, "log_lik_joint_collapsed", fixed = TRUE)
    expect_match(src, "quality_parameterization == 1 ? N_quality_coord : 0",
                 fixed = TRUE)
    expect_match(src, "vector[quality_parameterization == 1 ? N_quality : 0] b_quality_all",
                 fixed = TRUE)
    expect_false(grepl("matrix\\s*\\[\\s*N\\s*,\\s*N\\s*\\]", src,
                       perl = TRUE))
    expect_false(grepl("matrix\\s*\\[\\s*N_group\\s*,\\s*N_group\\s*-\\s*1\\s*\\]",
                       src, perl = TRUE))
    expect_false(grepl("multi_normal|cholesky", src, perl = TRUE,
                       ignore.case = TRUE))
  }

  for (src in code[3:4]) {
    expect_match(src, "real collapsed_share_ll = 0.0", fixed = TRUE)
    expect_match(src, "collapsed_joint_from_share_stats", fixed = TRUE)
    expect_match(src, "N_market == 1", fixed = TRUE)
    expect_match(src, "target += log_sum_exp(score)", fixed = TRUE)
    expect_match(src, "- collapsed_share_ll", fixed = TRUE)
    expect_match(src, "(quality_scope == 2 || N_market == 1)", fixed = TRUE)
    expect_match(src, "else for (pos in 1:market_obs_n[m])", fixed = TRUE)
    expect_match(src, "unit_margin_loglik", fixed = TRUE)
    expect_match(src, "posterior_component_prob", fixed = TRUE)
  }
})

test_that("production source keeps the common share score outside ordered-k mixtures", {
  stan_dir <- .ab_test_path("inst", "stan")
  src <- paste(readLines(file.path(stan_dir, "logit_stackelberg.stan"),
                         warn = FALSE), collapse = "\n")
  share_pos <- regexpr("collapsed_share_ll = centered_gaussian_from_stats", src,
                       fixed = TRUE)[[1L]]
  mixture_pos <- regexpr("target += log_sum_exp(score)", src, fixed = TRUE)[[1L]]
  joint_pos <- regexpr("score[kk] += collapsed_joint_from_share_stats", src,
                       fixed = TRUE)[[1L]]
  expect_true(share_pos > 0L && mixture_pos > share_pos && joint_pos > share_pos)
  expect_true(regexpr("- collapsed_share_ll", src, fixed = TRUE)[[1L]] > share_pos)
})
