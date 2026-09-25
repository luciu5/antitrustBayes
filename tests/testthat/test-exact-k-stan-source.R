test_that("exact ordered-k Stan sources use one market mixture", {
  stan_names <- c("logit_stackelberg.stan", "ces_stackelberg.stan")
  files <- file.path(.ab_test_path("inst", "stan"), stan_names)
  expect_true(all(nzchar(files) & file.exists(files)))
  code <- lapply(files, function(path) paste(readLines(path, warn = FALSE),
                                             collapse = "\n"))

  for (src in code) {
    expect_match(src, "int<lower=1, upper=2> structure", fixed = TRUE)
    expect_match(src, "array[N_market] int<lower=1> market_obs_start", fixed = TRUE)
    expect_match(src, "array[N_market_year] int<lower=1, upper=N_market_year> market_mky_order", fixed = TRUE)
    expect_match(src, "array[N_partition_component] int<lower=0> partition_component_k", fixed = TRUE)
    expect_match(src, "posterior_component_prob", fixed = TRUE)
    expect_match(src, "component_log_score", fixed = TRUE)
    expect_match(src, "posterior_k_mean", fixed = TRUE)
    expect_match(src, "posterior_k_median", fixed = TRUE)
    expect_match(src, "posterior_k_mode", fixed = TRUE)
    expect_match(src, "if (likelihood_mode != 2 && margin_obs[i] == 1 &&", fixed = TRUE)
    expect_match(src, "partition_component_k[cr] <= partition_component_k[cc]", fixed = TRUE)
    expect_match(src, "firm_market_core_prob", fixed = TRUE)
    expect_match(src, "current - previous", fixed = TRUE)
    expect_false(grepl("mu_b_fringe", src, fixed = TRUE))
    expect_false(grepl("cutoff_share", src, fixed = TRUE))
    expect_false(grepl("real<lower=0> cutoff_scale;", src, fixed = TRUE))
    expect_false(grepl("conduct_weight", src, fixed = TRUE))
  }

  logit <- code[[1L]]
  ces <- code[[2L]]
  expect_match(logit, "if (firm_mky_rank[me] > k) h = 1.0", fixed = TRUE)
  expect_match(ces, "if (firm_mky_rank[me] > k) m = 1 / gamma", fixed = TRUE)
  expect_match(logit, "target += log_sum_exp(score)", fixed = TRUE)
  expect_match(ces, "target += log_sum_exp(score)", fixed = TRUE)
})
