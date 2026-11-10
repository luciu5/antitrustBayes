# History-only diagnostic; snapshots and all production defaults stay unchanged.
list(profile = "pnb_history_inverse_assets_log_variance_linear_time_v1",
  asset_field = "total_assets", asset_timing = "contemporaneous_bank_year",
  asset_scaling = "unique_bank_years_across_all_banking_cases",
  fixed_variance_exponent = 1, variance_time_prior_sd = 1,
  trend_base_year = 1960L, trend_scale_years = 5,
  chains = 4L, warmup = 500L, retained_per_chain = 500L,
  parallel_chains = 4L, threads_per_chain = 1L, numerical_threads = 1L,
  policy_draws = 100L, budget_seconds = 5400L,
  bridge_sampling = FALSE, model_averaging = FALSE,
  adapt_delta = .999, max_treedepth = 13L, seed = 11210L,
  policy_seed = 10902L, target_year = "1960", history_years = 1960:1965,
  large_bank_tolerance_pp = c(2, 3), small_bank_accuracy_gate = FALSE)
