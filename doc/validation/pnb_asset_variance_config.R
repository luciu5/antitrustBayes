# Diagnostic only. Assets are institution totals, never local-market shares.
list(profile = "pnb_bertrand_institution_asset_variance_v1",
  asset_field = "total_assets", asset_timing = "contemporaneous_bank_year",
  asset_scaling = "unique_bank_years_across_all_banking_cases",
  gamma_prior_sd = .35, chains = 4L, warmup = 500L,
  retained_per_chain = 500L, parallel_chains = 4L, threads_per_chain = 1L,
  numerical_threads = 1L, policy_draws = 100L, budget_seconds = 5400L,
  bridge_sampling = FALSE, model_averaging = FALSE,
  adapt_delta = .999, max_treedepth = 13L, seed = 10910L,
  target_year = "1960", history_years = 1960:1965,
  large_bank_tolerance_pp = c(2, 3), small_bank_accuracy_gate = FALSE)
