# Diagnostic assumptions, not package defaults or manuscript selection.
list(profile = "pnb_bertrand_period_shock_fixed_measurement_v1",
  measurement_sd = .10, chains = 4L, warmup = 500L,
  retained_per_chain = 500L, parallel_chains = 4L, threads_per_chain = 1L,
  numerical_threads = 1L, policy_draws = 100L, budget_seconds = 5400L,
  bridge_sampling = FALSE, model_averaging = FALSE,
  adapt_delta = .999, max_treedepth = 13L, seed = 10710L,
  target_year = "1960", history_years = 1960:1965,
  large_bank_tolerance_pp = c(2, 3), small_bank_accuracy_gate = FALSE)
