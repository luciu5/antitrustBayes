# Six bounded known-truth fits; these are validation evidence, not paper estimates.
list(profile = "residual_model_integration_20261109", cases = c("input_logit_snapshot",
  "input_logit_history", "output_ces_panel"), models = c("combined", "demand_shock"),
  start_utc = "2026-11-09T13:51:24+00:00", deadline_utc = "2026-11-09T15:21:24+00:00",
  maximum_fits = 6L, chains = 4L, warmup = 500L, retained = 500L,
  cores = 4L, threads_per_chain = 1L, numerical_threads = 1L,
  adapt_delta = .999, max_treedepth = 13L, policy_draws = 100L,
  seed = 915113L, policy_seed = 915213L, quality_sd = .45,
  demand_shock_sd = .25, measurement_sd = .10, margin_sd = .12,
  rho_demand = .35, measurement_prior_median = .10,
  measurement_prior_log_sd = .5, bridge_sampling = FALSE, model_averaging = FALSE)
