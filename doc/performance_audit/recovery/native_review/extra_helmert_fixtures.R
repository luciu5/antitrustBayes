suppressPackageStartupMessages(library(antitrustBayes))

ns <- asNamespace("antitrustBayes")
resolve <- get("resolve_prior_manifest", ns)

multi <- list(
  market = c("M1", "M1", "M1", "M1", "M2", "M2", "M2", "M2"),
  year = c("Y1", "Y1", "Y2", "Y2", "Y1", "Y1", "Y2", "Y2"),
  group = c("G1", "G1", "G2", "G3", "G1", "G2", "G4", "G5"),
  prices = c(1.0, 1.15, 1.3, 1.45, 1.05, 1.2, 1.35, 1.5),
  quantities = c(1.2, 1.0, .9, .8, 1.1, .95, .85, .75),
  margins = c(.2, .22, .18, .21, .19, .23, .2, .24)
)

g1000_single <- list(
  market = rep("M1", 1000),
  year = rep("Y1", 1000),
  group = paste0("G", seq_len(1000)),
  prices = seq(1, 2, length.out = 1000),
  quantities = exp(seq(log(1.3), log(.35), length.out = 1000)),
  margins = seq(.18, .27, length.out = 1000)
)

g1000_panel <- list(
  market = rep(c("M1", "M2"), each = 1000),
  year = rep(c("Y1", "Y2"), each = 1000),
  group = rep(paste0("G", seq_len(1000)), 2),
  prices = rep(seq(1, 2, length.out = 1000), 2) +
    rep(c(0, .05), each = 1000),
  quantities = rep(exp(seq(log(1.3), log(.35), length.out = 1000)), 2),
  margins = rep(seq(.18, .27, length.out = 1000), 2)
)

make_extra_data <- function(demand, structure = "simultaneous",
                            quality = "sampled", base = multi,
                            supply_model = 1L, use_cutoff = FALSE,
                            use_hmt = FALSE, use_rho = FALSE,
                            use_market_resid_scale = FALSE,
                            use_alpha_market_pool = TRUE,
                            use_s0_market_year_pool = FALSE,
                            use_s0_market_pool = FALSE,
                            market_s0_parameterization = "legacy",
                            use_centered_alpha = FALSE,
                            use_centered_s0_pool = FALSE, X_s0 = NULL,
                            margins = base$margins, margin_observed = NULL,
                            outside_price = 1, likelihood_mode = 1L,
                            share_likelihood = "student_t",
                            margin_likelihood = NULL,
                            cutoff_scale = NULL,
                            partition_prior = "cutoff", k_support = NULL,
                            k_prior = NULL) {
  if (is.null(margin_likelihood)) margin_likelihood <- share_likelihood
  if (is.null(margin_observed)) margin_observed <- rep(TRUE, length(margins))
  args <- list(
    prices = base$prices, quantities = base$quantities, margins = margins,
    market = base$market, year = base$year, group = base$group,
    supply_model = as.integer(supply_model), use_cutoff = use_cutoff,
    use_rho = use_rho, use_market_resid_scale = use_market_resid_scale,
    use_hmt = use_hmt, use_alpha_market_pool = use_alpha_market_pool,
    use_s0_market_year_pool = use_s0_market_year_pool,
    use_s0_market_pool = use_s0_market_pool,
    use_centered_alpha = use_centered_alpha,
    use_centered_s0_pool = use_centered_s0_pool,
    market_s0_parameterization = market_s0_parameterization,
    X_s0 = X_s0, priors = resolve("baseline", list()), grainsize = 1L,
    margin_observed = margin_observed, outside_price = outside_price,
    likelihood_mode = as.integer(likelihood_mode), loo_cluster = base$market,
    estimate_df = FALSE, resid_df_fixed = 5, structure = structure,
    cutoff_scale = cutoff_scale, partition_prior = partition_prior,
    k_support = k_support, k_prior = k_prior, k_scope = "market",
    share_likelihood = share_likelihood,
    margin_likelihood = margin_likelihood,
    quality_parameterization = quality
  )
  fun <- if (identical(demand, "logit")) get("prepare_stan_data", ns)
    else get("prepare_ces_stan_data", ns)
  if (identical(demand, "ces")) args$market_side <- "output"
  dat <- do.call(fun, args)
  dat$structure <- switch(structure, simultaneous = 0L, fringe = 1L,
                          stackelberg = 2L)
  dat
}

## Four small fixtures exercise likelihood and conduct branches. The final two
## deliberately use all 999 sampled quality contrasts and expose both ends of
## the suffix recurrence; the panel repeats every group across markets.
extra_cases <- list(
  logit_sampled_share_only_panel = list(
    demand = "logit", structure = "simultaneous", quality = "sampled",
    base = multi, supply_model = 1L, use_rho = FALSE,
    likelihood_mode = 2L, share_likelihood = "student_t",
    margin_likelihood = "student_t", margins = rep(NA_real_, 8),
    margin_observed = rep(FALSE, 8)
  ),
  ces_sampled_margin_only_panel = list(
    demand = "ces", structure = "simultaneous", quality = "sampled",
    base = multi, supply_model = 3L, use_rho = FALSE,
    likelihood_mode = 3L, share_likelihood = "gaussian",
    margin_likelihood = "gaussian", margin_observed = rep(TRUE, 8)
  ),
  logit_smooth_cutoff_hmt = list(
    demand = "logit", structure = "fringe", quality = "sampled",
    base = multi, supply_model = 1L, use_cutoff = TRUE, use_hmt = TRUE,
    use_rho = FALSE, cutoff_scale = .17, partition_prior = "cutoff",
    share_likelihood = "gaussian", margin_likelihood = "gaussian"
  ),
  logit_bargaining = list(
    demand = "logit", structure = "simultaneous", quality = "sampled",
    base = multi, supply_model = 5L, use_rho = FALSE,
    share_likelihood = "student_t", margin_likelihood = "student_t"
  ),
  ces_second_score_auction = list(
    demand = "ces", structure = "simultaneous", quality = "sampled",
    base = multi, supply_model = 4L, use_rho = FALSE,
    share_likelihood = "gaussian", margin_likelihood = "gaussian"
  ),
  logit_g1000_single_product = list(
    demand = "logit", structure = "simultaneous", quality = "sampled",
    base = g1000_single, supply_model = 1L, use_rho = FALSE,
    share_likelihood = "gaussian", margin_likelihood = "gaussian"
  ),
  ces_g1000_repeated_panel = list(
    demand = "ces", structure = "simultaneous", quality = "sampled",
    base = g1000_panel, supply_model = 3L, use_rho = FALSE,
    share_likelihood = "gaussian", margin_likelihood = "gaussian"
  )
)
