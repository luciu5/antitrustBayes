suppressPackageStartupMessages({
  library(antitrustBayes)
})
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
single5 <- list(
  market = rep("M1", 5), year = rep("Y1", 5), group = paste0("G", 1:5),
  prices = seq(1, 1.5, length.out = 5), quantities = c(1.1, 1, .9, .8, .7),
  margins = rep(.2, 5)
)

make_data <- function(demand, structure = "simultaneous", quality = "sampled",
                      base = multi, supply_model = 1L,
                      share_likelihood = "student_t", margin_likelihood = NULL,
                      use_rho = FALSE, use_s0_market_pool = FALSE,
                      use_s0_market_year_pool = FALSE,
                      market_s0_parameterization = "legacy",
                      use_centered_alpha = FALSE,
                      use_centered_s0_pool = FALSE, X_s0 = NULL,
                      margins = base$margins, margin_observed = NULL,
                      k_support = NULL, k_prior = NULL) {
  if (is.null(margin_likelihood)) margin_likelihood <- share_likelihood
  args <- list(
    prices = base$prices, quantities = base$quantities, margins = margins,
    market = base$market, year = base$year, group = base$group,
    supply_model = supply_model, use_cutoff = FALSE,
    use_rho = use_rho, use_market_resid_scale = FALSE, use_hmt = FALSE,
    use_alpha_market_pool = TRUE,
    use_s0_market_year_pool = use_s0_market_year_pool,
    use_s0_market_pool = use_s0_market_pool,
    use_centered_alpha = use_centered_alpha,
    use_centered_s0_pool = use_centered_s0_pool,
    market_s0_parameterization = market_s0_parameterization,
    X_s0 = X_s0, priors = resolve("baseline", list()), grainsize = 1L,
    margin_observed = margin_observed, outside_price = 1,
    likelihood_mode = 1L, loo_cluster = base$market, estimate_df = FALSE,
    resid_df_fixed = 5, structure = structure, cutoff_scale = NULL,
    partition_prior = if (structure == "simultaneous") "cutoff" else "ordered_k",
    k_support = k_support, k_prior = k_prior, k_scope = "market",
    share_likelihood = share_likelihood,
    margin_likelihood = margin_likelihood,
    quality_parameterization = quality
  )
  fun <- if (demand == "logit") get("prepare_stan_data", ns)
    else get("prepare_ces_stan_data", ns)
  if (demand == "ces") args$market_side <- "output"
  dat <- do.call(fun, args)
  dat$structure <- switch(structure, simultaneous = 0L, fringe = 1L,
                          stackelberg = 2L)
  dat
}

cases <- list(
  sampled5_student_rho = list(
    demand = "logit", structure = "simultaneous", quality = "sampled",
    base = single5, share_likelihood = "student_t", margin_likelihood = "student_t",
    use_rho = TRUE, supply_model = 1L
  ),
  sampled_multi_gaussian = list(
    demand = "logit", structure = "simultaneous", quality = "sampled",
    base = multi, share_likelihood = "gaussian", margin_likelihood = "gaussian",
    use_rho = FALSE, supply_model = 1L
  ),
  collapsed_control = list(
    demand = "logit", structure = "simultaneous", quality = "collapsed",
    base = multi, share_likelihood = "gaussian", margin_likelihood = "gaussian",
    use_rho = FALSE, supply_model = 1L
  ),
  rho_missing = list(
    demand = "logit", structure = "simultaneous", quality = "sampled",
    base = multi, share_likelihood = "student_t", margin_likelihood = "student_t",
    use_rho = TRUE, supply_model = 1L,
    margins = c(NA, .22, NA, .21, .19, NA, .2, NA),
    margin_observed = c(FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE)
  ),
  pool_market_legacy_noncentered = list(
    demand = "logit", structure = "simultaneous", quality = "sampled",
    base = multi, share_likelihood = "student_t", margin_likelihood = "student_t",
    use_rho = TRUE, use_s0_market_pool = TRUE,
    market_s0_parameterization = "legacy", use_centered_s0_pool = FALSE,
    use_centered_alpha = FALSE, supply_model = 1L
  ),
  pool_market_legacy_centered = list(
    demand = "logit", structure = "simultaneous", quality = "sampled",
    base = multi, share_likelihood = "student_t", margin_likelihood = "student_t",
    use_rho = TRUE, use_s0_market_pool = TRUE,
    market_s0_parameterization = "legacy", use_centered_s0_pool = TRUE,
    use_centered_alpha = TRUE, supply_model = 1L
  ),
  pool_market_year_rotated_noncentered_X = list(
    demand = "logit", structure = "simultaneous", quality = "sampled",
    base = multi, share_likelihood = "student_t", margin_likelihood = "student_t",
    use_rho = TRUE, use_s0_market_year_pool = TRUE,
    market_s0_parameterization = "rotated", use_centered_s0_pool = FALSE,
    use_centered_alpha = FALSE, X_s0 = matrix(c(0, .1, -.1, .2), 4, 1),
    supply_model = 1L
  ),
  pool_market_year_rotated_centered = list(
    demand = "logit", structure = "simultaneous", quality = "sampled",
    base = multi, share_likelihood = "student_t", margin_likelihood = "student_t",
    use_rho = TRUE, use_s0_market_year_pool = TRUE,
    market_s0_parameterization = "rotated", use_centered_s0_pool = TRUE,
    use_centered_alpha = TRUE, supply_model = 1L
  ),
  exact_fringe_BF = list(
    demand = "logit", structure = "fringe", quality = "sampled", base = multi,
    share_likelihood = "student_t", margin_likelihood = "student_t", use_rho = FALSE,
    supply_model = 1L, k_support = 0:3, k_prior = rep(.25, 4)
  ),
  exact_fringe_CF = list(
    demand = "logit", structure = "fringe", quality = "sampled", base = multi,
    share_likelihood = "student_t", margin_likelihood = "student_t", use_rho = FALSE,
    supply_model = 3L, k_support = 0:3, k_prior = rep(.25, 4)
  ),
  exact_stack_BS = list(
    demand = "logit", structure = "stackelberg", quality = "sampled", base = multi,
    share_likelihood = "student_t", margin_likelihood = "student_t", use_rho = FALSE,
    supply_model = 1L, k_support = 1:3, k_prior = rep(1 / 3, 3)
  ),
  exact_stack_CS = list(
    demand = "logit", structure = "stackelberg", quality = "sampled", base = multi,
    share_likelihood = "student_t", margin_likelihood = "student_t", use_rho = FALSE,
    supply_model = 3L, k_support = 1:3, k_prior = rep(1 / 3, 3)
  ),
  ces_exact_fringe_BF = list(
    demand = "ces", structure = "fringe", quality = "sampled", base = multi,
    share_likelihood = "student_t", margin_likelihood = "student_t", use_rho = FALSE,
    supply_model = 1L, k_support = 0:3, k_prior = rep(.25, 4)
  ),
  ces_exact_fringe_CF = list(
    demand = "ces", structure = "fringe", quality = "sampled", base = multi,
    share_likelihood = "student_t", margin_likelihood = "student_t", use_rho = FALSE,
    supply_model = 3L, k_support = 0:3, k_prior = rep(.25, 4)
  ),
  ces_exact_stack_BS = list(
    demand = "ces", structure = "stackelberg", quality = "sampled", base = multi,
    share_likelihood = "student_t", margin_likelihood = "student_t", use_rho = FALSE,
    supply_model = 1L, k_support = 1:3, k_prior = rep(1 / 3, 3)
  ),
  ces_exact_stack_CS = list(
    demand = "ces", structure = "stackelberg", quality = "sampled", base = multi,
    share_likelihood = "student_t", margin_likelihood = "student_t", use_rho = FALSE,
    supply_model = 3L, k_support = 1:3, k_prior = rep(1 / 3, 3)
  ),
  ces_sampled_multi = list(
    demand = "ces", structure = "simultaneous", quality = "sampled", base = multi,
    share_likelihood = "student_t", margin_likelihood = "student_t", use_rho = TRUE,
    supply_model = 3L
  ),
  ces_collapsed = list(
    demand = "ces", structure = "simultaneous", quality = "collapsed", base = multi,
    share_likelihood = "gaussian", margin_likelihood = "gaussian", use_rho = FALSE,
    supply_model = 3L
  )
)

for (nm in names(cases)) {
  x <- cases[[nm]]
  ans <- tryCatch({
    dat <- do.call(make_data, x)
    c(status = "ok", N = dat$N, N_market = dat$N_market,
      N_market_year = dat$N_market_year, N_group = dat$N_group,
      K_s0 = dat$K_s0, rho_active = if ("rho_active" %in% names(dat)) dat$rho_active else NA)
  }, error = function(e) c(status = "error", message = conditionMessage(e)))
  cat(nm, paste(names(ans), ans, sep = "=", collapse = " | "), "\n")
}
