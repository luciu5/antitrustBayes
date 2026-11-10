# Prevent silently changing legacy likelihoods, dropping the new component on
# CES/public routes, and applying an unnormalized or inactive measurement prior.
.share_residual_fixture <- function(demand = "logit", ...) {
  prepare <- if (demand == "logit") prepare_stan_data else prepare_ces_stan_data
  args <- list(prices = c(1, 1.4, 1.8, 2.1), quantities = c(9, 6, 4, 2),
    margins = rep(.25, 4), market = rep(1L, 4), year = rep(1L, 4),
    group = 1:4, supply_model = 1L, use_cutoff = FALSE,
    use_rho = TRUE, use_hmt = FALSE, priors = list())
  do.call(prepare, utils::modifyList(args, list(...)))
}

test_that("residual decomposition defaults cannot silently change absolute or explicit combined targets", {
  for (demand in c("logit", "ces")) {
    d <- .share_residual_fixture(demand)
    expect_identical(d$share_residual_model, 2L)
    expect_identical(attr(d, "share_residual_model"), "demand_shock")
    combined <- .share_residual_fixture(demand, share_residual_model = "combined")
    expect_identical(combined$share_residual_model, 1L)
    expect_identical(attr(combined, "share_residual_model"), "combined")
    absolute <- .share_residual_fixture(demand, share_basis = "absolute")
    expect_identical(absolute$share_residual_model, 1L)
    expect_error(.share_residual_fixture(demand, share_basis = "absolute",
      share_residual_model = "demand_shock"), "requires relative Gaussian")
    expect_error(.share_residual_fixture(demand, share_basis = "absolute",
      share_likelihood = "student_t", quality_parameterization = "sampled",
      share_residual_model = "demand_shock"), "requires relative Gaussian")
  }
})

test_that("measurement prior is proper and identical across approved manifests", {
  for (manifest in c("baseline", "tight", "diffuse")) {
    p <- bayes_prior_manifest(manifest)
    expect_identical(p$prior_sigma_measurement, .10)
    expect_identical(p$prior_log_sigma_measurement_sd, .5)
  }
  for (demand in c("logit", "ces")) {
    d <- .share_residual_fixture(demand, priors = list(
      prior_sigma_measurement = .08, prior_log_sigma_measurement_sd = .3))
    expect_identical(d$prior_sigma_measurement, .08)
    expect_identical(d$prior_log_sigma_measurement_sd, .3)
    for (field in c("prior_sigma_measurement", "prior_log_sigma_measurement_sd"))
      for (bad in list(0, -1, NA_real_, Inf))
        expect_error(.share_residual_fixture(demand,
          priors = setNames(list(bad), field)), "Prior")
  }
})

test_that("public estimator forwards explicit residual choice and measurement priors before sampling", {
  captured <- NULL
  fields <- c("share_residual_model", "prior_sigma_measurement", "prior_log_sigma_measurement_sd")
  testthat::local_mocked_bindings(.cmdstan_model_for_kernel = function(kernel) list(
    variables = function() list(data = setNames(as.list(rep(1, length(fields))), fields)),
    sample = function(data, ...) { captured <<- data; stop("captured before HMC") }),
    .package = "antitrustBayes")
  for (demand in c("logit", "ces")) for (mode in c("auto", "combined", "demand_shock")) {
    expect_error(bayes_est(prices = c(1, 1.4, 1.8, 2.1), quantities = c(9, 6, 4, 2),
      margins = rep(.25, 4), market = rep(1L, 4), group = 1:4,
      demand_model = demand, share_residual_model = mode,
      priors = list(prior_sigma_measurement = .08, prior_log_sigma_measurement_sd = .3)),
      "captured before HMC")
    expect_identical(captured$share_residual_model, if (mode == "combined") 1L else 2L)
    expect_identical(captured$prior_sigma_measurement, .08)
    expect_identical(captured$prior_log_sigma_measurement_sd, .3)
  }
})

test_that("every native target gates measurement dimension and prior without an extra Jacobian", {
  for (kernel in c("logit_supply", "ces_supply", "logit_stackelberg", "ces_stackelberg")) {
    src <- paste(readLines(.ab_test_path("inst", "stan", paste0(kernel, ".stan")), warn = FALSE), collapse = "\n")
    expect_match(src, "vector[share_residual_model == 2 ? 1 : 0] log_sigma_measurement;", fixed = TRUE)
    expect_match(src, paste0("if (share_residual_model == 2)\n",
      "    target += normal_lpdf(log_sigma_measurement | log(prior_sigma_measurement),\n",
      "      prior_log_sigma_measurement_sd);"), fixed = TRUE)
    expect_match(src, "(prior_lkj - 1.0) * log1m(square(rho_demand_shock))", fixed = TRUE)
    # A second target contribution on this primitive would change normalized
    # evidence, even when posterior comparisons appear unchanged.
    contributions <- strsplit(src, "\n", fixed = TRUE)[[1]]
    contributions <- contributions[grepl("target +=", contributions, fixed = TRUE) &
      grepl("log_sigma_measurement", contributions, fixed = TRUE)]
    expect_length(contributions, 1L)
    expect_match(src, "rho_market[m] * sigma_share_market[m]", fixed = TRUE)
  }
})
