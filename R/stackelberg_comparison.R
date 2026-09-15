# A marginalized count belongs to the whole partition unit. Prefix likelihood
# accounting is not a leave-one-product density after integrating that count.
.bayes_is_stackelberg <- function(fit) {
  isTRUE(fit@stan_data$structure == 2L) ||
    identical(fit@metadata$structure, "stackelberg")
}
.bayes_comparison_units <- function(fits, loo_unit) {
  for (fit in fits) if (.bayes_is_stackelberg(fit)) {
    dat <- fit@stan_data
    if (loo_unit == "observation") {
      stop("Observation LOO is unsupported for whole-equilibrium Stackelberg mixtures.")
    }
    if (loo_unit == "period" && isTRUE(dat$k_scope_code == 1L)) {
      stop("Period LOO is unsupported when one Stackelberg k is shared across the market's years.")
    }
    if (loo_unit == "cluster") {
      unit <- dat$mky_to_partition_unit[dat$market_year_idx]
      counts <- tapply(dat$cluster, unit, function(x) length(unique(x)))
      if (any(counts != 1L)) stop("LOO clusters must contain whole Stackelberg partition units.")
    }
  }
  invisible(TRUE)
}
.bayes_validate_bridge <- function(fits) {
  for (fit in fits) {
    dat <- fit@stan_data
    normalization <- fit@metadata$prior_normalization
    if (!is.list(normalization) ||
        !identical(normalization$mode, "native_stan_normalized") ||
        !isTRUE(normalization$model_spec_version == 8L)) {
      stop("Marginal likelihood requires native normalized specification-8 fits. Refit older or unverified objects.")
    }
    if (isTRUE(dat$use_s0_market_pool == 1L) &&
        isTRUE(dat$use_s0_market_year_pool == 1L)) {
      stop("Both outside-share pooling flags declare coordinates without a prior; evidence is undefined.")
    }
    if (isTRUE(dat$use_hmt == 1L) && !isTRUE(dat$likelihood_mode == 2L)) {
      stop("Marginal likelihood is unsupported with the non-normalized HMT penalty. Use bridge = FALSE for predictive summaries.")
    }
  }
  invisible(TRUE)
}
