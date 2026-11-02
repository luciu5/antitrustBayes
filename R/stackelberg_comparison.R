# A marginalized count belongs to the whole partition unit. Prefix likelihood
# accounting is not a leave-one-product density after integrating that count.
.bayes_is_stackelberg <- function(fit) {
  isTRUE(fit@stan_data$structure == 2L) ||
    identical(fit@metadata$structure, "stackelberg")
}

.bayes_is_exact_ordered_k <- function(fit) {
  structure <- .stack_bridge_or(
    fit@metadata$structure,
    if (isTRUE(fit@stan_data$structure == 2L)) {
      "stackelberg"
    } else if (isTRUE(fit@stan_data$structure == 1L)) {
      "fringe"
    } else {
      "simultaneous"
    }
  )
  prior <- .stack_bridge_or(
    fit@metadata$partition_prior,
    if (isTRUE(fit@stan_data$partition_prior_code == 1L)) {
      "ordered_k"
    } else {
      "cutoff"
    }
  )
  structure %in% c("fringe", "stackelberg") &&
    identical(prior, "ordered_k")
}

.bayes_uses_cutoff_partition <- function(fit) {
  structure <- .stack_bridge_or(
    fit@metadata$structure,
    if (isTRUE(fit@stan_data$structure == 2L)) {
      "stackelberg"
    } else if (isTRUE(fit@stan_data$use_cutoff == 1L)) {
      "fringe"
    } else {
      "simultaneous"
    }
  )
  structure %in% c("fringe", "stackelberg") &&
    !.bayes_is_exact_ordered_k(fit) &&
    isTRUE(fit@stan_data$use_cutoff == 1L)
}

.bayes_is_smooth_fringe <- function(fit) {
  structure <- .stack_bridge_or(
    fit@metadata$structure,
    if (isTRUE(fit@stan_data$use_cutoff == 1L)) "fringe" else "simultaneous"
  )
  identical(structure, "fringe") &&
    .bayes_uses_cutoff_partition(fit)
}

.bayes_comparison_units <- function(fits, loo_unit) {
  for (fit in fits) {
    if (identical(.share_basis_from_fit(fit), "relative")) {
      if (loo_unit == "observation")
        stop("Observation LOO is unsupported for relative shares; use complete market-year cells.", call. = FALSE)
      if (loo_unit == "cluster") {
        counts <- tapply(fit@stan_data$cluster, fit@stan_data$market_year_idx,
                         function(x) length(unique(x)))
        if (any(counts != 1L))
          stop("Relative-share LOO clusters must contain complete market-year cells.", call. = FALSE)
      }
    }
  }
  for (fit in fits) if (.bayes_is_exact_ordered_k(fit)) {
    dat <- fit@stan_data
    if (loo_unit == "observation") {
      stop("Observation LOO is unsupported for whole-market ordered-k mixtures.")
    }
    if (loo_unit == "period" && isTRUE(dat$k_scope_code == 1L)) {
      stop("Period LOO is unsupported when one ordered k is shared across the market's years.")
    }
    if (loo_unit == "cluster") {
      unit <- dat$mky_to_partition_unit[dat$market_year_idx]
      counts <- tapply(dat$cluster, unit, function(x) length(unique(x)))
      if (any(counts != 1L)) {
        stop("LOO clusters must contain whole ordered-k partition units.")
      }
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
        !isTRUE(normalization$model_spec_version == .antitrust_bayes_model_spec_version)) {
      stop("Marginal likelihood requires native normalized current-specification fits. Refit older or unverified objects.")
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
