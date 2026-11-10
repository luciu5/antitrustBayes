## Central contract for extracting posterior quantities used by policy code.
## Keeping this in one place makes the amount of data pulled from CmdStanR
## explicit, and gives callers a common missing-field/alias failure.

.posterior_fit_draw_names <- function(fit) {
  if (methods::is(fit, "BayesSpecifiedFit") || is.null(fit@stanfit)) {
    return(character())
  }
  .bayes_backend_variable_names(fit@stanfit)
}

.posterior_contract_structure <- function(fit) {
  md <- fit@metadata
  sd <- if (methods::is(fit, "BayesSpecifiedFit")) NULL else fit@stan_data
  value <- md$structure %||% attr(sd, "structure", exact = TRUE)
  if (is.null(value) && !is.null(sd) && isTRUE(sd$use_cutoff == 1L)) value <- "fringe"
  if (is.null(value)) value <- "simultaneous"
  as.character(value[[1L]])
}

.posterior_contract_partition <- function(fit) {
  md <- fit@metadata
  sd <- if (methods::is(fit, "BayesSpecifiedFit")) NULL else fit@stan_data
  value <- md$partition_prior %||% attr(sd, "partition_prior", exact = TRUE)
  if (is.null(value)) value <- "cutoff"
  as.character(value[[1L]])
}

.posterior_contract_quality <- function(fit) {
  if (methods::is(fit, "BayesSpecifiedFit")) return("none")
  stan_value <- fit@stan_data$quality_parameterization
  value <- fit@metadata$quality_parameterization %||%
    attr(fit@stan_data, "quality_parameterization", exact = TRUE) %||%
    stan_value
  if (is.numeric(value) && length(value) == 1L) {
    value <- if (value == 1) "sampled" else if (value == 2) "collapsed" else NA_character_
  }
  if (is.character(value) && length(value) == 1L &&
      value %in% c("sampled", "collapsed")) return(value)
  if (isTRUE(fit@metadata$bayes_is_collapsed_quality)) return("collapsed")
  "none"
}

## Missing residual metadata denotes the frozen historical combined model.
.posterior_demand_shock_model <- function(fit) {
  if (methods::is(fit, "BayesSpecifiedFit")) return(FALSE)
  values <- list(fit@metadata$share_residual_model,
    attr(fit@stan_data, "share_residual_model", exact = TRUE),
    fit@stan_data$share_residual_model)
  values <- values[!vapply(values, is.null, logical(1))]
  if (!length(values)) return(FALSE)
  modes <- vapply(values, function(value) {
    if (length(value) != 1L || is.na(value))
      stop("Malformed share_residual_model metadata.", call. = FALSE)
    if (is.numeric(value)) value <- switch(as.character(value),
      `1` = "combined", `2` = "demand_shock", "invalid")
    if (!value %in% c("combined", "demand_shock"))
      stop("Malformed share_residual_model metadata.", call. = FALSE)
    value
  }, character(1))
  if (length(unique(modes)) != 1L)
    stop("Inconsistent share_residual_model metadata.", call. = FALSE)
  identical(modes[[1L]], "demand_shock")
}

.posterior_contract_rho_active <- function(fit) {
  if (methods::is(fit, "BayesSpecifiedFit")) return(FALSE)
  isTRUE(fit@stan_data$use_rho == 1L) &&
    isTRUE(fit@stan_data$likelihood_mode == 1L) &&
    length(fit@stan_data$margin_idx) > 0L
}

.posterior_contract_native <- function(fit) {
  structure <- .posterior_contract_structure(fit)
  identical(structure, "stackelberg") ||
    (identical(structure, "fringe") &&
       identical(.posterior_contract_partition(fit), "ordered_k"))
}

## Return canonical names and the actual Stan names available for this fit.
## `optional` fields are requested when present, preserving old fixtures that
## predate an optional leader-prior generated quantity.
.posterior_contract <- function(fit,
                                purpose = c("simulation", "policy_plan",
                                            "native_model", "posterior_quality",
                                            "cutoff_summary", "posterior_k",
                                            "posterior_predictive"),
                                conditioning = "auto") {
  purpose <- match.arg(purpose)
  if (methods::is(fit, "BayesSpecifiedFit")) {
    return(list(canonical = character(), actual = character(),
                aliases = character(), purpose = purpose))
  }
  structure <- .posterior_contract_structure(fit)
  partition <- .posterior_contract_partition(fit)
  native <- .posterior_contract_native(fit)
  quality <- .posterior_contract_quality(fit)
  intercept_collapsed <- identical(
    .market_intercept_parameterization_from_fit(fit), "collapsed"
  )
  relative <- identical(.share_basis_from_fit(fit), "relative")
  share_mean_field <- if (relative) {
    "relative_share_predictor"
  } else if (intercept_collapsed) {
    "share_mean_integrated"
  } else {
    "share_mean"
  }
  required <- character()
  optional <- character()
  if (purpose %in% c("simulation", "policy_plan", "native_model",
                     "posterior_predictive")) {
    required <- c(required, "a_market", "s0_prob")
    if (structure %in% c("fringe", "stackelberg") &&
        identical(partition, "cutoff")) required <- c(required, "cutoff_share")
    if (native) required <- c(required, "posterior_component_prob")
    if (length(fit@supply_model) == 1L && as.integer(fit@supply_model) == 5L)
      required <- c(required, "bargaining_power")
    if (purpose %in% c("policy_plan", "posterior_predictive")) {
      if (quality == "sampled") {
        required <- c(required, share_mean_field)
        ## One global group is deterministically centered at zero and the
        ## historical recovery path intentionally accepts no b_quality_all.
        n_quality <- fit@stan_data$N_quality %||% fit@stan_data$N_group %||%
          if (!is.null(fit@stan_data$group)) max(fit@stan_data$group) else NA_integer_
        if (length(n_quality) == 1L && is.finite(n_quality) && n_quality > 1L)
          required <- c(required, "b_quality_all")
      }
      if (quality == "collapsed") {
        required <- c(required, share_mean_field, "sigma_b_strat", "sigma_share")
        if (intercept_collapsed) {
          required <- c(required, "mu_b_market_vec", "sigma_b_market_vec")
          if (isTRUE(fit@stan_data$use_market_resid_scale == 1L))
            required <- c(required, "sigma_share_market")
        }
        scope <- if (methods::is(fit, "BayesLogit"))
          .quality_scope_from_fit(fit, strict = FALSE) else "global"
        if (.posterior_contract_rho_active(fit) &&
            identical(conditioning, "auto") &&
            !(native || identical(partition, "ordered_k"))) {
          required <- c(required, "rho", "sigma_margin", "pred_marginNormInv")
          if ((intercept_collapsed || relative) &&
              isTRUE(fit@stan_data$use_market_resid_scale == 1L))
            required <- c(required, "sigma_margin_market")
        }
        exact_quality <- .posterior_contract_rho_active(fit) &&
          identical(conditioning, "auto") &&
          (native || identical(partition, "ordered_k")) &&
          (identical(scope, "market") ||
           (identical(scope, "global") &&
            as.integer(fit@stan_data$N_market %||%
              if (!is.null(fit@stan_data$market)) max(fit@stan_data$market) else 1L) == 1L))
        if (exact_quality)
          required <- c(required, "rho", "sigma_margin", "posterior_component_prob",
                        "a_market", "s0_prob")
        if (exact_quality && (intercept_collapsed || relative) &&
            isTRUE(fit@stan_data$use_market_resid_scale == 1L))
          required <- c(required, "sigma_margin_market")
      }
    }
  }
  if (purpose == "posterior_predictive") {
    required <- c(required, "sigma_share", "sigma_margin")
    required <- c(required, share_mean_field)
    if (!native) required <- c(required, "pred_marginNormInv")
    if (.posterior_contract_rho_active(fit)) required <- c(required, "rho")
    if (isTRUE(fit@stan_data$use_market_resid_scale == 1L)) {
      required <- c(required, "sigma_share_market", "sigma_margin_market")
    }
    if (identical(fit@metadata$share_likelihood, "student_t") ||
        identical(fit@metadata$margin_likelihood, "student_t") ||
        isTRUE(fit@stan_data$share_likelihood == 1L) ||
        isTRUE(fit@stan_data$margin_likelihood == 1L)) {
      required <- c(required, "nu")
    }
  }
  if (purpose == "posterior_quality") {
    if (!identical(quality, "collapsed")) {
      stop("posterior_quality extraction requires a collapsed-quality fit.", call. = FALSE)
    }
    required <- c(share_mean_field, "sigma_b_strat", "sigma_share")
    if (intercept_collapsed) {
      required <- c(required, "mu_b_market_vec", "sigma_b_market_vec")
      if (isTRUE(fit@stan_data$use_market_resid_scale == 1L))
        required <- c(required, "sigma_share_market")
    }
    scope <- .quality_scope_from_fit(fit, strict = FALSE)
    if (.posterior_contract_rho_active(fit) && identical(conditioning, "auto") &&
        !(native || identical(partition, "ordered_k"))) {
      required <- c(required, "rho", "sigma_margin", "pred_marginNormInv")
      if ((intercept_collapsed || relative) &&
          isTRUE(fit@stan_data$use_market_resid_scale == 1L))
        required <- c(required, "sigma_margin_market")
    }
    exact_quality <- .posterior_contract_rho_active(fit) &&
      identical(conditioning, "auto") &&
      (native || identical(partition, "ordered_k")) &&
      (identical(scope, "market") ||
       (identical(scope, "global") &&
        as.integer(fit@stan_data$N_market %||%
          if (!is.null(fit@stan_data$market)) max(fit@stan_data$market) else 1L) == 1L))
    if (exact_quality)
      required <- c(required, "rho", "sigma_margin", "posterior_component_prob",
                    "a_market", "s0_prob")
    if (exact_quality && (intercept_collapsed || relative) &&
        isTRUE(fit@stan_data$use_market_resid_scale == 1L))
      required <- c(required, "sigma_margin_market")
  }
  if (purpose == "posterior_k") {
    if (!native || !identical(partition, "ordered_k")) {
      stop("posterior_k extraction requires an exact ordered-k fit.", call. = FALSE)
    }
    required <- c("posterior_component_prob", "prior_component_prob")
  }
  if (purpose == "cutoff_summary") {
    if (structure %in% c("stackelberg", "fringe") &&
        identical(partition, "ordered_k")) {
      required <- c("posterior_component_prob", "prior_component_prob",
                    "firm_leader_prob", "s0_prob")
      optional <- c(optional, "prior_firm_leader_prob")
      ## A direct ordered-k Stackelberg fit has no cutoff parameter.  The
      ## stackelberg summary only uses cutoff_share for a cutoff-prior fit.
    } else if (identical(structure, "stackelberg")) {
      ## Historical Stackelberg cutoff fits include both the smooth cutoff
      ## quantities and the component/leader probabilities used by the
      ## summary. Preserve that supported contract even though new fitting
      ## APIs prefer direct ordered-k Stackelberg specifications.
      required <- c("posterior_component_prob", "prior_component_prob",
                    "firm_leader_prob", "cutoff_share", "s0_prob")
    } else {
      required <- c("cutoff_share", "s0_prob")
    }
    if (identical(as.character(fit@demand_model), "ces"))
      required <- c(required, "R0_prob")
  }
  if (relative && purpose %in% c("policy_plan", "posterior_quality", "posterior_predictive")) {
    required <- c(required, "s0_prob", "relative_share_predictor")
    if (quality == "sampled") optional <- c(optional, "derived_market_year_intercept")
    if (quality == "collapsed" && isTRUE(fit@stan_data$use_market_resid_scale == 1L))
      required <- c(required, "sigma_share_market")
  }
  if (.posterior_demand_shock_model(fit) && purpose %in%
      c("policy_plan", "posterior_quality", "posterior_predictive")) {
    required <- c(required, "sigma_share", "sigma_demand_shock", "sigma_measurement",
                  "rho_demand_shock")
    if (isTRUE(fit@stan_data$use_market_resid_scale == 1L))
      required <- c(required, "sigma_demand_shock_market", "rho_market",
                    "sigma_share_market", "sigma_margin_market")
    if (!isTRUE(fit@stan_data$likelihood_mode == 2L))
      required <- c(required, "sigma_margin", if (native)
        "posterior_component_prob" else "pred_marginNormInv")
    ## Sampled q still needs the actual role's margin mean when recovering d.
    if (.posterior_contract_rho_active(fit) && quality == "sampled") {
      required <- c(required, "sigma_margin", if (native)
        "posterior_component_prob" else "pred_marginNormInv")
      if (isTRUE(fit@stan_data$use_market_resid_scale == 1L))
        required <- c(required, "sigma_margin_market")
    }
  }
  required <- unique(required)
  optional <- setdiff(unique(optional), required)
  available <- .posterior_fit_draw_names(fit)
  total_fraction_scale <- identical(
    fit@metadata$quality_scale_parameterization_effective %||%
      fit@metadata$quality_scale_parameterization,
    "total_fraction"
  ) || isTRUE(fit@stan_data$quality_scale_parameterization == 2L)
  alternatives <- list(
    s0_prob = c("s0_prob", "s0", "s0_logit"),
    share_mean_integrated = "share_mean_integrated",
    sigma_b_strat = if (total_fraction_scale) {
      "sigma_b_strat_economic"
    } else {
      c("sigma_b_strat_economic", "sigma_b_strat")
    }
  )
  resolve <- function(canonical, optional = FALSE) {
    candidates <- alternatives[[canonical]] %||% canonical
    hit <- candidates[candidates %in% available]
    if (length(hit)) return(hit[[1L]])
    if (optional) return(NA_character_)
    ## Test doubles and a few old serialized objects do not expose model_pars;
    ## request the canonical name and validate the returned list below.
    if (!length(available)) return(candidates[[1L]])
    stop("Estimated fit is missing required posterior draw field(s): ",
         paste(canonical, collapse = ", "), ".", call. = FALSE)
  }
  actual_required <- vapply(required, resolve, character(1), optional = FALSE)
  actual_optional <- if (length(optional)) {
    vapply(optional, resolve, character(1), optional = TRUE)
  } else character()
  keep_optional <- !is.na(actual_optional)
  actual <- c(actual_required, actual_optional[keep_optional])
  canonical <- c(required, optional[keep_optional])
  names(actual) <- canonical
  list(canonical = canonical, actual = unname(actual), aliases = actual_required,
       required = required, optional = optional[keep_optional], purpose = purpose,
       total_fraction_scale = total_fraction_scale,
       market_intercept_collapsed = intercept_collapsed)
}

.posterior_normalize_draw_names <- function(draws, contract) {
  if (!is.list(draws)) stop("Posterior extraction did not return a draw list.", call. = FALSE)
  out <- vector("list", length(contract$canonical))
  names(out) <- contract$canonical
  for (i in seq_along(contract$canonical)) {
    canonical <- contract$canonical[[i]]
    actual <- contract$actual[[i]]
    candidates <- switch(canonical,
      s0_prob = c("s0_prob", "s0", "s0_logit"),
      share_mean_integrated = "share_mean_integrated",
      sigma_b_strat = if (isTRUE(contract$total_fraction_scale)) {
        "sigma_b_strat_economic"
      } else {
        c("sigma_b_strat_economic", "sigma_b_strat")
      },
      canonical
    )
    value <- NULL
    source_name <- NULL
    ## Test doubles sometimes return a legacy alias even when the contract
    ## requested the canonical name. Prefer the requested retained name, then
    ## use a returned alias for the same compatibility behavior as simulate().
    returned <- unique(c(actual, candidates))
    returned <- returned[!is.na(returned) & returned %in% names(draws)]
    for (candidate in returned) {
      if (!is.null(draws[[candidate]])) {
        value <- draws[[candidate]]
        source_name <- candidate
        break
      }
    }
    if (is.null(value) || length(value) == 0L) {
      if (canonical %in% contract$optional) {
        out[i] <- list(NULL)
        next
      }
      stop("Posterior extraction is missing required field(s): ", canonical, ".",
           call. = FALSE)
    }
    ## `s0` and `s0_logit` were historically saved on the logit scale. Only
    ## the alias-to-canonical rename is transformed; native `s0_prob` remains
    ## untouched. plogis preserves vector/matrix dimensions.
    if (identical(canonical, "s0_prob") &&
        source_name %in% c("s0", "s0_logit")) {
      value <- stats::plogis(value)
    }
    out[[i]] <- value
  }
  if (isTRUE(contract$total_fraction_scale) &&
      "sigma_b_strat" %in% names(out)) {
    ## The canonical sigma_b_strat above is already the economic quality SD.
    ## Keep its explicit alias when forwarding this list to quality recovery:
    ## supplied total/fraction draws otherwise interpret that historical name
    ## as the raw total scale and require its fraction logit a second time.
    out$sigma_b_strat_economic <- out$sigma_b_strat
  }
  out
}

.posterior_extract <- function(fit, purpose = "simulation", conditioning = "auto") {
  if (methods::is(fit, "BayesSpecifiedFit")) return(NULL)
  contract <- .posterior_contract(fit, purpose, conditioning)
  raw <- .bayes_backend_extract(fit@stanfit, variables = contract$actual)
  .posterior_normalize_draw_names(raw, contract)
}

.posterior_draw_count_from_draws <- function(draws, fit = NULL) {
  if (!is.list(draws) || !length(draws)) {
    stop("The fitted object contains no posterior draws.", call. = FALSE)
  }
  ## The first dimension of retained `a_market` is the historical draw axis.
  ## Metadata can describe the original fit rather than a selected/subset
  ## draw list, so it must not override the actual array dimensions.
  fields <- unique(c("a_market", names(draws)))
  for (field in fields) {
    x <- draws[[field]]
    if (is.null(x)) next
    d <- dim(x)
    n <- if (is.null(d)) length(x) else d[[1L]]
    if (length(n) == 1L && is.finite(n) && n >= 1L && n == as.integer(n))
      return(as.integer(n))
  }
  stop("The fitted object contains no posterior draws.", call. = FALSE)
}

.posterior_canonical_draws <- function(draws, fit,
                                       purpose = "policy_plan") {
  contract <- .posterior_contract(fit, purpose)
  .posterior_normalize_draw_names(draws, contract)
}
