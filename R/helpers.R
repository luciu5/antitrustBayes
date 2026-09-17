# Internal helper used by both data-preparation routines.
safe_sd <- function(x, floor = 1e-8) {
  s <- stats::sd(x)
  if (!is.finite(s) || s < floor) 1.0 else s
}

## Resolve the statistical likelihood contract shared by bayes_est() and the
## two direct Stan-data preparation entry points.  Codes are deliberately
## stable and are passed as data to Stan: 1 = Student-t, 2 = Gaussian;
## 1 = sampled quality, 2 = collapsed quality.
.resolve_likelihood_modes <- function(share_likelihood = c("student_t", "gaussian"),
                                      quality_parameterization = c("sampled", "collapsed"),
                                      margin_likelihood = NULL,
                                      use_rho = TRUE,
                                      use_market_resid_scale = FALSE,
                                      estimate_df = FALSE,
                                      N_market = NULL,
                                      N_margin = NULL,
                                      likelihood_mode = 1L,
                                      exact_k = FALSE) {
  share_likelihood <- match.arg(share_likelihood)
  quality_parameterization <- match.arg(quality_parameterization)
  if (is.null(margin_likelihood)) margin_likelihood <- share_likelihood
  margin_likelihood <- match.arg(margin_likelihood, c("student_t", "gaussian"))
  if (quality_parameterization == "collapsed" && share_likelihood != "gaussian") {
    stop("quality_parameterization = 'collapsed' requires share_likelihood = 'gaussian'.")
  }
  if (quality_parameterization == "collapsed" && use_market_resid_scale) {
    stop("Collapsed Gaussian quality does not support use_market_resid_scale = TRUE.")
  }
  if (share_likelihood == "student_t" && margin_likelihood == "gaussian") {
    stop("Student-t shares cannot be combined with Gaussian margins.")
  }
  if (share_likelihood == "gaussian" && margin_likelihood == "student_t" &&
      !(quality_parameterization == "collapsed" && !isTRUE(use_rho))) {
    stop("The Gaussian-share/Student-t-margin hybrid requires collapsed quality and use_rho = FALSE.")
  }
  if (isTRUE(estimate_df) && share_likelihood == "gaussian") {
    stop("df = 'estimate' is only available for Student-t shares.")
  }
  if (isTRUE(estimate_df) && margin_likelihood == "gaussian") {
    stop("df = 'estimate' is unavailable with Gaussian margins.")
  }
  if (share_likelihood == "gaussian" && margin_likelihood == "student_t" &&
      isTRUE(estimate_df)) {
    stop("The collapsed Gaussian-share/Student-t-margin hybrid requires fixed df.")
  }
  if (share_likelihood == "gaussian" && margin_likelihood == "student_t" &&
      isTRUE(use_rho)) {
    stop("Gaussian-share/Student-t-margin hybrids require use_rho = FALSE.")
  }
  if (!is.null(N_market) && !is.null(N_margin) &&
      quality_parameterization == "collapsed" && isTRUE(use_rho) &&
      likelihood_mode == 1L && N_market > 1L && N_margin > 0L && exact_k) {
    stop("Collapsed Gaussian rho is currently supported only for one market.")
  }
  list(
    share_likelihood = share_likelihood,
    margin_likelihood = margin_likelihood,
    quality_parameterization = quality_parameterization,
    share_likelihood_code = if (share_likelihood == "student_t") 1L else 2L,
    margin_likelihood_code = if (margin_likelihood == "student_t") 1L else 2L,
    quality_parameterization_code = if (quality_parameterization == "sampled") 1L else 2L,
    mode = paste(share_likelihood, quality_parameterization, margin_likelihood,
                 sep = "/")
  )
}

## Resolve the outside-share pooling ladder.  `s0_pooling` is the public
## interface used by bayes_est(); direct preparation helpers continue to take
## the two legacy booleans and therefore resolve only their explicit mode.
.resolve_s0_pooling <- function(s0_pooling = "auto",
                                s0_pooling_supplied = FALSE,
                                use_s0_market_pool = FALSE,
                                use_s0_market_year_pool = FALSE,
                                legacy_supplied = FALSE,
                                N_market = NULL,
                                N_market_year = NULL) {
  if (!is.logical(s0_pooling_supplied) || length(s0_pooling_supplied) != 1L ||
      is.na(s0_pooling_supplied)) {
    stop("'s0_pooling_supplied' must be a non-missing logical scalar.")
  }
  if (!is.logical(legacy_supplied) || length(legacy_supplied) != 1L ||
      is.na(legacy_supplied)) {
    stop("'legacy_supplied' must be a non-missing logical scalar.")
  }
  if (!is.logical(use_s0_market_pool) || length(use_s0_market_pool) != 1L ||
      is.na(use_s0_market_pool) ||
      !is.logical(use_s0_market_year_pool) ||
      length(use_s0_market_year_pool) != 1L ||
      is.na(use_s0_market_year_pool)) {
    stop("The legacy outside-share pooling flags must be non-missing logical scalars.")
  }
  if (isTRUE(s0_pooling_supplied) && isTRUE(legacy_supplied)) {
    stop("'s0_pooling' cannot be supplied together with either legacy outside-share pooling flag.")
  }

  if (isTRUE(s0_pooling_supplied)) {
    requested <- match.arg(s0_pooling, c("auto", "none", "market", "market_year"))
    source <- "s0_pooling"
  } else if (isTRUE(legacy_supplied)) {
    if (use_s0_market_pool && use_s0_market_year_pool) {
      stop("Use at most one of 'use_s0_market_pool' (one outside share per ",
           "market) and 'use_s0_market_year_pool' (one per market-year).")
    }
    requested <- if (use_s0_market_year_pool) {
      "market_year"
    } else if (use_s0_market_pool) {
      "market"
    } else {
      "none"
    }
    source <- "legacy_flags"
  } else {
    requested <- "auto"
    source <- "auto_data"
  }

  if (identical(requested, "auto")) {
    if (is.null(N_market) || is.null(N_market_year)) {
      resolved <- "auto"
    } else if (N_market == 1L && N_market_year == 1L) {
      resolved <- "none"
    } else if (N_market == 1L) {
      resolved <- "market"
    } else {
      resolved <- "market"
    }
  } else {
    resolved <- requested
  }

  list(
    requested = requested,
    resolved = resolved,
    source = source,
    use_s0_market_pool = identical(resolved, "market"),
    use_s0_market_year_pool = identical(resolved, "market_year")
  )
}

## Resolve the market/outside-share coordinate convention.  The Stan data code
## is deliberately stable: 1 = legacy coordinates, 2 = rotated coordinates.
## The affine rotation is only active when there is more than one market;
## single-market fits retain their historical normalization and dimensions.
.resolve_market_s0_parameterization <- function(market_s0_parameterization =
                                                 c("legacy", "rotated"),
                                                 N_market = NULL) {
  market_s0_parameterization <- match.arg(
    market_s0_parameterization, c("legacy", "rotated")
  )
  if (!is.null(N_market) &&
      (!is.numeric(N_market) || length(N_market) != 1L ||
       !is.finite(N_market) || N_market < 1)) {
    stop("'N_market' must be a positive numeric scalar.")
  }
  active <- identical(market_s0_parameterization, "rotated") &&
    (is.null(N_market) || N_market > 1)
  list(
    requested = market_s0_parameterization,
    code = if (identical(market_s0_parameterization, "legacy")) 1L else 2L,
    effective = if (isTRUE(active)) "rotated" else "legacy",
    active = isTRUE(active),
    coordinate_convention = if (isTRUE(active)) {
      "rotated_market_outside_share"
    } else {
      "legacy_market_outside_share"
    }
  )
}

bounded01 <- function(x, eps = 1e-6) {
  pmin(pmax(x, eps), 1 - eps)
}

# Resolve a scalar or market-year vector of outside-good prices.
resolve_outside_price <- function(outside_price, N_market_year) {
  if (!is.numeric(outside_price) ||
      !length(outside_price) %in% c(1L, N_market_year) ||
      any(!is.finite(outside_price)) || any(outside_price <= 0)) {
    stop("'outside_price' must be a positive finite scalar or have one value per market-year.")
  }
  if (length(outside_price) == 1L) {
    rep(as.numeric(outside_price), N_market_year)
  } else {
    as.numeric(outside_price)
  }
}

validate_priors <- function(priors) {
  if (!is.list(priors) || is.null(names(priors))) {
    stop("'priors' must be a named list.")
  }
  for (nm in names(priors)) {
    value <- priors[[nm]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
      stop("Prior '", nm, "' must be a finite numeric scalar.")
    }
  }
  positive <- setdiff(
    names(priors),
    c("prior_alpha_mean", "s0_prior_mean", "bargaining_prior_lower")
  )
  bad <- positive[vapply(priors[positive], function(x) x <= 0, logical(1))]
  if (length(bad)) {
    stop("Prior scale/shape parameters must be strictly positive: ",
         paste(bad, collapse = ", "), ".")
  }
  if (all(c("bargaining_prior_lower", "bargaining_prior_upper") %in%
          names(priors))) {
    if (priors$bargaining_prior_lower < 0 ||
        priors$bargaining_prior_upper >= 1 ||
        priors$bargaining_prior_lower >= priors$bargaining_prior_upper) {
      stop("Bargaining prior bounds must define an increasing interval in [0, 1).")
    }
  }
  if ("resid_df_fixed" %in% names(priors) && priors$resid_df_fixed <= 2) {
    stop("'resid_df_fixed' must exceed 2.")
  }
  invisible(TRUE)
}

compositional_share_data <- function(shares, mky_int, N_market_year) {
  log_share <- log(pmax(shares, 1e-12))
  clr_share <- numeric(length(shares))
  for (m in seq_len(N_market_year)) {
    idx <- which(mky_int == m)
    clr_share[idx] <- log_share[idx] - mean(log_share[idx])
  }
  obs_order <- order(mky_int, seq_along(mky_int))
  mky_n <- tabulate(mky_int, nbins = N_market_year)
  mky_start <- cumsum(c(1L, utils::head(mky_n, -1L)))
  list(
    log_share = log_share,
    clr_share = clr_share,
    obs_order = as.integer(obs_order),
    mky_n = as.integer(mky_n),
    mky_start = as.integer(mky_start)
  )
}

#' Convert Positive Draws to Log-Normal Prior Parameters
#'
#' Computes prior mean and SD from FOC-implied parameter draws. By default uses
#' weighted mean and SD when weights are provided (emphasizing large products),
#' or robust quantiles (10th and 90th percentiles) when unweighted.
#'
#' @param x Numeric vector of positive draws (e.g., FOC-implied alphas).
#' @param weights Optional numeric vector of non-negative weights (e.g., market shares).
#'   If provided, uses weighted mean and SD directly. Must have same length as x
#'   after filtering to valid values - no vector recycling. Default NULL uses
#'   robust quantiles.
#' @param fallback_mean Fallback prior mean if draws are insufficient.
#' @param fallback_sd Fallback prior SD if draws are insufficient.
#' @return A list with elements \code{mean} and \code{sd} on log scale.
#' @keywords internal
log_prior_from_positive_draws <- function(x, weights = NULL,
                                          fallback_mean = 1.0,
                                          fallback_sd = 1.0) {
  # Filter to valid x first
  valid <- is.finite(x) & x > 0
  x <- x[valid]

  if (length(x) < 2L) {
    return(list(mean = log(fallback_mean), sd = fallback_sd))
  }

  # Handle and validate weights
  if (!is.null(weights)) {
    # Filter weights to match filtered x
    weights <- weights[valid]

    # Validation: length must match (no recycling)
    if (length(weights) != length(x)) {
      stop("length(weights) must equal length(x) after filtering. ",
           "Got ", length(weights), " weights and ", length(x), " valid x values.")
    }

    # Validation: must be non-negative
    if (any(weights < 0, na.rm = TRUE)) {
      stop("weights must be non-negative. Found negative values.")
    }

    # Normalize weights to sum to 1
    weights <- pmax(weights, 0)
    weight_sum <- sum(weights, na.rm = TRUE)
    if (weight_sum > 0) {
      weights <- weights / weight_sum
    } else {
      # All weights are zero - fall back to unweighted
      weights <- NULL
    }
  }

  # Compute weighted mean and SD if weights provided
  if (!is.null(weights)) {
    log_x <- log(x)
    mean_val <- sum(log_x * weights, na.rm = TRUE)

    # Weighted variance
    var_val <- sum(weights * (log_x - mean_val)^2, na.rm = TRUE)
    sd_val <- max(sqrt(var_val), 0.15)

    return(list(mean = mean_val, sd = sd_val))
  }

  # Unweighted fallback: use robust quantiles (original behavior)
  # Only used when weights explicitly NULL or all zero
  qs <- stats::quantile(log(x), c(0.10, 0.90), na.rm = TRUE, names = FALSE)
  if (!all(is.finite(qs)) || qs[2] <= qs[1]) {
    return(list(mean = mean(log(x), na.rm = TRUE), sd = fallback_sd))
  }

  list(
    mean = mean(qs),
    sd = max((qs[2] - qs[1]) / 2.56, 0.15)
  )
}

#' Prepare Stan Data for Logit Demand
#'
#' Transforms user-supplied vectors into the list expected by the Logit Stan
#' model. Computes conditional quantity shares and inverse normalized margins.
#'
#' @param prices Numeric vector of product prices.
#' @param quantities Numeric vector of product quantities.
#' @param margins Numeric vector of positive price-normalized margins. Values
#'   greater than 1 are valid for input-market markdowns.
#' @param market Factor or vector identifying markets.
#' @param year Factor or vector identifying time periods.
#' @param group Factor or vector identifying product groups (e.g., firms).
#' @param supply_model Integer 1--5.
#' @param market_side Character: \code{"input"} or \code{"output"}.
#' @param use_cutoff Logical.
#' @param use_rho Logical.
#' @param use_market_resid_scale Logical.
#' @param use_alpha_market_pool Logical.
#' @param use_s0_market_year_pool Logical.
#' @param use_s0_market_pool Logical. If true, use one latent outside quantity
#'   share per market before any market-year \code{X_s0} adjustment.
#'   The preparation helpers retain these explicit legacy controls; resolved
#'   pooling mode is recorded in the returned list attributes.
#' @param use_hmt Logical.
#' @param use_centered_alpha Logical.
#' @param use_centered_s0_pool Logical. Select the centered parameterization
#'   for either outside-share pooling mode.
#' @param market_s0_parameterization Character string selecting the legacy or
#'   market/outside-share rotated Stan coordinates. The rotation is a no-op for
#'   single-market data.
#' @param priors Named list of prior hyperparameters.
#' @param X_s0 Optional model matrix for s0 covariates.
#' @param outside_price Positive scalar or one value per market-year.
#' @param grainsize Positive integer reduce-sum chunk size.
#' @param margin_observed Optional logical margin-observation mask.
#' @param likelihood_mode Integer: 1 joint, 2 share only, 3 margin only.
#' @param loo_cluster Optional natural predictive-cluster identifier.
#' @param estimate_df Logical; estimate the Student-t degrees of freedom.
#' @param resid_df_fixed Fixed Student-t degrees of freedom when
#'   \code{estimate_df = FALSE}.
#' @param k_support Explicit unique integer leader counts. Fringe supports may
#'   include zero; Stackelberg supports must be positive.
#' @param k_prior Positional probabilities corresponding to \code{k_support};
#'   \code{NULL} gives a uniform prior.
#' @return A list suitable for passing to \code{rstan::sampling}.
#' @keywords internal
prepare_stan_data <- function(prices, quantities, margins, market, year, group,
                              supply_model, use_cutoff = FALSE,
                              market_side = "input",
                              use_rho = TRUE,
                              use_market_resid_scale = FALSE,
                              use_hmt = FALSE,
                              use_alpha_market_pool = TRUE,
                              use_s0_market_year_pool = FALSE,
                              use_s0_market_pool = FALSE,
                              use_centered_alpha = FALSE,
                              use_centered_s0_pool = FALSE,
                              market_s0_parameterization = c("legacy", "rotated"),
                              priors = list(), X_s0 = NULL,
                              grainsize = 1L,
                              margin_observed = NULL,
                              outside_price = 1,
                              likelihood_mode = 1L,
                              loo_cluster = market,
                              estimate_df = FALSE,
                              resid_df_fixed = NULL,
                              structure = NULL,
                              cutoff_scale = NULL,
                              partition_prior = "cutoff",
                              k_support = NULL,
                              k_prior = NULL,
                              k_scope = "market",
                              share_likelihood = c("student_t", "gaussian"),
                              quality_parameterization = c("sampled", "collapsed"),
                              margin_likelihood = NULL) {
  N <- length(prices)
  likelihood_modes <- .resolve_likelihood_modes(
    share_likelihood = share_likelihood,
    quality_parameterization = quality_parameterization,
    margin_likelihood = margin_likelihood,
    use_rho = use_rho,
    use_market_resid_scale = use_market_resid_scale,
    estimate_df = estimate_df,
    likelihood_mode = likelihood_mode
  )
  share_likelihood <- likelihood_modes$share_likelihood
  margin_likelihood <- likelihood_modes$margin_likelihood
  quality_parameterization <- likelihood_modes$quality_parameterization
  market_side <- match.arg(market_side, c("input", "output"))
  market_sign <- if (market_side == "input") 1L else -1L
  resolved_structure <- .stack_resolve_structure(
    structure = structure, use_cutoff = use_cutoff,
    use_cutoff_supplied = FALSE
  )
  partition_prior_supplied <- !missing(partition_prior)
  exact_structure <- resolved_structure$supplied &&
    resolved_structure$name %in% c("fringe", "stackelberg")
  if (exact_structure) {
    if (!partition_prior_supplied) partition_prior <- "ordered_k"
    partition_prior <- match.arg(partition_prior, c("cutoff", "ordered_k"))
    if (resolved_structure$name == "stackelberg" && partition_prior != "ordered_k") {
      stop("Explicit Stackelberg requires partition_prior = 'ordered_k'.")
    }
    if (partition_prior == "ordered_k") {
      if (is.null(k_support)) stop("Explicit ordered-k structures require 'k_support'.")
      if (!identical(k_scope, "market")) {
        stop("Exact ordered-k uses one partition per market; 'k_scope' must be 'market'.")
      }
      if (is.list(k_prior)) {
        stop("Exact ordered-k requires a positional numeric 'k_prior' vector.")
      }
    } else if (!is.null(k_support) || !is.null(k_prior)) {
      stop("'k_support' and 'k_prior' require partition_prior = 'ordered_k'.")
    }
  } else {
    ## Omitted structure retains legacy simultaneous/smooth-fringe behavior.
    partition_prior <- "cutoff"
    k_support <- NULL
    k_prior <- NULL
  }
  k_scope <- match.arg(k_scope, c("market", "period"))
  .stack_validate_public_k(k_support, k_prior, resolved_structure$name)
  if (resolved_structure$name %in% c("stackelberg", "fringe") &&
      partition_prior == "ordered_k" && !(supply_model %in% c(1L, 3L))) {
    stop("Exact ordered-k fringe/Stackelberg supports only Bertrand and Cournot supply models.")
  }
  cutoff_scale_supplied <- !missing(cutoff_scale) && !is.null(cutoff_scale)
  if (identical(resolved_structure$name, "simultaneous") && cutoff_scale_supplied) {
    stop("'cutoff_scale' is only defined for fringe and Stackelberg structures.")
  }
  if (identical(resolved_structure$name, "stackelberg") &&
      identical(partition_prior, "ordered_k") && cutoff_scale_supplied) {
    stop("'cutoff_scale' is not defined for partition_prior = 'ordered_k'.")
  }
  if (exact_structure && identical(partition_prior, "ordered_k") &&
      cutoff_scale_supplied) {
    stop("'cutoff_scale' is not defined for exact ordered-k structures.")
  }
  if (exact_structure && identical(partition_prior, "ordered_k")) {
    cutoff_scale <- NULL
  } else {
    cutoff_scale <- .stack_resolve_scale(resolved_structure$name, cutoff_scale)
  }

  if (is.null(margin_observed)) {
    margin_observed <- rep(TRUE, N)
  }
  margin_observed <- as.logical(margin_observed) & !is.na(margin_observed)

  # Encode factors
  market_f <- droplevels(as.factor(market))
  year_f <- droplevels(as.factor(year))
  group_f <- droplevels(as.factor(group))

  market_int <- as.integer(market_f)
  year_int <- as.integer(year_f)
  group_int <- as.integer(group_f)

  N_market <- nlevels(market_f)
  N_year <- nlevels(year_f)
  N_group <- nlevels(group_f)

  is_single_market <- as.integer(N_market == 1L)

  # Market-year indexing
  mky_labels <- paste(market_int, year_int, sep = "_")
  mky_f <- as.factor(mky_labels)
  mky_int <- as.integer(mky_f)
  N_market_year <- nlevels(mky_f)
  market_s0_resolution <- .resolve_market_s0_parameterization(
    market_s0_parameterization, N_market
  )
  s0_resolution <- .resolve_s0_pooling(
    s0_pooling_supplied = FALSE,
    use_s0_market_pool = use_s0_market_pool,
    use_s0_market_year_pool = use_s0_market_year_pool,
    legacy_supplied = TRUE,
    N_market = N_market,
    N_market_year = N_market_year
  )
  use_s0_market_pool <- s0_resolution$use_s0_market_pool
  use_s0_market_year_pool <- s0_resolution$use_s0_market_year_pool
  outside_price <- resolve_outside_price(outside_price, N_market_year)
  cluster_f <- droplevels(as.factor(loo_cluster))
  if (length(cluster_f) != N || anyNA(cluster_f)) {
    stop("'loo_cluster' must have one non-missing value per observation.")
  }
  cluster_int <- as.integer(cluster_f)
  N_cluster <- nlevels(cluster_f)
  margin_idx <- which(margin_observed)
  N_margin <- length(margin_idx)
  if (quality_parameterization == "collapsed" && isTRUE(use_rho) &&
      likelihood_mode == 1L && N_market > 1L && N_margin > 0L) {
    .resolve_likelihood_modes(
      share_likelihood = likelihood_modes$share_likelihood,
      quality_parameterization = likelihood_modes$quality_parameterization,
      margin_likelihood = likelihood_modes$margin_likelihood,
      use_rho = use_rho,
      use_market_resid_scale = use_market_resid_scale,
      estimate_df = estimate_df,
      N_market = N_market, N_margin = N_margin,
      likelihood_mode = likelihood_mode,
      exact_k = exact_structure && identical(partition_prior, "ordered_k")
    )
  }
  if (any(vapply(split(cluster_int, mky_int), function(x) length(unique(x)) > 1L,
                 logical(1)))) {
    stop("'loo_cluster' must be constant within each market-year equilibrium.")
  }

  # Build market_year -> market and market_year -> year maps
  mky_to_market <- integer(N_market_year)
  mky_to_year <- integer(N_market_year)
  for (i in seq_len(N)) {
    mky_to_market[mky_int[i]] <- market_int[i]
    mky_to_year[mky_int[i]] <- year_int[i]
  }

  # Compute quantity shares within each market-year
  shares <- numeric(N)
  for (m in seq_len(N_market_year)) {
    idx <- which(mky_int == m)
    total_q <- sum(quantities[idx])
    shares[idx] <- quantities[idx] / total_q
  }

  # Strategic supply FOCs use the within-market-year share of the firm,
  # not the share of an individual product. This is identical to shares for
  # single-product firms and preserves the multi-product ownership logic.
  firm_shareIn <- numeric(N)
  for (m in seq_len(N_market_year)) {
    idx <- which(mky_int == m)
    for (g in unique(group_int[idx])) {
      firm_idx <- idx[group_int[idx] == g]
      firm_shareIn[firm_idx] <- sum(shares[firm_idx])
    }
  }

  # Standardize prices
  price_sd <- safe_sd(prices)
  price_std <- (prices - mean(prices)) / price_sd

  # Put both demand systems on the same inverse normalized-margin response.
  # Masked rows
  # (margin_observed == FALSE) may carry NA/invalid margins -- Stan never
  # reads their margin likelihood term, but the vectors must still be finite,
  # so substitute a dummy placeholder margin of 1 before transforming it.
  margins_for_stan <- margins
  margins_for_stan[!margin_observed] <- 1.0
  marginNormInv <- 1.0 / margins_for_stan
  share_comp <- compositional_share_data(shares, mky_int, N_market_year)
  log_shareIn <- share_comp$log_share
  clr_shareIn <- share_comp$clr_share
  log_marginNormInv <- log(pmax(marginNormInv, 1e-12))

  # S0 covariates
  if (is.null(X_s0)) {
    K_s0 <- 0L
    X_s0_mat <- matrix(0, nrow = N_market_year, ncol = 0)
  } else {
    X_s0_mat <- as.matrix(X_s0)
    if (nrow(X_s0_mat) != N_market_year) {
      stop("'X_s0' must have one row per market-year.")
    }
    if (!is.numeric(X_s0_mat) || any(!is.finite(X_s0_mat))) {
      stop("'X_s0' must be a finite numeric matrix.")
    }
    K_s0 <- ncol(X_s0_mat)
  }

  # Default outside-share prior is logistic-normal with arithmetic mean 20%
  # and 5th percentile 10%: location = -1.4467864, scale = 0.4562340.
  default_priors <- bayes_prior_manifest("baseline")
  # Override with user priors
  for (nm in names(priors)) {
    if (nm %in% names(default_priors)) {
      default_priors[[nm]] <- priors[[nm]]
    } else {
      warning("Unknown prior parameter '", nm, "' ignored.")
    }
  }

  s0_prior_reference <- logit_normal_mean(
    default_priors$s0_prior_mean,
    default_priors$s0_prior_sd
  )
  if (!is.null(resid_df_fixed)) default_priors$resid_df_fixed <- resid_df_fixed
  validate_priors(default_priors)
  if (!likelihood_mode %in% 1:3) {
    stop("'likelihood_mode' must be 1 (joint), 2 (share), or 3 (margin).")
  }

  # HMT data use a dimensionless normalized margin, matching the elasticity.
  avg_price_hmt <- mean(prices)
  avg_margin_hmt <- if (any(margin_observed)) {
    mean(margins[margin_observed])
  } else {
    0.0
  }
  ssnip_hmt <- 0.05

  stan_data <- c(
    list(
      N = N,
      shareIn = array(shares, dim = N),
      firm_shareIn = array(firm_shareIn, dim = N),
      marginNormInv = array(marginNormInv, dim = N),
      log_shareIn = array(log_shareIn, dim = N),
      clr_shareIn = array(clr_shareIn, dim = N),
      log_marginNormInv = array(log_marginNormInv, dim = N),
      price = array(as.numeric(prices), dim = N),
      price_std = array(price_std, dim = N),
      price_sd = price_sd,
      outside_price = array(outside_price, dim = N_market_year),
      market_sign = market_sign,
      N_market = N_market,
      market = array(market_int, dim = N),
      N_group = N_group,
      group = array(group_int, dim = N),
      N_year = N_year,
      year = array(year_int, dim = N),
      supply_model = supply_model,
      use_cutoff = as.integer(use_cutoff),
      cutoff_scale = cutoff_scale,
      grainsize = as.integer(max(1L, grainsize)),
      N_market_year = N_market_year,
      market_year_idx = array(mky_int, dim = N),
      mky_n = array(share_comp$mky_n, dim = N_market_year),
      mky_start = array(share_comp$mky_start, dim = N_market_year),
      obs_order = array(share_comp$obs_order, dim = N),
      N_cluster = N_cluster,
      cluster = array(cluster_int, dim = N),
      likelihood_mode = as.integer(likelihood_mode),
      margin_obs = array(as.integer(margin_observed), dim = N),
      N_margin = as.integer(N_margin),
      margin_idx = array(as.integer(margin_idx), dim = N_margin),
      share_likelihood = likelihood_modes$share_likelihood_code,
      margin_likelihood = likelihood_modes$margin_likelihood_code,
      quality_parameterization = likelihood_modes$quality_parameterization_code,
      K_s0 = K_s0,
      X_s0 = X_s0_mat,
      mky_to_market = array(mky_to_market),
      mky_to_year = array(mky_to_year),
      is_single_market = is_single_market,
      use_rho = as.integer(use_rho),
      use_market_resid_scale = as.integer(use_market_resid_scale),
      use_alpha_market_pool = as.integer(use_alpha_market_pool),
      use_s0_market_year_pool = as.integer(use_s0_market_year_pool),
      use_s0_market_pool = as.integer(use_s0_market_pool),
      market_s0_parameterization = market_s0_resolution$code,
      use_centered_alpha = as.integer(use_centered_alpha),
      use_centered_s0_pool = as.integer(use_centered_s0_pool),
      use_hmt = as.integer(use_hmt),
      avg_price_hmt = avg_price_hmt,
      avg_margin_hmt = avg_margin_hmt,
      ssnip_hmt = ssnip_hmt,
      estimate_df = as.integer(estimate_df)
    ),
    default_priors
  )

  # Attach metadata for later use
  attr(stan_data, "market_levels") <- levels(market_f)
  attr(stan_data, "year_levels") <- levels(year_f)
  attr(stan_data, "group_levels") <- levels(group_f)
  attr(stan_data, "mky_levels") <- levels(mky_f)
  attr(stan_data, "cluster_levels") <- levels(cluster_f)
  attr(stan_data, "market_year_idx") <- mky_int
  attr(stan_data, "model_spec_version") <- .antitrust_bayes_model_spec_version
  attr(stan_data, "s0_prior_reference") <- s0_prior_reference
  attr(stan_data, "structure") <- resolved_structure$name
  attr(stan_data, "structure_code") <- resolved_structure$code
  attr(stan_data, "cutoff_scale") <- cutoff_scale
  attr(stan_data, "partition_prior") <- partition_prior
  attr(stan_data, "k_scope") <- k_scope
  attr(stan_data, "k_support") <- k_support
  attr(stan_data, "k_prior") <- k_prior
  attr(stan_data, "share_likelihood") <- likelihood_modes$share_likelihood
  attr(stan_data, "margin_likelihood") <- likelihood_modes$margin_likelihood
  attr(stan_data, "quality_parameterization") <- likelihood_modes$quality_parameterization
  attr(stan_data, "likelihood_mode_metadata") <- likelihood_modes
  attr(stan_data, "s0_pooling") <- s0_resolution$resolved
  attr(stan_data, "s0_pooling_requested") <- s0_resolution$requested
  attr(stan_data, "s0_pooling_resolved") <- s0_resolution$resolved
  attr(stan_data, "s0_pooling_source") <- s0_resolution$source
  attr(stan_data, "market_s0_parameterization") <-
    market_s0_resolution$requested
  attr(stan_data, "market_s0_parameterization_requested") <-
    market_s0_resolution$requested
  attr(stan_data, "market_s0_parameterization_effective") <-
    market_s0_resolution$effective
  attr(stan_data, "market_s0_rotation_active") <- market_s0_resolution$active
  attr(stan_data, "market_s0_coordinate_convention") <-
    market_s0_resolution$coordinate_convention
  if (resolved_structure$name == "stackelberg" ||
      (exact_structure && identical(resolved_structure$name, "fringe") &&
       identical(partition_prior, "ordered_k"))) {
    stack_data <- .stack_append_data(
      stan_data = stan_data, mky_int = mky_int, market_int = market_int,
      group_int = group_int, group_levels = levels(group_f), shares = shares,
      revShare = shares, structure = resolved_structure,
      demand_model = "logit", partition_prior = partition_prior,
      k_support = k_support, k_prior = k_prior, k_scope = k_scope,
      cutoff_scale = cutoff_scale
    )
    ranking_metadata <- attr(stack_data, "firm_market_ranking", exact = TRUE)
    partition_metadata <- attr(stack_data, "partition_metadata", exact = TRUE)
    stan_data <- c(stan_data, stack_data)
    if (!is.null(cutoff_scale)) stan_data$cutoff_scale <- cutoff_scale
    attr(stan_data, "market_levels") <- levels(market_f)
    attr(stan_data, "year_levels") <- levels(year_f)
    attr(stan_data, "group_levels") <- levels(group_f)
    attr(stan_data, "mky_levels") <- levels(mky_f)
    attr(stan_data, "cluster_levels") <- levels(cluster_f)
    attr(stan_data, "market_year_idx") <- mky_int
    attr(stan_data, "model_spec_version") <- .antitrust_bayes_model_spec_version
    attr(stan_data, "s0_prior_reference") <- s0_prior_reference
    attr(stan_data, "structure") <- resolved_structure$name
    attr(stan_data, "structure_code") <- resolved_structure$code
    attr(stan_data, "cutoff_scale") <- cutoff_scale
    attr(stan_data, "partition_prior") <- partition_prior
    attr(stan_data, "k_scope") <- k_scope
    attr(stan_data, "k_support") <- k_support
    attr(stan_data, "k_prior") <- k_prior
    attr(stan_data, "share_likelihood") <- likelihood_modes$share_likelihood
    attr(stan_data, "margin_likelihood") <- likelihood_modes$margin_likelihood
    attr(stan_data, "quality_parameterization") <- likelihood_modes$quality_parameterization
    attr(stan_data, "likelihood_mode_metadata") <- likelihood_modes
    attr(stan_data, "s0_pooling") <- s0_resolution$resolved
    attr(stan_data, "s0_pooling_requested") <- s0_resolution$requested
    attr(stan_data, "s0_pooling_resolved") <- s0_resolution$resolved
    attr(stan_data, "s0_pooling_source") <- s0_resolution$source
    attr(stan_data, "market_s0_parameterization") <-
      market_s0_resolution$requested
    attr(stan_data, "market_s0_parameterization_requested") <-
      market_s0_resolution$requested
    attr(stan_data, "market_s0_parameterization_effective") <-
      market_s0_resolution$effective
    attr(stan_data, "market_s0_rotation_active") <- market_s0_resolution$active
    attr(stan_data, "market_s0_coordinate_convention") <-
      market_s0_resolution$coordinate_convention
    attr(stan_data, "firm_market_ranking") <- ranking_metadata
    attr(stan_data, "partition_metadata") <- partition_metadata
  }

  stan_data
}

#' Prepare Stan Data for CES Demand
#'
#' Transforms user-supplied vectors into the list expected by the CES Stan
#' model. The likelihood uses the same conditional quantity-share and inverse
#' normalized-margin outcomes as the Logit model; revenue shares are retained
#' for the CES supply first-order conditions.
#'
#' @inheritParams prepare_stan_data
#' @return A list suitable for passing to \code{rstan::sampling}.
#' @keywords internal
prepare_ces_stan_data <- function(prices, quantities, margins, market, year,
                                  group, supply_model, use_cutoff, use_rho,
                                  use_hmt, priors, X_s0 = NULL,
                                  market_side = "output",
                                  grainsize = 1L,
                                  use_market_resid_scale = FALSE,
                                  use_alpha_market_pool = FALSE,
                                  use_s0_market_year_pool = FALSE,
                                  use_s0_market_pool = FALSE,
                                  use_centered_alpha = FALSE,
                                  use_centered_s0_pool = FALSE,
                                  market_s0_parameterization = c("legacy", "rotated"),
                                  margin_observed = NULL,
                                  outside_price = 1,
                                  likelihood_mode = 1L,
                                  loo_cluster = market,
                                  estimate_df = FALSE,
                                  resid_df_fixed = NULL,
                                  structure = NULL,
                                  cutoff_scale = NULL,
                                  partition_prior = "cutoff",
                                  k_support = NULL,
                                  k_prior = NULL,
                                  k_scope = "market",
                                  share_likelihood = c("student_t", "gaussian"),
                                  quality_parameterization = c("sampled", "collapsed"),
                                  margin_likelihood = NULL) {
  N <- length(prices)
  likelihood_modes <- .resolve_likelihood_modes(
    share_likelihood = share_likelihood,
    quality_parameterization = quality_parameterization,
    margin_likelihood = margin_likelihood,
    use_rho = use_rho,
    use_market_resid_scale = use_market_resid_scale,
    estimate_df = estimate_df,
    likelihood_mode = likelihood_mode
  )
  share_likelihood <- likelihood_modes$share_likelihood
  margin_likelihood <- likelihood_modes$margin_likelihood
  quality_parameterization <- likelihood_modes$quality_parameterization
  market_side <- match.arg(market_side, c("input", "output"))
  if (market_side != "output") {
    stop("CES estimation currently supports output markets only; use ",
         "demand_model = \"logit\" for input markets.")
  }
  market_sign <- -1L
  resolved_structure <- .stack_resolve_structure(
    structure = structure, use_cutoff = use_cutoff,
    use_cutoff_supplied = FALSE
  )
  partition_prior_supplied <- !missing(partition_prior)
  exact_structure <- resolved_structure$supplied &&
    resolved_structure$name %in% c("fringe", "stackelberg")
  if (exact_structure) {
    if (!partition_prior_supplied) partition_prior <- "ordered_k"
    partition_prior <- match.arg(partition_prior, c("cutoff", "ordered_k"))
    if (resolved_structure$name == "stackelberg" && partition_prior != "ordered_k") {
      stop("Explicit Stackelberg requires partition_prior = 'ordered_k'.")
    }
    if (partition_prior == "ordered_k") {
      if (is.null(k_support)) stop("Explicit ordered-k structures require 'k_support'.")
      if (!identical(k_scope, "market")) {
        stop("Exact ordered-k uses one partition per market; 'k_scope' must be 'market'.")
      }
      if (is.list(k_prior)) {
        stop("Exact ordered-k requires a positional numeric 'k_prior' vector.")
      }
    } else if (!is.null(k_support) || !is.null(k_prior)) {
      stop("'k_support' and 'k_prior' require partition_prior = 'ordered_k'.")
    }
  } else {
    ## Omitted structure retains legacy simultaneous/smooth-fringe behavior.
    partition_prior <- "cutoff"
    k_support <- NULL
    k_prior <- NULL
  }
  k_scope <- match.arg(k_scope, c("market", "period"))
  .stack_validate_public_k(k_support, k_prior, resolved_structure$name)
  if (resolved_structure$name %in% c("stackelberg", "fringe") &&
      partition_prior == "ordered_k" && !(supply_model %in% c(1L, 3L))) {
    stop("Exact ordered-k fringe/Stackelberg supports only Bertrand and Cournot supply models.")
  }
  cutoff_scale_supplied <- !missing(cutoff_scale) && !is.null(cutoff_scale)
  if (identical(resolved_structure$name, "simultaneous") && cutoff_scale_supplied) {
    stop("'cutoff_scale' is only defined for fringe and Stackelberg structures.")
  }
  if (identical(resolved_structure$name, "stackelberg") &&
      identical(partition_prior, "ordered_k") && cutoff_scale_supplied) {
    stop("'cutoff_scale' is not defined for partition_prior = 'ordered_k'.")
  }
  if (exact_structure && identical(partition_prior, "ordered_k") &&
      cutoff_scale_supplied) {
    stop("'cutoff_scale' is not defined for exact ordered-k structures.")
  }
  if (exact_structure && identical(partition_prior, "ordered_k")) {
    cutoff_scale <- NULL
  } else {
    cutoff_scale <- .stack_resolve_scale(resolved_structure$name, cutoff_scale)
  }

  if (is.null(margin_observed)) {
    margin_observed <- rep(TRUE, N)
  }
  margin_observed <- as.logical(margin_observed) & !is.na(margin_observed)

  # Masked rows may carry NA/invalid margins -- their likelihood term is
  # never read by Stan, so substitute a dummy placeholder before validating
  # and computing the inverse normalized margin.
  margins[!margin_observed] <- 0.5

  # CES FOCs require L < 1 (positive marginal costs), checked on observed rows.
  if (any(margins[margin_observed] >= 1)) {
    stop("CES demand requires Lerner indices strictly less than 1. ",
         "Found ", sum(margins[margin_observed] >= 1), " observation(s) with L >= 1. ",
         "Use demand_model = \"logit\" for input markets with L >= 1.")
  }

  # Encode factors
  market_f <- droplevels(as.factor(market))
  year_f <- droplevels(as.factor(year))
  group_f <- droplevels(as.factor(group))

  market_int <- as.integer(market_f)
  year_int <- as.integer(year_f)
  group_int <- as.integer(group_f)

  N_market <- nlevels(market_f)
  N_year <- nlevels(year_f)
  N_group <- nlevels(group_f)

  is_single_market <- as.integer(N_market == 1L)

  # Market-year indexing
  mky_labels <- paste(market_int, year_int, sep = "_")
  mky_f <- as.factor(mky_labels)
  mky_int <- as.integer(mky_f)
  N_market_year <- nlevels(mky_f)
  market_s0_resolution <- .resolve_market_s0_parameterization(
    market_s0_parameterization, N_market
  )
  s0_resolution <- .resolve_s0_pooling(
    s0_pooling_supplied = FALSE,
    use_s0_market_pool = use_s0_market_pool,
    use_s0_market_year_pool = use_s0_market_year_pool,
    legacy_supplied = TRUE,
    N_market = N_market,
    N_market_year = N_market_year
  )
  use_s0_market_pool <- s0_resolution$use_s0_market_pool
  use_s0_market_year_pool <- s0_resolution$use_s0_market_year_pool
  outside_price <- resolve_outside_price(outside_price, N_market_year)
  cluster_f <- droplevels(as.factor(loo_cluster))
  if (length(cluster_f) != N || anyNA(cluster_f)) {
    stop("'loo_cluster' must have one non-missing value per observation.")
  }
  cluster_int <- as.integer(cluster_f)
  N_cluster <- nlevels(cluster_f)
  margin_idx <- which(margin_observed)
  N_margin <- length(margin_idx)
  if (quality_parameterization == "collapsed" && isTRUE(use_rho) &&
      likelihood_mode == 1L && N_market > 1L && N_margin > 0L) {
    .resolve_likelihood_modes(
      share_likelihood = likelihood_modes$share_likelihood,
      quality_parameterization = likelihood_modes$quality_parameterization,
      margin_likelihood = likelihood_modes$margin_likelihood,
      use_rho = use_rho,
      use_market_resid_scale = use_market_resid_scale,
      estimate_df = estimate_df,
      N_market = N_market, N_margin = N_margin,
      likelihood_mode = likelihood_mode,
      exact_k = exact_structure && identical(partition_prior, "ordered_k")
    )
  }
  if (any(vapply(split(cluster_int, mky_int), function(x) length(unique(x)) > 1L,
                 logical(1)))) {
    stop("'loo_cluster' must be constant within each market-year equilibrium.")
  }

  # Build market_year -> market and market_year -> year maps
  mky_to_market <- integer(N_market_year)
  mky_to_year <- integer(N_market_year)
  for (i in seq_len(N)) {
    mky_to_market[mky_int[i]] <- market_int[i]
    mky_to_year[mky_int[i]] <- year_int[i]
  }

  # Compute the common conditional quantity-share response plus the
  # conditional revenue shares needed by CES supply FOCs.
  shareIn <- numeric(N)
  revenue <- prices * quantities
  revShare <- numeric(N)
  inside_price_index <- numeric(N_market_year)
  for (m in seq_len(N_market_year)) {
    idx <- which(mky_int == m)
    total_q <- sum(quantities[idx])
    total_rev <- sum(revenue[idx])
    shareIn[idx] <- quantities[idx] / total_q
    revShare[idx] <- revenue[idx] / total_rev
    inside_price_index[m] <- sum(prices[idx] * shareIn[idx])
  }

  firm_revShare <- numeric(N)
  for (m in seq_len(N_market_year)) {
    idx <- which(mky_int == m)
    for (g in unique(group_int[idx])) {
      firm_idx <- idx[group_int[idx] == g]
      firm_revShare[firm_idx] <- sum(revShare[firm_idx])
    }
  }

  # Standardize log-prices
  log_prices <- log(prices)
  log_price_sd <- safe_sd(log_prices)
  log_price_std <- (log_prices - mean(log_prices)) / log_price_sd

  # Common Logit/CES likelihood responses.
  marginNormInv <- 1.0 / margins
  share_comp <- compositional_share_data(shareIn, mky_int, N_market_year)
  log_shareIn <- share_comp$log_share
  clr_shareIn <- share_comp$clr_share
  log_marginNormInv <- log(pmax(marginNormInv, 1e-12))

  # S0 covariates
  if (is.null(X_s0)) {
    K_s0 <- 0L
    X_s0_mat <- matrix(0, nrow = N_market_year, ncol = 0)
  } else {
    X_s0_mat <- as.matrix(X_s0)
    if (nrow(X_s0_mat) != N_market_year) {
      stop("'X_s0' must have one row per market-year.")
    }
    if (!is.numeric(X_s0_mat) || any(!is.finite(X_s0_mat))) {
      stop("'X_s0' must be a finite numeric matrix.")
    }
    K_s0 <- ncol(X_s0_mat)
  }

  # Default outside-share prior is logistic-normal with arithmetic mean 20%
  # and 5th percentile 10%.
  default_priors <- bayes_prior_manifest("baseline")
  # Override with user priors
  for (nm in names(priors)) {
    if (nm %in% names(default_priors)) {
      default_priors[[nm]] <- priors[[nm]]
    } else {
      warning("Unknown prior parameter '", nm, "' ignored.")
    }
  }


  s0_prior_reference <- logit_normal_mean(
    default_priors$s0_prior_mean,
    default_priors$s0_prior_sd
  )
  R0_prior_reference <- stats::plogis(
    stats::qlogis(s0_prior_reference) + log(outside_price) -
      log(inside_price_index)
  )
  if (!is.null(resid_df_fixed)) default_priors$resid_df_fixed <- resid_df_fixed
  validate_priors(default_priors)
  if (!likelihood_mode %in% 1:3) {
    stop("'likelihood_mode' must be 1 (joint), 2 (share), or 3 (margin).")
  }

  # HMT data (uses Lerner index directly for CES)
  avg_price_hmt <- mean(prices)
  avg_margin_hmt <- if (any(margin_observed)) {
    mean(margins[margin_observed])
  } else {
    0.0
  }
  ssnip_hmt <- 0.05

  stan_data <- c(
    list(
      N = N,
      shareIn = array(shareIn, dim = N),
      revShare = array(revShare, dim = N),
      firm_revShare = array(firm_revShare, dim = N),
      marginNormInv = array(marginNormInv, dim = N),
      log_shareIn = array(log_shareIn, dim = N),
      clr_shareIn = array(clr_shareIn, dim = N),
      log_marginNormInv = array(log_marginNormInv, dim = N),
      price = array(as.numeric(prices), dim = N),
      log_price_std = array(log_price_std, dim = N),
      log_price_sd = log_price_sd,
      outside_price = array(outside_price, dim = N_market_year),
      inside_price_index = array(inside_price_index, dim = N_market_year),
      market_sign = market_sign,
      N_market = N_market,
      market = array(market_int, dim = N),
      N_group = N_group,
      group = array(group_int, dim = N),
      N_year = N_year,
      year = array(year_int, dim = N),
      supply_model = supply_model,
      use_cutoff = as.integer(use_cutoff),
      cutoff_scale = cutoff_scale,
      grainsize = as.integer(max(1L, grainsize)),
      N_market_year = N_market_year,
      market_year_idx = array(mky_int, dim = N),
      mky_n = array(share_comp$mky_n, dim = N_market_year),
      mky_start = array(share_comp$mky_start, dim = N_market_year),
      obs_order = array(share_comp$obs_order, dim = N),
      N_cluster = N_cluster,
      cluster = array(cluster_int, dim = N),
      likelihood_mode = as.integer(likelihood_mode),
      margin_obs = array(as.integer(margin_observed), dim = N),
      N_margin = as.integer(N_margin),
      margin_idx = array(as.integer(margin_idx), dim = N_margin),
      share_likelihood = likelihood_modes$share_likelihood_code,
      margin_likelihood = likelihood_modes$margin_likelihood_code,
      quality_parameterization = likelihood_modes$quality_parameterization_code,
      K_s0 = K_s0,
      X_s0 = X_s0_mat,
      mky_to_market = array(mky_to_market),
      mky_to_year = array(mky_to_year),
      is_single_market = is_single_market,
      use_rho = as.integer(use_rho),
      use_market_resid_scale = as.integer(use_market_resid_scale),
      use_alpha_market_pool = as.integer(use_alpha_market_pool),
      use_s0_market_year_pool = as.integer(use_s0_market_year_pool),
      use_s0_market_pool = as.integer(use_s0_market_pool),
      market_s0_parameterization = market_s0_resolution$code,
      use_centered_alpha = as.integer(use_centered_alpha),
      use_centered_s0_pool = as.integer(use_centered_s0_pool),
      use_hmt = as.integer(use_hmt),
      avg_price_hmt = avg_price_hmt,
      avg_margin_hmt = avg_margin_hmt,
      ssnip_hmt = ssnip_hmt,
      estimate_df = as.integer(estimate_df)
    ),
    default_priors
  )

  # Attach metadata for later use
  attr(stan_data, "market_levels") <- levels(market_f)
  attr(stan_data, "year_levels") <- levels(year_f)
  attr(stan_data, "group_levels") <- levels(group_f)
  attr(stan_data, "mky_levels") <- levels(mky_f)
  attr(stan_data, "cluster_levels") <- levels(cluster_f)
  attr(stan_data, "market_year_idx") <- mky_int
  attr(stan_data, "model_spec_version") <- .antitrust_bayes_model_spec_version
  attr(stan_data, "s0_prior_reference") <- s0_prior_reference
  attr(stan_data, "R0_prior_reference") <- R0_prior_reference
  attr(stan_data, "structure") <- resolved_structure$name
  attr(stan_data, "structure_code") <- resolved_structure$code
  attr(stan_data, "cutoff_scale") <- cutoff_scale
  attr(stan_data, "partition_prior") <- partition_prior
  attr(stan_data, "k_scope") <- k_scope
  attr(stan_data, "k_support") <- k_support
  attr(stan_data, "k_prior") <- k_prior
  attr(stan_data, "share_likelihood") <- likelihood_modes$share_likelihood
  attr(stan_data, "margin_likelihood") <- likelihood_modes$margin_likelihood
  attr(stan_data, "quality_parameterization") <- likelihood_modes$quality_parameterization
  attr(stan_data, "likelihood_mode_metadata") <- likelihood_modes
  attr(stan_data, "s0_pooling") <- s0_resolution$resolved
  attr(stan_data, "s0_pooling_requested") <- s0_resolution$requested
  attr(stan_data, "s0_pooling_resolved") <- s0_resolution$resolved
  attr(stan_data, "s0_pooling_source") <- s0_resolution$source
  attr(stan_data, "market_s0_parameterization") <-
    market_s0_resolution$requested
  attr(stan_data, "market_s0_parameterization_requested") <-
    market_s0_resolution$requested
  attr(stan_data, "market_s0_parameterization_effective") <-
    market_s0_resolution$effective
  attr(stan_data, "market_s0_rotation_active") <- market_s0_resolution$active
  attr(stan_data, "market_s0_coordinate_convention") <-
    market_s0_resolution$coordinate_convention
  if (resolved_structure$name == "stackelberg" ||
      (exact_structure && identical(resolved_structure$name, "fringe") &&
       identical(partition_prior, "ordered_k"))) {
    stack_data <- .stack_append_data(
      stan_data = stan_data, mky_int = mky_int, market_int = market_int,
      group_int = group_int, group_levels = levels(group_f), shares = shareIn,
      revShare = revShare, structure = resolved_structure,
      demand_model = "ces", partition_prior = partition_prior,
      k_support = k_support, k_prior = k_prior, k_scope = k_scope,
      cutoff_scale = cutoff_scale
    )
    ranking_metadata <- attr(stack_data, "firm_market_ranking", exact = TRUE)
    partition_metadata <- attr(stack_data, "partition_metadata", exact = TRUE)
    stan_data <- c(stan_data, stack_data)
    if (!is.null(cutoff_scale)) stan_data$cutoff_scale <- cutoff_scale
    attr(stan_data, "market_levels") <- levels(market_f)
    attr(stan_data, "year_levels") <- levels(year_f)
    attr(stan_data, "group_levels") <- levels(group_f)
    attr(stan_data, "mky_levels") <- levels(mky_f)
    attr(stan_data, "cluster_levels") <- levels(cluster_f)
    attr(stan_data, "market_year_idx") <- mky_int
    attr(stan_data, "model_spec_version") <- .antitrust_bayes_model_spec_version
    attr(stan_data, "s0_prior_reference") <- s0_prior_reference
    attr(stan_data, "R0_prior_reference") <- R0_prior_reference
    attr(stan_data, "structure") <- resolved_structure$name
    attr(stan_data, "structure_code") <- resolved_structure$code
    attr(stan_data, "cutoff_scale") <- cutoff_scale
    attr(stan_data, "partition_prior") <- partition_prior
    attr(stan_data, "k_scope") <- k_scope
    attr(stan_data, "k_support") <- k_support
    attr(stan_data, "k_prior") <- k_prior
    attr(stan_data, "share_likelihood") <- likelihood_modes$share_likelihood
    attr(stan_data, "margin_likelihood") <- likelihood_modes$margin_likelihood
    attr(stan_data, "quality_parameterization") <- likelihood_modes$quality_parameterization
    attr(stan_data, "likelihood_mode_metadata") <- likelihood_modes
    attr(stan_data, "s0_pooling") <- s0_resolution$resolved
    attr(stan_data, "s0_pooling_requested") <- s0_resolution$requested
    attr(stan_data, "s0_pooling_resolved") <- s0_resolution$resolved
    attr(stan_data, "s0_pooling_source") <- s0_resolution$source
    attr(stan_data, "market_s0_parameterization") <-
      market_s0_resolution$requested
    attr(stan_data, "market_s0_parameterization_requested") <-
      market_s0_resolution$requested
    attr(stan_data, "market_s0_parameterization_effective") <-
      market_s0_resolution$effective
    attr(stan_data, "market_s0_rotation_active") <- market_s0_resolution$active
    attr(stan_data, "market_s0_coordinate_convention") <-
      market_s0_resolution$coordinate_convention
    attr(stan_data, "firm_market_ranking") <- ranking_metadata
    attr(stan_data, "partition_metadata") <- partition_metadata
  }

  stan_data
}

#' Compute Alpha Prior from Data and Supply FOC
#'
#' Uses **share-weighted mean and SD** of observation-level FOC-implied
#' alphas to derive a log-scale prior for the price coefficient. Large products
#' receive more weight in calibration, reducing sensitivity to small-product outliers.
#'
#' @param shares Numeric vector of within-market-year quantity shares.
#' @param margins_abs Numeric vector of absolute price-cost margins (not Lerner indices).
#' @param s0_reference Probability-scale outside share used to calibrate the
#'   FOC prior. Conditional inside shares do not identify this quantity.
#' @param mky_int Integer vector mapping observations to market-years.
#' @param supply_model Integer 1--5.
#' @param price_sd Standard deviation of prices used by the standardized Stan
#'   coefficient. The FOC-implied structural alpha is multiplied by this
#'   value before taking logs.
#' @return A list with elements \code{mean} and \code{sd} on log scale.
#' @keywords internal
compute_alpha_prior <- function(shares, margins_abs, s0_reference, mky_int,
                                supply_model, price_sd = 1.0) {
  shares <- bounded01(shares)
  margins_abs <- pmax(margins_abs, 1e-12)
  s0 <- if (length(s0_reference) == 1L) {
    rep(bounded01(s0_reference), length(shares))
  } else {
    bounded01(s0_reference[mky_int])
  }
  margin_inv <- 1.0 / margins_abs

  conduct <- if (supply_model == 1L) {
    1.0 - shares * (1.0 - s0)
  } else if (supply_model == 2L) {
    rep(1.0, length(shares))
  } else if (supply_model == 3L) {
    s0 / (s0 + shares)
  } else {
    denom <- pmin(0.99, shares * (1.0 - s0))
    denom / (-log1p(-denom))
  }

  alpha_draws <- margin_inv / pmax(conduct, 1e-12)
  # Stan estimates a = alpha * sd(price), not structural alpha.
  alpha_draws <- alpha_draws * max(price_sd, 1e-8)

  # Pass shares as weights for share-weighted calibration
  log_prior_from_positive_draws(alpha_draws, weights = shares)
}

#' Compute CES Gamma Prior from Data and Supply FOC
#'
#' Uses **revenue-share-weighted mean and SD** of observation-level
#' FOC-implied CES coefficients to derive a log-scale prior for the standardized
#' elasticity coefficient.
#'
#' @param revShare Numeric vector of revenue shares.
#' @param margins Numeric vector of Lerner indices (p-c)/p.
#' @param s0_reference Probability-scale outside revenue share used to
#'   convert conditional inside shares to unconditional shares.
#' @param mky_int Integer vector mapping observations to market-years.
#' @param supply_model Integer 1--5.
#' @param log_price_sd Scale used to convert structural gamma to the
#'   standardized coefficient estimated by the Stan model.
#' @return A list with elements \code{mean} and \code{sd} on log scale.
#' @keywords internal
compute_gamma_prior <- function(revShare, margins, s0_reference, mky_int,
                                supply_model, log_price_sd = 1.0) {
  revShare <- bounded01(revShare)
  margins <- bounded01(margins)
  s0 <- if (length(s0_reference) == 1L) {
    rep(bounded01(s0_reference), length(revShare))
  } else {
    bounded01(s0_reference[mky_int])
  }
  rev_share_uncond <- revShare * (1.0 - s0)

  gamma_draws <- if (supply_model == 1L) {
    (1 / margins - rev_share_uncond) / (1 - rev_share_uncond)
  } else if (supply_model == 2L) {
    1 / margins
  } else if (supply_model == 3L) {
    denom <- margins - rev_share_uncond
    out <- rep(NA_real_, length(margins))
    ## Solve the full CES Marshallian Cournot FOC for gamma.  With q the
    ## unconditional firm revenue share and s0 the outside revenue share,
    ## L = [1 + (gamma - 1)(s0 + q)] /
    ##     [gamma {1 + (gamma - 1)s0}].
    A <- margins * s0
    B <- margins * (1 - s0) - (s0 + rev_share_uncond)
    C <- -(1 - s0 - rev_share_uncond)
    disc <- B^2 - 4 * A * C
    ok <- is.finite(A) & is.finite(B) & is.finite(C) &
      is.finite(disc) & disc >= 0 & A > 0
    out[ok] <- (-B[ok] + sqrt(disc[ok])) / (2 * A[ok])
    out
  } else {
    1.0 + log1p(-rev_share_uncond) / log1p(-margins)
  }

  # Filter to finite, valid gamma values
  valid <- is.finite(gamma_draws) & gamma_draws > 1.0
  gamma_draws_valid <- gamma_draws[valid]
  revShare_valid <- revShare[valid]

  a_draws <- (gamma_draws_valid - 1.0) * max(log_price_sd, 1e-8)

  # Pass revenue shares as weights for share-weighted calibration
  log_prior_from_positive_draws(a_draws, weights = revShare_valid)
}

#' Validate Inputs for bayes_est
#' @keywords internal
validate_inputs <- function(prices, quantities, margins, market, year, group,
                             margin_observed = NULL) {
  N <- length(prices)

  if (N == 0L) {
    stop("At least one observation is required.")
  }
  if (!is.numeric(prices) || !is.numeric(quantities) || !is.numeric(margins)) {
    stop("'prices', 'quantities', and 'margins' must be numeric vectors.")
  }

  if (length(quantities) != N) {
    stop("'quantities' must have the same length as 'prices'.")
  }
  if (length(margins) != N) {
    stop("'margins' must have the same length as 'prices'.")
  }
  if (is.null(margin_observed)) {
    if (any(margins <= 0 | !is.finite(margins))) {
      stop("'margins' must be positive and finite.")
    }
  } else {
    if (length(margin_observed) != N) {
      stop("'margin_observed' must have the same length as 'prices'.")
    }
    observed <- as.logical(margin_observed) & !is.na(margin_observed)
    observed_margins <- margins[observed]
    if (any(observed_margins <= 0 | !is.finite(observed_margins))) {
      stop("'margins' must be positive and finite for observed (margin_observed == TRUE) rows.")
    }
  }
  if (any(quantities <= 0 | !is.finite(quantities))) {
    stop("'quantities' must be positive and finite.")
  }
  if (any(prices <= 0 | !is.finite(prices))) {
    stop("'prices' must be positive and finite.")
  }
  if (length(market) != N) {
    stop("'market' must have the same length as 'prices'.")
  }
  if (length(year) != N) {
    stop("'year' must have the same length as 'prices'.")
  }
  if (length(group) != N) {
    stop("'group' must have the same length as 'prices'.")
  }
  if (anyNA(market) || anyNA(year) || anyNA(group)) {
    stop("'market', 'year', and 'group' identifiers may not be missing.")
  }

  invisible(TRUE)
}
# Current fitted-object/model contract. Increment whenever the Stan likelihood
# changes in a way that makes serialized fits unsuitable for comparison or
# simulation under the current package.
.antitrust_bayes_model_spec_version <- 9L

# Arithmetic mean of a logistic-normal random variable. The outside-share
# prior is parameterized on the logit scale in Stan, while FOC prior
# calibration needs a probability-scale reference.
logit_normal_mean <- function(location, scale) {
  if (!is.finite(location) || !is.finite(scale) || scale < 0) {
    stop("Logistic-normal location and scale must be finite, with scale >= 0.")
  }
  if (scale == 0) return(stats::plogis(location))

  stats::integrate(
    function(z) stats::plogis(location + scale * z) * stats::dnorm(z),
    lower = -Inf,
    upper = Inf,
    rel.tol = 1e-10
  )$value
}

validate_fit_spec <- function(fit) {
  fit_version <- attr(fit@stan_data, "model_spec_version", exact = TRUE)
  if (is.null(fit_version) ||
      !identical(as.integer(fit_version), .antitrust_bayes_model_spec_version)) {
    found <- if (is.null(fit_version)) "pre-versioned" else as.character(fit_version)
    stop(
      "This fit uses antitrustBayes model specification ", found,
      "; the current specification is ", .antitrust_bayes_model_spec_version,
      ". Refit it with bayes_est() before comparison or simulation.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}
