# Internal helper used by both data-preparation routines.
safe_sd <- function(x, floor = 1e-8) {
  s <- stats::sd(x)
  if (!is.finite(s) || s < floor) 1.0 else s
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
#' @param use_hmt Logical.
#' @param use_centered_alpha Logical.
#' @param use_centered_s0_pool Logical.
#' @param priors Named list of prior hyperparameters.
#' @param X_s0 Optional model matrix for s0 covariates.
#' @param outside_price Positive scalar or one value per market-year.
#' @param grainsize Positive integer reduce-sum chunk size.
#' @param margin_observed Optional logical margin-observation mask.
#' @param likelihood_mode Integer: 1 joint, 2 share only, 3 margin only.
#' @param loo_cluster Optional natural predictive-cluster identifier.
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
                              priors = list(), X_s0 = NULL,
                              grainsize = 1L,
                              margin_observed = NULL,
                              outside_price = 1,
                              likelihood_mode = 1L,
                              loo_cluster = market,
                              estimate_df = FALSE,
                              resid_df_fixed = NULL) {
  N <- length(prices)
  market_side <- match.arg(market_side, c("input", "output"))
  market_sign <- if (market_side == "input") 1L else -1L

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
  outside_price <- resolve_outside_price(outside_price, N_market_year)
  cluster_f <- droplevels(as.factor(loo_cluster))
  if (length(cluster_f) != N || anyNA(cluster_f)) {
    stop("'loo_cluster' must have one non-missing value per observation.")
  }
  cluster_int <- as.integer(cluster_f)
  N_cluster <- nlevels(cluster_f)
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
      shareIn = shares,
      firm_shareIn = firm_shareIn,
      marginNormInv = marginNormInv,
      log_shareIn = log_shareIn,
      clr_shareIn = clr_shareIn,
      log_marginNormInv = log_marginNormInv,
      price = as.numeric(prices),
      price_std = price_std,
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
                                  use_centered_alpha = FALSE,
                                  use_centered_s0_pool = FALSE,
                                  margin_observed = NULL,
                                  outside_price = 1,
                                  likelihood_mode = 1L,
                                  loo_cluster = market,
                                  estimate_df = FALSE,
                                  resid_df_fixed = NULL) {
  N <- length(prices)
  market_side <- match.arg(market_side, c("input", "output"))
  if (market_side != "output") {
    stop("CES estimation currently supports output markets only; use ",
         "demand_model = \"logit\" for input markets.")
  }
  market_sign <- -1L

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
  outside_price <- resolve_outside_price(outside_price, N_market_year)
  cluster_f <- droplevels(as.factor(loo_cluster))
  if (length(cluster_f) != N || anyNA(cluster_f)) {
    stop("'loo_cluster' must have one non-missing value per observation.")
  }
  cluster_int <- as.integer(cluster_f)
  N_cluster <- nlevels(cluster_f)
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
      shareIn = shareIn,
      revShare = revShare,
      firm_revShare = firm_revShare,
      marginNormInv = marginNormInv,
      log_shareIn = log_shareIn,
      clr_shareIn = clr_shareIn,
      log_marginNormInv = log_marginNormInv,
      price = as.numeric(prices),
      log_price_std = log_price_std,
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
      K_s0 = K_s0,
      X_s0 = X_s0_mat,
      mky_to_market = array(mky_to_market),
      mky_to_year = array(mky_to_year),
      is_single_market = is_single_market,
      use_rho = as.integer(use_rho),
      use_market_resid_scale = as.integer(use_market_resid_scale),
      use_alpha_market_pool = as.integer(use_alpha_market_pool),
      use_s0_market_year_pool = as.integer(use_s0_market_year_pool),
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
.antitrust_bayes_model_spec_version <- 6L

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
