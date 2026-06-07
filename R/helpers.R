#' Prepare Stan Data for logit Demand
#'
#' Transforms user-supplied vectors into the list expected by the logit Stan
#' model. Computes quantity shares from raw quantities and converts positive
#' price-normalized margins to inverse absolute margins and log residual inputs.
#'
#' @param prices Numeric vector of product prices.
#' @param quantities Numeric vector of product quantities.
#' @param margins Numeric vector of positive price-normalized margins. Values
#'   greater than 1 are valid for input-market markdowns.
#' @param market Factor or vector identifying markets.
#' @param year Factor or vector identifying time periods.
#' @param group Factor or vector identifying product groups (e.g., firms).
#' @param supply_model Integer 1--4.
#' @param use_cutoff Logical.
#' @param use_rho Logical.
#' @param use_market_resid_scale Logical.
#' @param use_alpha_market_pool Logical.
#' @param use_s0_market_year_pool Logical.
#' @param use_hmt Logical.
#' @param priors Named list of prior hyperparameters.
#' @param X_s0 Optional model matrix for s0 covariates.
#' @param grainsize Positive integer reduce-sum chunk size.
#' @return A list suitable for passing to \code{rstan::sampling}.
#' @keywords internal
safe_sd <- function(x, floor = 1e-8) {
  s <- stats::sd(x)
  if (!is.finite(s) || s < floor) 1.0 else s
}

bounded01 <- function(x, eps = 1e-6) {
  pmin(pmax(x, eps), 1 - eps)
}

log_prior_from_positive_draws <- function(x, fallback_mean = 1.0,
                                          fallback_sd = 1.0) {
  x <- x[is.finite(x) & x > 0]
  if (length(x) < 2L) {
    return(list(mean = log(fallback_mean), sd = fallback_sd))
  }

  qs <- stats::quantile(log(x), c(0.10, 0.90), na.rm = TRUE, names = FALSE)
  if (!all(is.finite(qs)) || qs[2] <= qs[1]) {
    return(list(mean = mean(log(x), na.rm = TRUE), sd = fallback_sd))
  }

  list(
    mean = mean(qs),
    sd = max((qs[2] - qs[1]) / 2.56, 0.15)
  )
}

prepare_stan_data <- function(prices, quantities, margins, market, year, group,
                              supply_model, use_cutoff = FALSE,
                              use_rho = TRUE,
                              use_market_resid_scale = FALSE,
                              use_hmt = FALSE,
                              use_alpha_market_pool = TRUE,
                              use_s0_market_year_pool = FALSE,
                              priors = list(), X_s0 = NULL,
                              grainsize = 1L) {
  N <- length(prices)

  # Encode factors
  market_f <- as.factor(market)
  year_f <- as.factor(year)
  group_f <- as.factor(group)

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

  # Standardize prices
  price_sd <- safe_sd(prices)
  price_std <- (prices - mean(prices)) / price_sd

  # Convert price-normalized margin to inverse absolute margin.
  marginInv <- 1.0 / (margins * prices)
  log_shareIn <- log(pmax(shares, 1e-12))
  log_marginInv <- log(pmax(marginInv, 1e-12))

  # Implied outside share per market-year (quantity-based)
  implied_s0 <- numeric(N_market_year)
  for (m in seq_len(N_market_year)) {
    idx <- which(mky_int == m)
    implied_s0[m] <- max(0.01, min(0.99, 1.0 - sum(shares[idx])))
  }

  # S0 covariates
  if (is.null(X_s0)) {
    K_s0 <- 0L
    X_s0_mat <- matrix(0, nrow = N_market_year, ncol = 0)
  } else {
    X_s0_mat <- as.matrix(X_s0)
    if (nrow(X_s0_mat) != N_market_year) {
      stop("'X_s0' must have one row per market-year.")
    }
    K_s0 <- ncol(X_s0_mat)
  }

  # Compute model-specific alpha prior from data and FOC
  # Convert price-normalized margins to absolute margins for FOC inversion
  margins_abs <- margins * prices
  alpha_prior <- compute_alpha_prior(
    shares, margins_abs, implied_s0,
    mky_int, supply_model
  )

  # Default priors (matched to mergerBayes defaults)
  default_priors <- list(
    prior_sigma_share = 0.5,
    prior_sigma_margin = 0.75,
    prior_sigma_meanval_strat = 1.0,
    prior_sigma_meanval_fringe = 1.0,
    prior_alpha_mean = alpha_prior$mean,
    prior_alpha_sd = alpha_prior$sd,
    s0_prior_mean = -3.0,
    s0_prior_sd = 1.2,
    prior_sigma_alpha = 0.5,
    prior_sigma_s0_pool = 0.5,
    prior_sigma_beta_s0 = 1.0,
    prior_sigma_resid_market = 0.35,
    prior_lkj = 2.0,
    prior_cutoff_alpha = 3.0,
    prior_cutoff_beta = 100.0
  )
  # Override with user priors
  for (nm in names(priors)) {
    if (nm %in% names(default_priors)) {
      default_priors[[nm]] <- priors[[nm]]
    } else {
      warning("Unknown prior parameter '", nm, "' ignored.")
    }
  }

  # HMT data (uses absolute margins)
  avg_price_hmt <- mean(prices)
  avg_margin_hmt <- mean(margins_abs)
  ssnip_hmt <- 0.05

  stan_data <- c(
    list(
      N = N,
      shareIn = shares,
      marginInv = marginInv,
      log_shareIn = log_shareIn,
      log_marginInv = log_marginInv,
      price_std = price_std,
      price_sd = price_sd,
      N_market = N_market,
      market = market_int,
      N_group = N_group,
      group = group_int,
      N_year = N_year,
      year = year_int,
      supply_model = supply_model,
      use_cutoff = as.integer(use_cutoff),
      grainsize = as.integer(max(1L, grainsize)),
      N_market_year = N_market_year,
      market_year_idx = mky_int,
      K_s0 = K_s0,
      X_s0 = X_s0_mat,
      mky_to_market = array(mky_to_market),
      mky_to_year = array(mky_to_year),
      is_single_market = is_single_market,
      use_rho = as.integer(use_rho),
      use_market_resid_scale = as.integer(use_market_resid_scale),
      use_alpha_market_pool = as.integer(use_alpha_market_pool),
      use_s0_market_year_pool = as.integer(use_s0_market_year_pool),
      use_hmt = as.integer(use_hmt),
      avg_price_hmt = avg_price_hmt,
      avg_margin_hmt = avg_margin_hmt,
      ssnip_hmt = ssnip_hmt,
      implied_s0 = array(implied_s0)
    ),
    default_priors
  )

  # Attach metadata for later use
  attr(stan_data, "market_levels") <- levels(market_f)
  attr(stan_data, "year_levels") <- levels(year_f)
  attr(stan_data, "group_levels") <- levels(group_f)
  attr(stan_data, "mky_levels") <- levels(mky_f)
  attr(stan_data, "market_year_idx") <- mky_int

  stan_data
}

#' Prepare Stan Data for CES Demand
#'
#' Transforms user-supplied vectors into the list expected by the CES Stan
#' model. Computes revenue shares from prices and quantities, and converts
#' bounded Lerner margins to inverse Lerner indices and log residual inputs.
#'
#' @inheritParams prepare_stan_data
#' @return A list suitable for passing to \code{rstan::sampling}.
#' @keywords internal
prepare_ces_stan_data <- function(prices, quantities, margins, market, year,
                                  group, supply_model, use_cutoff, use_rho,
                                  use_hmt, priors, X_s0 = NULL,
                                  grainsize = 1L) {
  N <- length(prices)

  # CES FOCs require L < 1 (positive marginal costs)
  if (any(margins >= 1)) {
    stop("CES demand requires Lerner indices strictly less than 1. ",
         "Found ", sum(margins >= 1), " observation(s) with L >= 1. ",
         "Use demand_model = \"logit\" for input markets with L >= 1.")
  }

  # Encode factors
  market_f <- as.factor(market)
  year_f <- as.factor(year)
  group_f <- as.factor(group)

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

  # Build market_year -> market and market_year -> year maps
  mky_to_market <- integer(N_market_year)
  mky_to_year <- integer(N_market_year)
  for (i in seq_len(N)) {
    mky_to_market[mky_int[i]] <- market_int[i]
    mky_to_year[mky_int[i]] <- year_int[i]
  }

  # Compute revenue shares within each market-year
  revenue <- prices * quantities
  revShare <- numeric(N)
  for (m in seq_len(N_market_year)) {
    idx <- which(mky_int == m)
    total_rev <- sum(revenue[idx])
    revShare[idx] <- revenue[idx] / total_rev
  }

  # Standardize log-prices
  log_prices <- log(prices)
  log_price_sd <- safe_sd(log_prices)
  log_price_std <- (log_prices - mean(log_prices)) / log_price_sd

  # Inverse Lerner index: 1/L = p/(p-c)
  lernerInv <- 1.0 / margins
  log_revShare <- log(pmax(revShare, 1e-12))
  log_lernerInv <- log(pmax(lernerInv, 1e-12))

  # Implied outside revenue share per market-year
  implied_s0 <- numeric(N_market_year)
  for (m in seq_len(N_market_year)) {
    idx <- which(mky_int == m)
    implied_s0[m] <- max(0.01, min(0.99, 1.0 - sum(revShare[idx])))
  }

  # S0 covariates
  if (is.null(X_s0)) {
    K_s0 <- 0L
    X_s0_mat <- matrix(0, nrow = N_market_year, ncol = 0)
  } else {
    X_s0_mat <- as.matrix(X_s0)
    if (nrow(X_s0_mat) != N_market_year) {
      stop("'X_s0' must have one row per market-year.")
    }
    K_s0 <- ncol(X_s0_mat)
  }

  # Compute CES-specific gamma prior from data and FOC
  gamma_prior <- compute_gamma_prior(
    revShare, margins, implied_s0,
    mky_int, supply_model, log_price_sd
  )

  # Default priors
  default_priors <- list(
    prior_sigma_share = 0.5,
    prior_sigma_margin = 0.75,
    prior_sigma_meanval_strat = 1.0,
    prior_sigma_meanval_fringe = 1.0,
    prior_alpha_mean = gamma_prior$mean,
    prior_alpha_sd = gamma_prior$sd,
    s0_prior_mean = -3.0,
    s0_prior_sd = 1.2,
    prior_sigma_alpha = 1.0,
    prior_sigma_beta_s0 = 1.0,
    prior_lkj = 2.0,
    prior_cutoff_alpha = 3.0,
    prior_cutoff_beta = 100.0
  )
  # Override with user priors
  for (nm in names(priors)) {
    if (nm %in% names(default_priors)) {
      default_priors[[nm]] <- priors[[nm]]
    } else {
      warning("Unknown prior parameter '", nm, "' ignored.")
    }
  }

  # HMT data (uses Lerner index directly for CES)
  avg_price_hmt <- mean(prices)
  avg_margin_hmt <- mean(margins)
  ssnip_hmt <- 0.05

  stan_data <- c(
    list(
      N = N,
      revShare = revShare,
      lernerInv = lernerInv,
      log_revShare = log_revShare,
      log_lernerInv = log_lernerInv,
      log_price_std = log_price_std,
      log_price_sd = log_price_sd,
      N_market = N_market,
      market = market_int,
      N_group = N_group,
      group = group_int,
      N_year = N_year,
      year = year_int,
      supply_model = supply_model,
      use_cutoff = as.integer(use_cutoff),
      grainsize = as.integer(max(1L, grainsize)),
      N_market_year = N_market_year,
      market_year_idx = mky_int,
      K_s0 = K_s0,
      X_s0 = X_s0_mat,
      mky_to_market = array(mky_to_market),
      mky_to_year = array(mky_to_year),
      is_single_market = is_single_market,
      use_rho = as.integer(use_rho),
      use_hmt = as.integer(use_hmt),
      avg_price_hmt = avg_price_hmt,
      avg_margin_hmt = avg_margin_hmt,
      ssnip_hmt = ssnip_hmt,
      implied_s0 = array(implied_s0)
    ),
    default_priors
  )

  # Attach metadata for later use
  attr(stan_data, "market_levels") <- levels(market_f)
  attr(stan_data, "year_levels") <- levels(year_f)
  attr(stan_data, "group_levels") <- levels(group_f)
  attr(stan_data, "mky_levels") <- levels(mky_f)
  attr(stan_data, "market_year_idx") <- mky_int

  stan_data
}

#' Compute Model-Specific Alpha Prior from Data and Supply FOC
#'
#' Uses robust quantiles of observation-level FOC-implied coefficients to
#' derive a log-scale prior for the price coefficient alpha.
#'
#' @param shares Numeric vector of product shares.
#' @param margins_abs Numeric vector of absolute price-cost margins (p-c).
#' @param implied_s0 Numeric vector of implied outside shares per market-year.
#' @param mky_int Integer vector mapping observations to market-years.
#' @param supply_model Integer 1--4.
#' @return A list with elements \code{mean} and \code{sd}.
#' @keywords internal
compute_alpha_prior <- function(shares, margins_abs, implied_s0, mky_int,
                                supply_model) {
  shares <- bounded01(shares)
  margins_abs <- pmax(margins_abs, 1e-12)
  s0 <- bounded01(implied_s0[mky_int])
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
  log_prior_from_positive_draws(alpha_draws)
}

#' Compute CES Gamma Prior from Data and Supply FOC
#'
#' Uses robust quantiles of observation-level FOC-implied CES coefficients to
#' derive a log-scale prior for the standardized elasticity coefficient.
#'
#' @param revShare Numeric vector of revenue shares.
#' @param margins Numeric vector of Lerner indices (p-c)/p.
#' @param implied_s0 Numeric vector of implied outside revenue shares.
#' @param mky_int Integer vector mapping observations to market-years.
#' @param supply_model Integer 1--4.
#' @param log_price_sd Scale used to convert structural gamma to the
#'   standardized coefficient estimated by the Stan model.
#' @return A list with elements \code{mean} and \code{sd}.
#' @keywords internal
compute_gamma_prior <- function(revShare, margins, implied_s0, mky_int,
                                supply_model, log_price_sd = 1.0) {
  revShare <- bounded01(revShare)
  margins <- bounded01(margins)

  gamma_draws <- if (supply_model == 1L) {
    (1 / margins - revShare) / (1 - revShare)
  } else if (supply_model == 2L) {
    1 / margins
  } else if (supply_model == 3L) {
    denom <- margins - revShare
    out <- rep(NA_real_, length(margins))
    ok <- is.finite(denom) & denom > 1e-8
    out[ok] <- (1 - revShare[ok]) / denom[ok]
    out
  } else {
    1.0 - log1p(-revShare) / log1p(-margins)
  }

  gamma_draws <- gamma_draws[is.finite(gamma_draws) & gamma_draws > 1.0]
  a_draws <- (gamma_draws - 1.0) * max(log_price_sd, 1e-8)
  log_prior_from_positive_draws(a_draws)
}

#' Validate Inputs for bayes_est
#' @keywords internal
validate_inputs <- function(prices, quantities, margins, market, year, group) {
  N <- length(prices)

  if (length(quantities) != N) {
    stop("'quantities' must have the same length as 'prices'.")
  }
  if (length(margins) != N) {
    stop("'margins' must have the same length as 'prices'.")
  }
  if (any(margins <= 0 | !is.finite(margins))) {
    stop("'margins' must be positive and finite.")
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

  invisible(TRUE)
}
