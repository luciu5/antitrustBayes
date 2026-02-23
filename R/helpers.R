#' Prepare Stan Data for logit Demand
#'
#' Transforms user-supplied vectors into the list expected by the logit Stan
#' model. Computes quantity shares from raw quantities and converts Lerner
#' margins to inverse absolute margins.
#'
#' @param prices Numeric vector of product prices.
#' @param quantities Numeric vector of product quantities.
#' @param margins Numeric vector of Lerner indices (p-c)/p, between 0 and 1.
#' @param market Factor or vector identifying markets.
#' @param year Factor or vector identifying time periods.
#' @param group Factor or vector identifying product groups (e.g., firms).
#' @param supply_model Integer 1--4.
#' @param use_cutoff Logical.
#' @param use_rho Logical.
#' @param use_hmt Logical.
#' @param priors Named list of prior hyperparameters.
#' @param X_s0 Optional model matrix for s0 covariates.
#' @return A list suitable for passing to \code{rstan::sampling}.
#' @keywords internal
prepare_stan_data <- function(prices, quantities, margins, market, year, group,
                              supply_model, use_cutoff, use_rho, use_hmt,
                              priors, X_s0 = NULL) {
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
  price_sd <- stats::sd(prices)
  price_std <- (prices - mean(prices)) / price_sd

  # Convert Lerner index to inverse absolute margin: 1/(L*p) = 1/(p-c)
  marginInv <- 1.0 / (margins * prices)

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
    K_s0 <- ncol(X_s0_mat)
  }

  # Compute model-specific alpha prior from data and FOC
  # Convert Lerner to absolute margins for FOC inversion
  margins_abs <- margins * prices
  alpha_prior <- compute_alpha_prior(
    shares, margins_abs, implied_s0,
    mky_int, supply_model
  )

  # Default priors (matched to mergerBayes defaults)
  default_priors <- list(
    prior_sigma_share = 0.01,
    prior_sigma_margin = 1.0,
    prior_sigma_meanval_strat = 1.0,
    prior_sigma_meanval_fringe = 1.0,
    prior_alpha_mean = alpha_prior$mean,
    prior_alpha_sd = alpha_prior$sd,
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

  # HMT data (uses absolute margins)
  avg_price_hmt <- mean(prices)
  avg_margin_hmt <- mean(margins_abs)
  ssnip_hmt <- 0.05

  stan_data <- c(
    list(
      N = N,
      shareIn = shares,
      marginInv = marginInv,
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
      grainsize = 1L,
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

#' Prepare Stan Data for CES Demand
#'
#' Transforms user-supplied vectors into the list expected by the CES Stan
#' model. Computes revenue shares from prices and quantities, and converts
#' Lerner margins to inverse Lerner indices.
#'
#' @inheritParams prepare_stan_data
#' @return A list suitable for passing to \code{rstan::sampling}.
#' @keywords internal
prepare_ces_stan_data <- function(prices, quantities, margins, market, year,
                                  group, supply_model, use_cutoff, use_rho,
                                  use_hmt, priors, X_s0 = NULL) {
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
  log_price_sd <- stats::sd(log_prices)
  log_price_std <- (log_prices - mean(log_prices)) / log_price_sd

  # Inverse Lerner index: 1/L = p/(p-c)
  lernerInv <- 1.0 / margins

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
    K_s0 <- ncol(X_s0_mat)
  }

  # Compute CES-specific gamma prior from data and FOC
  gamma_prior <- compute_gamma_prior(
    revShare, margins, implied_s0,
    mky_int, supply_model
  )

  # Default priors
  default_priors <- list(
    prior_sigma_share = 0.01,
    prior_sigma_margin = 1.0,
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
      grainsize = 1L,
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
#' Uses 5th/95th percentile of shares and margins to derive a prior range
#' for the price coefficient alpha that is consistent with the chosen supply
#' model's first-order conditions.
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
  s_low <- stats::quantile(shares, 0.10)
  s_high <- stats::quantile(shares, 0.90)
  m_low <- stats::quantile(margins_abs, 0.10)
  m_high <- stats::quantile(margins_abs, 0.90)

  if (supply_model == 1L) {
    # Bertrand: alpha ~ 1 / ((p-c) * (1 - s))
    a_low <- 1 / (m_high * (1 - s_low))
    a_high <- 1 / (m_low * (1 - s_high))
  } else if (supply_model == 2L) {
    # Monopolistic competition: alpha = 1 / (p-c)
    a_low <- 1 / m_high
    a_high <- 1 / m_low
  } else if (supply_model == 3L) {
    # Cournot: alpha ~ (s0 + s) / (m * s0)
    s0_proxy <- 0.01
    a_low <- (s0_proxy + s_low) / (m_high * s0_proxy)
    a_high <- (s0_proxy + s_high) / (m_low * s0_proxy)
  } else {
    # Auction: alpha ~ -log(1 - s) / (m * s)
    a_low <- -log(1 - s_low) / (m_high * s_low)
    a_high <- -log(1 - s_high) / (m_low * s_high)
  }

  # Cover [a_low, a_high] with ~95% mass
  p_mean <- (a_high + a_low) / 2
  p_sd <- (a_high - a_low) / 4

  list(mean = p_mean, sd = p_sd)
}

#' Compute CES Gamma Prior from Data and Supply FOC
#'
#' Uses 5th/95th percentile of revenue shares and Lerner indices to derive
#' a prior range for the CES elasticity of substitution gamma (= sigma)
#' consistent with the chosen supply model's first-order conditions.
#'
#' @param revShare Numeric vector of revenue shares.
#' @param margins Numeric vector of Lerner indices (p-c)/p.
#' @param implied_s0 Numeric vector of implied outside revenue shares.
#' @param mky_int Integer vector mapping observations to market-years.
#' @param supply_model Integer 1--4.
#' @return A list with elements \code{mean} and \code{sd}.
#' @keywords internal
compute_gamma_prior <- function(revShare, margins, implied_s0, mky_int,
                                supply_model) {
  r_low <- stats::quantile(revShare, 0.10)
  r_high <- stats::quantile(revShare, 0.90)
  L_low <- stats::quantile(margins, 0.10)
  L_high <- stats::quantile(margins, 0.90)

  if (supply_model == 1L) {
    # Bertrand CES: L = (gamma-1)*r / (gamma - (gamma-1)*r)
    # => gamma = (1/L - r) / (1 - r) ... rearranging lernerInv = gamma - (gamma-1)*r
    # => gamma = (lernerInv - r) / (1 - r) ... but lernerInv = 1/L
    g_low <- (1 / L_high - r_high) / (1 - r_high)
    g_high <- (1 / L_low - r_low) / (1 - r_low)
  } else if (supply_model == 2L) {
    # Monopolistic competition: L = 1/gamma => gamma = 1/L
    g_low <- 1 / L_high
    g_high <- 1 / L_low
  } else if (supply_model == 3L) {
    # Cournot CES: L = (1 + (gamma-1)*r) / gamma
    # => gamma * L = 1 + (gamma-1)*r => gamma*(L - r) = 1 - r
    # => gamma = (1 - r) / (L - r)
    g_low <- (1 - r_low) / (L_high - r_low)
    g_high <- (1 - r_high) / (L_low - r_high)
  } else {
    # Auction CES: (1-L)^(gamma-1) = (1-r)
    # => gamma = 1 - log(1-r) / log(1-L)
    g_low <- 1.0 - log(1 - r_low) / log(1 - L_high)
    g_high <- 1.0 - log(1 - r_high) / log(1 - L_low)
  }

  # Ensure gamma > 1
  g_low <- max(g_low, 1.01)
  g_high <- max(g_high, g_low + 0.1)

  # Cover [g_low, g_high] with ~95% mass
  p_mean <- (g_high + g_low) / 2
  p_sd <- (g_high - g_low) / 4

  list(mean = p_mean, sd = p_sd)
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
    stop("'margins' (Lerner index (p-c)/p) must be positive and finite.")
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
