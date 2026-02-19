#' Prepare Stan Data from User Inputs
#'
#' Transforms user-supplied vectors into the list expected by the Stan model.
#'
#' @param prices Numeric vector of product prices.
#' @param shares Numeric vector of product shares (within-market conditional).
#' @param margins Numeric vector of product price-cost margins (price/MC).
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
prepare_stan_data <- function(prices, shares, margins, market, year, group,
                              supply_model, use_cutoff, use_rho, use_hmt,
                              priors, X_s0 = NULL) {

  N <- length(prices)

  # Encode factors
  market_f <- as.factor(market)
  year_f   <- as.factor(year)
  group_f  <- as.factor(group)

  market_int <- as.integer(market_f)
  year_int   <- as.integer(year_f)
  group_int  <- as.integer(group_f)

  N_market <- nlevels(market_f)
  N_year   <- nlevels(year_f)
  N_group  <- nlevels(group_f)

  is_single_market <- as.integer(N_market == 1L)

  # Normalize shares within each market-year
  mky_labels <- paste(market_int, year_int, sep = "_")
  mky_f <- as.factor(mky_labels)
  mky_int <- as.integer(mky_f)
  N_market_year <- nlevels(mky_f)

  # Build market_year -> market and market_year -> year maps
  mky_to_market <- integer(N_market_year)
  mky_to_year   <- integer(N_market_year)
  for (i in seq_len(N)) {
    mky_to_market[mky_int[i]] <- market_int[i]
    mky_to_year[mky_int[i]]   <- year_int[i]
  }

  # Standardize prices
  price_sd <- stats::sd(prices)
  price_std <- (prices - mean(prices)) / price_sd

  # Inverse margin
  marginInv <- 1.0 / margins

  # Implied outside share per market-year
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
  alpha_prior <- compute_alpha_prior(shares, margins, implied_s0,
                                     mky_int, supply_model)

  # Default priors (matched to mergerBayes defaults)
  default_priors <- list(
    prior_sigma_share         = 0.01,
    prior_sigma_margin        = 1.0,
    prior_sigma_meanval_strat = 1.0,
    prior_sigma_meanval_fringe = 1.0,
    prior_alpha_mean          = alpha_prior$mean,
    prior_alpha_sd            = alpha_prior$sd,
    s0_prior_mean             = -3.0,
    s0_prior_sd               = 1.2,
    prior_sigma_alpha         = 1.0,
    prior_sigma_beta_s0       = 1.0,
    prior_lkj                 = 2.0
  )
  # Override with user priors
  for (nm in names(priors)) {
    if (nm %in% names(default_priors)) {
      default_priors[[nm]] <- priors[[nm]]
    } else {
      warning("Unknown prior parameter '", nm, "' ignored.")
    }
  }

  # HMT data
  avg_price_hmt  <- mean(prices)
  avg_margin_hmt <- mean(margins)
  ssnip_hmt      <- 0.05

  stan_data <- c(
    list(
      N              = N,
      shareIn        = shares,
      marginInv      = marginInv,
      price_std      = price_std,
      price_sd       = price_sd,
      N_market       = N_market,
      market         = market_int,
      N_group        = N_group,
      group          = group_int,
      N_year         = N_year,
      year           = year_int,
      supply_model   = supply_model,
      use_cutoff     = as.integer(use_cutoff),
      grainsize      = 1L,
      N_market_year  = N_market_year,
      market_year_idx = mky_int,
      K_s0           = K_s0,
      X_s0           = X_s0_mat,
      mky_to_market  = array(mky_to_market),
      mky_to_year    = array(mky_to_year),
      is_single_market = is_single_market,
      use_rho        = as.integer(use_rho),
      use_hmt        = as.integer(use_hmt),
      avg_price_hmt  = avg_price_hmt,
      avg_margin_hmt = avg_margin_hmt,
      ssnip_hmt      = ssnip_hmt,
      implied_s0     = array(implied_s0)
    ),
    default_priors
  )

  # Attach metadata for later use
  attr(stan_data, "market_levels")  <- levels(market_f)
  attr(stan_data, "year_levels")    <- levels(year_f)
  attr(stan_data, "group_levels")   <- levels(group_f)
  attr(stan_data, "mky_levels")     <- levels(mky_f)
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
#' @param margins Numeric vector of price-cost margins.
#' @param implied_s0 Numeric vector of implied outside shares per market-year.
#' @param mky_int Integer vector mapping observations to market-years.
#' @param supply_model Integer 1--4.
#' @return A list with elements \code{mean} and \code{sd}.
#' @keywords internal
compute_alpha_prior <- function(shares, margins, implied_s0, mky_int,
                                supply_model) {
  s_low  <- stats::quantile(shares, 0.05)
  s_high <- stats::quantile(shares, 0.95)
  m_low  <- stats::quantile(margins, 0.05)
  m_high <- stats::quantile(margins, 0.95)

  if (supply_model == 1L) {
    # Bertrand: alpha ~ 1 / ((p-c) * (1 - s))
    a_low  <- 1 / (m_high * (1 - s_low))
    a_high <- 1 / (m_low  * (1 - s_high))
  } else if (supply_model == 2L) {
    # Monopolistic competition: alpha = 1 / (p-c)
    a_low  <- 1 / m_high
    a_high <- 1 / m_low
  } else if (supply_model == 3L) {
    # Cournot: alpha ~ (s0 + s) / (m * s0)
    s0_proxy <- max(0.01, stats::median(implied_s0))
    a_low  <- (s0_proxy + s_low)  / (m_high * s0_proxy)
    a_high <- (s0_proxy + s_high) / (m_low  * s0_proxy)
  } else {
    # Auction: alpha ~ -log(1 - s) / (m * s)
    a_low  <- -log(1 - s_low)  / (m_high * s_low)
    a_high <- -log(1 - s_high) / (m_low  * s_high)
  }

  # Cover [a_low, a_high] with ~95% mass
  p_mean <- (a_high + a_low) / 2
  p_sd   <- (a_high - a_low) / 4

  list(mean = p_mean, sd = p_sd)
}

#' Validate Inputs for bayes_est
#' @keywords internal
validate_inputs <- function(prices, shares, margins, market, year, group) {
  N <- length(prices)

  if (length(shares) != N) {
    stop("'shares' must have the same length as 'prices'.")
  }
  if (length(margins) != N) {
    stop("'margins' must have the same length as 'prices'.")
  }
  if (any(margins <= 0 | !is.finite(margins))) {
    stop("'margins' must be positive and finite.")
  }
  if (any(shares <= 0 | shares >= 1)) {
    stop("'shares' must be strictly between 0 and 1.")
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
