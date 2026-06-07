#' Simulate Merger Effects from Bayesian Posterior
#'
#' Takes a fitted \code{\linkS4class{BayesLogit}} object and simulates the
#' price effects of a merger using posterior draws of demand parameters.
#' For each draw, the function constructs an \pkg{antitrust} class object
#' matching the estimated supply model, calibrates marginal costs via
#' \code{\link[antitrust]{calcMC}}, and solves for post-merger equilibrium
#' prices via \code{\link[antitrust]{calcPrices}}, yielding posterior
#' distributions of price changes, share changes, and other competition
#' metrics.
#'
#' @param fit A \code{\linkS4class{BayesLogit}} object from
#'   \code{\link{bayes_est}}.
#' @param ownerPost Post-merger ownership: either a vector of firm IDs
#'   (same length as products in each market) or a named list of ownership
#'   matrices, one per market.
#' @param n_draws Integer number of posterior draws to use for simulation.
#'   Default \code{NULL}, which uses \code{min(200, total_draws)}.
#' @param mcDelta Numeric vector of proportional marginal cost changes.
#'   Default \code{NULL} (no cost changes).
#' @param subset Logical vector indicating which products to include in
#'   post-merger simulation. Default \code{NULL} (all products).
#' @param insideSize Total inside market size. Default 1.
#' @param priceOutside Price of the outside good. Default 0.
#' @param labels Character vector of product labels. If \code{NULL},
#'   uses group identifiers from the estimation data.
#' @param market Which market(s) to simulate. Default \code{NULL} (all).
#' @param year Which year(s) to simulate. Default \code{NULL}, preserving
#'   market-level behavior for existing callers. When supplied, simulations
#'   are run for each selected market-year cell.
#' @param ... Additional arguments passed to
#'   \code{\link[antitrust]{calcPrices}}.
#'
#' @details
#' For each posterior draw \eqn{i} and each market \eqn{m}, the function:
#' \enumerate{
#'   \item Extracts the price coefficient \eqn{\alpha_m^{(i)}} and outside
#'     share \eqn{s_{0,m}^{(i)}} from the posterior.
#'   \item Computes unconditional shares:
#'     \eqn{s_j = s_j^{cond} \times (1 - s_0)}.
#'   \item Inverts for mean valuations:
#'     \eqn{\delta_j = \log(s_j) - \log(s_0) - \alpha \times p_j}.
#'   \item Constructs the appropriate \pkg{antitrust} class object
#'     (\code{Logit}, \code{LogitCournot}, or \code{Auction2ndLogit})
#'     matching the estimated supply model, with the draw's demand
#'     parameters set as \code{slopes}.
#'   \item Calibrates pre- and post-merger marginal costs using
#'     \code{\link[antitrust]{calcMC}}.
#'   \item Sets observed prices as the pre-merger equilibrium and solves
#'     for post-merger equilibrium prices using
#'     \code{\link[antitrust]{calcPrices}}.
#' }
#'
#' The supply model determines which \pkg{antitrust} class is used:
#' \describe{
#'   \item{\code{"bertrand"}}{Creates a \code{Logit} object (Nash-Bertrand).}
#'   \item{\code{"monopolistic"}}{Creates a \code{Logit} object
#'     (monopolistic competition).}
#'   \item{\code{"cournot"}}{Creates a \code{LogitCournot} object.}
#'   \item{\code{"auction"}}{Creates an \code{Auction2ndLogit} object.}
#' }
#'
#' @return A \code{\linkS4class{BayesSim}} object.
#'
#' @seealso \code{\link{bayes_est}},
#'   \code{\link[antitrust]{calcMC}},
#'   \code{\link[antitrust]{calcPrices}},
#'   \code{\link[antitrust]{calcShares}}
#'
#' @examples
#' \donttest{
#' # See bayes_est() for a complete PNB merger example
#' data(pnb)
#' # fit <- bayes_est(prices = pnb$prices, quantities = pnb$quantities,
#' #   margins = pnb$margins, market = rep("PNB", nrow(pnb)),
#' #   group = pnb$group, supply_model = "bertrand")
#' #
#' # Merge acquirer (298816) and target (889018)
#' # owner_post <- pnb$group
#' # owner_post[owner_post == "889018"] <- "298816"
#' # result <- bayes_sim(fit, ownerPost = owner_post, n_draws = 200)
#' # summary(result)
#' }
#'
#' @export
bayes_sim <- function(fit, ownerPost,
                      n_draws = NULL,
                      mcDelta = NULL,
                      subset = NULL,
                      insideSize = 1,
                      priceOutside = 0,
                      labels = NULL,
                      market = NULL,
                      year = NULL,
                      ...) {

  cl <- match.call()

  if (!is(fit, "BayesLogit")) {
    stop("'fit' must be a BayesLogit object from bayes_est().")
  }

  draws <- rstan::extract(fit@stanfit)
  df    <- fit@data
  price_sd <- fit@price_sd

  # Map supply_model integer to antitrust class name
  supply_class <- switch(as.character(fit@supply_model),
    "1" = "Logit",
    "2" = "Logit",
    "3" = "LogitCournot",
    "4" = "Auction2ndLogit"
  )
  is_auction <- fit@supply_model == 4L

  # Determine markets/cells to simulate
  market_f <- as.factor(df$market)
  market_levels <- levels(market_f)
  if (!is.null(market)) {
    market_levels <- intersect(as.character(market), market_levels)
    if (length(market_levels) == 0L) {
      stop("None of the specified markets found in the fitted data.")
    }
  }

  year_f <- as.factor(df$year)
  year_levels <- levels(year_f)
  if (!is.null(year)) {
    year_levels <- intersect(as.character(year), year_levels)
    if (length(year_levels) == 0L) {
      stop("None of the specified years found in the fitted data.")
    }
  }

  n_total <- nrow(draws$a_market)
  if (is.null(n_draws)) {
    n_draws <- min(200L, n_total)
  }
  n_draws <- min(as.integer(n_draws), n_total)
  draw_idx <- sample(seq_len(n_total), n_draws)
  n_draws_actual <- length(draw_idx)

  # Prepare market-year index mapping
  mky_labels <- paste(as.integer(as.factor(df$market)),
                      as.integer(as.factor(df$year)), sep = "_")
  mky_f <- as.factor(mky_labels)

  if (is.null(year)) {
    sim_units <- data.frame(
      unit_id = market_levels,
      market = market_levels,
      year = NA_character_,
      stringsAsFactors = FALSE
    )
  } else {
    sim_units <- unique(data.frame(
      market = as.character(df$market),
      year = as.character(df$year),
      stringsAsFactors = FALSE
    ))
    sim_units <- sim_units[
      sim_units$market %in% market_levels & sim_units$year %in% year_levels,
      ,
      drop = FALSE
    ]
    sim_units$unit_id <- paste(sim_units$market, sim_units$year, sep = "_")
    sim_units <- sim_units[, c("unit_id", "market", "year")]
    rownames(sim_units) <- NULL
    if (nrow(sim_units) == 0L) {
      stop("No market-year cells match the requested market/year filters.")
    }
  }

  all_sims     <- list()
  price_pre    <- list()
  price_post   <- list()
  price_change <- list()
  shares_pre   <- list()
  shares_post  <- list()
  prod_labels  <- list()
  n_valid_vec  <- integer(nrow(sim_units))
  names(n_valid_vec) <- sim_units$unit_id

  select_cell_arg <- function(x, idx, unit_id, mk) {
    if (is.list(x)) {
      if (!is.null(x[[unit_id]])) return(x[[unit_id]])
      return(x[[mk]])
    }
    if (is.matrix(x)) return(x)
    x[idx]
  }

  for (u in seq_len(nrow(sim_units))) {
    mk <- sim_units$market[u]
    yr <- sim_units$year[u]
    unit_id <- sim_units$unit_id[u]
    mkt_idx <- if (is.null(year)) {
      which(df$market == mk)
    } else {
      which(df$market == mk & as.character(df$year) == yr)
    }
    if (length(mkt_idx) == 0L) next

    mdata <- df[mkt_idx, , drop = FALSE]
    n_prod <- nrow(mdata)

    # Market index in the factor encoding
    market_int_val <- as.integer(as.factor(df$market))[mkt_idx[1]]

    # Market-year index (take the first one for this market)
    mky_int_val <- as.integer(mky_f)[mkt_idx[1]]

    # Ownership: pass vectors or matrices; ownerToMatrix handles conversion
    ownerPre_m <- as.character(mdata$group)
    ownerPost_m <- select_cell_arg(ownerPost, mkt_idx, unit_id, mk)
    if (!is.matrix(ownerPost_m)) {
      ownerPost_m <- as.character(ownerPost_m)
    }

    prices_m <- as.numeric(mdata$prices)
    quantities_m <- as.numeric(mdata$quantities)
    margins_m <- as.numeric(mdata$margins)

    # Compute quantity shares within this market
    shares_cond <- quantities_m / sum(quantities_m)

    if (is.null(labels)) {
      labels_m <- as.character(mdata$group)
    } else {
      labels_m <- select_cell_arg(labels, mkt_idx, unit_id, mk)
    }

    # mcDelta
    if (is.null(mcDelta)) {
      mcDelta_m <- rep(0, n_prod)
    } else {
      mcDelta_m <- select_cell_arg(mcDelta, mkt_idx, unit_id, mk)
    }

    # subset
    if (is.null(subset)) {
      subset_m <- rep(TRUE, n_prod)
    } else {
      subset_m <- select_cell_arg(subset, mkt_idx, unit_id, mk)
    }

    # Pre-allocate result matrices
    ppre  <- matrix(NA_real_, nrow = n_draws_actual, ncol = n_prod)
    ppost <- matrix(NA_real_, nrow = n_draws_actual, ncol = n_prod)
    spre  <- matrix(NA_real_, nrow = n_draws_actual, ncol = n_prod)
    spost <- matrix(NA_real_, nrow = n_draws_actual, ncol = n_prod)
    colnames(ppre) <- colnames(ppost) <- labels_m
    colnames(spre) <- colnames(spost) <- labels_m
    mkt_sims <- vector("list", n_draws_actual)
    n_ok <- 0L

    for (d in seq_along(draw_idx)) {
      i <- draw_idx[d]

      # Extract posterior draw parameters
      alpha_raw <- draws$a_market[i, market_int_val]
      alpha <- alpha_raw / price_sd

      s0 <- if ("s0_prob" %in% names(draws)) {
        draws$s0_prob[i, mky_int_val]
      } else if ("s0" %in% names(draws)) {
        plogis(draws$s0[i, mky_int_val])
      } else {
        plogis(draws$s0_logit[i, mky_int_val])
      }

      # Compute unconditional shares and mean valuations
      shares_uncond <- shares_cond * (1 - s0)
      meanval <- log(shares_uncond) - log(s0) - alpha * prices_m

      res <- tryCatch({

        # Construct the appropriate antitrust class object
        # (following the structure of antitrust::sim)
        obj <- new(supply_class,
          prices       = prices_m,
          shares       = shares_uncond,
          margins      = margins_m,
          weights      = rep(1, n_prod),
          ownerPre     = ownerPre_m,
          ownerPost    = ownerPost_m,
          mcDelta      = mcDelta_m,
          subset       = subset_m,
          insideSize   = insideSize,
          priceOutside = priceOutside,
          priceStart   = prices_m,
          shareInside  = 1 - s0,
          normIndex    = NA,
          labels       = labels_m
        )

        # Set demand slopes from posterior draw
        obj@slopes <- list(alpha = alpha, meanval = meanval)

        # Convert ownership vectors to matrices
        obj@ownerPre  <- ownerToMatrix(obj, TRUE)
        obj@ownerPost <- ownerToMatrix(obj, FALSE)

        # Calibrate marginal costs from FOC + observed prices
        obj@mcPre  <- calcMC(obj, preMerger = TRUE)
        obj@mcPost <- calcMC(obj, preMerger = FALSE)

        # Use observed prices as pre-merger equilibrium
        obj@pricePre <- prices_m

        # Solve for post-merger equilibrium prices
        # (auction models don't pass subset to calcPrices)
        if (is_auction) {
          obj@pricePost <- calcPrices(obj, preMerger = FALSE, ...)
        } else {
          obj@pricePost <- calcPrices(obj, preMerger = FALSE,
                                      subset = subset_m, ...)
        }

        # Set market size
        obj@mktSize <- insideSize / sum(calcShares(obj))

        obj
      },
        error = function(e) NULL
      )

      if (!is.null(res)) {
        n_ok <- n_ok + 1L
        mkt_sims[[d]] <- res
        ppre[d, ]  <- res@pricePre
        ppost[d, ] <- res@pricePost
        spre[d, ]  <- calcShares(res, preMerger = TRUE)
        spost[d, ] <- calcShares(res, preMerger = FALSE)
      }
    }

    all_sims[[unit_id]]     <- mkt_sims
    price_pre[[unit_id]]    <- ppre
    price_post[[unit_id]]   <- ppost
    price_change[[unit_id]] <- ppost / ppre - 1
    shares_pre[[unit_id]]   <- spre
    shares_post[[unit_id]]  <- spost
    prod_labels[[unit_id]]  <- labels_m
    n_valid_vec[unit_id]    <- n_ok
  }

  new("BayesSim",
    sims           = all_sims,
    price_pre      = price_pre,
    price_post     = price_post,
    price_change   = price_change,
    shares_pre     = shares_pre,
    shares_post    = shares_post,
    market_ids     = sim_units$unit_id,
    product_labels = prod_labels,
    n_draws        = as.integer(n_draws_actual),
    n_valid        = n_valid_vec,
    call           = cl
  )
}
