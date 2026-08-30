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
#' @param insideSize Total inside quantity. For CES, this is converted to
#'   inside expenditure using the observed quantity-weighted price. Default 1.
#' @param priceOutside Optional positive price of the outside good. By default,
#'   uses the market-year outside price stored in the fitted model.
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
#'   \item Extracts the Logit price coefficient or CES elasticity and the
#'     outside quantity share from the posterior.
#'   \item Computes unconditional shares:
#'     \eqn{s_j = s_j^{cond} \times (1 - s_0)}.
#'   \item Inverts either the Logit or CES demand equation for mean values,
#'     converting CES quantity shares to revenue shares where required by
#'     \pkg{antitrust}.
#'   \item Constructs an internal conduct-aware subclass of \pkg{antitrust}'s
#'     \code{Logit} or \code{CES} class with the draw's demand and cutoff
#'     parameters.
#'   \item Calibrates pre- and post-merger marginal costs using
#'     \code{\link[antitrust]{calcMC}}.
#'   \item Sets observed prices as the pre-merger equilibrium and solves
#'     for post-merger equilibrium prices using
#'     \code{\link[antitrust]{calcPrices}}.
#' }
#'
#' Internal subclasses of \pkg{antitrust}'s \code{Logit} and \code{CES}
#' classes override only the margin equation. This lets all five estimated
#' conduct models (including the posterior fringe cutoff) use the upstream
#' equilibrium, elasticity, diversion, welfare, and HMT machinery.
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
                      priceOutside = NULL,
                      labels = NULL,
                      market = NULL,
                      year = NULL,
                      ...) {

  cl <- match.call()

  if (!is(fit, "BayesLogit")) {
    stop("'fit' must be a BayesLogit object from bayes_est().")
  }
  validate_fit_spec(fit)

  draws <- rstan::extract(fit@stanfit)
  df    <- fit@data
  price_sd <- fit@price_sd
  is_output <- identical(fit@market_side, "output")
  is_ces <- identical(fit@demand_model, "ces")

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

  if (is.null(year)) {
    year_counts <- tapply(as.character(df$year), df$market,
                          function(x) length(unique(x)))
    if (any(year_counts[market_levels] > 1L, na.rm = TRUE)) {
      stop("year must be supplied when simulating a panel fit; otherwise ",
           "different years would be treated as simultaneous products.",
           call. = FALSE)
    }
  }

  n_total <- dim(draws$a_market)[1]
  if (is.null(n_total)) n_total <- length(draws$a_market)
  if (length(n_total) != 1L || !is.finite(n_total) || n_total < 1L) {
    stop("The fitted object contains no posterior draws.", call. = FALSE)
  }
  if (is.null(n_draws)) {
    n_draws <- min(200L, n_total)
  }
  if (length(n_draws) != 1L || !is.finite(n_draws) || n_draws < 1 ||
      n_draws != as.integer(n_draws)) {
    stop("'n_draws' must be a positive integer.")
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
  failures     <- list()
  n_valid_vec  <- integer(nrow(sim_units))
  names(n_valid_vec) <- sim_units$unit_id

  select_cell_arg <- function(x, idx, unit_id, mk, scalar = FALSE) {
    if (is.list(x)) {
      if (!is.null(x[[unit_id]])) return(x[[unit_id]])
      return(x[[mk]])
    }
    if (is.matrix(x)) return(x)
    if (scalar && length(x) == 1L) return(x)
    if (length(x) == nrow(df)) return(x[idx])
    x
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
    if (is.null(ownerPost_m)) {
      stop("No post-merger ownership supplied for simulation unit '",
           unit_id, "'.")
    }
    if (!is.matrix(ownerPost_m)) {
      ownerPost_m <- as.character(ownerPost_m)
    }

    prices_m <- as.numeric(mdata$prices)
    quantities_m <- as.numeric(mdata$quantities)
    margins_m <- as.numeric(mdata$margins)
    margins_m[!is.finite(margins_m)] <- 0.5

    # Compute quantity shares within this market
    shares_cond <- quantities_m / sum(quantities_m)

    if (is.null(labels)) {
      labels_m <- as.character(mdata$group)
    } else {
      labels_m <- select_cell_arg(labels, mkt_idx, unit_id, mk)
    }
    if (length(labels_m) != n_prod) {
      stop("'labels' must have one value per product in '", unit_id, "'.")
    }

    # mcDelta
    if (is.null(mcDelta)) {
      mcDelta_m <- rep(0, n_prod)
    } else {
      mcDelta_m <- select_cell_arg(mcDelta, mkt_idx, unit_id, mk)
    }
    if (length(mcDelta_m) != n_prod || any(!is.finite(mcDelta_m)) ||
        any(mcDelta_m <= -1)) {
      stop("'mcDelta' must contain finite proportional changes greater than -1 ",
           "for every product in '", unit_id, "'.")
    }

    # subset
    if (is.null(subset)) {
      subset_m <- rep(TRUE, n_prod)
    } else {
      subset_m <- select_cell_arg(subset, mkt_idx, unit_id, mk)
    }
    subset_m <- as.logical(subset_m)
    if (length(subset_m) != n_prod || anyNA(subset_m) || !any(subset_m)) {
      stop("'subset' must select at least one product in '", unit_id, "'.")
    }

    outside_price_m <- if (is.null(priceOutside)) {
      unique(mdata$outside_price)
    } else {
      select_cell_arg(priceOutside, mkt_idx, unit_id, mk, scalar = TRUE)
    }
    if (length(outside_price_m) != 1L || !is.finite(outside_price_m) ||
        outside_price_m <= 0) {
      stop("'priceOutside' must resolve to one positive finite value for '",
           unit_id, "'.")
    }
    inside_size_m <- select_cell_arg(insideSize, mkt_idx, unit_id, mk,
                                     scalar = TRUE)
    if (length(inside_size_m) != 1L || !is.finite(inside_size_m) ||
        inside_size_m <= 0) {
      stop("'insideSize' must resolve to one positive finite value for '",
           unit_id, "'.")
    }

    if (is.matrix(ownerPost_m)) {
      if (!identical(dim(ownerPost_m), c(n_prod, n_prod))) {
        stop("An ownership matrix for '", unit_id, "' must be ", n_prod,
             " by ", n_prod, ".")
      }
    } else if (length(ownerPost_m) != n_prod) {
      stop("'ownerPost' must have one owner per product in '", unit_id, "'.")
    }

    # Pre-allocate result matrices
    ppre  <- matrix(NA_real_, nrow = n_draws_actual, ncol = n_prod)
    ppost <- matrix(NA_real_, nrow = n_draws_actual, ncol = n_prod)
    spre  <- matrix(NA_real_, nrow = n_draws_actual, ncol = n_prod)
    spost <- matrix(NA_real_, nrow = n_draws_actual, ncol = n_prod)
    colnames(ppre) <- colnames(ppost) <- labels_m
    colnames(spre) <- colnames(spost) <- labels_m
    mkt_sims <- vector("list", n_draws_actual)
    mkt_failures <- rep(NA_character_, n_draws_actual)
    n_ok <- 0L

    for (d in seq_along(draw_idx)) {
      i <- draw_idx[d]

      # Extract posterior draw parameters
      a_raw <- if (is.null(dim(draws$a_market))) {
        draws$a_market[i]
      } else {
        draws$a_market[i, market_int_val]
      }

      s0 <- if ("s0_prob" %in% names(draws)) {
        draws$s0_prob[i, mky_int_val]
      } else if ("s0" %in% names(draws)) {
        plogis(draws$s0[i, mky_int_val])
      } else {
        plogis(draws$s0_logit[i, mky_int_val])
      }

      bargaining_draw <- if (fit@supply_model == 5L) {
        as.numeric(draws$bargaining_power[i])
      } else {
        0
      }

      # Compute unconditional quantity shares and invert demand at observed
      # prices. These parameterizations reproduce the fitted shares exactly
      # under antitrust::calcShares().
      shares_uncond <- shares_cond * (1 - s0)

      res <- tryCatch({
        if (is_ces) {
          gamma <- 1 + a_raw / price_sd
          meanval <- (shares_uncond / s0) *
            (prices_m / outside_price_m)^gamma
          total_expenditure_index <- s0 * outside_price_m +
            sum(prices_m * shares_uncond)
          R0 <- s0 * outside_price_m / total_expenditure_index
          shares_for_object <- prices_m * shares_uncond /
            total_expenditure_index
          share_inside_for_object <- 1 - R0
          inside_size_for_object <- inside_size_m *
            sum(prices_m * shares_cond)
          class_name <- "BayesConductCES"
          slopes <- list(gamma = gamma, meanval = meanval,
                         alpha = R0 / (1 - R0))
        } else {
          alpha <- a_raw / price_sd
          if (is_output) alpha <- -alpha
          meanval <- log(shares_uncond / s0) -
            alpha * (prices_m - outside_price_m)
          shares_for_object <- shares_uncond
          share_inside_for_object <- 1 - s0
          inside_size_for_object <- inside_size_m
          class_name <- "BayesConductLogit"
          slopes <- list(alpha = alpha, meanval = meanval)
        }

        cutoff_draw <- if (isTRUE(fit@stan_data$use_cutoff == 1L)) {
          if (is.null(dim(draws$cutoff_share))) {
            draws$cutoff_share[i]
          } else {
            draws$cutoff_share[i, market_int_val]
          }
        } else {
          0
        }

        obj <- new(class_name,
          prices       = prices_m,
          shares       = shares_for_object,
          margins      = margins_m,
          weights      = rep(1, n_prod),
          ownerPre     = ownerPre_m,
          ownerPost    = ownerPost_m,
          mcDelta      = mcDelta_m,
          subset       = subset_m,
          insideSize   = inside_size_for_object,
          priceOutside = outside_price_m,
          priceStart   = prices_m,
          pricePre     = prices_m,
          pricePost    = prices_m,
          shareInside  = share_inside_for_object,
          normIndex    = NA,
          output       = is_output,
          labels       = labels_m,
          supply_model = fit@supply_model,
          use_cutoff   = isTRUE(fit@stan_data$use_cutoff == 1L),
          cutoff_share = as.numeric(cutoff_draw),
          cutoff_scale = 0.2,
          bargaining_power = bargaining_draw
        )

        obj@slopes <- slopes

        # Convert ownership vectors to matrices
        obj@ownerPre  <- ownerToMatrix(obj, TRUE)
        obj@ownerPost <- ownerToMatrix(obj, FALSE)

        # Calibrate marginal costs from FOC + observed prices
        obj@mcPre  <- calcMC(obj, preMerger = TRUE)
        obj@mcPost <- calcMC(obj, preMerger = FALSE)

        # Use observed prices as pre-merger equilibrium
        obj@pricePre <- prices_m

        # Under monopolistic competition, ownership alone cannot change the
        # equilibrium. Preserve this implication exactly rather than asking a
        # nonlinear solver to rediscover the observed-price root.
        ownership_only_moncom <- fit@supply_model == 2L &&
          all(mcDelta_m == 0) && all(subset_m)
        if (ownership_only_moncom) {
          obj@pricePost <- prices_m
        } else {
          obj@pricePost <- calcPrices(obj, preMerger = FALSE, ...)
        }

        obj@mktSize <- inside_size_for_object / share_inside_for_object

        obj
      },
        error = function(e) {
          mkt_failures[d] <<- conditionMessage(e)
          NULL
        }
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
    failures[[unit_id]]     <- mkt_failures
    n_valid_vec[unit_id]    <- n_ok
  }

  n_failed <- sum(vapply(failures, function(x) sum(!is.na(x)), integer(1)))
  if (n_failed > 0L) {
    warning(n_failed, " of ", n_draws_actual * nrow(sim_units),
            " market-draw simulations failed. Inspect result@failures for details.",
            call. = FALSE)
  }

  sim_metadata <- fit@metadata
  sim_metadata$source_supply_model <- fit@supply_model
  sim_metadata$posterior_draws_requested <- n_draws_actual

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
    failures       = failures,
    call           = cl,
    metadata       = sim_metadata
  )
}
