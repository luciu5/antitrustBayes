#' Simulate Merger Effects from Bayesian Posterior
#'
#' Takes a fitted \code{\linkS4class{BayesLogit}} object and simulates the
#' price effects of a merger using posterior draws of demand parameters.
#' Each draw produces a complete merger simulation via
#' \code{\link[antitrust]{sim}}, yielding posterior distributions of price
#' changes, share changes, and other competition metrics.
#'
#' @param fit A \code{\linkS4class{BayesLogit}} object from
#'   \code{\link{bayes_est}}.
#' @param ownerPost Post-merger ownership: either a vector of firm IDs
#'   (same length as products in each market) or a named list of ownership
#'   matrices, one per market.
#' @param supply Character string for the \pkg{antitrust} supply model.
#'   Default \code{"bertrand"}.
#' @param demand Character string for the \pkg{antitrust} demand model.
#'   Default \code{"Logit"}.
#' @param n_draws Integer number of posterior draws to use. Default 200.
#' @param mcDelta Numeric vector of proportional marginal cost changes.
#'   Default \code{NULL} (no cost changes).
#' @param subset Logical vector indicating which products to include in
#'   post-merger simulation. Default \code{NULL} (all products).
#' @param insideSize Total inside market size. Default 1.
#' @param priceOutside Price of the outside good. Default 0.
#' @param labels Character vector of product labels. If \code{NULL},
#'   uses group identifiers from the estimation data.
#' @param market Which market(s) to simulate. Default \code{NULL} (all).
#' @param ... Additional arguments passed to \code{\link[antitrust]{sim}}.
#'
#' @details
#' For each posterior draw \eqn{i} and each market \eqn{m}, the function:
#' \enumerate{
#'   \item Extracts the price coefficient \eqn{\alpha^{(i)}_m} and outside
#'     share \eqn{s_0^{(i)}_m} from the posterior.
#'   \item Computes unconditional shares:
#'     \eqn{s_j = s_j^{cond} \times (1 - s_0)}.
#'   \item Inverts for mean valuations:
#'     \eqn{\delta_j = \log(s_j) - \log(s_0) - \alpha \times p_j}.
#'   \item Calls \code{\link[antitrust]{sim}} with these demand parameters.
#'   \item Extracts pre- and post-merger prices from the simulation object.
#' }
#'
#' @return A \code{\linkS4class{BayesSim}} object.
#'
#' @seealso \code{\link{bayes_est}}, \code{\link[antitrust]{sim}}
#'
#' @examples
#' \donttest{
#' # See bayes_est() for a complete PNB merger example
#' data(pnb)
#' # fit <- bayes_est(prices = pnb$prices, shares = pnb$shares,
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
                      supply = "bertrand",
                      demand = "Logit",
                      n_draws = 200L,
                      mcDelta = NULL,
                      subset = NULL,
                      insideSize = 1,
                      priceOutside = 0,
                      labels = NULL,
                      market = NULL,
                      ...) {

  cl <- match.call()

  if (!is(fit, "BayesLogit")) {
    stop("'fit' must be a BayesLogit object from bayes_est().")
  }

  draws <- rstan::extract(fit@stanfit)
  df    <- fit@data
  sd_data <- fit@stan_data
  price_sd <- fit@price_sd

  # Determine markets to simulate
  market_f <- as.factor(df$market)
  market_levels <- levels(market_f)
  if (!is.null(market)) {
    market_levels <- intersect(as.character(market), market_levels)
    if (length(market_levels) == 0L) {
      stop("None of the specified markets found in the fitted data.")
    }
  }

  n_total <- nrow(draws$a_market)
  draw_idx <- sample(seq_len(n_total), min(n_draws, n_total))
  n_draws_actual <- length(draw_idx)

  # Prepare market-year index mapping
  mky_labels <- paste(as.integer(as.factor(df$market)),
                      as.integer(as.factor(df$year)), sep = "_")
  mky_f <- as.factor(mky_labels)

  all_sims     <- list()
  price_pre    <- list()
  price_post   <- list()
  price_change <- list()
  shares_pre   <- list()
  shares_post  <- list()
  prod_labels  <- list()
  n_valid_vec  <- integer(length(market_levels))
  names(n_valid_vec) <- market_levels

  for (mk in market_levels) {
    mkt_idx <- which(df$market == mk)
    if (length(mkt_idx) == 0L) next

    mdata <- df[mkt_idx, , drop = FALSE]
    n_prod <- nrow(mdata)

    # Market index in the factor encoding
    market_int_val <- as.integer(as.factor(df$market))[mkt_idx[1]]

    # Market-year index (take the first one for this market)
    mky_int_val <- as.integer(mky_f)[mkt_idx[1]]

    # Ownership
    ownerPre_m <- diag(n_prod)
    if (is.list(ownerPost)) {
      ownerPost_m <- ownerPost[[mk]]
    } else if (is.matrix(ownerPost)) {
      ownerPost_m <- ownerPost
    } else {
      # Vector of firm IDs: convert to matrix
      owner_pre_vec <- as.character(mdata$group)
      if (length(ownerPost) == n_prod) {
        owner_post_vec <- as.character(ownerPost)
      } else {
        stop("'ownerPost' vector length does not match products in market '",
             mk, "'.")
      }
      ownerPost_m <- outer(owner_post_vec, owner_post_vec, "==") * 1.0
      ownerPre_m  <- outer(owner_pre_vec, owner_pre_vec, "==") * 1.0
    }

    prices_m <- as.numeric(mdata$prices)
    shares_cond <- as.numeric(mdata$shares)

    if (is.null(labels)) {
      labels_m <- as.character(mdata$group)
    } else if (is.list(labels)) {
      labels_m <- labels[[mk]]
    } else {
      labels_m <- labels[mkt_idx]
    }

    # mcDelta
    if (is.null(mcDelta)) {
      mcDelta_m <- rep(0, n_prod)
    } else if (is.list(mcDelta)) {
      mcDelta_m <- mcDelta[[mk]]
    } else {
      mcDelta_m <- mcDelta[mkt_idx]
    }

    # subset
    if (is.null(subset)) {
      subset_m <- rep(TRUE, n_prod)
    } else if (is.list(subset)) {
      subset_m <- subset[[mk]]
    } else {
      subset_m <- subset[mkt_idx]
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

      alpha_raw <- draws$a_market[i, market_int_val]
      alpha <- alpha_raw / price_sd

      s0 <- plogis(draws$s0_logit[i, mky_int_val])

      shares_uncond <- shares_cond * (1 - s0)

      meanval <- log(shares_uncond) - log(s0) - alpha * prices_m

      res <- tryCatch(
        antitrust::sim(
          prices       = prices_m,
          supply       = supply,
          demand       = demand,
          demand.param = list(alpha = alpha, meanval = meanval),
          ownerPre     = ownerPre_m,
          ownerPost    = ownerPost_m,
          mcDelta      = mcDelta_m,
          subset       = subset_m,
          insideSize   = insideSize,
          priceOutside = priceOutside,
          labels       = labels_m,
          ...
        ),
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

    all_sims[[mk]]     <- mkt_sims
    price_pre[[mk]]    <- ppre
    price_post[[mk]]   <- ppost
    price_change[[mk]] <- ppost / ppre - 1
    shares_pre[[mk]]   <- spre
    shares_post[[mk]]  <- spost
    prod_labels[[mk]]  <- labels_m
    n_valid_vec[mk]    <- n_ok
  }

  new("BayesSim",
    sims           = all_sims,
    price_pre      = price_pre,
    price_post     = price_post,
    price_change   = price_change,
    shares_pre     = shares_pre,
    shares_post    = shares_post,
    market_ids     = market_levels,
    product_labels = prod_labels,
    n_draws        = as.integer(n_draws_actual),
    n_valid        = n_valid_vec,
    call           = cl
  )
}
