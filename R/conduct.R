# Conduct-aware margin equations used by antitrust's inherited equilibrium
# solver.  Keeping these as subclasses of Logit/CES preserves the upstream
# package's demand, elasticity, diversion, welfare, and HMT implementations.

conduct_weight <- function(object, firm_share) {
  if (!isTRUE(object@use_cutoff)) {
    return(rep(0, length(firm_share)))
  }
  stats::plogis(
    -(firm_share - object@cutoff_share) / object@cutoff_scale
  )
}

conduct_state <- function(object, preMerger, revenue = FALSE) {
  nprod <- length(object@shares)
  if (preMerger) {
    owner <- object@ownerPre
    prices <- object@pricePre
    subset <- rep(TRUE, nprod)
  } else {
    owner <- object@ownerPost
    prices <- object@pricePost
    subset <- object@subset
  }

  shares <- antitrust::calcShares(
    object,
    preMerger = preMerger,
    revenue = revenue
  )
  inside <- sum(shares, na.rm = TRUE)
  if (!is.finite(inside) || inside <= 0 || inside >= 1) {
    stop("Demand shares must imply a strictly positive inside and outside share.")
  }
  firm_share <- as.vector(owner %*% shares)
  list(
    shares = shares,
    firm_share = firm_share,
    firm_share_cond = firm_share / inside,
    outside_share = 1 - inside,
    prices = prices,
    subset = subset
  )
}

#' @importFrom antitrust calcMargins
setMethod(
  "calcMargins",
  signature(object = "BayesConductLogit"),
  function(object, preMerger = TRUE, exAnte = FALSE, level = FALSE) {
    state <- conduct_state(object, preMerger, revenue = FALSE)
    alpha <- abs(object@slopes$alpha)
    if (!is.finite(alpha) || alpha <= 0) {
      stop("The Logit price coefficient must be finite and nonzero.")
    }

    w <- conduct_weight(object, state$firm_share_cond)
    firm_share <- pmin(pmax(state$firm_share, 0), 1 - 1e-12)
    s0 <- state$outside_share

    conduct <- switch(
      as.character(object@supply_model),
      "1" = w + (1 - w) * (1 - firm_share),
      "2" = rep(1, length(firm_share)),
      "3" = w + (1 - w) * s0 / (s0 + firm_share),
      "4" = {
        auction_conduct <- rep(1, length(firm_share))
        positive <- firm_share > 1e-10
        auction_conduct[positive] <- firm_share[positive] /
          -log1p(-firm_share[positive])
        w + (1 - w) * auction_conduct
      },
      "5" = {
        beta <- object@bargaining_power
        if (length(beta) != 1L || !is.finite(beta) || beta < 0 || beta >= 1) {
          stop("Bargaining power must lie in [0, 1).")
        }
        kappa <- beta / (1 - beta)
        barg_conduct <- (
          kappa * firm_share - (1 - firm_share) * log1p(-firm_share)
        ) / -log1p(-firm_share)
        small <- firm_share <= 1e-10
        barg_conduct[small] <- 1
        w + (1 - w) * barg_conduct
      },
      stop("Unknown supply model in BayesConductLogit.")
    )

    margin <- 1 / (alpha * pmax(conduct, 1e-12))
    # Match antitrust's Auction2ndLogit API: ex-ante auction margins are
    # weighted by the product's quantity share.
    if (isTRUE(exAnte) && object@supply_model == 4L) {
      margin <- margin * state$shares
    }
    if (!level) margin <- margin / state$prices
    margin[!state$subset] <- NA_real_
    names(margin) <- object@labels
    margin
  }
)

setMethod(
  "calcMargins",
  signature(object = "BayesConductCES"),
  function(object, preMerger = TRUE, exAnte = FALSE, level = FALSE) {
    state <- conduct_state(object, preMerger, revenue = TRUE)
    gamma <- object@slopes$gamma
    if (!is.finite(gamma) || gamma <= 1) {
      stop("The CES elasticity parameter gamma must exceed one.")
    }

    w <- conduct_weight(object, state$firm_share_cond)
    firm_share <- pmin(pmax(state$firm_share, 0), 1 - 1e-12)

    margin_inv <- switch(
      as.character(object@supply_model),
      "1" = w * gamma + (1 - w) *
        (gamma - (gamma - 1) * firm_share),
      "2" = rep(gamma, length(firm_share)),
      "3" = w * gamma + (1 - w) * gamma /
        (1 + (gamma - 1) * firm_share),
      "4" = {
        auction_inv <- 1 /
          (1 - (1 - firm_share)^(1 / (gamma - 1)))
        w * gamma + (1 - w) * auction_inv
      },
      "5" = {
        beta <- object@bargaining_power
        if (length(beta) != 1L || !is.finite(beta) || beta < 0 || beta >= 1) {
          stop("Bargaining power must lie in [0, 1).")
        }
        bertrand_inv <- gamma - (gamma - 1) * firm_share
        w * gamma + (1 - w) * bertrand_inv / (1 - beta)
      },
      stop("Unknown supply model in BayesConductCES.")
    )

    margin <- 1 / pmax(margin_inv, 1e-12)
    # Match antitrust's Auction2ndCES API: ex-ante auction margins are
    # weighted by the product's revenue share.
    if (isTRUE(exAnte) && object@supply_model == 4L) {
      margin <- margin * state$shares
    }
    if (level) margin <- margin * state$prices
    margin[!state$subset] <- NA_real_
    names(margin) <- object@labels
    margin
  }
)
