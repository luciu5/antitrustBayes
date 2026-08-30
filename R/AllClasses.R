#' @title S4 Classes for antitrustBayes
#' @name BayesLogit-class
#' @aliases BayesLogit
#'
#' @description S4 class holding the results of Bayesian estimation of a
#' logit demand and supply model via Stan.
#'
#' @slot stanfit A \code{stanfit} object from \pkg{rstan}.
#' @slot data A \code{data.frame} of the original input data.
#' @slot stan_data A \code{list} of the data passed to Stan.
#' @slot supply_model Integer (1--4) identifying the supply model.
#' @slot demand_model Character: \code{"logit"} or \code{"ces"}.
#' @slot market_side Character: \code{"input"} or \code{"output"}.
#' @slot price_sd Numeric scalar: the standard deviation of prices (logit) or
#'   log-prices (CES) used for standardization.
#' @slot call The matched call to \code{\link{bayes_est}}.
#'
#' @export
setClass(
  "BayesLogit",
  representation(
    stanfit      = "ANY",
    data         = "data.frame",
    stan_data    = "list",
    supply_model = "integer",
    demand_model = "character",
    market_side  = "character",
    price_sd     = "numeric",
    call         = "call"
  )
)

# Internal antitrust-compatible classes.  They inherit antitrust's demand,
# elasticity, diversion, welfare, and equilibrium methods, while replacing
# only the supply-margin equation with the conduct model estimated in Stan.
setClass(
  "BayesConductLogit",
  contains = "Logit",
  slots = c(
    supply_model = "integer",
    use_cutoff = "logical",
    cutoff_share = "numeric",
    cutoff_scale = "numeric"
  )
)

setClass(
  "BayesConductCES",
  contains = "CES",
  slots = c(
    supply_model = "integer",
    use_cutoff = "logical",
    cutoff_share = "numeric",
    cutoff_scale = "numeric"
  )
)

#' @title S4 Class for Bayesian Merger Simulation Results
#' @name BayesSim-class
#' @aliases BayesSim
#'
#' @description S4 class holding the results of merger simulations computed
#' over posterior draws from a \code{\linkS4class{BayesLogit}} object.
#'
#' @slot sims A \code{list} of \code{antitrust} simulation objects, one per
#'   valid posterior draw and market.
#' @slot price_pre A \code{list} of matrices (draws x products) of pre-merger
#'   prices, one element per market.
#' @slot price_post A \code{list} of matrices (draws x products) of
#'   post-merger prices, one element per market.
#' @slot price_change A \code{list} of matrices (draws x products) of
#'   proportional price changes, one element per market.
#' @slot shares_pre A \code{list} of matrices (draws x products) of
#'   pre-merger shares.
#' @slot shares_post A \code{list} of matrices (draws x products) of
#'   post-merger shares.
#' @slot market_ids Character vector of market identifiers.
#' @slot product_labels A \code{list} of character vectors of product labels,
#'   one per market.
#' @slot n_draws Integer: number of posterior draws used.
#' @slot n_valid A named integer vector: number of successful simulations per
#'   market.
#' @slot failures A named list containing simulation failure messages by
#'   market and draw.
#' @slot call The matched call to \code{\link{bayes_sim}}.
#'
#' @export
setClass(
  "BayesSim",
  representation(
    sims           = "list",
    price_pre      = "list",
    price_post     = "list",
    price_change   = "list",
    shares_pre     = "list",
    shares_post    = "list",
    market_ids     = "character",
    product_labels = "list",
    n_draws        = "integer",
    n_valid        = "integer",
    failures       = "list",
    call           = "call"
  )
)
