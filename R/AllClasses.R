setOldClass(c("bayes_params", "list"))
setOldClass(c("bayes_market", "list"))

#' @title Bayesian structural fit classes
#' @name BayesFit-class
#' @aliases BayesFit BayesLogit BayesSpecifiedFit
#'
#' @description \code{BayesFit} is the shared virtual representation for a
#' Bayesian distribution over structural models. It extends
#' \code{antitrust::StructuralFit}; \code{BayesLogit} stores a Stan-estimated
#' posterior and \code{BayesSpecifiedFit} stores directly supplied draws.
#'
#' @slot data A \code{data.frame} containing the observed market data.
#' @slot supply_model Integer identifying the supply model.
#' @slot demand_model Character demand-system identifier.
#' @slot market_side Character market-side identifier.
#' @slot call The matched construction call.
#' @slot metadata Bayesian provenance and prior metadata.
#' @slot stanfit A \code{stanfit} object from \pkg{rstan}.
#' @slot stan_data A \code{list} of the data passed to Stan.
#' @slot price_sd Numeric scalar: the standard deviation of prices (logit) or
#'   log-prices (CES) used for standardization.
#' @slot params A \code{bayes_params} list for a specified fit.
#' @slot market A \code{bayes_market} list for a specified fit.
#'
#' @export
setClass(
  "BayesFit",
  representation(
    data         = "data.frame",
    supply_model = "integer",
    demand_model = "character",
    market_side  = "character",
    call         = "call",
    metadata     = "list"
  ),
  contains = c("StructuralFit", "VIRTUAL"),
  prototype = list(
    model = NULL,
    data = data.frame(),
    supply_model = integer(),
    demand_model = character(),
    market_side = character(),
    metadata = list()
  )
)

#' @rdname BayesFit-class
#' @export
setClass(
  "BayesLogit",
  representation(
    stanfit = "ANY",
    stan_data = "list",
    price_sd = "numeric"
  ),
  contains = "BayesFit"
)

# Directly specified parameter draws. The class is intentionally internal;
# users construct it with bayes_specify() and interact with it through the
# exported BayesFit contract.
setClass(
  "BayesSpecifiedFit",
  representation(
    params = "list",
    market = "list"
  ),
  contains = "BayesFit"
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
    cutoff_scale = "numeric",
    bargaining_power = "numeric"
  )
)

setClass(
  "BayesConductCES",
  contains = "CES",
  slots = c(
    supply_model = "integer",
    use_cutoff = "logical",
    cutoff_share = "numeric",
    cutoff_scale = "numeric",
    bargaining_power = "numeric"
  )
)

#' @title S4 Class for Bayesian Counterfactual Results
#' @name BayesSim-class
#' @aliases BayesSim
#'
#' @description S4 class holding counterfactual results computed over draws
#' from a \code{\linkS4class{BayesFit}} object.
#'
#' @slot sims A \code{list} of \code{antitrust} simulation objects, one per
#'   valid posterior draw and market.
#' @slot price_pre A \code{list} of matrices (draws x products) of baseline
#'   prices, one element per market.
#' @slot price_post A \code{list} of matrices (draws x products) of
#'   post-counterfactual prices, one element per market.
#' @slot price_change A \code{list} of matrices (draws x products) of
#'   proportional price changes, one element per market.
#' @slot shares_pre A \code{list} of matrices (draws x products) of
#'   pre-merger shares.
#' @slot shares_post A \code{list} of matrices (draws x products) of
#'   post-counterfactual shares.
#' @slot market_ids Character vector of market identifiers.
#' @slot product_labels A \code{list} of character vectors of product labels,
#'   one per market.
#' @slot n_draws Integer: number of posterior draws used.
#' @slot n_valid A named integer vector: number of successful simulations per
#'   market.
#' @slot failures A named list containing simulation failure messages by
#'   market and draw.
#' @slot call The matched call to \code{simulate()} or \code{\link{bayes_sim}}.
#' @slot metadata Provenance and source-model metadata.
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
    call           = "call",
    metadata       = "list"
  ),
  prototype = list(metadata = list())
)
