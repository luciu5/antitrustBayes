#' Bayesian Estimation of Logit Demand and Supply
#'
#' Estimates a hierarchical Bayesian model of differentiated products demand
#' (Logit) and supply using Stan MCMC. The model jointly estimates demand
#' parameters (price sensitivity, mean valuations, outside share) and supply
#' parameters under one of four supply models.
#'
#' @param prices Numeric vector of product prices.
#' @param shares Numeric vector of within-market product shares. Must be
#'   strictly between 0 and 1.
#' @param margins Numeric vector of product price-cost margins
#'   (price / marginal cost). Must be positive.
#' @param market Factor or vector identifying markets. Each unique value
#'   defines a separate geographic or product market.
#' @param year Factor or vector identifying time periods. If \code{NULL}
#'   (default), all observations are treated as a single period.
#' @param group Factor or vector identifying product groups (e.g., parent
#'   firms or holding companies). Used for hierarchical mean-value effects.
#' @param supply_model Character string specifying the supply model:
#'   \code{"bertrand"} (default), \code{"monopolistic"}, \code{"cournot"},
#'   or \code{"competitive"}.
#' @param use_cutoff Logical. If \code{TRUE}, estimate a fringe-firm share
#'   cutoff. Default \code{FALSE}.
#' @param use_rho Logical. If \code{TRUE}, estimate share-margin error
#'   correlation. Default \code{TRUE}.
#' @param use_hmt Logical. If \code{TRUE}, apply a soft Hypothetical
#'   Monopolist Test penalty. Default \code{FALSE}.
#' @param priors Named list of prior hyperparameters. See Details.
#' @param stan_options Named list of arguments passed to
#'   \code{\link[rstan]{sampling}} (e.g., \code{chains}, \code{iter},
#'   \code{warmup}, \code{cores}).
#' @param ... Additional arguments (currently unused).
#'
#' @details
#' The model estimates a Logit demand system jointly with a supply-side
#' first-order condition. Four supply models are supported:
#' \describe{
#'   \item{1. \code{"bertrand"}}{Nash-Bertrand differentiated products pricing.}
#'   \item{2. \code{"monopolistic"}}{Monopolistic competition.}
#'   \item{3. \code{"cournot"}}{Cournot quantity competition.}
#'   \item{4. \code{"competitive"}}{Perfectly competitive (no markup).}
#' }
#'
#' The \code{priors} argument accepts a named list with any of:
#' \describe{
#'   \item{\code{prior_sigma_share}}{Half-normal scale for share error SD
#'     (default 0.5).}
#'   \item{\code{prior_sigma_margin}}{Half-normal scale for margin error SD
#'     (default 0.5).}
#'   \item{\code{prior_sigma_meanval_strat}}{Half-normal scale for strategic
#'     group mean-value SD (default 1.0).}
#'   \item{\code{prior_alpha_mean}}{Prior mean for log(alpha) on price
#'     sensitivity (default 1.0).}
#'   \item{\code{prior_alpha_sd}}{Prior SD for log(alpha) (default 1.0).}
#'   \item{\code{prior_lkj}}{LKJ prior parameter for share-margin correlation
#'     (default 2.0).}
#' }
#'
#' @return A \code{\linkS4class{BayesLogit}} object.
#'
#' @seealso \code{\link{bayes_sim}} to simulate merger effects from the
#'   fitted model.
#'
#' @examples
#' \donttest{
#' # Philadelphia National Bank example (single market, 39 banks)
#' data(pnb)
#'
#' # Estimate Bayesian Logit demand + Bertrand supply
#' fit <- bayes_est(
#'   prices  = pnb$prices,
#'   shares  = pnb$shares,
#'   margins = pnb$margins,
#'   market  = rep("PNB", nrow(pnb)),
#'   group   = pnb$group,
#'   supply_model = "bertrand",
#'   stan_options = list(chains = 2, iter = 1000, warmup = 500)
#' )
#' fit
#'
#' # Simulate the PNB-Girard Trust merger
#' # Build post-merger ownership: merge acquirer (298816) and target (889018)
#' owner_post <- pnb$group
#' owner_post[owner_post == "889018"] <- "298816"
#'
#' result <- bayes_sim(fit, ownerPost = owner_post, n_draws = 200)
#' summary(result)
#' }
#'
#' @export
bayes_est <- function(prices, shares, margins,
                      market, year = NULL, group,
                      supply_model = c("bertrand", "monopolistic",
                                       "cournot", "competitive"),
                      use_cutoff = FALSE,
                      use_rho = TRUE,
                      use_hmt = FALSE,
                      priors = list(),
                      stan_options = list(),
                      ...) {

  cl <- match.call()
  supply_model <- match.arg(supply_model)
  supply_int <- switch(supply_model,
    bertrand      = 1L,
    monopolistic  = 2L,
    cournot       = 3L,
    competitive   = 4L
  )

  # Default year if not supplied
  if (is.null(year)) {
    year <- rep(1L, length(prices))
  }

  # Validate

  validate_inputs(prices, shares, margins, market, year, group)

  # Build data frame
  df <- data.frame(
    prices  = prices,
    shares  = shares,
    margins = margins,
    market  = market,
    year    = year,
    group   = group,
    stringsAsFactors = FALSE
  )

  # Prepare Stan data
  stan_data <- prepare_stan_data(
    prices       = prices,
    shares       = shares,
    margins      = margins,
    market       = market,
    year         = year,
    group        = group,
    supply_model = supply_int,
    use_cutoff   = use_cutoff,
    use_rho      = use_rho,
    use_hmt      = use_hmt,
    priors       = priors
  )

  # Get pre-compiled model
  stan_mod <- stanmodels[["logit_supply"]]

  # Default Stan options
  default_stan <- list(
    chains  = 4L,
    iter    = 2000L,
    warmup  = 1000L,
    cores   = 1L,
    refresh = 200L
  )
  for (nm in names(stan_options)) {
    default_stan[[nm]] <- stan_options[[nm]]
  }

  # Run sampling
  stanfit <- do.call(rstan::sampling, c(
    list(object = stan_mod, data = stan_data),
    default_stan
  ))

  new("BayesLogit",
    stanfit      = stanfit,
    data         = df,
    stan_data    = stan_data,
    supply_model = supply_int,
    price_sd     = stan_data$price_sd,
    call         = cl
  )
}
