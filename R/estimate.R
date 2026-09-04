#' Bayesian Estimation of Demand and Supply
#'
#' Estimates a hierarchical Bayesian model of differentiated products demand
#' (logit or CES) and supply using Stan MCMC. The model jointly estimates
#' demand parameters (price sensitivity, mean valuations, outside share) and
#' supply parameters under one of five supply models.
#'
#' @param prices Numeric vector of product prices. Must be positive.
#' @param quantities Numeric vector of product quantities (e.g., units sold,
#'   deposits). Must be positive. Within each market-year, the package
#'   normalizes these into conditional inside quantity shares for both Logit
#'   and CES. The outside-share level is a
#'   model parameter identified jointly with the supply FOCs.
#' @param margins Numeric vector of positive price-normalized margins. For
#'   output markets this is the usual Lerner index \eqn{(p - c) / p}; for
#'   input markets it may be an input markdown ratio \eqn{(v - w) / w}.
#'   Values greater than 1 are valid for input markets such as deposits.
#' @param market Factor or vector identifying markets. Each unique value
#'   defines a separate geographic or product market.
#' @param year Factor or vector identifying time periods. If \code{NULL}
#'   (default), all observations are treated as a single period.
#' @param group Factor or vector identifying product groups (e.g., parent
#'   firms or holding companies). Used for hierarchical mean-value effects.
#' @param demand_model Character string specifying the demand model:
#'   \code{"logit"} (default) or \code{"ces"} (Constant Elasticity of
#'   Substitution).
#' @param market_side Character string specifying whether prices are input
#'   prices (\code{"input"}, default for Logit; positive price coefficient)
#'   or output prices (\code{"output"}, negative price coefficient). CES
#'   currently supports output markets only; if omitted, its default is
#'   \code{"output"}.
#' @param supply_model Character string specifying the supply model:
#'   \code{"bertrand"} (default), \code{"monopolistic"}, \code{"cournot"},
#'   \code{"auction"}, or \code{"bargaining"}.
#' @param use_cutoff Logical. If \code{TRUE}, estimate a fringe-firm share
#'   cutoff. Default \code{FALSE}.
#' @param use_rho Logical. If \code{TRUE} (default), couple the share and
#'   margin likelihoods per observation through a free cross-equation
#'   correlation \code{rho}, using a bivariate Student-t factored as
#'   \code{p(share) * p(margin | share)}. \code{rho} is only ever
#'   instantiated when \code{likelihood = "joint"} and at least one
#'   observation has an observed margin; otherwise it is structurally fixed
#'   at zero (never sampled). Note that \code{likelihood = "share"} omits the
#'   margin equation and therefore the equilibrium-margin information that
#'   principally identifies the outside share's level -- under share-only
#'   fits, the outside share may be weakly identified or prior-dominated,
#'   depending on the market/time/intercept structure of the data. This is an
#'   expected property of the model, not a defect to be tuned away.
#' @param df Degrees of freedom for the joint bivariate (or, when \code{rho}
#'   is inactive, marginal) Student-t likelihood. Either a single finite
#'   numeric greater than 2 (fixed, default \code{5}, preferred for sparse
#'   information sets) or the character sentinel \code{"estimate"}, which
#'   places a \code{gamma(prior_nu_shape, prior_nu_rate)} prior on
#'   \code{nu - 2} and estimates \code{nu} jointly with the rest of the
#'   model.
#' @param use_market_resid_scale Logical. If \code{TRUE}, allow market-level
#'   share and inverse-margin residual scales with hierarchical shrinkage.
#'   Default \code{FALSE}.
#' @param use_alpha_market_pool Logical. If \code{TRUE}, allow market-level
#'   logit price coefficients with hierarchical shrinkage. Default \code{TRUE}.
#' @param use_s0_market_year_pool Logical. If \code{TRUE}, allow market-year
#'   outside shares with hierarchical shrinkage (one outside share per
#'   market-year). Default \code{FALSE}. Mutually exclusive with
#'   \code{use_s0_market_pool}. Note that market-year pooling introduces one
#'   weakly-identified outside share per market-year cell; when a market spans
#'   several periods with sparse per-cell information, its hyper-scale can be
#'   poorly identified (low BFMI). Prefer \code{use_s0_market_pool} in that
#'   regime.
#' @param use_s0_market_pool Logical. If \code{TRUE}, pool the outside share at
#'   the \emph{market} level: one outside share per market, shared across all of
#'   that market's periods, with hierarchical shrinkage across markets
#'   (structurally parallel to \code{use_alpha_market_pool}). Better identified
#'   than \code{use_s0_market_year_pool} when markets are observed over multiple
#'   periods, because each market-level effect is disciplined by all of the
#'   market's observations. Default \code{FALSE}. Mutually exclusive with
#'   \code{use_s0_market_year_pool}.
#' @param use_centered_alpha Logical. If \code{TRUE}, use the centered
#'   parameterization for market-level price coefficients. Default \code{FALSE}.
#' @param use_centered_s0_pool Logical. If \code{TRUE}, use the centered
#'   parameterization for market-year outside-share pooling. Default \code{FALSE}.
#' @param use_hmt Logical. If \code{TRUE}, apply a soft Hypothetical
#'   Monopolist Test penalty. Default \code{FALSE}.
#' @param X_s0 Optional model matrix of market-year covariates for the outside
#'   share equation. Rows must match the market-year combinations generated by
#'   the input data after factor encoding.
#' @param outside_price Positive outside-good price used to put CES outside
#'   quantity and revenue shares on a common scale. Supply a scalar (default
#'   \code{1}) or one value per market-year. Stored for Logit fits as part of
#'   the common Logit/CES comparison contract, but otherwise unused by Logit.
#' @param prior_manifest Prespecified proper prior manifest: \code{"baseline"},
#'   \code{"tight"}, or \code{"diffuse"}. User-supplied \code{priors}
#'   override individual entries in the selected manifest.
#' @param priors Named list of prior hyperparameters. See Details.
#' @param likelihood Which likelihood block to fit: \code{"joint"} (default),
#'   \code{"share"}, or \code{"margin"}. Share-only fits permit an exact
#'   marginal-likelihood decomposition.
#' @param loo_cluster Optional vector identifying the natural predictive
#'   cluster (for example, a Supreme Court case). Defaults to \code{market}.
#' @param model_prior Optional scalar recording the prespecified prior
#'   probability assigned to this conduct model. It is provenance metadata;
#'   posterior model odds are computed by \code{bayes_compare}.
#' @param stan_options Named list of arguments passed to
#'   \code{\link[rstan]{sampling}} (e.g., \code{chains}, \code{iter},
#'   \code{warmup}, \code{cores}).
#' @param margin_observed Optional logical vector (length \code{N}, one per
#'   row) indicating which observations have a usable margin. \code{NA} is
#'   treated as \code{FALSE}. Rows with \code{margin_observed == FALSE} are
#'   excluded from the margin/conduct likelihood term (and \code{margins} may
#'   be \code{NA}/invalid for those rows) but still fully participate in the
#'   share equation and the market share denominator -- shares are never
#'   renormalized to exclude masked rows. Default \code{NULL} treats every
#'   row as observed (current behavior).
#' @param ... Additional arguments (currently unused).
#'
#' @details
#' Two demand models are supported:
#' \describe{
#'   \item{\code{"logit"}}{Standard multinomial logit. Uses quantity shares
#'     and price levels. The absolute price coefficient alpha measures price
#'     sensitivity; its sign is determined by \code{market_side}.}
#'   \item{\code{"ces"}}{Constant Elasticity of Substitution. Uses revenue
#'     shares and log-prices in output markets. The elasticity of substitution
#'     gamma (= sigma) controls substitution patterns.}
#' }
#'
#' Five supply models are supported:
#' \describe{
#'   \item{1. \code{"bertrand"}}{Nash-Bertrand differentiated products pricing.}
#'   \item{2. \code{"monopolistic"}}{Monopolistic competition
#'     (single-product Bertrand).}
#'   \item{3. \code{"cournot"}}{Cournot quantity competition.}
#'   \item{4. \code{"auction"}}{Second-score auction.}
#'   \item{5. \code{"bargaining"}}{Nash bargaining with an estimated
#'     buyer-side bargaining weight.}
#' }
#'
#' The default prior on the standardized price coefficient is fixed at
#' \eqn{N(0,1)} on the log scale and does not use the estimation sample.
#'
#' The \code{priors} argument accepts a named list with any of:
#' \describe{
#'   \item{\code{prior_sigma_share}}{Lognormal center for log-share error SD
#'     (default 0.5).}
#'   \item{\code{prior_sigma_margin}}{Lognormal center for log inverse-margin
#'     error SD (default 0.75).}
#'   \item{\code{prior_sigma_meanval_strat}}{Half-normal scale for strategic
#'     group mean-value SD (default 1.0).}
#'   \item{\code{prior_alpha_mean}}{Prior mean for the log standardized
#'     price coefficient (default 0).}
#'   \item{\code{prior_alpha_sd}}{Prior SD for the log standardized price
#'     coefficient (default 1).}
#'   \item{\code{s0_prior_mean}}{Prior location for the outside share on the
#'     logit scale (default -1.4467864). Together with the default scale, this
#'     gives an arithmetic prior mean of 20 percent.}
#'   \item{\code{s0_prior_sd}}{Prior SD for the outside share on the logit
#'     scale (default 0.4562340). The default prior's 5th percentile is
#'     10 percent.}
#'   \item{\code{prior_log_sigma_share_sd}}{Prior SD for the log share-error
#'     scale (default 0.5).}
#'   \item{\code{prior_log_sigma_margin_sd}}{Prior SD for the log
#'     inverse-margin-error scale (default 1).}
#'   \item{\code{resid_df_fixed}}{Fixed degrees of freedom for the joint
#'     bivariate (or marginal) Student-t likelihood when \code{df} is
#'     numeric (default 5).}
#'   \item{\code{prior_nu_shape}, \code{prior_nu_rate}}{Gamma shape/rate for
#'     the prior on \code{nu - 2} when \code{df = "estimate"} (defaults 2 and
#'     0.25, i.e. mean 8, mode 4; provisional, pending prior-predictive
#'     checks).}
#'   \item{\code{prior_lkj}}{LKJ-style concentration parameter for the
#'     share-margin correlation \code{rho} (default 2).}
#'   \item{\code{bargaining_prior_lower}, \code{bargaining_prior_upper}}{Bounds
#'     for the bargaining-power prior (default 0 and 0.9999). The baseline is
#'     uniform on this interval; beta shape parameters provide a high-power
#'     sensitivity specification.}
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
#' # Estimate Bayesian logit demand + Bertrand supply
#' fit <- bayes_est(
#'   prices = pnb$prices,
#'   quantities = pnb$quantities,
#'   margins = pnb$margins,
#'   market = rep("PNB", nrow(pnb)),
#'   group = pnb$group,
#'   demand_model = "logit",
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
bayes_est <- function(prices, quantities, margins,
                      market, year = NULL, group,
                      demand_model = c("logit", "ces"),
                      market_side = NULL,
                      supply_model = c(
                        "bertrand", "monopolistic",
                        "cournot", "auction", "bargaining"
                      ),
                      use_cutoff = FALSE,
                      use_rho = TRUE,
                      df = 5,
                      use_market_resid_scale = FALSE,
                      use_alpha_market_pool = TRUE,
                      use_s0_market_year_pool = FALSE,
                      use_s0_market_pool = FALSE,
                      use_centered_alpha = FALSE,
                      use_centered_s0_pool = FALSE,
                      use_hmt = FALSE,
                      X_s0 = NULL,
                      prior_manifest = c("baseline", "tight", "diffuse"),
                      priors = list(),
                      likelihood = c("joint", "share", "margin"),
                      loo_cluster = NULL,
                      model_prior = NA_real_,
                      stan_options = list(),
                      margin_observed = NULL,
                      outside_price = 1,
                      ...) {
  cl <- match.call()
  demand_model <- match.arg(demand_model)
  prior_manifest <- match.arg(prior_manifest)
  likelihood <- match.arg(likelihood)
  if (is.null(market_side)) {
    market_side <- if (demand_model == "ces") "output" else "input"
  } else {
    market_side <- match.arg(market_side, c("input", "output"))
  }
  if (demand_model == "ces" && market_side != "output") {
    stop("CES estimation currently supports output markets only; use ",
         "demand_model = \"logit\" for input markets.")
  }
  if (!is.null(margin_observed)) {
    if (length(margin_observed) != length(prices)) {
      stop("'margin_observed' must have the same length as 'prices'.")
    }
    margin_observed <- as.logical(margin_observed) & !is.na(margin_observed)
  }
  supply_model <- match.arg(supply_model)
  supply_int <- switch(supply_model,
    bertrand      = 1L,
    monopolistic  = 2L,
    cournot       = 3L,
    auction       = 4L,
    bargaining    = 5L
  )

  # Default year if not supplied
  if (is.null(year)) {
    year <- rep(1L, length(prices))
  }
  if (is.null(loo_cluster)) loo_cluster <- market
  if (length(loo_cluster) != length(prices) || anyNA(loo_cluster)) {
    stop("'loo_cluster' must have one non-missing value per observation.")
  }
  if (length(model_prior) != 1L ||
      (!is.na(model_prior) && (!is.finite(model_prior) || model_prior < 0 ||
                               model_prior > 1))) {
    stop("'model_prior' must be NA or a probability in [0, 1].")
  }

  # Cutoff is meaningless under monopolistic competition (all firms are
  # single-product price-setters, so there is no strategic/fringe distinction)
  if (use_cutoff && supply_model == "monopolistic") {
    warning("use_cutoff ignored for monopolistic competition (no strategic/fringe distinction).")
    use_cutoff <- FALSE
  }
  # Validate
  validate_inputs(prices, quantities, margins, market, year, group,
                  margin_observed = margin_observed)
  observed_margins <- if (is.null(margin_observed)) margins else margins[margin_observed]
  if (market_side == "output" && any(observed_margins >= 1)) {
    stop("Output-market Lerner indices must be strictly less than 1. ",
         "Use market_side = \"input\" for positive markdown ratios.")
  }
  estimate_df <- identical(df, "estimate")
  resid_df_fixed <- NULL
  if (!estimate_df) {
    if (!is.numeric(df) || length(df) != 1L || !is.finite(df) || df <= 2) {
      stop("'df' must be a single finite numeric greater than 2, or the ",
           "character sentinel \"estimate\".")
    }
    resid_df_fixed <- df
  }
  flags <- list(
    use_cutoff = use_cutoff,
    use_rho = use_rho,
    use_market_resid_scale = use_market_resid_scale,
    use_alpha_market_pool = use_alpha_market_pool,
    use_s0_market_year_pool = use_s0_market_year_pool,
    use_s0_market_pool = use_s0_market_pool,
    use_centered_alpha = use_centered_alpha,
    use_centered_s0_pool = use_centered_s0_pool,
    use_hmt = use_hmt
  )
  bad_flags <- names(flags)[!vapply(
    flags,
    function(x) is.logical(x) && length(x) == 1L && !is.na(x),
    logical(1)
  )]
  if (length(bad_flags)) {
    stop("The following arguments must be non-missing logical scalars: ",
         paste(bad_flags, collapse = ", "), ".")
  }
  if (use_s0_market_pool && use_s0_market_year_pool) {
    stop("Use at most one of 'use_s0_market_pool' (one outside share per ",
         "market) and 'use_s0_market_year_pool' (one per market-year).")
  }
  priors <- resolve_prior_manifest(prior_manifest, priors)
  if (use_hmt && length(observed_margins) == 0L) {
    stop("use_hmt = TRUE requires at least one observed margin.")
  }

  # Build data frame (store raw inputs for later use in simulation)
  df <- data.frame(
    prices = prices,
    quantities = quantities,
    margins = margins,
    market = as.character(market),
    year = as.character(year),
    group = as.character(group),
    loo_cluster = as.character(loo_cluster),
    stringsAsFactors = FALSE
  )
  df$margin_observed <- if (is.null(margin_observed)) TRUE else margin_observed

  # Default Stan options. Defaults favor stable geometry; users can relax these
  # when a particular application is well behaved.
  default_stan <- list(
    chains  = 4L,
    iter    = 2000L,
    warmup  = 1000L,
    cores   = 1L,
    refresh = 200L,
    control = list(adapt_delta = 0.95, max_treedepth = 13L)
  )
  for (nm in names(stan_options)) {
    default_stan[[nm]] <- stan_options[[nm]]
  }

  # Convenience: allow adapt_delta and max_treedepth at top level
  control_params <- c("adapt_delta", "max_treedepth")
  found <- intersect(names(default_stan), control_params)
  if (length(found) > 0L) {
    if (is.null(default_stan$control)) default_stan$control <- list()
    for (p in found) {
      default_stan$control[[p]] <- default_stan[[p]]
      default_stan[[p]] <- NULL
    }
  }

  threads_per_chain <- default_stan$threads_per_chain
  if (is.null(threads_per_chain)) {
    threads_per_chain <- getOption("rstan.threads_per_chain", 1L)
  }
  target_chunks <- max(1L, as.integer(default_stan$chains) * as.integer(threads_per_chain))
  grainsize <- if (!is.null(default_stan$grainsize)) {
    as.integer(default_stan$grainsize)
  } else {
    max(1L, ceiling(length(prices) / target_chunks))
  }
  default_stan$threads_per_chain <- NULL
  default_stan$grainsize <- NULL

  # Prepare Stan data based on demand model
  if (demand_model == "logit") {
    stan_data <- prepare_stan_data(
      prices       = prices,
      quantities   = quantities,
      margins      = margins,
      market       = market,
      year         = year,
      group        = group,
      supply_model = supply_int,
      use_cutoff   = use_cutoff,
      market_side  = market_side,
      use_rho      = use_rho,
      use_market_resid_scale = use_market_resid_scale,
      use_alpha_market_pool = use_alpha_market_pool,
      use_s0_market_year_pool = use_s0_market_year_pool,
      use_s0_market_pool = use_s0_market_pool,
      use_centered_alpha = use_centered_alpha,
      use_centered_s0_pool = use_centered_s0_pool,
      use_hmt      = use_hmt,
      priors       = priors,
      likelihood_mode = match(likelihood, c("joint", "share", "margin")),
      loo_cluster = loo_cluster,
      X_s0         = X_s0,
      outside_price = outside_price,
      grainsize    = grainsize,
      margin_observed = margin_observed,
      estimate_df  = estimate_df,
      resid_df_fixed = resid_df_fixed
    )
    stan_mod <- stanmodels[["logit_supply"]]
    price_sd_val <- stan_data$price_sd
  } else {
    stan_data <- prepare_ces_stan_data(
      prices       = prices,
      quantities   = quantities,
      margins      = margins,
      market       = market,
      year         = year,
      group        = group,
      supply_model = supply_int,
      market_side  = market_side,
      use_cutoff   = use_cutoff,
      use_rho      = use_rho,
      use_hmt      = use_hmt,
      priors       = priors,
      likelihood_mode = match(likelihood, c("joint", "share", "margin")),
      loo_cluster = loo_cluster,
      X_s0         = X_s0,
      outside_price = outside_price,
      grainsize    = grainsize,
      use_market_resid_scale = use_market_resid_scale,
      use_alpha_market_pool = use_alpha_market_pool,
      use_s0_market_year_pool = use_s0_market_year_pool,
      use_centered_alpha = use_centered_alpha,
      use_centered_s0_pool = use_centered_s0_pool,
      margin_observed = margin_observed,
      estimate_df  = estimate_df,
      resid_df_fixed = resid_df_fixed
    )
    stan_mod <- stanmodels[["ces_supply"]]
    price_sd_val <- stan_data$log_price_sd
  }

  # Store the resolved market-year outside price on every observation so that
  # simulation and cross-demand comparison use the fitted normalization.
  df$outside_price <- as.numeric(stan_data$outside_price[
    stan_data$market_year_idx
  ])

  # Run sampling
  stanfit <- do.call(rstan::sampling, c(
    list(object = stan_mod, data = stan_data),
    default_stan
  ))

  metadata <- package_provenance(df, prior_manifest, model_prior)
  metadata$likelihood <- likelihood
  metadata$supply_model <- supply_model
  metadata$demand_model <- demand_model
  metadata$resolved_priors <- priors
  metadata$df_mode <- if (estimate_df) "estimated" else "fixed"
  metadata$resid_df_fixed <- if (estimate_df) NA_real_ else resid_df_fixed
  metadata$nu_prior <- list(
    shape = priors$prior_nu_shape,
    rate = priors$prior_nu_rate
  )

  new("BayesLogit",
    stanfit      = stanfit,
    data         = df,
    stan_data    = stan_data,
    supply_model = supply_int,
    demand_model = demand_model,
    market_side  = market_side,
    price_sd     = price_sd_val,
    call         = cl,
    metadata     = metadata
  )
}
