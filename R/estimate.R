#' Resolve the number of independent Stan chain processes.
#'
#' This policy deliberately uses the user's \code{mc.cores} option as the
#' default ceiling. It does not inspect host CPU topology or alter
#' within-chain threading. An explicit \code{stan_options$cores} value is
#' validated and retained unchanged.
#'
#' @keywords internal
.resolve_stan_cores <- function(stan_options = list(), chains = 1L) {
  validate_count <- function(value, label) {
    valid <- is.numeric(value) && length(value) == 1L &&
      !is.na(value) && is.finite(value) && value >= 1 &&
      value <= .Machine$integer.max && value == floor(value)
    if (!isTRUE(valid)) {
      stop("'", label, "' must be a positive finite integer.",
           call. = FALSE)
    }
    as.integer(value)
  }

  resolved_chains <- validate_count(chains, "stan_options$chains")
  explicit <- "cores" %in% names(stan_options) &&
    !is.null(stan_options$cores)
  if (explicit) return(validate_count(stan_options$cores,
                                      "stan_options$cores"))

  option_cores <- getOption("mc.cores", 1L)
  if (is.null(option_cores)) option_cores <- 1L
  option_cores <- validate_count(option_cores, "getOption('mc.cores')")
  min(resolved_chains, option_cores)
}

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
#' @param share_likelihood Share residual family: \code{"student_t"}
#'   (default) or \code{"gaussian"}.
#' @param quality_parameterization Quality effects are \code{"sampled"}
#'   (default) or \code{"collapsed"}; collapsed quality integrates Gaussian
#'   group effects and disallows market residual scales.
#' @param margin_likelihood Margin residual family, \code{"student_t"} or
#'   \code{"gaussian"}; \code{NULL} (default) follows
#'   \code{share_likelihood}. Gaussian-share Student-t margins require
#'   collapsed quality, fixed \code{df}, and \code{use_rho = FALSE}.
#' @param use_market_resid_scale Logical. If \code{TRUE}, allow market-level
#'   share and inverse-margin residual scales with hierarchical shrinkage.
#'   Default \code{FALSE}.
#' @param use_alpha_market_pool Logical. If \code{TRUE}, allow the standardized
#'   price effect to vary by market with hierarchical shrinkage. For Logit this
#'   gives market-level price coefficients; for CES it gives market-level
#'   structural elasticities \code{gamma}. Default \code{TRUE}.
#' @param use_s0_market_year_pool Logical legacy control. If \code{TRUE}, allow
#'   market-year outside shares with hierarchical shrinkage (one outside share
#'   per market-year). It is mutually exclusive with
#'   \code{use_s0_market_pool}; omit both legacy controls to use
#'   \code{s0_pooling = "auto"}. Note that market-year pooling introduces one
#'   weakly-identified outside share per market-year cell; when a market spans
#'   several periods with sparse per-cell information, its hyper-scale can be
#'   poorly identified (low BFMI). Prefer \code{use_s0_market_pool} in that
#'   regime.
#' @param use_s0_market_pool Logical legacy control. If \code{TRUE}, pool the outside quantity
#'   share at the \emph{market} level: one latent share per market, shared
#'   across all of that market's periods before any market-year \code{X_s0}
#'   adjustment, with hierarchical shrinkage across markets
#'   (structurally parallel to \code{use_alpha_market_pool}). Better identified
#'   than \code{use_s0_market_year_pool} when markets are observed over multiple
#'   periods, because each market-level effect is disciplined by all of the
#'   market's observations. It is mutually exclusive with
#'   \code{use_s0_market_year_pool}; omit both legacy controls to use
#'   \code{s0_pooling = "auto"}.
#' @param s0_pooling Character string selecting the outside-share pooling
#'   ladder: \code{"auto"} (default), \code{"none"}, \code{"market"}, or
#'   \code{"market_year"}. In automatic mode, a single market-year uses one
#'   direct outside share, while panels use one base outside share per market.
#'   The explicit legacy flags \code{use_s0_market_pool} and
#'   \code{use_s0_market_year_pool} remain supported when \code{s0_pooling} is
#'   omitted; supplying both interfaces together is ambiguous and errors.
#' @param market_s0_parameterization Character string selecting the legacy or
#'   market/outside-share rotated parameterization. The default is
#'   \code{"legacy"}; \code{"rotated"} is active only for multi-market fits
#'   and is a documented no-op for single-market fits.
#' @param use_centered_alpha Logical. If \code{TRUE}, use the centered
#'   parameterization for the selected market-level price-effect pool (Logit
#'   price coefficients or CES structural elasticities). Default \code{FALSE}.
#' @param use_centered_s0_pool Logical. If \code{TRUE}, use the centered
#'   parameterization for either outside-share pooling mode. Default \code{FALSE}.
#' @param use_hmt Logical. If \code{TRUE}, apply a soft Hypothetical
#'   Monopolist Test penalty. Default \code{FALSE}.
#' @param X_s0 Optional model matrix of market-year covariates for the outside
#'   share equation. Rows must match the market-year combinations generated by
#'   the input data after factor encoding.
#' @param outside_price Positive outside-good price used to put CES outside
#'   quantity and revenue shares on a common scale. Supply a scalar (default
#'   \code{1}) or one value per market-year. For Logit it does not enter the
#'   fitted likelihood, but is retained for posterior structural-model
#'   reconstruction and counterfactual routing.
#' @param structure Bayesian timing structure.  One of \code{"simultaneous"},
#'   \code{"fringe"}, or \code{"stackelberg"}; when omitted, the legacy
#'   \code{use_cutoff} flag selects simultaneous or fringe.
#' @param cutoff_scale Prespecified positive transition scale.  Legacy fringe
#'   defaults to \code{0.2}; Stackelberg defaults to \code{0.01}.
#' @param partition_prior Explicit fringe may use \code{"cutoff"} for the
#'   smooth legacy robustness mode; exact fringe and Stackelberg use
#'   \code{"ordered_k"}.
#' @param k_support Explicit vector of unique integer leader counts. Fringe
#'   supports may include zero; Stackelberg supports must be positive and
#'   feasible in every selected market.
#' @param k_prior Positional probabilities corresponding to \code{k_support}.
#'   Values are finite, non-negative, and must sum to one. \code{NULL} gives
#'   a uniform prior.
#' @param k_scope Legacy cutoff scope selector. Exact ordered-k uses one
#'   partition and one firm ranking for the complete market panel.
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
#'   \code{warmup}, \code{cores}). If \code{cores} is omitted or
#'   \code{NULL}, the default is \code{min(chains, getOption("mc.cores", 1L))};
#'   an explicit \code{cores} value is validated and used unchanged.
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
#' The default is the historical Student-t share and margin likelihood with
#' sampled quality. Gaussian sampled and Gaussian collapsed quality are
#' selected with \code{share_likelihood} and
#' \code{quality_parameterization}. Collapsed fits retain normalized bridge
#' evidence; their observation-level predictive arrays are unavailable for
#' LOO/WAIC.
#'
#' Two demand models are supported:
#' \describe{
#'   \item{\code{"logit"}}{Standard multinomial logit. Uses quantity shares
#'     and price levels. The absolute price coefficient alpha measures price
#'     sensitivity; its sign is determined by \code{market_side}.}
#'   \item{\code{"ces"}}{Constant Elasticity of Substitution. Its demand
#'     likelihood uses conditional quantity shares and log-prices in output
#'     markets; revenue shares enter the CES supply first-order conditions.
#'     The elasticity of substitution gamma (= sigma) controls substitution
#'     patterns.}
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
#'   \item{\code{prior_sigma_meanval_strat}}{Scale parameter for the
#'     positive-truncated Normal strategic group mean-value SD, centered at
#'     1 (default 1.0).}
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
#' merger <- antitrust::counterfactual(ownership = owner_post)
#' result <- simulate(fit, merger, n_draws = 200)
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
                      share_likelihood = c("student_t", "gaussian"),
                      quality_parameterization = c("sampled", "collapsed"),
                      margin_likelihood = NULL,
                      use_market_resid_scale = FALSE,
                      use_alpha_market_pool = TRUE,
                      use_s0_market_year_pool = FALSE,
                      use_s0_market_pool = FALSE,
                      s0_pooling = c("auto", "none", "market", "market_year"),
                      market_s0_parameterization = c("legacy", "rotated"),
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
                      structure = NULL,
                      cutoff_scale = NULL,
                      partition_prior = "cutoff",
                      k_support = NULL,
                      k_prior = NULL,
                      k_scope = "market",
                      ...) {
  cl <- match.call()
  s0_pooling_supplied <- !missing(s0_pooling)
  legacy_s0_supplied <- !missing(use_s0_market_pool) ||
    !missing(use_s0_market_year_pool)
  use_cutoff_supplied <- !missing(use_cutoff)
  resolved_structure <- .stack_resolve_structure(
    structure = structure, use_cutoff = use_cutoff,
    use_cutoff_supplied = use_cutoff_supplied
  )
  structure <- resolved_structure$name
  use_cutoff <- resolved_structure$cutoff
  cutoff_scale_supplied <- !missing(cutoff_scale) && !is.null(cutoff_scale)
  if (identical(structure, "simultaneous") && cutoff_scale_supplied) {
    stop("'cutoff_scale' is only defined for fringe and Stackelberg structures.")
  }
  demand_model <- match.arg(demand_model)
  prior_manifest <- match.arg(prior_manifest)
  likelihood <- match.arg(likelihood)
  share_likelihood <- match.arg(share_likelihood)
  quality_parameterization <- match.arg(quality_parameterization)
  if (!is.null(margin_likelihood)) {
    margin_likelihood <- match.arg(margin_likelihood, c("student_t", "gaussian"))
  }
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
  if (structure == "stackelberg" && !(supply_model %in% c("bertrand", "cournot"))) {
    stop("Bayesian Stackelberg supports only Bertrand and Cournot supply models.")
  }
  partition_prior_supplied <- !missing(partition_prior)
  exact_structure <- resolved_structure$supplied &&
    structure %in% c("fringe", "stackelberg")
  if (exact_structure) {
    if (!partition_prior_supplied) partition_prior <- "ordered_k"
    partition_prior <- match.arg(partition_prior, c("cutoff", "ordered_k"))
    if (structure == "stackelberg" && partition_prior != "ordered_k") {
      stop("Explicit Stackelberg requires partition_prior = 'ordered_k'.")
    }
    if (partition_prior == "ordered_k") {
      if (is.null(k_support)) {
        stop("Explicit ordered-k structures require 'k_support'.")
      }
      if (cutoff_scale_supplied) {
        stop("'cutoff_scale' is not defined for exact ordered-k structures.")
      }
      if (!isTRUE(k_scope == "market")) {
        stop("Exact ordered-k uses one partition per market; 'k_scope' must be 'market'.")
      }
      if (is.list(k_prior)) {
        stop("Exact ordered-k requires a positional numeric 'k_prior' vector.")
      }
    } else if (!is.null(k_support) || !is.null(k_prior)) {
      stop("'k_support' and 'k_prior' require partition_prior = 'ordered_k'.")
    }
  } else {
    ## No explicit structure means the historical simultaneous or smooth
    ## fringe path. Its k arguments are deliberately ignored.
    partition_prior <- "cutoff"
    k_support <- NULL
    k_prior <- NULL
  }
  k_scope <- match.arg(k_scope, c("market", "period"))
  .stack_validate_public_k(k_support, k_prior, structure)
  cutoff_scale <- if (exact_structure && partition_prior == "ordered_k") {
    NULL
  } else {
    .stack_resolve_scale(structure, cutoff_scale)
  }
  if (exact_structure && structure == "fringe" && partition_prior == "ordered_k" &&
      !(supply_model %in% c("bertrand", "cournot"))) {
    stop("Exact ordered-k fringe supports only Bertrand and Cournot supply models.")
  }
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
  s0_resolution <- .resolve_s0_pooling(
    s0_pooling = if (s0_pooling_supplied) s0_pooling else "auto",
    s0_pooling_supplied = s0_pooling_supplied,
    use_s0_market_pool = use_s0_market_pool,
    use_s0_market_year_pool = use_s0_market_year_pool,
    legacy_supplied = legacy_s0_supplied,
    N_market = nlevels(droplevels(as.factor(market))),
    N_market_year = {
      market_resolve <- droplevels(as.factor(market))
      year_resolve <- droplevels(as.factor(year))
      length(unique(paste(as.integer(market_resolve),
                          as.integer(year_resolve), sep = "_")))
    }
  )
  s0_pooling_requested <- s0_resolution$requested
  s0_pooling <- s0_resolution$resolved
  use_s0_market_pool <- s0_resolution$use_s0_market_pool
  use_s0_market_year_pool <- s0_resolution$use_s0_market_year_pool
  market_s0_resolution <- .resolve_market_s0_parameterization(
    market_s0_parameterization,
    nlevels(droplevels(as.factor(market)))
  )
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
  likelihood_modes <- .resolve_likelihood_modes(
    share_likelihood = share_likelihood,
    quality_parameterization = quality_parameterization,
    margin_likelihood = margin_likelihood,
    use_rho = use_rho,
    use_market_resid_scale = use_market_resid_scale,
    estimate_df = estimate_df,
    likelihood_mode = match(likelihood, c("joint", "share", "margin"))
  )
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
  ## Keep the chain count needed for the default core policy even when a
  ## caller explicitly supplies chains = NULL.  A non-NULL cores entry is
  ## validated and takes precedence over this policy.
  if (is.null(default_stan$chains)) default_stan$chains <- 4L
  default_stan$cores <- .resolve_stan_cores(
    stan_options = stan_options,
    chains = default_stan$chains
  )

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
      market_s0_parameterization = market_s0_resolution$requested,
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
      resid_df_fixed = resid_df_fixed,
      structure = structure,
      cutoff_scale = if (identical(structure, "simultaneous") ||
                         (identical(structure, "stackelberg") &&
                          identical(partition_prior, "ordered_k"))) NULL else cutoff_scale,
      partition_prior = partition_prior,
      k_support = k_support,
      k_prior = k_prior,
      k_scope = k_scope,
      share_likelihood = likelihood_modes$share_likelihood,
      quality_parameterization = likelihood_modes$quality_parameterization,
      margin_likelihood = likelihood_modes$margin_likelihood
    )
    stan_mod <- stanmodels[[if (structure == "stackelberg" ||
                                (exact_structure && structure == "fringe" &&
                                 partition_prior == "ordered_k")) {
      "logit_stackelberg"
    } else "logit_supply"]]
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
      use_s0_market_pool = use_s0_market_pool,
      market_s0_parameterization = market_s0_resolution$requested,
      use_centered_alpha = use_centered_alpha,
      use_centered_s0_pool = use_centered_s0_pool,
      margin_observed = margin_observed,
      estimate_df  = estimate_df,
      resid_df_fixed = resid_df_fixed,
      structure = structure,
      cutoff_scale = if (identical(structure, "simultaneous") ||
                         (identical(structure, "stackelberg") &&
                          identical(partition_prior, "ordered_k"))) NULL else cutoff_scale,
      partition_prior = partition_prior,
      k_support = k_support,
      k_prior = k_prior,
      k_scope = k_scope,
      share_likelihood = likelihood_modes$share_likelihood,
      quality_parameterization = likelihood_modes$quality_parameterization,
      margin_likelihood = likelihood_modes$margin_likelihood
    )
    stan_mod <- stanmodels[[if (structure == "stackelberg" ||
                                (exact_structure && structure == "fringe" &&
                                 partition_prior == "ordered_k")) {
      "ces_stackelberg"
    } else "ces_supply"]]
    price_sd_val <- stan_data$log_price_sd
  }

  attr(stan_data, "s0_pooling") <- s0_pooling
  attr(stan_data, "s0_pooling_requested") <- s0_pooling_requested
  attr(stan_data, "s0_pooling_resolved") <- s0_pooling
  attr(stan_data, "s0_pooling_source") <- s0_resolution$source
  attr(stan_data, "market_s0_parameterization") <-
    market_s0_resolution$requested
  attr(stan_data, "market_s0_parameterization_requested") <-
    market_s0_resolution$requested
  attr(stan_data, "market_s0_parameterization_effective") <-
    market_s0_resolution$effective
  attr(stan_data, "market_s0_rotation_active") <- market_s0_resolution$active
  attr(stan_data, "market_s0_coordinate_convention") <-
    market_s0_resolution$coordinate_convention

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
  metadata$stan_chains <- as.integer(default_stan$chains)
  metadata$stan_cores <- as.integer(default_stan$cores)
  metadata$likelihood <- likelihood
  metadata$share_likelihood <- likelihood_modes$share_likelihood
  metadata$margin_likelihood <- likelihood_modes$margin_likelihood
  metadata$quality_parameterization <- likelihood_modes$quality_parameterization
  metadata$likelihood_mode <- likelihood_modes$mode
  metadata$bayes_is_collapsed_quality <-
    identical(likelihood_modes$quality_parameterization, "collapsed")
  metadata$supply_model <- supply_model
  metadata$demand_model <- demand_model
  metadata$s0_pooling <- s0_pooling
  metadata$s0_pooling_requested <- s0_pooling_requested
  metadata$s0_pooling_resolved <- s0_pooling
  metadata$s0_pooling_source <- s0_resolution$source
  metadata$market_s0_parameterization <- market_s0_resolution$requested
  metadata$market_s0_parameterization_requested <- market_s0_resolution$requested
  metadata$market_s0_parameterization_effective <- market_s0_resolution$effective
  metadata$market_s0_rotation_active <- market_s0_resolution$active
  metadata$market_s0_coordinate_convention <-
    market_s0_resolution$coordinate_convention
  metadata$structure <- structure
  metadata$structure_code <- as.integer(.stack_structure_code[[structure]])
  metadata$cutoff_scale <- cutoff_scale
  metadata$partition_prior <- partition_prior
  metadata$k_support <- k_support
  metadata$k_scope <- k_scope
  metadata$k_prior <- k_prior
  partition_metadata <- attr(stan_data, "partition_metadata", exact = TRUE)
  metadata$k_support_resolved <- if (is.null(partition_metadata)) NULL else
    partition_metadata$partition_supports
  metadata$k_prior_resolved <- if (is.null(partition_metadata)) NULL else
    partition_metadata$k_prior_by_unit
  metadata$firm_market_ranking <- attr(stan_data, "firm_market_ranking", exact = TRUE)
  metadata$resolved_priors <- priors
  metadata$df_mode <- if (estimate_df) "estimated" else "fixed"
  metadata$resid_df_fixed <- if (estimate_df) NA_real_ else resid_df_fixed
  metadata$nu_prior <- list(
    shape = priors$prior_nu_shape,
    rate = priors$prior_nu_rate
  )
  metadata$prior_normalization <- list(
    mode = "native_stan_normalized",
    model_spec_version = 9L,
    share_likelihood = likelihood_modes$share_likelihood,
    margin_likelihood = likelihood_modes$margin_likelihood,
    quality_parameterization = likelihood_modes$quality_parameterization,
    market_s0_parameterization = market_s0_resolution$requested,
    market_s0_parameterization_requested = market_s0_resolution$requested,
    market_s0_parameterization_effective = market_s0_resolution$effective,
    market_s0_rotation_active = market_s0_resolution$active,
    market_s0_coordinate_convention = market_s0_resolution$coordinate_convention,
    s0_pooling = s0_pooling,
    s0_pooling_requested = s0_pooling_requested,
    s0_pooling_resolved = s0_pooling,
    s0_pooling_source = s0_resolution$source,
    partition_support = k_support,
    k_support = k_support,
    k_prior = k_prior,
    k_support_resolved = metadata$k_support_resolved,
    k_prior_resolved = metadata$k_prior_resolved,
    firm_market_ranking = metadata$firm_market_ranking
  )

  spec <- .bayes_model_spec(list(
    demand_model = demand_model,
    supply_model = supply_int
  ))
  observed <- list(
    prices = as.numeric(df$prices),
    quantities = as.numeric(df$quantities),
    margins = as.numeric(df$margins),
    market = as.character(df$market),
    year = as.character(df$year),
    group = as.character(df$group),
    margin_observed = as.logical(df$margin_observed),
    outside_price = as.numeric(df$outside_price)
  )
  parameters <- list(
    source = "stanfit",
    demand_parameter = if (identical(demand_model, "ces")) "gamma" else "alpha",
    supply_model = supply_model
  )

  new("BayesLogit",
    spec         = spec,
    model        = NULL,
    parameters   = parameters,
    observed     = observed,
    diagnostics  = list(
      status = "completed",
      source = "estimated",
      route = "bayes_est",
      model_spec_version = attr(stan_data, "model_spec_version", exact = TRUE),
      provenance = metadata
    ),
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
