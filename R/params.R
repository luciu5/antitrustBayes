# Parameter-driven path to a simulate()-ready antitrust fit.
#
# This file extracts the market/parameter assembly that used to live only
# inside bayes_sim()'s Stan-draw loop (simulate.R) into Stan-free building
# blocks. bayes_params()/bayes_market() let a caller hand-specify demand and
# conduct parameters -- a single point estimate or a matrix of draws -- and
# build_conduct_object()/as_antitrust_fit() turn those into the same internal
# BayesConductLogit/BayesConductCES realization bayes_sim() has always built,
# then route exact no-cutoff variants to antitrust's registered concrete model
# classes. The result is an antitrust::AntitrustFit, so antitrust's
# counterfactual()/simulate() machinery can run without ever touching Stan.

#' Bayesian demand and conduct parameters for one simulation unit
#'
#' A validated container for the demand and conduct parameters that
#' \code{\link{bayes_sim}} otherwise reads out of posterior draws. Supply
#' either a single point estimate (one simulation) or a matrix of draws
#' (one row per draw, preserving the posterior-distribution workflow).
#'
#' @param alpha Logit price coefficient. A single numeric, or a numeric
#'   vector with one value per draw. Exactly one of \code{alpha}/\code{gamma}
#'   must be supplied.
#' @param gamma CES elasticity parameter (> 1). A single numeric, or a
#'   numeric vector with one value per draw.
#' @param s0 Outside-good probability/share in \code{(0, 1)}. A single
#'   numeric, or a numeric vector with one value per draw.
#' @param bargaining_power Numeric in \code{[0, 1)}; only used when
#'   \code{supply_model = 5}. Default \code{0}.
#' @param cutoff_share Posterior fringe-cutoff share. Default \code{0}.
#' @param cutoff_scale Fringe-cutoff logistic scale. Default \code{0.2}.
#' @param use_cutoff Whether the posterior fringe cutoff applies. Default
#'   \code{FALSE}; set \code{TRUE} together with a nonzero \code{cutoff_share}
#'   to activate it.
#' @param supply_model Integer 1-5 identifying the conduct model (see
#'   \code{\link{bayes_sim}}).
#' @param structure Supply structure: \code{"simultaneous"},
#'   \code{"fringe"}, or \code{"stackelberg"}.
#' @param partition_prior For explicit fringe, \code{"cutoff"} retains the
#'   smooth legacy robustness mode; \code{"ordered_k"} selects the exact
#'   finite support. Explicit Stackelberg uses \code{"ordered_k"}.
#' @param k_support Ordered integer counts. Fringe supports may include zero;
#'   Stackelberg supports must be positive. The vector must be feasible in
#'   every selected market.
#' @param k_prior Positional probabilities corresponding to \code{k_support};
#'   \code{NULL} gives a uniform prior. Values must be finite, non-negative,
#'   and sum to one.
#' @param k_scope Retained for legacy cutoff preparation; exact ordered-k is
#'   always common across the complete market panel.
#'
#' @return A list of class \code{bayes_params} with one element per draw
#'   (length 1 for a point estimate).
#' @export
bayes_params <- function(alpha = NULL, gamma = NULL, s0,
                          bargaining_power = 0, cutoff_share = 0,
                          cutoff_scale = 0.2, use_cutoff = FALSE,
                          supply_model, structure = NULL,
                          partition_prior = "cutoff", k_support = NULL,
                          k_prior = NULL,
                          k_scope = "market") {
  if (is.null(alpha) == is.null(gamma)) {
    stop("Supply exactly one of 'alpha' (logit) or 'gamma' (ces).")
  }
  is_ces <- !is.null(gamma)
  slope <- if (is_ces) gamma else alpha
  if (!is.numeric(slope) || length(slope) < 1L || anyNA(slope)) {
    stop("'", if (is_ces) "gamma" else "alpha", "' must be a non-missing numeric vector.")
  }
  if (is_ces && any(slope <= 1)) {
    stop("'gamma' must exceed one for every draw.")
  }
  if (!is_ces && any(slope == 0)) {
    stop("'alpha' must be nonzero for every draw.")
  }
  if (missing(supply_model) || length(supply_model) != 1L ||
      !supply_model %in% 1:5) {
    stop("'supply_model' must be a single integer in 1:5.")
  }

  resolved_structure <- .stack_resolve_structure(
    structure = structure, use_cutoff = use_cutoff,
    use_cutoff_supplied = !missing(use_cutoff)
  )
  k_scope <- match.arg(k_scope, c("market", "period"))
  if (resolved_structure$name == "stackelberg" && !(supply_model %in% c(1L, 3L))) {
    stop("Bayesian Stackelberg supports only Bertrand (1) and Cournot (3) supply models.")
  }
  cutoff_scale_supplied <- !missing(cutoff_scale) && !is.null(cutoff_scale)
  partition_prior_supplied <- !missing(partition_prior)
  exact_structure <- resolved_structure$supplied &&
    resolved_structure$name %in% c("fringe", "stackelberg")
  if (exact_structure) {
    if (!partition_prior_supplied) partition_prior <- "ordered_k"
    partition_prior <- match.arg(partition_prior, c("cutoff", "ordered_k"))
    if (identical(resolved_structure$name, "stackelberg") &&
        partition_prior != "ordered_k") {
      stop("Explicit Stackelberg requires partition_prior = 'ordered_k'.")
    }
    if (identical(partition_prior, "cutoff") &&
        (!is.null(k_support) || !is.null(k_prior))) {
      stop("'k_support' and 'k_prior' require partition_prior = 'ordered_k'.")
    }
    if (identical(partition_prior, "ordered_k")) {
      if (is.null(k_support)) {
        stop("Explicit ordered-k structures require 'k_support'.")
      }
      if (cutoff_scale_supplied) {
        stop("'cutoff_scale' is not defined for exact ordered-k structures.")
      }
      if (!identical(k_scope, "market")) {
        stop("Exact ordered-k uses one partition per market; 'k_scope' must be 'market'.")
      }
      if (is.list(k_prior)) {
        stop("Exact ordered-k requires a positional numeric 'k_prior' vector.")
      }
    }
  } else {
    ## Calls without an explicit structure retain the historical simultaneous
    ## or smooth-fringe interpretation, including ignoring k arguments.
    partition_prior <- "cutoff"
    k_scope <- "market"
    k_support <- NULL
    k_prior <- NULL
  }
  .stack_validate_public_k(k_support, k_prior, resolved_structure$name)
  if (exact_structure && resolved_structure$name == "fringe" &&
      identical(partition_prior, "ordered_k") &&
      !(supply_model %in% c(1L, 3L))) {
    stop("Exact ordered-k fringe supports only Bertrand and Cournot supply models.")
  }
  n <- length(slope)
  recycle <- function(x, name) {
    if (length(x) == 1L) return(rep(x, n))
    if (length(x) != n) {
      stop("'", name, "' must have length 1 or the same length as '",
           if (is_ces) "gamma" else "alpha", "'.")
    }
    x
  }
  s0 <- recycle(s0, "s0")
  if (any(!is.finite(s0)) || any(s0 <= 0) || any(s0 >= 1)) {
    stop("'s0' must lie strictly between 0 and 1 for every draw.")
  }
  bargaining_power <- recycle(bargaining_power, "bargaining_power")
  if (supply_model == 5L &&
      (any(!is.finite(bargaining_power)) || any(bargaining_power < 0) ||
       any(bargaining_power >= 1))) {
    stop("'bargaining_power' must lie in [0, 1) for every draw.")
  }
  cutoff_share <- recycle(cutoff_share, "cutoff_share")
  if (is.null(cutoff_scale) || missing(cutoff_scale)) {
    cutoff_scale <- if (resolved_structure$name == "stackelberg") 0.01 else 0.2
  }
  cutoff_scale <- recycle(cutoff_scale, "cutoff_scale")
  if (any(!is.finite(cutoff_scale)) || any(cutoff_scale <= 0)) {
    stop("'cutoff_scale' must be positive and finite for every draw.")
  }

  .bayes_params_raw(
    alpha = if (is_ces) NULL else slope,
    gamma = if (is_ces) slope else NULL,
    s0 = s0,
    bargaining_power = bargaining_power,
    cutoff_share = cutoff_share,
    cutoff_scale = cutoff_scale,
    use_cutoff = isTRUE(use_cutoff),
    supply_model = as.integer(supply_model),
    structure = resolved_structure$name,
    structure_code = as.integer(resolved_structure$code),
    partition_prior = partition_prior,
    k_support = k_support,
    k_prior = k_prior,
    k_scope = k_scope
  )
}

## Build a bayes_params container without validation, for internal reuse in
## hot per-draw loops (bayes_sim()'s market x draw loop) that already know
## their inputs are well-formed and rely on tryCatch, not upfront stop(), to
## isolate a single bad draw.
.bayes_params_raw <- function(alpha = NULL, gamma = NULL, s0,
                               bargaining_power, cutoff_share,
                               cutoff_scale, use_cutoff, supply_model,
                               structure = NULL, structure_code = NULL,
                               partition_prior = "cutoff", k_support = NULL,
                               k_prior = NULL,
                               k_scope = "market") {
  is_ces <- !is.null(gamma)
  resolved_structure <- .stack_resolve_structure(
    structure = structure, use_cutoff = use_cutoff,
    use_cutoff_supplied = FALSE
  )
  structure(
    list(
      alpha = alpha,
      gamma = gamma,
      s0 = s0,
      bargaining_power = bargaining_power,
      cutoff_share = cutoff_share,
      cutoff_scale = cutoff_scale,
      use_cutoff = use_cutoff,
      supply_model = as.integer(supply_model),
      structure = resolved_structure$name,
      structure_code = as.integer(if (is.null(structure_code)) resolved_structure$code else structure_code),
      partition_prior = partition_prior,
      k_support = k_support,
      k_prior = k_prior,
      k_scope = k_scope,
      demand_model = if (is_ces) "ces" else "logit",
      n_draws = length(if (is_ces) gamma else alpha)
    ),
    class = "bayes_params"
  )
}

#' Market geometry for a parameter-driven Bayesian simulation
#'
#' The observed market data that \code{\link{bayes_sim}} otherwise reads out
#' of \code{fit@data}/\code{fit@stan_data}/\code{fit@price_sd}.
#'
#' @param prices,quantities,margins Numeric vectors, one value per product.
#' @param group Character or factor firm identifier, one per product.
#' @param outside_price Positive scalar price of the outside good.
#' @param market_side \code{"input"} or \code{"output"}.
#' @param insideSize Total inside quantity (logit) or a CES scale factor.
#'   Default \code{1}.
#' @param labels Character vector of product labels. Defaults to
#'   \code{as.character(group)}.
#'
#' @return A list of class \code{bayes_market}.
#' @export
bayes_market <- function(prices, quantities, margins, group,
                          outside_price, market_side = c("output", "input"),
                          insideSize = 1, labels = NULL) {
  market_side <- match.arg(market_side)
  n_prod <- length(prices)
  if (length(quantities) != n_prod || length(group) != n_prod) {
    stop("'prices', 'quantities', and 'group' must have the same length.")
  }
  if (missing(margins) || is.null(margins)) {
    margins <- rep(NA_real_, n_prod)
  }
  margins <- as.numeric(margins)
  margin_observed <- is.finite(margins)
  margins[!is.finite(margins)] <- 0.5
  if (length(margins) != n_prod) {
    stop("'margins' must have one value per product.")
  }
  if (!is.numeric(outside_price) || length(outside_price) != 1L ||
      !is.finite(outside_price) || outside_price <= 0) {
    stop("'outside_price' must be a single positive finite number.")
  }
  if (!is.numeric(insideSize) || length(insideSize) != 1L ||
      !is.finite(insideSize) || insideSize <= 0) {
    stop("'insideSize' must be a single positive finite number.")
  }
  if (is.null(labels)) labels <- as.character(group)
  if (length(labels) != n_prod) {
    stop("'labels' must have one value per product.")
  }

  structure(
    list(
      prices = as.numeric(prices),
      quantities = as.numeric(quantities),
      margins = margins,
      margin_observed = margin_observed,
      group = as.character(group),
      outside_price = outside_price,
      market_side = market_side,
      insideSize = insideSize,
      labels = labels,
      n_prod = n_prod
    ),
    class = "bayes_market"
  )
}
#' Specify a Bayesian structural fit from parameter draws
#'
#' Combines directly supplied demand/conduct parameters and an observed market
#' into a \code{\linkS4class{BayesFit}}. The resulting object uses the same
#' \code{simulate()} S4 method as a fit returned by \code{\link{bayes_est}},
#' without running Stan.
#' No CmdStan installation or C++ toolchain is needed for this supplied-state
#' route. Its calibrated parameters are distinct from an estimated quality posterior.
#'
#' @param params A \code{\link{bayes_params}} object containing one or more
#'   parameter draws.
#' @param market A \code{\link{bayes_market}} object describing one market.
#'
#' @return A \code{BayesSpecifiedFit}, inheriting from
#'   \code{\linkS4class{BayesFit}} and \code{antitrust::StructuralFit}.
#' @export
#'
#' @examples
#' data(pnb)
#' params <- bayes_params(alpha = 1, s0 = 0.015, supply_model = 1L)
#' market <- bayes_market(
#'   prices = pnb$prices, quantities = pnb$quantities, margins = pnb$margins,
#'   group = pnb$group, outside_price = 1, market_side = "input",
#'   insideSize = attr(pnb, "insideSize"), labels = pnb$bank
#' )
#' fit <- bayes_specify(params, market)
#' owner_post <- pnb$group
#' owner_post[owner_post == attr(pnb, "target")] <- attr(pnb, "acquirer")
#' result <- simulate(fit, antitrust::counterfactual(ownership = owner_post))
bayes_specify <- function(params, market) {
  if (!inherits(params, "bayes_params")) {
    stop("'params' must be a bayes_params object; see bayes_params().")
  }
  if (!inherits(market, "bayes_market")) {
    stop("'market' must be a bayes_market object; see bayes_market().")
  }

  cl <- match.call()
  spec <- .bayes_model_spec(params)
  specified_metadata <- .stack_specified_metadata(params, market)
  resolved_partition <- specified_metadata$partition
  market_year_size <- .market_year_quantity_size(
    quantities = market$quantities,
    observation_to_market_year = rep.int(1L, market$n_prod),
    market_year_labels = "unit_1_1"
  )
  data <- data.frame(
    prices = market$prices,
    quantities = market$quantities,
    margins = market$margins,
    market = rep("unit_1", market$n_prod),
    year = rep("1", market$n_prod),
    group = market$group,
    margin_observed = market$margin_observed,
    outside_price = rep(market$outside_price, market$n_prod),
    market_year_idx = rep.int(1L, market$n_prod),
    market_year = rep("unit_1_1", market$n_prod),
    N_mt = rep(market_year_size$N_mt, market$n_prod),
    stringsAsFactors = FALSE
  )
  observed <- list(
    prices = market$prices,
    quantities = market$quantities,
    margins = market$margins,
    ownerPre = market$group,
    labels = market$labels,
    priceOutside = market$outside_price,
    insideSize = market$insideSize,
    market_year_idx = market_year_size$observation_to_market_year,
    N_mt = rep(market_year_size$N_mt,
               length(market_year_size$observation_to_market_year))
  )
  metadata <- list(
    source = "bayes_specify",
    demand_model = params$demand_model,
    supply_model = params$supply_model,
    posterior_draws = params$n_draws,
    structure = params$structure %||% if (isTRUE(params$use_cutoff)) "fringe" else "simultaneous",
    structure_code = params$structure_code %||% if (isTRUE(params$use_cutoff)) 1L else 0L,
    cutoff_scale = params$cutoff_scale,
    partition_prior = params$partition_prior %||% "cutoff",
    k_support = .stack_whole_or(params$k_support, NULL),
    k_prior = .stack_whole_or(params$k_prior, NULL),
    k_scope = params$k_scope %||% "market",
    quality_scope = "global",
    quality_scope_code = 1L,
    prior_normalization = list(
      mode = "supplied_draws",
      model_spec_version = .antitrust_bayes_model_spec_version,
      quality_scope = "global",
      quality_scope_code = 1L,
      partition_support = if (is.null(resolved_partition)) NULL else
        resolved_partition$partition_supports,
      k_support = .stack_whole_or(params$k_support, NULL),
      k_prior = .stack_whole_or(params$k_prior, NULL),
      k_support_resolved = if (is.null(resolved_partition)) NULL else
        resolved_partition$partition_supports,
      k_prior_resolved = if (is.null(resolved_partition)) NULL else
        resolved_partition$k_prior_by_unit,
      firm_market_ranking = specified_metadata$ranking,
      source = "bayes_specify"
    ),
    k_support_resolved = if (is.null(resolved_partition)) NULL else
      resolved_partition$partition_supports,
    k_prior_resolved = if (is.null(resolved_partition)) NULL else
      resolved_partition$k_prior_by_unit,
    firm_market_ranking = specified_metadata$ranking,
    market_year_size = market_year_size,
    N_mt = market_year_size$N_mt,
    market_year_idx = market_year_size$observation_to_market_year,
    market_year_labels = market_year_size$market_year_labels,
    market_year_size_units = market_year_size$units,
    market_year_size_source = market_year_size$source
  )

  methods::new(
    "BayesSpecifiedFit",
    spec = spec,
    model = NULL,
    parameters = params,
    observed = observed,
    diagnostics = list(
      status = "completed",
      source = "specified",
      route = "bayes_specify",
      posterior_draws = params$n_draws
    ),
    data = data,
    supply_model = as.integer(params$supply_model),
    demand_model = params$demand_model,
    market_side = market$market_side,
    call = cl,
    metadata = metadata,
    params = params,
    market = market
  )
}

## Invert demand for one draw and build the antitrust-compatible conduct
## object. This is the same math bayes_sim() has always run
## (simulate.R:344-407); moving it here lets both the calibration-based and
## the manual-parameter paths share one implementation.
build_conduct_object <- function(params, draw_idx, market, ownerPre, ownerPost,
                                  mcDelta = NULL, subset = NULL) {
  n_prod <- market$n_prod
  if (is.null(mcDelta)) mcDelta <- rep(0, n_prod)
  if (is.null(subset)) subset <- rep(TRUE, n_prod)

  is_ces <- identical(params$demand_model, "ces")
  is_output <- identical(market$market_side, "output")

  prices_m <- market$prices
  quantities_m <- market$quantities
  outside_price_m <- market$outside_price
  shares_cond <- quantities_m / sum(quantities_m)

  s0 <- params$s0[draw_idx]
  shares_uncond <- shares_cond * (1 - s0)

  if (is_ces) {
    gamma <- params$gamma[draw_idx]
    meanval <- (shares_uncond / s0) * (prices_m / outside_price_m)^gamma
    total_expenditure_index <- s0 * outside_price_m +
      sum(prices_m * shares_uncond)
    R0 <- s0 * outside_price_m / total_expenditure_index
    shares_for_object <- prices_m * shares_uncond / total_expenditure_index
    share_inside_for_object <- 1 - R0
    inside_size_for_object <- market$insideSize * sum(prices_m * shares_cond)
    class_name <- "BayesConductCES"
    slopes <- list(gamma = gamma, meanval = meanval, alpha = R0 / (1 - R0))
  } else {
    alpha <- params$alpha[draw_idx]
    if (is_output) alpha <- -alpha
    meanval <- log(shares_uncond / s0) - alpha * (prices_m - outside_price_m)
    shares_for_object <- shares_uncond
    share_inside_for_object <- 1 - s0
    inside_size_for_object <- market$insideSize
    class_name <- "BayesConductLogit"
    slopes <- list(alpha = alpha, meanval = meanval)
  }

  bargaining_draw <- if (params$supply_model == 5L) {
    params$bargaining_power[draw_idx]
  } else {
    0
  }
  use_cutoff <- isTRUE(params$use_cutoff)
  cutoff_draw <- params$cutoff_share[draw_idx]

  obj <- withCallingHandlers(methods::new(
    class_name,
    prices = prices_m,
    shares = shares_for_object,
    margins = market$margins,
    weights = rep(1, n_prod),
    ownerPre = ownerPre,
    ownerPost = ownerPost,
    mcDelta = mcDelta,
    subset = subset,
    insideSize = inside_size_for_object,
    priceOutside = outside_price_m,
    priceStart = prices_m,
    pricePre = prices_m,
    pricePost = prices_m,
    shareInside = share_inside_for_object,
    normIndex = NA,
    output = is_output,
    labels = market$labels,
    supply_model = params$supply_model,
    use_cutoff = use_cutoff,
    cutoff_share = as.numeric(cutoff_draw),
    cutoff_scale = params$cutoff_scale[draw_idx],
    bargaining_power = bargaining_draw
  ), warning = function(condition) {
    ## A conduct realization needs both ownership slots before a scenario
    ## exists. Silence only this exact no-change notice for that baseline;
    ## validation, domain and numerical warnings remain visible.
    if (identical(ownerPre, ownerPost) &&
        identical(conditionMessage(condition),
                  "'ownerPost' and 'ownerPre' are the same")) {
      tryInvokeRestart("muffleWarning")
    }
  })
  obj@slopes <- slopes
  obj@ownerPre <- antitrust::ownerToMatrix(obj, TRUE)
  obj@ownerPost <- antitrust::ownerToMatrix(obj, FALSE)
  obj
}

## Map each no-cutoff Bayesian supply model to the complete antitrust
## registry entry with the same conduct equation: 1 Bertrand, 2 MonCom,
## 3 Cournot, 4 second-price auction, and 5 bargaining. Cutoff draws and
## the legacy multiproduct Logit bargaining equation remain on the custom
## BayesConduct class, while their specification still records the intended
## registered conduct.
.bayes_model_spec <- function(params) {
  if (!is.list(params) || is.null(params$demand_model) ||
      is.null(params$supply_model)) {
    stop("'params' must provide demand_model and supply_model.")
  }
  demand <- params$demand_model
  supply_model <- as.integer(params$supply_model)
  if (length(demand) != 1L || !demand %in% c("logit", "ces") ||
      length(supply_model) != 1L || is.na(supply_model) ||
      !supply_model %in% 1:5) {
    stop("'params' has an unsupported Bayesian demand or supply model.")
  }
  conduct <- switch(as.character(supply_model),
    "1" = "bertrand",
    "2" = "moncom",
    "3" = "cournot",
    "4" = "auction2nd",
    "5" = "bargaining",
    stop("'params' has an unsupported Bayesian supply model.")
  )
  antitrust::model_spec(demand, conduct)
}

.bayes_core_conduct <- function(object) {
  if (isTRUE(object@use_cutoff)) return(FALSE)
  ## Bayes Logit bargaining uses a firm-share scalar, while the registered
  ## antitrust BargainingLogit uses a product-level ownership system. Keep
  ## the Bayes class until those equations are made identical.
  if (methods::is(object, "BayesConductLogit") &&
      object@supply_model == 5L) return(FALSE)
  TRUE
}

## Build a registered antitrust model from the already specified posterior
## demand parameters. antitrust::specify() is the public parameterized
## constructor; it does not recalibrate demand. Supplying ownerPre for both
## sides makes this realization scenario-neutral: all post-merger changes are
## applied later by antitrust::simulate().
.bayes_core_model <- function(object, spec, ...) {
  core_parameters <- object@slopes
  margins <- antitrust::calcMargins(object, preMerger = TRUE, level = FALSE)

  ## Auction2ndLogit's legacy demand state treats meanval as log shares and
  ## keeps those shares fixed while solving the auction. Translate the
  ## posterior mean values to that registered representation while retaining
  ## the posterior alpha and observed shares exactly.
  if (spec$conduct == "auction2nd" &&
      methods::is(object, "BayesConductLogit")) {
    outside_share <- 1 - object@shareInside
    if (!is.finite(outside_share) || outside_share <= 0) {
      stop("Auction posterior realization requires a positive outside share.")
    }
    core_parameters$meanval <- log(object@shares / outside_share)
  }

  args <- list(
    demand = spec$demand,
    conduct = spec$conduct,
    prices = object@prices,
    parameters = core_parameters,
    ownerPre = object@ownerPre,
    shares = object@shares,
    margins = margins,
    insideSize = object@insideSize,
    priceOutside = if (spec$conduct == "auction2nd" &&
                       methods::is(object, "BayesConductLogit")) {
      0
    } else {
      object@priceOutside
    },
    priceStart = object@priceStart,
    labels = object@labels,
    output = object@output
  )
  if (spec$conduct == "bargaining") {
    beta <- as.numeric(object@bargaining_power)
    args$bargpowerPre <- rep(beta, length(object@prices))
    args$bargpowerPost <- rep(beta, length(object@prices))
  }

  ## Keep public constructor controls available to callers. Scenario fields
  ## are rejected by specify(), where they belong in simulate().
  dots <- list(...)
  duplicate <- intersect(names(dots), names(args))
  if (length(duplicate)) {
    stop("argument(s) supplied more than once: ",
         paste(duplicate, collapse = ", "))
  }
  args <- c(args, dots)

  ## Newer antitrust constructors can bind an observed structural baseline
  ## without solving an unchanged post state. Keep the public legacy route
  ## usable with older antitrust installations; explicit lifecycle controls
  ## supplied through ... still take precedence.
  if ("baseline" %in% names(formals(antitrust::specify)) &&
      !"baseline" %in% names(args)) {
    args$baseline <- "observed"
  }

  ## specify() returns the registered model with its baseline cost state and
  ## baseline prices initialized by antitrust. No private state or slots are
  ## rewritten here, and no counterfactual equilibrium is solved locally.
  ## The core constructor still uses these shares for market-mass
  ## normalization, despite its non-BLP notice about demand inversion. Keep
  ## that state and silence only this exact, known notice in this adapter.
  shares_notice <- paste0(
    "'shares' argument is only used for BLP demand (for BLP contraction mapping). For '",
    if (spec$demand == "ces") "CES" else "Logit", "', shares will be ignored.")
  model <- withCallingHandlers(do.call(antitrust::specify, args),
    warning = function(condition) {
      if (identical(conditionMessage(condition), shares_notice)) {
        tryInvokeRestart("muffleWarning")
      }
    })@model

  model
}

#' Wrap a Bayesian conduct object as an antitrust fit
#'
#' Specifies one posterior realization, recovers marginal costs from the
#' conduct model FOC (the same \code{calcMC()} fallback
#' \code{\link{bayes_sim}} relies on), and wraps the result as an
#' \code{antitrust::AntitrustFit}, so it can be passed to
#' \code{antitrust::counterfactual()} and \code{antitrust::simulate()}
#' without demand recalibration.
#'
#' @param object A \code{BayesConductLogit} or \code{BayesConductCES} built
#'   by the package's posterior-realization helper.
#' @param spec An \code{antitrust_model_spec}. Defaults to a mapping from
#'   \code{object@supply_model} (see \code{\link{model_spec}} in
#'   \pkg{antitrust}).
#' @param ... Additional arguments passed to the public
#'   \code{antitrust::specify()} constructor for the baseline realization.
#'   When supported, its observed-baseline lifecycle avoids an unchanged
#'   equilibrium solve; \code{baseline = "solve"} selects the legacy path.
#'   Counterfactual scenario options belong to \code{antitrust::simulate()}.
#'
#' @return An \code{antitrust::AntitrustFit}.
#' @export
as_antitrust_fit <- function(object, spec = NULL, ...) {
  if (!methods::is(object, "BayesConductLogit") &&
      !methods::is(object, "BayesConductCES")) {
    stop("'object' must be a BayesConductLogit or BayesConductCES.")
  }
  expected_spec <- .bayes_model_spec(list(
    demand_model = if (methods::is(object, "BayesConductCES")) "ces" else "logit",
    supply_model = object@supply_model
  ))
  if (is.null(spec)) {
    spec <- expected_spec
  } else if (!inherits(spec, "antitrust_model_spec") ||
             !identical(spec, expected_spec)) {
    stop("'spec' must match the Bayesian demand and supply model.")
  }

  use_core <- .bayes_core_conduct(object)
  model <- if (use_core) {
    .bayes_core_model(object, spec, ...)
  } else {
    ## Custom cutoff and multiproduct Logit bargaining equations remain
    ## Bayes-owned, but the returned fit still represents the observed
    ## baseline. Scenario state carried by a manually constructed conduct
    ## object belongs in simulate(), just as it does for registered models.
    model <- object
    nprod <- length(model@prices)
    model@ownerPost <- model@ownerPre
    model@mcDelta <- rep(0, nprod)
    model@subset <- rep(TRUE, nprod)
    model@priceStart <- model@prices
    model@pricePre <- model@prices
    model@pricePost <- model@prices
    model@mcPre <- antitrust::calcMC(model, preMerger = TRUE)
    model@mcPost <- antitrust::calcMC(model, preMerger = FALSE)
    antitrust::initialize_cost_state(model)
  }
  model@mktSize <- model@insideSize / model@shareInside

  custom_reason <- if (isTRUE(object@use_cutoff)) {
    "cutoff_conduct_blend"
  } else if (methods::is(object, "BayesConductLogit") &&
             object@supply_model == 5L) {
    "multiproduct_logit_bargaining_equation"
  } else {
    NULL
  }

  methods::new(
    "AntitrustFit",
    spec = spec,
    model = model,
    parameters = object@slopes,
    observed = list(
      prices = object@prices,
      margins = object@margins,
      ownerPre = object@ownerPre,
      labels = object@labels,
      priceOutside = object@priceOutside
    ),
    diagnostics = list(
      status = "completed",
      source = "posterior_realization",
      route = "as_antitrust_fit",
      model_class = class(model)[[1]],
      specification = "one posterior realization",
      demand_recalibrated = FALSE,
      cost_recovery = "model_foc",
      core_model = use_core,
      custom_economics = custom_reason
    )
  )
}
