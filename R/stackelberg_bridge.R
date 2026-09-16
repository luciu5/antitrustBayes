# Public reconstruction and posterior summaries for the native coordination
# Stackelberg classes.  The bridge deliberately calls coordination directly;
# it never registers a Stackelberg realization as an ordinary antitrust model.

.stack_bridge_or <- function(x, y) if (is.null(x) || !length(x)) y else x[[1L]]

.stack_component_role_probability <- function(k, probability, n_firm) {
  k <- as.integer(k)
  probability <- as.numeric(probability)
  if (length(k) != length(probability) || !length(k) || anyNA(k) ||
      any(!is.finite(probability)) || any(probability < 0) ||
      sum(probability) <= 0 || length(n_firm) != 1L ||
      !is.finite(n_firm) || n_firm != as.integer(n_firm) || n_firm < 0L) {
    stop("The ordered-k component distribution is malformed.", call. = FALSE)
  }
  probability <- probability / sum(probability)
  vapply(seq_len(as.integer(n_firm)), function(rank) {
    sum(probability[k >= rank])
  }, numeric(1))
}

.stack_discrete_quantile <- function(k, probability, prob) {
  k <- as.integer(k)
  probability <- as.numeric(probability)
  if (length(k) != length(probability) || !length(k) || anyNA(k) ||
      any(!is.finite(probability)) || any(probability < 0) ||
      sum(probability) <= 0 || length(prob) != 1L || !is.finite(prob) ||
      prob < 0 || prob > 1) {
    stop("The ordered-k quantile inputs are malformed.", call. = FALSE)
  }
  ord <- order(k)
  k <- k[ord]
  probability <- probability[ord] / sum(probability)
  k[[which(cumsum(probability) >= prob - 1e-12)[1L]]]
}

.stack_bridge_units <- function(fit, market = NULL, year = NULL) {
  if (methods::is(fit, "BayesSpecifiedFit")) {
    d <- fit@data
    return(list(data = d, idx = seq_len(nrow(d)), market = "unit_1", year = "1",
                market_index = 1L, market_year_index = 1L))
  }
  d <- fit@data
  if (is.null(market)) market <- unique(as.character(d$market))
  if (is.null(year)) year <- unique(as.character(d$year))
  keep <- as.character(d$market) %in% as.character(market) &
    as.character(d$year) %in% as.character(year)
  ids <- unique(paste(d$market[keep], d$year[keep], sep = "\r"))
  if (length(ids) != 1L) {
    stop("bayes_stackelberg_model() requires one unambiguous market-year cell.")
  }
  idx <- which(keep)
  mm <- match(as.character(d$market[idx[1L]]),
              attr(fit@stan_data, "market_levels", exact = TRUE))
  mky <- fit@stan_data$market_year_idx[idx[1L]]
  list(data = d[idx, , drop = FALSE], idx = idx,
       market = as.character(d$market[idx[1L]]),
       year = as.character(d$year[idx[1L]]), market_index = mm,
       market_year_index = as.integer(mky))
}

.stack_bridge_sizes <- function(df, demand_model) {
  q <- as.numeric(df$quantities)
  g <- as.character(df$group)
  product_share <- if (identical(demand_model, "ces")) {
    p <- as.numeric(df$prices)
    (p * q) / sum(p * q)
  } else q / sum(q)
  ## Match .stack_unique_firm_mky: normalize products first, aggregate by
  ## firm, then rank by the declared character firm ID as the tie breaker.
  firm_ids <- sort(unique(g))
  out <- vapply(firm_ids, function(id) sum(product_share[g == id]), numeric(1))
  names(out) <- firm_ids
  out[order(-out, names(out))]
}

## Ordered-k estimation fixes one ranking for the complete market panel.
## Preserve that ranking during reconstruction while retaining the selected
## period's shares as the size values used by diagnostics. Firms absent from
## the selected period keep their market rank and receive size zero.
.stack_bridge_ranked_sizes <- function(fit, unit, df, demand_model) {
  local <- .stack_bridge_sizes(df, demand_model)
  prior <- .stack_bridge_or(fit@metadata$partition_prior, "cutoff")
  if (!identical(prior, "ordered_k")) return(local)

  ranking <- fit@metadata$firm_market_ranking
  required <- c(
    "firm_market_market", "firm_market_group_label", "firm_market_rank"
  )
  if (!is.list(ranking) || !all(required %in% names(ranking))) {
    stop("Exact ordered-k reconstruction requires the stored firm-market ranking.",
         call. = FALSE)
  }
  market_index <- if (methods::is(fit, "BayesSpecifiedFit")) {
    1L
  } else {
    as.integer(unit$market_index)
  }
  if (length(market_index) != 1L || is.na(market_index)) {
    stop("Exact ordered-k reconstruction has no valid market index.",
         call. = FALSE)
  }
  keep <- as.integer(ranking$firm_market_market) == market_index
  ids <- as.character(ranking$firm_market_group_label[keep])
  ranks <- as.integer(ranking$firm_market_rank[keep])
  n <- length(ids)
  if (!n || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids) ||
      anyNA(ranks) || !identical(sort(ranks), seq_len(n))) {
    stop("The stored firm-market ranking is malformed.", call. = FALSE)
  }
  if (any(!names(local) %in% ids)) {
    stop("A simulation firm is absent from the stored firm-market ranking.",
         call. = FALSE)
  }
  out <- stats::setNames(numeric(n), ids)
  out[names(local)] <- local
  out[order(ranks)]
}

.stack_bridge_params <- function(fit, unit, draw, draws = NULL) {
  if (methods::is(fit, "BayesSpecifiedFit")) {
    p <- fit@params
    n <- p$n_draws
    if (draw < 1L || draw > n) stop("'draw' is outside the supplied posterior draws.")
    return(list(
      alpha = if (identical(p$demand_model, "logit")) p$alpha[draw] else NULL,
      gamma = if (identical(p$demand_model, "ces")) p$gamma[draw] else NULL,
      s0 = p$s0[draw], structure = p$structure,
      cutoff = p$cutoff_share[draw], cutoff_scale = p$cutoff_scale[draw],
      partition_prior = p$partition_prior, k_support = p$k_support,
      k_prior = p$k_prior,
      k_scope = p$k_scope, supply_model = p$supply_model
    ))
  }
  dr <- if (is.null(draws)) rstan::extract(fit@stanfit) else draws
  p <- .estimated_draw_params(fit, dr, draw, unit)
  list(
    alpha = p$alpha[[1L]], gamma = p$gamma[[1L]], s0 = p$s0[[1L]],
    structure = p$structure, cutoff = p$cutoff_share[[1L]],
    cutoff_scale = p$cutoff_scale[[1L]], partition_prior = p$partition_prior,
    k_support = p$k_support, k_prior = p$k_prior,
    k_scope = p$k_scope, supply_model = p$supply_model
  )
}

## rstan drops singleton vector dimensions during extraction.  Normalize the
## few generated-quantity slices used by the native bridge while rejecting an
## absent or malformed estimated-fit quantity instead of silently reverting to
## a prior-only role draw.
.stack_bridge_draw_slice <- function(x, draw, n_draw, start, count, name) {
  if (is.null(x)) stop("Estimated fit is missing generated quantity '", name, "'.",
                       call. = FALSE)
  xdim <- dim(x)
  if (is.null(xdim)) {
    if (count == 1L && length(x) == n_draw) return(as.numeric(x[[draw]]))
    if (n_draw == 1L && length(x) >= start + count - 1L) {
      return(as.numeric(x[start:(start + count - 1L)]))
    }
    stop("Estimated generated quantity '", name,
         "' has malformed singleton dimensions.", call. = FALSE)
  }
  if (length(xdim) != 2L || xdim[[1L]] != n_draw ||
      xdim[[2L]] < start + count - 1L) {
    stop("Estimated generated quantity '", name,
         "' has malformed dimensions.", call. = FALSE)
  }
  as.numeric(x[draw, start:(start + count - 1L), drop = TRUE])
}

.stack_bridge_prior <- function(params, sizes, group, market_label = NULL) {
  n <- length(sizes)
  if (identical(params$partition_prior, "ordered_k")) {
    supports <- list(as.integer(params$k_support))
    units <- data.frame(market = 1L, market_label = .stack_bridge_or(market_label, "unit_1"))
    pr <- .stack_validate_k_prior(params$k_prior, supports, units, "market",
                                  k_support = params$k_support)[[1L]]
    return(list(k = supports[[1L]], probability = as.numeric(pr), cutoff = NULL))
  }
  lp <- stackelberg_partition_logprob(sizes, params$cutoff,
                                     params$cutoff_scale)
  ## Stackelberg estimation conditions on at least one leader.  The raw
  ## boundary helper retains k = 0 for deterministic nesting checks, while
  ## the production prior renormalizes the positive-count components.
  lp <- lp[-1L]
  lp <- lp - .stack_log_inv_logit((sizes[1L] - params$cutoff) /
                                  params$cutoff_scale)
  list(k = seq_len(n), probability = exp(lp), cutoff = params$cutoff)
}

.stack_bridge_hard_k <- function(sizes, cutoff) {
  above <- sum(sizes > cutoff)
  if (above > 0L) return(as.integer(above))
  as.integer(sum(sizes == max(sizes)))
}

.stack_bridge_component_prob <- function(fit, unit, params, sizes, draw,
                                         draws = NULL) {
  if (!methods::is(fit, "BayesSpecifiedFit")) {
    d <- fit@stan_data
    pu <- if (identical(params$k_scope, "market")) {
      d$mky_to_partition_unit[unit$market_year_index]
    } else d$mky_to_partition_unit[unit$market_year_index]
    cs <- d$partition_unit_start[pu]
    ce <- cs + d$partition_unit_count[pu] - 1L
    dr <- if (is.null(draws)) rstan::extract(fit@stanfit) else draws
    if (is.null(dr$a_market)) {
      stop("Estimated fit is missing the demand-draw index needed for Stackelberg reconstruction.",
           call. = FALSE)
    }
    n_draw <- if (is.null(dim(dr$a_market))) length(dr$a_market) else dim(dr$a_market)[[1L]]
    if (length(n_draw) != 1L || !is.finite(n_draw) || n_draw < 1L) {
      stop("Estimated fit has malformed posterior draw dimensions.", call. = FALSE)
    }
    z <- .stack_bridge_draw_slice(dr$posterior_component_prob, draw, n_draw,
                                  cs, ce - cs + 1L,
                                  "posterior_component_prob")
    if (length(z) != ce - cs + 1L || any(!is.finite(z)) || sum(z) <= 0) {
      stop("Estimated generated quantity 'posterior_component_prob' is malformed.",
           call. = FALSE)
    }
    return(list(k = d$partition_component_k[cs:ce],
                probability = as.numeric(z / sum(z)),
                cutoff = if (identical(params$partition_prior, "cutoff")) params$cutoff else NULL))
  }
  out <- .stack_bridge_prior(params, sizes, group = NULL,
                             market_label = unit$market)
  if (identical(params$partition_prior, "ordered_k")) out$cutoff <- NULL
  out
}

.stack_bridge_select_k <- function(component, params, sizes, role_mapping,
                                   fixed_k, role_offset, support_min = NULL,
                                   support_max = NULL) {
  if (is.null(support_min)) support_min <- min(component$k)
  if (is.null(support_max)) support_max <- max(component$k)
  if (length(support_min) != 1L || !is.finite(support_min) ||
      support_min != as.integer(support_min) || support_min < 0L ||
      length(support_max) != 1L || !is.finite(support_max) ||
      support_max != as.integer(support_max) || support_max < support_min ||
      support_max > length(sizes)) {
    stop("The supplied ordered-k partition support is malformed.")
  }
  if (!is.null(fixed_k)) {
    if (length(fixed_k) != 1L || !is.finite(fixed_k) || fixed_k != as.integer(fixed_k) ||
        !as.integer(fixed_k) %in% component$k) {
      stop("'fixed_k' must be an integer in the supplied support {",
           paste(component$k, collapse = ", "), "}.")
    }
    return(as.integer(fixed_k))
  }
  if (length(role_offset) != 1L || !is.finite(role_offset) || role_offset != as.integer(role_offset)) {
    stop("'role_offset' must be an integer scalar.")
  }
  role_mapping <- match.arg(role_mapping, c("mixture", "map", "threshold"))
  base_k <- if (role_mapping == "threshold") {
    if (is.null(component$cutoff)) stop("role_mapping = 'threshold' requires partition_prior = 'cutoff'.")
    .stack_bridge_hard_k(sizes, component$cutoff)
  } else if (role_mapping == "map") {
    component$k[which.max(component$probability)]
  } else {
    sample(component$k, 1L, prob = component$probability)
  }
  k <- base_k + as.integer(role_offset)
  if (!k %in% component$k) {
    stop("'role_offset' selects a role count outside the supplied support {",
         paste(component$k, collapse = ", "), "}.")
  }
  ## A sensitivity offset changes the rank boundary.  Reject only a boundary
  ## that actually splits an observed tie; proximity of the continuous
  ## cutoff to a firm size is not itself evidence of an unsupported split.
  if (role_offset != 0L && !is.null(component$cutoff) && k > 0L &&
      k < length(sizes) &&
      isTRUE(sizes[[k]] == sizes[[k + 1L]])) {
    stop("'role_offset' cannot split a tied firm-size boundary.")
  }
  as.integer(k)
}

#' Reconstruct a native coordination Stackelberg model from a Bayesian draw.
#' @param fit A \code{BayesFit} or \code{BayesSpecifiedFit} with
#'   \code{structure = "stackelberg"}.
#' @param draw Posterior draw index.
#' @param market,year Optional selectors; exactly one market-year is required.
#' @param role_mapping One of \code{"mixture"}, \code{"map"}, or
#'   \code{"threshold"}.
#' @param fixed_k Optional explicit leader count.
#' @param role_offset Integer adjacent-partition sensitivity offset.
#' @param insideSize Bayesian inside quantity convention.
#' @param priceOutside Outside-good price override.
#' @param labels Product labels override.
#' @param ... Passed to \code{coordination::stackelberg()}.
#' @return A native \code{coordination::StackelbergLogit} or
#'   \code{coordination::StackelbergCES} object.
#' The reconstruction anchors mean utilities/values to observed conditional
#' product shares for each demand draw; this anchoring is recorded in the
#' returned diagnostics. CES `insideSize` is converted once to inside
#' expenditure before entering the native coordination model.
#' @export
bayes_stackelberg_model <- function(fit, draw = 1L, market = NULL, year = NULL,
                                    role_mapping = c("mixture", "map", "threshold"),
                                    fixed_k = NULL, role_offset = 0L,
                                    insideSize = NULL, priceOutside = NULL,
                                    labels = NULL, ...) {
  if (!methods::is(fit, "BayesFit")) stop("'fit' must be a BayesFit object.")
  structure_name <- .stack_bridge_or(fit@metadata$structure,
    if (isTRUE(fit@stan_data$use_cutoff == 1L)) "fringe" else "simultaneous")
  if (!identical(structure_name, "stackelberg")) stop("fit@metadata$structure must be 'stackelberg'.")
  role_mapping <- match.arg(role_mapping)
  if (!requireNamespace("coordination", quietly = TRUE) ||
      utils::packageVersion("coordination") < "0.0.0.9003") {
    stop("native Bayesian reconstruction requires coordination >= 0.0.0.9003.")
  }
  unit <- .stack_bridge_units(fit, market, year)
  draws <- if (methods::is(fit, "BayesSpecifiedFit")) NULL else {
    rstan::extract(fit@stanfit)
  }
  .stack_bridge_model_draw(
    fit = fit, unit = unit, draw = as.integer(draw),
    role_mapping = role_mapping, fixed_k = fixed_k, role_offset = role_offset,
    insideSize = insideSize, priceOutside = priceOutside, labels = labels,
    draws = draws, dots = list(...)
  )
}

## Internal draw-aware reconstruction.  Simulation passes its already
## extracted draw arrays here so each draw is indexed without a second full
## rstan::extract() call.  The public wrapper above remains responsible for
## validating selectors and performing one extraction when needed.
.stack_bridge_model_draw <- function(fit, unit, draw,
                                     role_mapping, fixed_k, role_offset,
                                     insideSize, priceOutside, labels,
                                     draws = NULL, dots = list()) {
  structure_name <- .stack_bridge_or(fit@metadata$structure,
    if (isTRUE(fit@stan_data$use_cutoff == 1L)) "fringe" else "simultaneous")
  if (!identical(structure_name, "stackelberg")) {
    stop("fit@metadata$structure must be 'stackelberg'.", call. = FALSE)
  }
  params <- .stack_bridge_params(fit, unit, as.integer(draw), draws = draws)
  df <- unit$data
  demand <- fit@demand_model
  sizes <- .stack_bridge_ranked_sizes(fit, unit, df, demand)
  comp <- .stack_bridge_component_prob(
    fit, unit, params, sizes, as.integer(draw), draws = draws
  )
  k <- .stack_bridge_select_k(
    comp, params, sizes, role_mapping, fixed_k, role_offset,
    support_max = max(comp$k)
  )
  ranked <- names(sizes)
  core_ids <- if (k) ranked[seq_len(k)] else character()
  ## A market-level core firm may be absent in one panel period. Its rank is
  ## retained, but only active owners are passed to the deterministic game.
  leaders <- intersect(core_ids, unique(as.character(df$group)))
  p <- as.numeric(df$prices)
  q <- as.numeric(df$quantities)
  t <- q / sum(q)
  p0 <- if (is.null(priceOutside)) {
    if ("outside_price" %in% names(df)) unique(df$outside_price)[1L] else .stack_bridge_or(fit@observed$priceOutside, 1)
  } else as.numeric(priceOutside)
  if (length(p0) != 1L || !is.finite(p0) || p0 <= 0) stop("priceOutside must be positive and finite.")
  inside <- if (is.null(insideSize)) {
    if (methods::is(fit, "BayesSpecifiedFit")) fit@market$insideSize else 1
  } else as.numeric(insideSize)
  if (length(inside) != 1L || !is.finite(inside) || inside <= 0) stop("insideSize must be positive and finite.")
  if (is.null(labels)) labels <- if (methods::is(fit, "BayesSpecifiedFit")) fit@market$labels else paste0("product_", seq_along(p))
  if (length(labels) != length(p)) stop("labels must have one value per product.")
  if (identical(demand, "ces")) {
    avgp <- sum(p * t)
    R0 <- params$s0 * p0 / (params$s0 * p0 + (1 - params$s0) * avgp)
    shares <- p * t * (1 - R0) / avgp
    coord_inside <- inside * avgp
    args <- list(prices = p, shares = shares, margins = as.numeric(df$margins),
                 ownerPre = as.character(df$group), leadersPre = leaders,
                 demand = "ces", conduct = if (params$supply_model == 1L) "bertrand" else "cournot",
                 output = TRUE, insideSize = coord_inside, normIndex = NA,
                 priceOutside = p0, gamma = params$gamma, labels = labels,
                 control.equ = list(implicitCheck = FALSE))
  } else {
    shares <- t * (1 - params$s0)
    args <- list(prices = p, shares = shares, margins = as.numeric(df$margins),
                 ownerPre = as.character(df$group), leadersPre = leaders,
                 demand = "logit", conduct = if (params$supply_model == 1L) "bertrand" else "cournot",
                 output = identical(fit@market_side, "output"), insideSize = inside,
                 normIndex = NA, priceOutside = p0, alpha = params$alpha,
                 labels = labels, control.equ = list(implicitCheck = FALSE))
  }
  if (length(dots)) {
    dup <- intersect(names(dots), names(args)); if (length(dup)) stop("argument(s) supplied more than once: ", paste(dup, collapse = ", "))
    args <- c(args, dots)
  }
  out <- do.call(coordination::stackelberg, args)
  residual_diag <- out@diagnostics$residualsPre
  residual_vector <- if (is.list(residual_diag)) residual_diag$all else residual_diag
  normalized_diag <- if (is.list(residual_diag)) residual_diag$normalizedAll else {
    residual_vector / pmax(abs(p), 1e-12)
  }
  out@diagnostics$bayes_reconstruction <- list(
    draw = as.integer(draw), market = unit$market, year = unit$year,
    structure = structure_name, role_mapping = role_mapping, selected_k = k,
    leaders = leaders, market_core_ids = core_ids,
    component_k = comp$k, component_probability = comp$probability,
    partition_prior = params$partition_prior, k_scope = params$k_scope,
    normalization = .stack_whole_or(fit@metadata$prior_normalization, list(
      mode = "native_stan_normalized", model_spec_version = 9L,
      partition_support = "1:N"
    )),
    demand_anchoring = list(
      method = "observed_conditional_shares",
      description = paste(
        "Mean utilities/values are reconstructed from the observed",
        "conditional product shares for each posterior demand draw."
      ),
      conditional_share_residual = NA_real_
    ),
    residuals = residual_vector,
    normalized_residuals = normalized_diag
  )
  out
}

#' Reconstruct a native coordination exact core-fringe model
#'
#' Reconstructs one hard ordered-k realization from a Bayesian fringe fit.
#' The top k firms in the stored whole-market ranking form the simultaneous
#' strategic core; all remaining active firms are monopolistically competitive.
#'
#' @param fit A \code{BayesFit} with \code{structure = "fringe"} and
#'   \code{partition_prior = "ordered_k"}.
#' @param draw Posterior draw index.
#' @param market,year Optional selectors; exactly one market-year is required.
#' @param role_mapping One of \code{"mixture"} or \code{"map"}.
#' @param fixed_k Optional explicit strategic-core count, including zero.
#' @param role_offset Integer adjacent-partition sensitivity offset.
#' @param insideSize Bayesian inside quantity convention.
#' @param priceOutside Outside-good price override.
#' @param labels Product labels override.
#' @param ... Passed to \code{coordination::core_fringe()}.
#' @return A native \code{coordination::CoreFringeLogit} or
#'   \code{coordination::CoreFringeCES} object.
#' @export
bayes_core_fringe_model <- function(fit, draw = 1L, market = NULL, year = NULL,
                                    role_mapping = c("mixture", "map"),
                                    fixed_k = NULL, role_offset = 0L,
                                    insideSize = NULL, priceOutside = NULL,
                                    labels = NULL, ...) {
  if (!methods::is(fit, "BayesFit")) stop("'fit' must be a BayesFit object.")
  structure_name <- .stack_bridge_or(
    fit@metadata$structure,
    if (isTRUE(fit@stan_data$use_cutoff == 1L)) "fringe" else "simultaneous"
  )
  partition_prior <- .stack_bridge_or(fit@metadata$partition_prior, "cutoff")
  if (!identical(structure_name, "fringe") ||
      !identical(partition_prior, "ordered_k")) {
    stop("fit must be an exact ordered-k fringe model.", call. = FALSE)
  }
  role_mapping <- match.arg(role_mapping)
  if (!requireNamespace("coordination", quietly = TRUE) ||
      utils::packageVersion("coordination") < "0.0.0.9003") {
    stop("native Bayesian reconstruction requires coordination >= 0.0.0.9003.")
  }
  unit <- .stack_bridge_units(fit, market, year)
  draws <- if (methods::is(fit, "BayesSpecifiedFit")) NULL else {
    rstan::extract(fit@stanfit)
  }
  .core_fringe_bridge_model_draw(
    fit = fit, unit = unit, draw = as.integer(draw),
    role_mapping = role_mapping, fixed_k = fixed_k, role_offset = role_offset,
    insideSize = insideSize, priceOutside = priceOutside, labels = labels,
    draws = draws, dots = list(...)
  )
}

.core_fringe_bridge_model_draw <- function(fit, unit, draw,
                                           role_mapping, fixed_k, role_offset,
                                           insideSize, priceOutside, labels,
                                           draws = NULL, dots = list()) {
  structure_name <- .stack_bridge_or(fit@metadata$structure, "simultaneous")
  partition_prior <- .stack_bridge_or(fit@metadata$partition_prior, "cutoff")
  if (!identical(structure_name, "fringe") ||
      !identical(partition_prior, "ordered_k")) {
    stop("fit must be an exact ordered-k fringe model.", call. = FALSE)
  }
  params <- .stack_bridge_params(fit, unit, as.integer(draw), draws = draws)
  df <- unit$data
  demand <- fit@demand_model
  sizes <- .stack_bridge_ranked_sizes(fit, unit, df, demand)
  comp <- .stack_bridge_component_prob(
    fit, unit, params, sizes, as.integer(draw), draws = draws
  )
  k <- .stack_bridge_select_k(
    comp, params, sizes, role_mapping, fixed_k, role_offset,
    support_min = min(comp$k), support_max = max(comp$k)
  )
  ranked <- names(sizes)
  core_ids <- if (k > 0L) ranked[seq_len(k)] else character()
  active_core <- intersect(core_ids, unique(as.character(df$group)))
  p <- as.numeric(df$prices)
  q <- as.numeric(df$quantities)
  t <- q / sum(q)
  p0 <- if (is.null(priceOutside)) {
    if ("outside_price" %in% names(df)) unique(df$outside_price)[1L] else 1
  } else as.numeric(priceOutside)
  if (length(p0) != 1L || !is.finite(p0) || p0 <= 0) {
    stop("priceOutside must be positive and finite.")
  }
  inside <- if (is.null(insideSize)) {
    if (methods::is(fit, "BayesSpecifiedFit")) fit@market$insideSize else 1
  } else as.numeric(insideSize)
  if (length(inside) != 1L || !is.finite(inside) || inside <= 0) {
    stop("insideSize must be positive and finite.")
  }
  if (is.null(labels)) {
    labels <- if (methods::is(fit, "BayesSpecifiedFit")) {
      fit@market$labels
    } else {
      paste0("product_", seq_along(p))
    }
  }
  if (length(labels) != length(p)) stop("labels must have one value per product.")
  conduct <- if (params$supply_model == 1L) "bertrand" else "cournot"
  if (identical(demand, "ces")) {
    avgp <- sum(p * t)
    R0 <- params$s0 * p0 / (params$s0 * p0 + (1 - params$s0) * avgp)
    shares <- p * t * (1 - R0) / avgp
    args <- list(
      prices = p, shares = shares, margins = as.numeric(df$margins),
      ownerPre = as.character(df$group), corePre = active_core,
      demand = "ces", conduct = conduct, output = TRUE,
      insideSize = inside * avgp, normIndex = NA, priceOutside = p0,
      gamma = params$gamma, labels = labels
    )
  } else {
    shares <- t * (1 - params$s0)
    args <- list(
      prices = p, shares = shares, margins = as.numeric(df$margins),
      ownerPre = as.character(df$group), corePre = active_core,
      demand = "logit", conduct = conduct,
      output = identical(fit@market_side, "output"), insideSize = inside,
      normIndex = NA, priceOutside = p0, alpha = params$alpha,
      labels = labels
    )
  }
  if (length(dots)) {
    dup <- intersect(names(dots), names(args))
    if (length(dup)) {
      stop("argument(s) supplied more than once: ", paste(dup, collapse = ", "))
    }
    args <- c(args, dots)
  }
  out <- do.call(coordination::core_fringe, args)
  residual_diag <- out@diagnostics$residualsPre
  residual_vector <- if (is.list(residual_diag)) residual_diag$all else residual_diag
  normalized_diag <- if (is.list(residual_diag)) residual_diag$normalizedAll else {
    residual_vector / pmax(abs(p), 1e-12)
  }
  out@diagnostics$bayes_reconstruction <- list(
    draw = as.integer(draw), market = unit$market, year = unit$year,
    structure = structure_name, role_mapping = role_mapping, selected_k = k,
    strategic_core = active_core, market_core_ids = core_ids,
    component_k = comp$k, component_probability = comp$probability,
    partition_prior = params$partition_prior, k_scope = params$k_scope,
    normalization = .stack_whole_or(fit@metadata$prior_normalization, list(
      mode = "native_stan_normalized", model_spec_version = 9L,
      partition_support = "0:N"
    )),
    demand_anchoring = list(
      method = "observed_conditional_shares",
      description = paste(
        "Mean utilities/values are reconstructed from the observed",
        "conditional product shares for each posterior demand draw."
      ),
      conditional_share_residual = NA_real_
    ),
    residuals = residual_vector,
    normalized_residuals = normalized_diag
  )
  out
}

#' Summarize exact hard Stackelberg partitions.
#' @param fit A Stackelberg Bayesian fit.
#' @param level Confidence interval level.
#' @return A list containing one row per market-year and partition metadata.
#' @export
bayes_cutoff_summary <- function(fit, level = 0.95) {
  if (!methods::is(fit, "BayesFit")) stop("'fit' must be a BayesFit object.")
  structure_name <- .stack_bridge_or(fit@metadata$structure, "simultaneous")
  if (!structure_name %in% c("stackelberg", "fringe")) {
    stop("cutoff summaries require structure = 'stackelberg' or 'fringe'.")
  }
  if (length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
    stop("'level' must lie in (0, 1).")
  }

  df <- fit@data
  if (!nrow(df)) stop("The fitted data contain no observations.")
  cells <- split(seq_len(nrow(df)),
                 interaction(df$market, df$year, drop = TRUE, lex.order = TRUE))
  stack <- identical(structure_name, "stackelberg")
  direct <- identical(
    .stack_bridge_or(fit@metadata$partition_prior, "cutoff"),
    "ordered_k"
  )
  estimated <- !methods::is(fit, "BayesSpecifiedFit")

  ## Extract each vector-valued generated quantity once.  rstan drops the
  ## matrix dimension when a vector has length one, so restore it here.
  draws <- NULL
  ndraw <- NULL
  if (estimated) {
    gq_pars <- if (direct || stack) {
      c("posterior_component_prob", "prior_component_prob",
        "firm_leader_prob", "prior_firm_leader_prob", "s0_prob")
    } else {
      c("cutoff_share", "s0_prob")
    }
    if (stack && !direct) gq_pars <- c(gq_pars, "cutoff_share")
    if (identical(fit@demand_model, "ces")) gq_pars <- c(gq_pars, "R0_prob")
    draws <- rstan::extract(fit@stanfit, pars = gq_pars)
    a_market <- rstan::extract(fit@stanfit, pars = "a_market")$a_market
    ndraw <- if (is.null(dim(a_market))) length(a_market) else dim(a_market)[1L]
    if (!is.finite(ndraw) || ndraw < 1L) stop("The fitted object contains no posterior draws.")
  }

  as_draw_matrix <- function(x, ncol, name) {
    if (is.null(x)) stop("Generated quantity '", name, "' is unavailable.")
    if (!is.null(dim(x))) {
      d <- dim(x)
      if (length(d) != 2L || d[1L] != ndraw || d[2L] != ncol) {
        stop("Generated quantity '", name, "' has unexpected dimensions.")
      }
      return(matrix(as.numeric(x), nrow = ndraw, ncol = ncol))
    }
    if (length(x) != ndraw * ncol) {
      stop("Generated quantity '", name, "' has unexpected length.")
    }
    matrix(as.numeric(x), nrow = ndraw, ncol = ncol)
  }
  posterior <- prior <- leader <- cutoff_draws <- NULL
  if (estimated) {
    d <- fit@stan_data
    if (direct || stack) {
      ncomp <- as.integer(d$N_partition_component)
      nfirm <- as.integer(d$N_firm_mky)
      nmarket <- as.integer(d$N_market)
      posterior <- as_draw_matrix(draws$posterior_component_prob, ncomp,
                                  "posterior_component_prob")
      prior <- as_draw_matrix(draws$prior_component_prob, ncomp,
                              "prior_component_prob")
      leader <- as_draw_matrix(draws$firm_leader_prob, nfirm,
                               "firm_leader_prob")
      prior_leader <- if (!is.null(draws$prior_firm_leader_prob)) {
        as_draw_matrix(draws$prior_firm_leader_prob, nfirm,
                       "prior_firm_leader_prob")
      } else NULL
      if (!direct) {
        cutoff_draws <- as_draw_matrix(draws$cutoff_share, nmarket, "cutoff_share")
      }
    } else {
      cutoff_draws <- as_draw_matrix(draws$cutoff_share,
                                     as.integer(d$N_market), "cutoff_share")
    }
  }

  stan_data <- if (estimated) fit@stan_data else NULL
  params <- if (!estimated) fit@params else NULL
  requested_scope <- .stack_bridge_or(fit@metadata$k_scope, "market")
  prior_name <- .stack_bridge_or(fit@metadata$partition_prior, "cutoff")

  quantile_scalar <- function(x, prob) as.numeric(stats::quantile(
    x, prob = prob, names = FALSE, type = 8
  ))
  weighted_quantile <- function(x, w, prob) {
    keep <- is.finite(x) & is.finite(w) & w >= 0
    if (!any(keep) || sum(w[keep]) <= 0) return(NA_real_)
    x <- as.numeric(x[keep]); w <- as.numeric(w[keep])
    ord <- order(x); x <- x[ord]; w <- w[ord] / sum(w)
    x[[which(cumsum(w) >= prob - 1e-12)[1L]]]
  }
  component_for_cell <- function(ii, sizes, market_label, year_label) {
    if (!direct && !stack) return(NULL)
    if (estimated) {
      ## mky indexes are the canonical mapping from original observations to
      ## unique firm-market-year rows and then to the flattened unit array.
      mky <- stan_data$market_year_idx[ii[[1L]]]
      pu <- stan_data$mky_to_partition_unit[mky]
      cstart <- stan_data$partition_unit_start[pu]
      ccount <- stan_data$partition_unit_count[pu]
      cc <- seq.int(cstart, length.out = ccount)
      k <- as.integer(stan_data$partition_component_k[cc])
      post <- posterior[, cc, drop = FALSE]
      pri <- prior[, cc, drop = FALSE]
      post <- post / rowSums(post)
      pri <- pri / rowSums(pri)
      fstart <- stan_data$firm_mky_start[mky]
      fcount <- stan_data$firm_mky_n[mky]
      ff <- seq.int(fstart, length.out = fcount)
      list(k = k, posterior = colMeans(post), prior = colMeans(pri),
           posterior_draws = post, prior_draws = pri, unit = pu,
           mky = mky, leader_probability = colMeans(leader[, ff, drop = FALSE]),
           prior_leader_probability = if (!is.null(prior_leader)) {
             colMeans(prior_leader[, ff, drop = FALSE])
           } else vapply(seq_along(ff), function(j) {
             sum(colMeans(pri)[k >= j])
           }, numeric(1)))
    } else {
      if (direct) {
        comp <- .stack_bridge_prior(
          list(partition_prior = "ordered_k", k_support = params$k_support,
               k_prior = params$k_prior),
          sizes, group = NULL, market_label = "unit_1"
        )
        list(k = comp$k, posterior = comp$probability, prior = comp$probability,
             posterior_draws = matrix(comp$probability, nrow = params$n_draws,
                                      ncol = length(comp$k), byrow = TRUE),
             prior_draws = matrix(comp$probability, nrow = params$n_draws,
                                  ncol = length(comp$k), byrow = TRUE),
             unit = 1L, mky = 1L)
      } else {
        cuts <- as.numeric(params$cutoff_share)
        k <- seq_len(length(sizes))
        pd <- t(vapply(seq_along(cuts), function(draw_idx) {
          .stack_bridge_prior(
            list(partition_prior = "cutoff", cutoff = cuts[[draw_idx]],
                 cutoff_scale = params$cutoff_scale[draw_idx]),
            sizes, group = NULL, market_label = "unit_1"
          )$probability
        }, numeric(length(k))))
        list(k = k, posterior = colMeans(pd), prior = colMeans(pd),
             posterior_draws = pd, prior_draws = pd, unit = 1L, mky = 1L)
      }
    }
  }

  rows <- lapply(cells, function(ii) {
    dcell <- df[ii, , drop = FALSE]
    sizes <- if (direct) {
      unit <- list(
        market_index = if (estimated) {
          as.integer(stan_data$market[ii[[1L]]])
        } else {
          1L
        }
      )
      .stack_bridge_ranked_sizes(fit, unit, dcell, fit@demand_model)
    } else {
      .stack_bridge_sizes(dcell, fit@demand_model)
    }
    firm_ids <- names(sizes)
    if (is.null(firm_ids)) firm_ids <- paste0("firm_", seq_along(sizes))
    base <- list(
      market = as.character(dcell$market[[1L]]),
      year = as.character(dcell$year[[1L]]),
      firm_ids = I(list(as.character(firm_ids))),
      firm_leader_probability = I(list(NULL)),
      firm_strategic_probability = I(list(NULL)),
      firm_core_probability = I(list(NULL)),
      prior_firm_leader_probability = I(list(NULL)),
      realized_k_posterior = I(list(NULL)),
      realized_k_mean = NA_real_, realized_k_median = NA_real_,
      realized_k_lower = NA_real_, realized_k_upper = NA_real_,
      prior_k_pmf = I(list(NULL)),
      prior_k_pmf_type = if (direct) "supplied_market_prior" else if (stack) {
        "conditional_on_draw_averaged"
      } else NA_character_,
      core_share_conditional = NA_real_, core_share_conditional_median = NA_real_,
      core_share_conditional_lower = NA_real_, core_share_conditional_upper = NA_real_,
      core_share_unconditional = NA_real_, core_share_unconditional_median = NA_real_,
      core_share_unconditional_lower = NA_real_, core_share_unconditional_upper = NA_real_,
      requested_k_scope = requested_scope,
      cutoff_mean = NA_real_, cutoff_median = NA_real_,
      cutoff_lower = NA_real_, cutoff_upper = NA_real_,
      raw_k_pmf = I(list(NULL)), raw_k_mean = NA_real_, raw_k_median = NA_real_,
      raw_k_draws = I(list(NULL)),
      hard_k_pmf = I(list(NULL)), hard_k_mean = NA_real_, hard_k_median = NA_real_,
      hard_k_lower = NA_real_, hard_k_upper = NA_real_, hard_k_draws = I(list(NULL)),
      hard_core_share_conditional = NA_real_, hard_core_share_unconditional = NA_real_,
      smooth_role_probability = I(list(NULL)),
      smooth_core_share_conditional = NA_real_, smooth_core_share_unconditional = NA_real_
    )
    if (direct || stack) {
      comp <- component_for_cell(ii, sizes, as.character(dcell$market[[1L]]),
                                 as.character(dcell$year[[1L]]))
      k <- comp$k
      post_p <- as.numeric(comp$posterior); post_p <- post_p / sum(post_p)
      prior_p <- as.numeric(comp$prior); prior_p <- prior_p / sum(prior_p)
      ## Component probabilities and the fixed whole-market rank determine
      ## role probabilities for every ranked firm, including firms absent in
      ## this market-year. The Stan firm GQ only has active firm-year rows.
      role <- .stack_component_role_probability(k, post_p, length(sizes))
      prior_role <- .stack_component_role_probability(k, prior_p, length(sizes))
      names(role) <- firm_ids; names(prior_role) <- firm_ids
      if (stack) {
        base$firm_leader_probability <- I(list(role))
      } else {
        base$firm_strategic_probability <- I(list(role))
      }
      base$firm_core_probability <- I(list(role))
      base$prior_firm_leader_probability <- I(list(prior_role))
      base$realized_k_posterior <- I(list(stats::setNames(post_p, as.character(k))))
      base$realized_k_mean <- sum(k * post_p)
      base$realized_k_median <- .stack_discrete_quantile(k, post_p, .5)
      base$realized_k_lower <- .stack_discrete_quantile(
        k, post_p, (1 - level) / 2
      )
      base$realized_k_upper <- .stack_discrete_quantile(
        k, post_p, 1 - (1 - level) / 2
      )
      ## For cutoff mixtures this is E_theta[pi_k(theta)], while direct
      ## ordered-k fits report the supplied market prior itself.
      base$prior_k_pmf <- I(list(stats::setNames(prior_p, as.character(k))))
      component_core <- vapply(k, function(kk) sum(sizes[seq_len(kk)]), numeric(1))
      s0_draws <- if (estimated) {
        s0 <- as_draw_matrix(draws$s0_prob, as.integer(stan_data$N_market_year), "s0_prob")
        as.numeric(s0[, comp$mky])
      } else as.numeric(params$s0)
      outside_core <- if (identical(fit@demand_model, "ces")) {
        if (estimated) {
          r0 <- as_draw_matrix(draws$R0_prob, as.integer(stan_data$N_market_year), "R0_prob")
          as.numeric(r0[, comp$mky])
        } else {
          pp <- as.numeric(dcell$prices); qq <- as.numeric(dcell$quantities)
          avgp <- sum(pp * qq) / sum(qq)
          p0 <- if ("outside_price" %in% names(dcell)) as.numeric(dcell$outside_price[[1L]]) else fit@observed$priceOutside
          as.numeric(params$s0) * p0 / (as.numeric(params$s0) * p0 +
                                         (1 - as.numeric(params$s0)) * avgp)
        }
      } else s0_draws
      nd <- nrow(comp$posterior_draws)
      core_values <- rep(component_core, times = nd)
      core_weights <- as.vector(t(comp$posterior_draws)) / nd
      uncond_values <- core_values * rep(1 - outside_core, each = length(k))
      base$core_share_conditional <- sum(core_values * core_weights)
      base$core_share_conditional_median <- weighted_quantile(core_values, core_weights, .5)
      base$core_share_conditional_lower <- weighted_quantile(core_values, core_weights, (1 - level) / 2)
      base$core_share_conditional_upper <- weighted_quantile(core_values, core_weights, 1 - (1 - level) / 2)
      base$core_share_unconditional <- sum(uncond_values * core_weights)
      base$core_share_unconditional_median <- weighted_quantile(uncond_values, core_weights, .5)
      base$core_share_unconditional_lower <- weighted_quantile(uncond_values, core_weights, (1 - level) / 2)
      base$core_share_unconditional_upper <- weighted_quantile(uncond_values, core_weights, 1 - (1 - level) / 2)
      if (!direct) {
        cuts <- if (estimated) cutoff_draws[, stan_data$mky_to_market[comp$mky]] else {
          as.numeric(params$cutoff_share)
        }
        raw <- vapply(cuts, function(cut) sum(sizes > cut), integer(1))
        hard <- vapply(cuts, function(cut) .stack_bridge_hard_k(sizes, cut), integer(1))
        hard_core <- vapply(hard, function(kk) sum(sizes[seq_len(kk)]), numeric(1))
        hard_pmf <- tabulate(hard, nbins = length(sizes)) / length(hard)
        raw_pmf <- tabulate(raw + 1L, nbins = length(sizes) + 1L) / length(raw)
        base$cutoff_mean <- mean(cuts)
        base$cutoff_median <- stats::median(cuts)
        base$cutoff_lower <- quantile_scalar(cuts, (1 - level) / 2)
        base$cutoff_upper <- quantile_scalar(cuts, 1 - (1 - level) / 2)
        base$raw_k_pmf <- I(list(stats::setNames(raw_pmf, as.character(0:length(sizes)))))
        base$raw_k_mean <- mean(raw); base$raw_k_median <- stats::median(raw)
        base$raw_k_draws <- I(list(raw))
        base$hard_k_pmf <- I(list(stats::setNames(hard_pmf, as.character(seq_along(sizes)))))
        base$hard_k_mean <- mean(hard); base$hard_k_median <- stats::median(hard)
        base$hard_k_lower <- quantile_scalar(hard, (1 - level) / 2)
        base$hard_k_upper <- quantile_scalar(hard, 1 - (1 - level) / 2)
        base$hard_k_draws <- I(list(hard))
        base$hard_core_share_conditional <- mean(hard_core)
        base$hard_core_share_unconditional <- mean(hard_core * (1 - outside_core))
      }
    } else {
      cuts <- if (estimated) cutoff_draws[, stan_data$mky_to_market[
        stan_data$market_year_idx[ii[[1L]]]
      ]] else as.numeric(params$cutoff_share)
      outside_core <- if (identical(fit@demand_model, "ces")) {
        if (estimated) {
          r0 <- as_draw_matrix(draws$R0_prob, as.integer(stan_data$N_market_year), "R0_prob")
          as.numeric(r0[, stan_data$market_year_idx[ii[[1L]]]])
        } else {
          pp <- as.numeric(dcell$prices); qq <- as.numeric(dcell$quantities)
          avgp <- sum(pp * qq) / sum(qq)
          p0 <- if ("outside_price" %in% names(dcell)) as.numeric(dcell$outside_price[[1L]]) else fit@observed$priceOutside
          as.numeric(params$s0) * p0 / (as.numeric(params$s0) * p0 +
                                         (1 - as.numeric(params$s0)) * avgp)
        }
      } else if (estimated) {
        s0 <- as_draw_matrix(draws$s0_prob, as.integer(stan_data$N_market_year), "s0_prob")
        as.numeric(s0[, stan_data$market_year_idx[ii[[1L]]]])
      } else as.numeric(params$s0)
      raw <- vapply(cuts, function(cut) sum(sizes > cut), integer(1))
      hard <- raw
      hard_core <- vapply(hard, function(kk) {
        if (!kk) 0 else sum(sizes[seq_len(kk)])
      }, numeric(1))
      hard_pmf <- tabulate(hard + 1L, nbins = length(sizes) + 1L) / length(hard)
      tau_source <- if (estimated) stan_data$cutoff_scale else params$cutoff_scale
      if (is.null(tau_source) || !length(tau_source)) tau_source <- .2
      tau_draws <- rep(as.numeric(tau_source), length.out = length(cuts))
      if (any(!is.finite(tau_draws)) || any(tau_draws <= 0)) {
        stop("The fitted cutoff scale must be positive and finite.")
      }
      role_draws <- vapply(seq_along(sizes), function(j) mean(sizes[[j]] > cuts), numeric(1))
      smooth_draws <- vapply(seq_along(sizes), function(j) {
        mean(stats::plogis((sizes[[j]] - cuts) / tau_draws))
      }, numeric(1))
      names(role_draws) <- names(sizes); names(smooth_draws) <- names(sizes)
      base$firm_strategic_probability <- I(list(role_draws))
      base$firm_core_probability <- I(list(role_draws))
      base$prior_firm_leader_probability <- I(list(NULL))
      base$cutoff_mean <- mean(cuts); base$cutoff_median <- stats::median(cuts)
      base$cutoff_lower <- quantile_scalar(cuts, (1 - level) / 2)
      base$cutoff_upper <- quantile_scalar(cuts, 1 - (1 - level) / 2)
      base$raw_k_pmf <- I(list(stats::setNames(hard_pmf, as.character(0:length(sizes)))))
      base$raw_k_mean <- mean(raw); base$raw_k_median <- stats::median(raw)
      base$raw_k_draws <- I(list(raw)); base$hard_k_pmf <- I(list(stats::setNames(hard_pmf, as.character(0:length(sizes)))))
      base$hard_k_mean <- mean(hard); base$hard_k_median <- stats::median(hard)
      base$hard_k_lower <- quantile_scalar(hard, (1 - level) / 2)
      base$hard_k_upper <- quantile_scalar(hard, 1 - (1 - level) / 2)
      base$hard_k_draws <- I(list(hard))
      base$hard_core_share_conditional <- mean(hard_core)
      base$hard_core_share_unconditional <- mean(hard_core * (1 - outside_core))
      base$core_share_conditional <- mean(hard_core)
      base$core_share_conditional_median <- quantile_scalar(hard_core, .5)
      base$core_share_conditional_lower <- quantile_scalar(hard_core, (1 - level) / 2)
      base$core_share_conditional_upper <- quantile_scalar(hard_core, 1 - (1 - level) / 2)
      uncond_core <- hard_core * (1 - outside_core)
      base$core_share_unconditional <- mean(uncond_core)
      base$core_share_unconditional_median <- quantile_scalar(uncond_core, .5)
      base$core_share_unconditional_lower <- quantile_scalar(uncond_core, (1 - level) / 2)
      base$core_share_unconditional_upper <- quantile_scalar(uncond_core, 1 - (1 - level) / 2)
      base$smooth_role_probability <- I(list(smooth_draws))
      smooth_core_draws <- vapply(seq_along(cuts), function(draw_idx) {
        sum(sizes * stats::plogis((sizes - cuts[[draw_idx]]) /
                                  tau_draws[[draw_idx]]))
      }, numeric(1))
      base$smooth_core_share_conditional <- mean(smooth_core_draws)
      ## Keep cutoff and outside-share posterior dependence by multiplying
      ## within draw.  Products of posterior means lose this dependence.
      base$smooth_core_share_unconditional <- mean(
        smooth_core_draws * (1 - outside_core)
      )
    }
    as.data.frame(base, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  base::structure(
    list(
      by_market_year = out,
      cutoff_posterior = if (direct) NULL else TRUE,
      structure = structure_name,
      partition_prior = prior_name,
      k_scope = requested_scope
    ),
    class = "bayes_cutoff_summary"
  )
}
