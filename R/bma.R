validate_comparison_contract <- function(fits) {
  reference <- fits[[1L]]
  columns <- c(
    "prices", "quantities", "margins", "market", "year", "group",
    "margin_observed", "outside_price"
  )
  missing <- setdiff(columns, names(reference@data))
  if (length(missing)) {
    stop("The reference fit lacks comparison data: ",
         paste(missing, collapse = ", "), ". Refit all models.")
  }

  for (i in seq_along(fits)[-1L]) {
    candidate <- fits[[i]]
    if (!identical(candidate@market_side, reference@market_side)) {
      stop("Models must use the same market_side for comparison.")
    }
    if (!all(columns %in% names(candidate@data)) ||
        !isTRUE(all.equal(
          candidate@data[columns], reference@data[columns],
          check.attributes = FALSE
        ))) {
      stop("Models must be fitted to the same ordered observations, margin ",
           "mask, identifiers, and outside prices.")
    }
    common_outcomes <- c("shareIn", "log_shareIn", "marginNormInv",
                         "log_marginNormInv", "margin_obs")
    if (!all(common_outcomes %in% names(candidate@stan_data)) ||
        !all(common_outcomes %in% names(reference@stan_data)) ||
        !isTRUE(all.equal(
          candidate@stan_data[common_outcomes],
          reference@stan_data[common_outcomes],
          check.attributes = FALSE
        ))) {
      stop("Models do not share the current quantity-space likelihood contract.")
    }
  }
  invisible(TRUE)
}

#' Compare Bayesian Models Using Information Criteria and Bridge Sampling
#'
#' Computes LOO-IC, WAIC, and optionally log marginal likelihoods via bridge
#' sampling for one or more fitted \code{\linkS4class{BayesLogit}} objects.
#' Returns a comparison table with model weights.
#'
#' @param ... One or more \code{BayesLogit} objects, optionally named.
#' @param bridge Logical. If \code{TRUE} (default), compute log marginal
#'   likelihoods via \code{\link[bridgesampling]{bridge_sampler}} for Bayes
#'   factor model comparison. If \code{FALSE}, only LOO/WAIC are computed.
#' @param bridge_maxiter Maximum iterations for bridge sampling.
#'   Default 5000.
#' @param prior_weights Optional numeric vector of prior model probabilities.
#'   Must be the same length as the number of models and sum to 1. If
#'   \code{NULL} (default), equal prior probabilities are used.
#'
#' @return A data frame with columns for model name, supply model, LOOIC,
#'   WAIC, stacking weights, and (if \code{bridge = TRUE}) log marginal
#'   likelihoods, Bayes factors (vs. best model), and posterior model
#'   probabilities.
#'
#' @details
#' Three model comparison approaches are provided:
#' \describe{
#'   \item{LOO-IC}{Leave-one-out cross-validation via \pkg{loo}. Stacking
#'     weights reflect predictive performance.}
#'   \item{WAIC}{Widely Applicable Information Criterion.}
#'   \item{Bridge Sampling}{Computes the log marginal likelihood
#'     \eqn{p(y | M_k)} via \code{\link[bridgesampling]{bridge_sampler}}.
#'     Posterior model probabilities assume equal prior model probabilities.
#'     Bayes factors are computed relative to the best model.}
#' }
#'
#' Requires that the Stan model includes a \code{log_lik} generated quantity
#' (for LOO/WAIC). All fits must use the same ordered observations, market
#' side, outside prices, margin mask, and current quantity-space likelihood.
#' This permits Logit/CES comparison without comparing likelihoods defined on
#' different observed outcomes.
#'
#' @examples
#' \donttest{
#' # fit1 <- bayes_est(..., supply_model = "bertrand")
#' # fit2 <- bayes_est(..., supply_model = "cournot")
#' # bayes_compare(Bertrand = fit1, Cournot = fit2)
#' # bayes_compare(Bertrand = fit1, Cournot = fit2, bridge = FALSE)
#' }
#'
#' @export
bayes_compare <- function(..., bridge = TRUE, bridge_maxiter = 5000L,
                          prior_weights = NULL) {
  args <- list(...)

  # Separate BayesLogit objects from named arguments
  fits <- args[vapply(args, is, logical(1), "BayesLogit")]
  if (length(fits) == 0L) {
    stop("At least one BayesLogit fit is required.")
  }
  invisible(lapply(fits, validate_fit_spec))
  validate_comparison_contract(fits)

  # Handle named and unnamed arguments
  nms <- names(fits)
  if (is.null(nms)) {
    nms <- paste0("Model_", seq_along(fits))
  } else {
    empty <- nms == ""
    nms[empty] <- paste0("Model_", which(empty))
  }

  # Prior model probabilities
  if (is.null(prior_weights)) {
    prior_weights <- rep(1 / length(fits), length(fits))
  } else {
    if (length(prior_weights) != length(fits)) {
      stop("'prior_weights' must have the same length as ",
           "the number of models.")
    }
    if (!is.numeric(prior_weights) || any(!is.finite(prior_weights)) ||
        any(prior_weights < 0) || sum(prior_weights) <= 0) {
      stop("'prior_weights' must be finite, non-negative probabilities.")
    }
    if (abs(sum(prior_weights) - 1) > 1e-6) {
      stop("'prior_weights' must sum to 1.")
    }
  }

  supply_names <- c("Bertrand", "Monopolistic", "Cournot", "Auction")

  loo_list  <- list()
  waic_list <- list()
  bridge_list <- list()

  for (i in seq_along(fits)) {
    fit <- fits[[i]]

    log_lik <- loo::extract_log_lik(fit@stanfit, parameter_name = "log_lik",
                                     merge_chains = FALSE)
    r_eff <- loo::relative_eff(exp(log_lik))

    loo_list[[i]]  <- loo::loo(log_lik, r_eff = r_eff)
    waic_list[[i]] <- loo::waic(log_lik)

    # Bridge sampling for marginal likelihood
    if (bridge) {
      bridge_list[[i]] <- tryCatch(
        bridgesampling::bridge_sampler(fit@stanfit, silent = TRUE,
                                        maxiter = bridge_maxiter),
        error = function(e) {
          warning("Bridge sampling failed for model ", nms[i], ": ",
                  conditionMessage(e))
          NULL
        }
      )
    }
  }

  names(loo_list) <- nms
  comp <- loo::loo_compare(loo_list)

  # Stacking weights from LOO
  weights_stacking <- loo::loo_model_weights(loo_list, method = "stacking")

  result <- data.frame(
    model         = nms,
    supply_model  = vapply(fits, function(f)
      supply_names[f@supply_model], character(1)),
    demand_model  = vapply(fits, slot, character(1), "demand_model"),
    market_side   = vapply(fits, slot, character(1), "market_side"),
    looic         = vapply(loo_list, function(x) x$estimates["looic", "Estimate"],
                           numeric(1)),
    looic_se      = vapply(loo_list, function(x) x$estimates["looic", "SE"],
                           numeric(1)),
    p_loo         = vapply(loo_list, function(x) x$estimates["p_loo", "Estimate"],
                           numeric(1)),
    waic          = vapply(waic_list, function(x) x$estimates["waic", "Estimate"],
                           numeric(1)),
    weight_stacking = as.numeric(weights_stacking),
    stringsAsFactors = FALSE
  )

  # Add bridge sampling results
  if (bridge && length(bridge_list) > 0L) {
    logml <- vapply(bridge_list, function(b) {
      if (!is.null(b)) b$logml else NA_real_
    }, numeric(1))

    finite <- is.finite(logml)
    diff_logml <- bf_vs_best <- post_prob <- rep(NA_real_, length(logml))
    if (any(finite)) {
      max_logml <- max(logml[finite])
      diff_logml[finite] <- logml[finite] - max_logml
      bf_vs_best[finite] <- exp(diff_logml[finite])
      unnorm <- exp(diff_logml[finite]) * prior_weights[finite]
      if (sum(unnorm) > 0) {
        post_prob[finite] <- unnorm / sum(unnorm)
      }
    }

    result$logml        <- logml
    result$bf_vs_best   <- bf_vs_best
    result$post_prob    <- post_prob
  }

  rownames(result) <- NULL

  attr(result, "loo_compare") <- comp
  attr(result, "loo_list") <- loo_list
  if (bridge) attr(result, "bridge_list") <- bridge_list
  result
}

#' Extract Model Weights from a Comparison Table
#'
#' Convenience function to extract weights from a \code{\link{bayes_compare}}
#' result using either bridge sampling posterior probabilities or LOO stacking
#' weights.
#'
#' @param comp A data frame returned by \code{\link{bayes_compare}}.
#' @param method Character: \code{"bridge"} (default) for posterior model
#'   probabilities from bridge sampling, or \code{"stacking"} for LOO
#'   stacking weights. Falls back to stacking if bridge results are
#'   unavailable.
#'
#' @return A named numeric vector of weights summing to 1.
#'
#' @examples
#' \donttest{
#' # comp <- bayes_compare(Bertrand = fit1, Cournot = fit2)
#' # bayes_weights(comp, method = "bridge")
#' # bayes_weights(comp, method = "stacking")
#' }
#'
#' @export
bayes_weights <- function(comp, method = c("bridge", "stacking")) {
  method <- match.arg(method)

  bridge_ok <- "post_prob" %in% names(comp) &&
    any(is.finite(comp$post_prob) & comp$post_prob > 0)
  if (method == "bridge" && bridge_ok) {
    w <- comp$post_prob
    w[!is.finite(w)] <- 0
  } else {
    if (method == "bridge" && !bridge_ok) {
      message("Bridge sampling results not available; using stacking weights.")
    }
    w <- comp$weight_stacking
  }

  names(w) <- comp$model
  if (any(!is.finite(w)) || any(w < 0) || sum(w) <= 0) {
    stop("The requested model weights are not finite non-negative values.")
  }
  w / sum(w)
}


#' Bayesian Model Averaging of Merger Simulations
#'
#' Combines merger simulation results from multiple models using
#' model weights (e.g., from \code{\link{bayes_compare}}). Returns a
#' \code{\linkS4class{BayesSim}} object containing a weighted pool of
#' simulation draws, so that all \code{BayesSim} methods
#' (\code{summary}, \code{elast}, \code{diversion}, \code{CV},
#' \code{HypoMonTest}) work on the model-averaged result.
#'
#' @param sims A list of \code{\linkS4class{BayesSim}} objects.
#' @param weights Numeric vector of model weights (must sum to 1).
#'   Typically \code{post_prob} (from bridge sampling) or
#'   \code{weight_stacking} (from LOO) in the output of
#'   \code{\link{bayes_compare}}.
#' @param n_draws Total number of draws in the pooled result.
#'   Default 1000.
#'
#' @return A \code{\linkS4class{BayesSim}} object with draws resampled
#'   across models proportional to their weights. All \code{BayesSim}
#'   methods (e.g., \code{summary}, \code{elast}, \code{CV}) can be
#'   called on the result.
#'
#' @examples
#' \donttest{
#' # comp <- bayes_compare(Bertrand = fit1, Cournot = fit2)
#' # sim1 <- bayes_sim(fit1, ownerPost = ...)
#' # sim2 <- bayes_sim(fit2, ownerPost = ...)
#' #
#' # # BMA using bridge sampling posterior probabilities:
#' # avg_bridge <- bayes_avg(list(sim1, sim2),
#' #                         weights = bayes_weights(comp, "bridge"))
#' #
#' # # BMA using LOO stacking weights:
#' # avg_loo <- bayes_avg(list(sim1, sim2),
#' #                      weights = bayes_weights(comp, "stacking"))
#' #
#' # summary(avg_bridge)
#' }
#'
#' @export
bayes_avg <- function(sims, weights, n_draws = 1000L) {
  if (!is.list(sims) || length(sims) < 2L) {
    stop("'sims' must be a list of at least 2 BayesSim objects.")
  }
  if (!all(vapply(sims, is, logical(1), "BayesSim"))) {
    stop("Every element of 'sims' must be a BayesSim object.")
  }
  if (length(sims) != length(weights)) {
    stop("'sims' and 'weights' must have the same length.")
  }
  if (!is.numeric(weights) || any(!is.finite(weights)) || any(weights < 0) ||
      sum(weights) <= 0) {
    stop("'weights' must be finite and non-negative, with a positive sum.")
  }
  if (abs(sum(weights) - 1) > 1e-6) {
    stop("'weights' must sum to 1.")
  }
  weights <- weights / sum(weights)
  if (length(n_draws) != 1L || !is.finite(n_draws) || n_draws < 1 ||
      n_draws != as.integer(n_draws)) {
    stop("'n_draws' must be a positive integer.")
  }
  n_draws <- as.integer(n_draws)

  # Get common markets
  all_markets <- Reduce(intersect, lapply(sims, function(s) s@market_ids))
  if (length(all_markets) == 0L) {
    stop("No common markets across simulation objects.")
  }

  # Largest-remainder allocation preserves zero weights and returns exactly
  # n_draws, unlike independent rounding with a minimum of one per model.
  raw_counts <- weights * n_draws
  n_per_model <- floor(raw_counts)
  remainder <- n_draws - sum(n_per_model)
  if (remainder > 0L) {
    take_extra <- order(raw_counts - n_per_model, decreasing = TRUE)[
      seq_len(remainder)
    ]
    n_per_model[take_extra] <- n_per_model[take_extra] + 1L
  }
  n_per_model <- as.integer(n_per_model)

  # Draw indices are selected once per source model and reused in every
  # market, preserving the joint posterior draw across markets.
  selected <- vector("list", length(sims))
  for (i in seq_along(sims)) {
    if (n_per_model[i] == 0L) next
    common_valid <- Reduce(intersect, lapply(all_markets, function(mk) {
      which(!vapply(sims[[i]]@sims[[mk]], is.null, logical(1)))
    }))
    if (length(common_valid) == 0L) {
      stop("Simulation model ", i,
           " has no draw that is valid in every common market.")
    }
    selected[[i]] <- sample(common_valid, n_per_model[i], replace = TRUE)
  }
  source_order <- rep(seq_along(sims), n_per_model)
  source_position <- ave(source_order, source_order, FUN = seq_along)
  permutation <- sample.int(n_draws)
  source_order <- source_order[permutation]
  source_position <- source_position[permutation]

  pooled_sims     <- list()
  pooled_ppre     <- list()
  pooled_ppost    <- list()
  pooled_pchg     <- list()
  pooled_spre     <- list()
  pooled_spost    <- list()
  pooled_labels   <- list()
  pooled_n_valid  <- integer(length(all_markets))
  names(pooled_n_valid) <- all_markets

  for (mk in all_markets) {
    mk_sims  <- list()
    mk_ppre  <- NULL
    mk_ppost <- NULL
    mk_spre  <- NULL
    mk_spost <- NULL

    reference_labels <- sims[[1L]]@product_labels[[mk]]
    for (i in seq_along(sims)) {
      if (!identical(sims[[i]]@product_labels[[mk]], reference_labels)) {
        stop("Product labels or ordering differ across models in market '",
             mk, "'.")
      }
    }

    for (row in seq_len(n_draws)) {
      i <- source_order[row]
      j <- selected[[i]][source_position[row]]
      s <- sims[[i]]
      mk_sims[[row]] <- s@sims[[mk]][[j]]
      mk_ppre <- rbind(mk_ppre, s@price_pre[[mk]][j, , drop = FALSE])
      mk_ppost <- rbind(mk_ppost, s@price_post[[mk]][j, , drop = FALSE])
      mk_spre <- rbind(mk_spre, s@shares_pre[[mk]][j, , drop = FALSE])
      mk_spost <- rbind(mk_spost, s@shares_post[[mk]][j, , drop = FALSE])
    }

    pooled_sims[[mk]]    <- mk_sims
    pooled_ppre[[mk]]    <- mk_ppre
    pooled_ppost[[mk]]   <- mk_ppost
    pooled_pchg[[mk]]    <- mk_ppost / mk_ppre - 1
    pooled_spre[[mk]]    <- mk_spre
    pooled_spost[[mk]]   <- mk_spost
    pooled_labels[[mk]]  <- reference_labels
    pooled_n_valid[mk]   <- n_draws
  }

  new("BayesSim",
    sims           = pooled_sims,
    price_pre      = pooled_ppre,
    price_post     = pooled_ppost,
    price_change   = pooled_pchg,
    shares_pre     = pooled_spre,
    shares_post    = pooled_spost,
    market_ids     = all_markets,
    product_labels = pooled_labels,
    n_draws        = n_draws,
    n_valid        = pooled_n_valid,
    failures       = stats::setNames(
      replicate(length(all_markets), rep(NA_character_, n_draws),
                simplify = FALSE),
      all_markets
    ),
    call           = match.call()
  )
}
