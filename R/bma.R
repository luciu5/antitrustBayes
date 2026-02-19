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
#'
#' @return A data frame with columns for model name, supply model, LOOIC,
#'   WAIC, stacking weights, and (if \code{bridge = TRUE}) log marginal
#'   likelihoods, Bayes factors (vs. best model), and posterior model
#'   probabilities (assuming equal prior model probabilities).
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
#' (for LOO/WAIC).
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
bayes_compare <- function(..., bridge = TRUE, bridge_maxiter = 5000L) {
  args <- list(...)

  # Separate BayesLogit objects from named arguments
  fits <- args[vapply(args, is, logical(1), "BayesLogit")]

  # Handle named and unnamed arguments
  nms <- names(fits)
  if (is.null(nms)) {
    nms <- paste0("Model_", seq_along(fits))
  } else {
    empty <- nms == ""
    nms[empty] <- paste0("Model_", which(empty))
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

    max_logml <- max(logml, na.rm = TRUE)
    diff_logml <- logml - max_logml
    bf_vs_best <- exp(diff_logml)

    # Posterior model probabilities (equal prior weights)
    unnorm <- exp(diff_logml)
    unnorm[is.na(unnorm)] <- 0
    post_prob <- unnorm / sum(unnorm)

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

  if (method == "bridge" && "post_prob" %in% names(comp)) {
    w <- comp$post_prob
  } else {
    if (method == "bridge" && !"post_prob" %in% names(comp)) {
      message("Bridge sampling results not available; using stacking weights.")
    }
    w <- comp$weight_stacking
  }

  names(w) <- comp$model
  w
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
  if (length(sims) != length(weights)) {
    stop("'sims' and 'weights' must have the same length.")
  }
  if (abs(sum(weights) - 1) > 1e-6) {
    stop("'weights' must sum to 1.")
  }

  n_draws <- as.integer(n_draws)

  # Get common markets
  all_markets <- Reduce(intersect, lapply(sims, function(s) s@market_ids))
  if (length(all_markets) == 0L) {
    stop("No common markets across simulation objects.")
  }

  # Number of draws to take from each model
  n_per_model <- pmax(1L, round(weights * n_draws))

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

    for (i in seq_along(sims)) {
      s <- sims[[i]]
      sims_mk <- s@sims[[mk]]
      valid <- which(!vapply(sims_mk, is.null, logical(1)))
      if (length(valid) == 0L) next

      # Resample from this model's valid draws
      n_take <- min(n_per_model[i], length(valid))
      sel <- sample(valid, n_take, replace = TRUE)

      for (j in sel) {
        mk_sims <- c(mk_sims, list(sims_mk[[j]]))
      }

      mk_ppre  <- rbind(mk_ppre,  s@price_pre[[mk]][sel, , drop = FALSE])
      mk_ppost <- rbind(mk_ppost, s@price_post[[mk]][sel, , drop = FALSE])
      mk_spre  <- rbind(mk_spre,  s@shares_pre[[mk]][sel, , drop = FALSE])
      mk_spost <- rbind(mk_spost, s@shares_post[[mk]][sel, , drop = FALSE])
    }

    n_total <- length(mk_sims)
    pooled_sims[[mk]]    <- mk_sims
    pooled_ppre[[mk]]    <- mk_ppre
    pooled_ppost[[mk]]   <- mk_ppost
    pooled_pchg[[mk]]    <- if (n_total > 0L) mk_ppost / mk_ppre - 1 else NULL
    pooled_spre[[mk]]    <- mk_spre
    pooled_spost[[mk]]   <- mk_spost
    pooled_labels[[mk]]  <- sims[[1]]@product_labels[[mk]]
    pooled_n_valid[mk]   <- n_total
  }

  n_actual <- max(vapply(pooled_sims, length, integer(1)))

  new("BayesSim",
    sims           = pooled_sims,
    price_pre      = pooled_ppre,
    price_post     = pooled_ppost,
    price_change   = pooled_pchg,
    shares_pre     = pooled_spre,
    shares_post    = pooled_spost,
    market_ids     = all_markets,
    product_labels = pooled_labels,
    n_draws        = n_actual,
    n_valid        = pooled_n_valid,
    call           = match.call()
  )
}
