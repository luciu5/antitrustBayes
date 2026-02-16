#' Compare Bayesian Models Using Information Criteria
#'
#' Computes WAIC and LOO-IC for one or more fitted
#' \code{\linkS4class{BayesLogit}} objects and returns a comparison table.
#'
#' @param ... One or more \code{BayesLogit} objects, optionally named.
#'
#' @return A data frame with columns for model name, WAIC, LOO-IC,
#'   effective number of parameters, and model weights.
#'
#' @details
#' Requires that the Stan model includes a \code{log_lik} generated quantity.
#' Uses \code{\link[loo]{loo}} and \code{\link[loo]{waic}} for computation
#' and \code{\link[loo]{loo_model_weights}} (stacking) for weights.
#'
#' @examples
#' \donttest{
#' # fit1 <- bayes_est(..., supply_model = "bertrand")
#' # fit2 <- bayes_est(..., supply_model = "cournot")
#' # bayes_compare(Bertrand = fit1, Cournot = fit2)
#' }
#'
#' @export
bayes_compare <- function(...) {
  fits <- list(...)

  # Handle named and unnamed arguments
  nms <- names(fits)
  if (is.null(nms)) {
    nms <- paste0("Model_", seq_along(fits))
  } else {
    empty <- nms == ""
    nms[empty] <- paste0("Model_", which(empty))
  }

  supply_names <- c("Bertrand", "Monopolistic", "Cournot", "Competitive")

  loo_list  <- list()
  waic_list <- list()

  for (i in seq_along(fits)) {
    fit <- fits[[i]]
    if (!is(fit, "BayesLogit")) {
      stop("Argument ", i, " is not a BayesLogit object.")
    }

    log_lik <- loo::extract_log_lik(fit@stanfit, parameter_name = "log_lik",
                                     merge_chains = FALSE)
    r_eff <- loo::relative_eff(exp(log_lik))

    loo_list[[i]]  <- loo::loo(log_lik, r_eff = r_eff)
    waic_list[[i]] <- loo::waic(log_lik)
  }

  names(loo_list) <- nms
  comp <- loo::loo_compare(loo_list)

  # Model weights via stacking
  weights <- loo::loo_model_weights(loo_list, method = "stacking")

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
    weight        = as.numeric(weights),
    stringsAsFactors = FALSE
  )
  rownames(result) <- NULL

  attr(result, "loo_compare") <- comp
  attr(result, "loo_list") <- loo_list
  result
}

#' Bayesian Model Averaging of Merger Simulations
#'
#' Combines merger simulation results from multiple models using
#' model weights (e.g., from \code{\link{bayes_compare}}).
#'
#' @param sims A list of \code{\linkS4class{BayesSim}} objects.
#' @param weights Numeric vector of model weights (must sum to 1).
#'   Typically from \code{\link{bayes_compare}}.
#'
#' @return A data frame of model-averaged price change summaries with
#'   columns for market, product, weighted mean, weighted SD, and quantiles.
#'
#' @examples
#' \donttest{
#' # comp <- bayes_compare(Bertrand = fit1, Cournot = fit2)
#' # sim1 <- bayes_sim(fit1, ownerPost = ...)
#' # sim2 <- bayes_sim(fit2, ownerPost = ...)
#' # bayes_avg(list(sim1, sim2), weights = comp$weight)
#' }
#'
#' @export
bayes_avg <- function(sims, weights) {
  if (!is.list(sims) || length(sims) < 2L) {
    stop("'sims' must be a list of at least 2 BayesSim objects.")
  }
  if (length(sims) != length(weights)) {
    stop("'sims' and 'weights' must have the same length.")
  }
  if (abs(sum(weights) - 1) > 1e-6) {
    stop("'weights' must sum to 1.")
  }

  # Get common markets
  all_markets <- Reduce(intersect, lapply(sims, function(s) s@market_ids))
  if (length(all_markets) == 0L) {
    stop("No common markets across simulation objects.")
  }

  results <- list()

  for (mk in all_markets) {
    # Pool weighted draws
    pooled_pchg <- NULL

    for (i in seq_along(sims)) {
      s <- sims[[i]]
      pchg <- s@price_change[[mk]]
      valid_rows <- which(rowSums(is.na(pchg)) == 0)
      if (length(valid_rows) == 0L) next

      pchg_valid <- pchg[valid_rows, , drop = FALSE]

      # Resample proportional to model weight
      n_resample <- max(1L, round(weights[i] * 1000))
      resample_idx <- sample(nrow(pchg_valid), n_resample, replace = TRUE)
      pooled_pchg <- rbind(pooled_pchg, pchg_valid[resample_idx, , drop = FALSE])
    }

    if (is.null(pooled_pchg) || nrow(pooled_pchg) == 0L) next

    labs <- sims[[1]]@product_labels[[mk]]

    tbl <- data.frame(
      market  = mk,
      product = labs,
      mean    = colMeans(pooled_pchg),
      sd      = apply(pooled_pchg, 2, sd),
      q025    = apply(pooled_pchg, 2, quantile, probs = 0.025),
      median  = apply(pooled_pchg, 2, quantile, probs = 0.5),
      q975    = apply(pooled_pchg, 2, quantile, probs = 0.975),
      prob_positive = colMeans(pooled_pchg > 0),
      stringsAsFactors = FALSE
    )
    rownames(tbl) <- NULL
    results[[mk]] <- tbl
  }

  do.call(rbind, results)
}
