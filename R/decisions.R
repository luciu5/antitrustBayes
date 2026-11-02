# Decision-focused summaries of posterior merger effects.

#' Posterior Merger Decision
#'
#' Implements the Bayes action for a binary harm threshold under asymmetric
#' false-positive and false-negative losses.
#'
#' @param draws Numeric posterior draws of a merger harm measure, with larger
#'   values indicating greater harm.
#' @param threshold Harm threshold.
#' @param false_positive_cost Cost of challenging a non-harmful merger.
#' @param false_negative_cost Cost of clearing a harmful merger.
#' @return A list containing the posterior harm probability, decision cutoff,
#'   and binary challenge decision.
#' @export
merger_decision <- function(draws, threshold = 0,
                            false_positive_cost = 1,
                            false_negative_cost = 1) {
  if (!is.numeric(draws) || !length(draws) || any(!is.finite(draws))) {
    stop("'draws' must be a non-empty finite numeric vector.")
  }
  costs <- c(false_positive_cost, false_negative_cost)
  if (any(!is.finite(costs)) || any(costs <= 0)) {
    stop("Decision costs must be positive finite scalars.")
  }
  probability <- mean(draws > threshold)
  probability_cutoff <- false_positive_cost /
    (false_positive_cost + false_negative_cost)
  structure(list(
    probability_harm = probability,
    probability_cutoff = probability_cutoff,
    challenge = probability > probability_cutoff,
    threshold = threshold,
    false_positive_cost = false_positive_cost,
    false_negative_cost = false_negative_cost
  ), class = "merger_decision")
}

#' Merger Decision Loss and Regret
#'
#' @param decision Logical vector; \code{TRUE} means challenge.
#' @param truth Numeric vector of realized merger effects.
#' @param threshold Harm threshold.
#' @param false_positive_cost,false_negative_cost Positive decision costs.
#' @param oracle_decision Optional logical oracle decisions. Defaults to the
#'   no-error decision implied by \code{truth > threshold}.
#' @return A data frame containing loss, oracle loss, and regret.
#' @export
decision_regret <- function(decision, truth, threshold = 0,
                            false_positive_cost = 1,
                            false_negative_cost = 1,
                            oracle_decision = truth > threshold) {
  decision <- as.logical(decision)
  oracle_decision <- as.logical(oracle_decision)
  if (!is.numeric(truth) || any(!is.finite(truth)) ||
      length(decision) != length(truth) ||
      length(oracle_decision) != length(truth) ||
      anyNA(decision) || anyNA(oracle_decision)) {
    stop("Decisions and finite truth values must have conformable lengths.")
  }
  harmful <- truth > threshold
  loss_for <- function(action) {
    false_positive_cost * (action & !harmful) +
      false_negative_cost * (!action & harmful)
  }
  loss <- loss_for(decision)
  oracle_loss <- loss_for(oracle_decision)
  data.frame(
    harmful = harmful,
    decision = decision,
    loss = loss,
    oracle_loss = oracle_loss,
    regret = loss - oracle_loss
  )
}
