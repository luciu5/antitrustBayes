# One-market conditional-k computations use ordinary BayesFit components.
# No alternate likelihood, conduct class, or posterior policy state is introduced.

.conditional_k_prior <- function(k_support, k_prior) {
  labels <- as.character(k_support)
  if (is.null(k_prior)) k_prior <- rep(1 / length(labels), length(labels))
  if (!is.numeric(k_prior) || length(k_prior) != length(labels) ||
      any(!is.finite(k_prior)) || any(k_prior < 0) || abs(sum(k_prior) - 1) > 1e-8)
    stop("'k_prior' must give finite nonnegative probabilities summing to one.", call. = FALSE)
  if (!is.null(names(k_prior)) && !identical(names(k_prior), labels))
    stop("Named 'k_prior' must match k_support in order.", call. = FALSE)
  stats::setNames(as.numeric(k_prior), labels)
}

#' Prepare Separate Conditional-k Fits for One Market
#'
#' Constructs ordinary \code{\link{bayes_est}} calls for each ordered-k state
#' in a single-market fringe or Stackelberg hypothesis. This supports either
#' Bertrand or Cournot supply. Every component retains the same observations,
#' continuous priors and likelihood, with singleton \code{k_support} and
#' \code{k_prior = 1}. The original k prior is restored when combining evidence.
#'
#' @param ... Named arguments to \code{bayes_est}, including \code{market},
#'   \code{group}, \code{structure}, and \code{supply_model}.
#' @param k_support Explicit ordered-k support, feasible in the selected market.
#' @param k_prior Original probabilities in support order; defaults to uniform.
#' @param component_stan_options Optional list, named by k in support order,
#'   of per-component \code{stan_options} overrides (e.g., seeds/output paths).
#' @return A named list of argument lists, with \code{k_support} and
#'   \code{k_prior} attributes. No estimation occurs.
#' @details Multiple periods are supported, with market-scoped k. More than
#'   one market is rejected: shared parameters prevent separate per-market
#'   marginal likelihoods from being multiplied to recover the joint family.
#' @export
bayes_conditional_k_arguments <- function(..., k_support, k_prior = NULL,
                                          component_stan_options = NULL) {
  args <- list(...)
  if (is.null(names(args)) || any(!nzchar(names(args))) || anyDuplicated(names(args)))
    stop("Supply unique named bayes_est arguments.", call. = FALSE)
  if (!all(c("market", "group", "structure", "supply_model") %in% names(args)) ||
      !args$structure %in% c("fringe", "stackelberg") ||
      !args$supply_model %in% c("bertrand", "cournot"))
    stop("Supply a Bertrand/Cournot fringe or Stackelberg hypothesis.", call. = FALSE)
  if (!length(args$market) || anyNA(args$market) || length(unique(args$market)) != 1L)
    stop("Conditional-k fitting requires exactly one market.", call. = FALSE)
  if (length(args$group) != length(args$market) || anyNA(args$group))
    stop("Firm identities must match the observations.", call. = FALSE)
  if (!is.null(args$k_scope) && !identical(args$k_scope, "market"))
    stop("Conditional-k fitting requires market-scoped k.", call. = FALSE)
  if (!is.null(args$partition_prior) && !identical(args$partition_prior, "ordered_k"))
    stop("Conditional-k fitting requires partition_prior = 'ordered_k'.", call. = FALSE)
  support <- .stack_validate_k_support(k_support, args$structure, length(unique(args$group)))
  if (is.null(support)) stop("Supply explicit k_support.", call. = FALSE)
  prior <- .conditional_k_prior(support, k_prior)
  labels <- as.character(support)
  if (!is.null(component_stan_options) &&
      (!is.list(component_stan_options) || !identical(names(component_stan_options), labels) ||
       any(!vapply(component_stan_options, is.list, logical(1)))))
    stop("component_stan_options must contain one named list per k in support order.", call. = FALSE)
  out <- lapply(seq_along(support), function(i) {
    component <- args
    component$partition_prior <- "ordered_k"
    component$k_support <- as.integer(support[i]); component$k_prior <- 1
    if (!is.null(component_stan_options)) component$stan_options <- utils::modifyList(
      component$stan_options %||% list(), component_stan_options[[i]])
    component
  })
  names(out) <- labels
  attr(out, "k_support") <- as.integer(support); attr(out, "k_prior") <- prior
  out
}

#' Estimate Separate Conditional-k Components
#'
#' @inheritParams bayes_conditional_k_arguments
#' @return A list with named \code{fits}, original \code{k_support}, and
#'   \code{k_prior}. Each successful component is an ordinary \code{BayesFit}.
#'   Failed components retain their error condition and do not cancel later fits.
#' @details Fits execute sequentially; independent components may instead be
#'   scheduled using \code{bayes_conditional_k_arguments}. HMC budgets and seeds
#'   remain explicit caller inputs. Review every component's diagnostics before
#'   combining evidence; missing components cannot be dropped.
#' @export
bayes_est_conditional_k <- function(..., k_support, k_prior = NULL,
                                   component_stan_options = NULL) {
  calls <- bayes_conditional_k_arguments(..., k_support = k_support,
    k_prior = k_prior, component_stan_options = component_stan_options)
  fits <- lapply(calls, function(args) tryCatch(do.call(bayes_est, args), error = identity))
  list(fits = fits, k_support = attr(calls, "k_support"), k_prior = attr(calls, "k_prior"))
}

#' Combine Conditional-k Marginal Likelihoods
#'
#' Restores the original k prior exactly once: \eqn{Z = \sum_k \pi_k Z_k}.
#' Posterior k weights are \eqn{\pi_k Z_k / Z}. These weights can mix posterior
#' predictions or policy draws from the same conditional components.
#' @param logml Named finite log marginal likelihoods, one per original k state,
#'   from normalized singleton-prior fits (e.g., \code{\link{bayes_compare}}).
#' @param k_support Complete original support, in order, supplied explicitly so
#'   a missing component cannot redefine the hypothesis.
#' @param k_prior Original k probabilities; defaults to uniform.
#' @param repetitions Optional finite matrix with repetitions in rows and k
#'   states in columns, named in support order. The reported SD is bridge repeat
#'   variability, not a certificate of evidence accuracy.
#' @return A list with \code{logml}, \code{posterior_k}, original \code{k_prior},
#'   and combined repetitions/\code{logml_sd} when supplied. This represents one
#'   scientific conduct hypothesis, not additional top-level model candidates.
#' @export
bayes_combine_conditional_k <- function(logml, k_support,
                                       k_prior = NULL, repetitions = NULL) {
  if (!is.numeric(k_support) || !length(k_support) || any(!is.finite(k_support)) ||
      any(k_support != floor(k_support)) || anyDuplicated(k_support) || any(k_support < 0))
    stop("Supply the complete original integer k_support.", call. = FALSE)
  labels <- as.character(k_support)
  if (!is.numeric(logml) || !identical(names(logml), labels) ||
      length(logml) != length(labels) || any(!is.finite(logml)))
    stop("Every original k state requires finite named conditional evidence.", call. = FALSE)
  prior <- .conditional_k_prior(k_support, k_prior)
  log_sum <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }
  log_weight <- logml + log(prior)
  total <- log_sum(log_weight)
  out <- list(logml = total, k_prior = prior,
    posterior_k = data.frame(k = as.integer(k_support), probability = exp(log_weight - total)),
    component_logml = logml, repetitions = NULL, logml_sd = NA_real_)
  if (!is.null(repetitions)) {
    if (!is.matrix(repetitions) || !is.numeric(repetitions) || nrow(repetitions) < 2L ||
        !identical(colnames(repetitions), labels) || any(!is.finite(repetitions)))
      stop("repetitions must have at least two finite rows and all original k columns.", call. = FALSE)
    out$repetitions <- apply(repetitions, 1L, function(x) log_sum(x + log(prior)))
    out$logml_sd <- stats::sd(out$repetitions)
  }
  out
}
