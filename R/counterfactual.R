# Non-calibration entry point: run counterfactuals directly from supplied
# demand/conduct parameters, using antitrust's refactor-branch lifecycle
# (counterfactual()/simulate()) as the engine. No Stan fit is required.
#
#' @include params.R
NULL

## Default parallel worker count for bayes_counterfactual()'s per-draw loop:
## all available cores (parallel::detectCores()), capped at 4 to avoid
## surprising resource use on a shared machine, or 1 on Windows, where
## parallel::mclapply() cannot fork and would otherwise silently ignore
## mc.cores > 1 and run single-threaded anyway.
.default_mc_cores <- function() {
  if (identical(.Platform$OS.type, "windows")) return(1L)
  cores <- tryCatch(parallel::detectCores(logical = TRUE),
                    error = function(e) NA_integer_)
  if (is.na(cores) || cores < 1L) cores <- 1L
  min(4L, cores)
}

#' Run a counterfactual from directly-supplied Bayesian parameters
#'
#' Builds the same internal \code{BayesConductLogit}/\code{BayesConductCES}
#' realization \code{\link{bayes_sim}} builds from posterior draws, but from
#' hand-specified parameters (a point estimate or a matrix of draws). Exact
#' no-cutoff conduct variants are routed to antitrust's registered model
#' classes; cutoff blends and the legacy multiproduct Logit bargaining equation
#' retain their Bayesian conduct class. Each realization is wrapped as an
#' \code{antitrust::AntitrustFit}, and the counterfactual/merger solve is
#' delegated to \code{antitrust::simulate()}. No Bayesian estimation
#' (\code{\link{bayes_est}}) is required.
#'
#' @param params A \code{\link{bayes_params}} object.
#' @param market A \code{\link{bayes_market}} object.
#' @param ownerPost Post-merger ownership (a vector of firm IDs, or an
#'   ownership matrix), or an \pkg{antitrust} \code{Counterfactual} object
#'   built with \code{\link[antitrust]{counterfactual}}/
#'   \code{\link[antitrust]{add_step}} for cost shocks, entry, exit,
#'   quality changes, or multi-step paths.
#' @param mcDelta,subset As in \code{\link{bayes_sim}}; ignored when
#'   \code{ownerPost} is a \code{Counterfactual}.
#' @param n_draws Integer number of draws to use when \code{params} holds
#'   more than one draw. Default \code{NULL} uses all draws (or, for large
#'   posteriors, \code{min(200, n_draws_available)}).
#' @param mc.cores Integer number of parallel workers used to solve draws via
#'   \code{parallel::mclapply()}. Default is all available cores capped at
#'   4, or \code{1} on Windows (where \code{mclapply()} cannot fork and would
#'   otherwise run single-threaded regardless of this argument). Ignored
#'   when there is only one draw to solve.
#' @param ... Additional arguments passed to \code{antitrust::simulate()}.
#'
#' @return A single \pkg{antitrust} result object (or \code{CounterfactualPath})
#'   when \code{params} holds one draw; a \code{\linkS4class{BayesSim}} when
#'   it holds more than one.
#' @seealso \code{\link{bayes_sim}}, \code{\link{bayes_params}},
#'   \code{\link{bayes_market}}, \code{\link[antitrust]{counterfactual}},
#'   \code{\link[antitrust]{simulate}}
#'
#' @examples
#' # PNB (single market, single year): simulate the merger from a
#' # hand-specified demand/conduct parameterization, with no MCMC estimation.
#' data(pnb)
#'
#' market <- bayes_market(
#'   prices = pnb$prices, quantities = pnb$quantities, margins = pnb$margins,
#'   group = pnb$group, outside_price = 1, market_side = "output",
#'   insideSize = attr(pnb, "insideSize"), labels = pnb$bank
#' )
#' params <- bayes_params(alpha = 1, s0 = 0.015, supply_model = 1L)
#'
#' owner_post <- pnb$group
#' owner_post[owner_post == attr(pnb, "target")] <- attr(pnb, "acquirer")
#'
#' result <- bayes_counterfactual(params, market, ownerPost = owner_post)
#' merging <- pnb$group %in% c(attr(pnb, "acquirer"), attr(pnb, "target"))
#' data.frame(bank = pnb$bank[merging],
#'            pricePre = result@pricePre[merging],
#'            pricePost = result@pricePost[merging])
#' @export
bayes_counterfactual <- function(params, market, ownerPost,
                                  mcDelta = NULL, subset = NULL,
                                  n_draws = NULL,
                                  mc.cores = .default_mc_cores(), ...) {
  if (!inherits(params, "bayes_params")) {
    stop("'params' must be a bayes_params object; see bayes_params().")
  }
  if (!inherits(market, "bayes_market")) {
    stop("'market' must be a bayes_market object; see bayes_market().")
  }
  if (length(mc.cores) != 1L || !is.finite(mc.cores) || mc.cores < 1 ||
      mc.cores != as.integer(mc.cores)) {
    stop("'mc.cores' must be a positive integer.")
  }
  mc.cores <- as.integer(mc.cores)

  ownerPre <- market$group
  n_total <- params$n_draws
  if (is.null(n_draws)) {
    n_draws <- if (n_total > 1L) min(200L, n_total) else 1L
  }
  if (length(n_draws) != 1L || !is.finite(n_draws) || n_draws < 1 ||
      n_draws != as.integer(n_draws)) {
    stop("'n_draws' must be a positive integer.")
  }
  n_draws <- min(as.integer(n_draws), n_total)
  draw_idx <- if (n_total == 1L) 1L else sample(seq_len(n_total), n_draws)

  is_cf <- methods::is(ownerPost, "Counterfactual")

  run_one <- function(i) {
    if (is_cf) {
      ## The Counterfactual's own step(s) carry the ownership/cost/entry/exit
      ## change; build the object at the pre-merger baseline and let
      ## antitrust::simulate() apply the step(s).
      obj <- build_conduct_object(params, i, market, ownerPre,
                                  ownerPost = ownerPre)
      fit <- as_antitrust_fit(obj)
      return(antitrust::simulate(fit, ownerPost, ...))
    }
    ## No Counterfactual: bake ownerPost into construction and let
    ## as_antitrust_fit() fully solve pre/post prices directly (the same
    ## single solve bayes_sim() has always performed for this case).
    obj <- build_conduct_object(params, i, market, ownerPre, ownerPost,
                                mcDelta = mcDelta, subset = subset)
    as_antitrust_fit(obj, ...)@model
  }

  if (length(draw_idx) == 1L) {
    return(run_one(draw_idx[[1]]))
  }

  ## Each draw is solved independently (see build_conduct_object()), so the
  ## per-draw loop parallelizes with no shared state. mclapply() forks
  ## worker processes; a failure inside a forked call is captured here and
  ## turned into a plain list rather than left to propagate as a try-error,
  ## so the caller's error message matches the single-core path exactly.
  run_one_safe <- function(i) {
    tryCatch(
      list(result = run_one(i), error = NA_character_),
      error = function(e) list(result = NULL, error = conditionMessage(e))
    )
  }
  raw <- parallel::mclapply(draw_idx, run_one_safe, mc.cores = mc.cores)

  ## A crashed fork (segfault, OOM kill) surfaces from mclapply() as
  ## something other than the list(result=, error=) shape run_one_safe()
  ## always returns on a normal condition; treat that as a failed draw too.
  malformed <- !vapply(raw, function(x) {
    is.list(x) && identical(names(x), c("result", "error"))
  }, logical(1))
  raw[malformed] <- list(list(result = NULL,
                              error = "worker process did not return a result"))

  results <- lapply(raw, `[[`, "result")
  failures <- vapply(raw, `[[`, character(1), "error")
  .assemble_bayes_sim(results, failures, market, params, match.call())
}

## Re-exports of antitrust's Counterfactual grammar, so callers can build
## multi-step / entry / exit / quality counterfactuals without attaching
## antitrust themselves. Each is a bound name in this package's namespace
## (not just an imported symbol), so @export produces a valid NAMESPACE entry.

#' @name counterfactual
#' @rdname bayes_counterfactual
#' @export
counterfactual <- antitrust::counterfactual

#' @name add_step
#' @rdname bayes_counterfactual
#' @export
add_step <- antitrust::add_step

#' @name combine_counterfactuals
#' @rdname bayes_counterfactual
#' @export
combine_counterfactuals <- antitrust::combine_counterfactuals

#' @name entrant
#' @rdname bayes_counterfactual
#' @export
entrant <- antitrust::entrant

#' @name final_result
#' @rdname bayes_counterfactual
#' @export
final_result <- antitrust::final_result

#' @name result_at
#' @rdname bayes_counterfactual
#' @export
result_at <- antitrust::result_at
