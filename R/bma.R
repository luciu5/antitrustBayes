.bayes_is_collapsed_quality <- function(fit) {
  metadata <- fit@metadata
  isTRUE(metadata$bayes_is_collapsed_quality) ||
    identical(metadata$quality_parameterization, "collapsed") ||
    isTRUE(fit@stan_data$quality_parameterization == 2L)
}

validate_comparison_contract <- function(fits) {
  reference <- fits[[1L]]
  columns <- c(
    "prices", "quantities", "margins", "market", "year", "group",
    "loo_cluster", "margin_observed", "outside_price"
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
    common_outcomes <- c(
      "shareIn", "clr_shareIn", "marginNormInv", "log_marginNormInv",
      "margin_obs", "market_year_idx", "cluster", "likelihood_mode"
    )
    if (!all(common_outcomes %in% names(candidate@stan_data)) ||
        !all(common_outcomes %in% names(reference@stan_data)) ||
        !isTRUE(all.equal(
          candidate@stan_data[common_outcomes],
          reference@stan_data[common_outcomes],
          check.attributes = FALSE
        ))) {
      stop("Models do not share the current quantity-space likelihood contract.")
    }
    shared_prior_names <- setdiff(
      names(reference@metadata$resolved_priors),
      grep("^bargaining_", names(reference@metadata$resolved_priors), value = TRUE)
    )
    if (!(.bayes_uses_cutoff_partition(reference) &&
          .bayes_uses_cutoff_partition(candidate))) {
      shared_prior_names <- setdiff(shared_prior_names,
        c("prior_cutoff_alpha", "prior_cutoff_beta"))
    }
    if (!(.bayes_is_smooth_fringe(reference) &&
          .bayes_is_smooth_fringe(candidate))) {
      shared_prior_names <- setdiff(shared_prior_names, "prior_sigma_meanval_fringe")
    }
    if (!identical(candidate@metadata$prior_manifest,
                   reference@metadata$prior_manifest) ||
        !isTRUE(all.equal(
          candidate@metadata$resolved_priors[shared_prior_names],
          reference@metadata$resolved_priors[shared_prior_names]
        ))) {
      stop("Models must use the same prior manifest on shared parameters.")
    }
  }
  invisible(TRUE)
}

.bayes_bridge_kernel_names <- c(
  "ces_stackelberg", "ces_supply", "logit_stackelberg", "logit_supply"
)

.bayes_stan_source_text <- function(code) {
  if (!is.character(code) || !length(code) || anyNA(code) ||
      !nzchar(paste(code, collapse = ""))) {
    stop("Stan source code is unavailable.", call. = FALSE)
  }
  code <- paste(code, collapse = "\n")
  lines <- strsplit(code, "\r?\n", perl = TRUE)[[1L]]
  paste0(paste(lines, collapse = "\n"), "\n")
}

.bayes_cmdstan_source_contract <- function(metadata, source_code,
                                           require_hmc = TRUE) {
  if (!is.list(metadata)) {
    stop("Fit provenance metadata is unavailable.", call. = FALSE)
  }
  if (!identical(metadata$backend, "cmdstanr")) {
    stop("Evidence requires a fit produced by the CmdStanR backend.",
         call. = FALSE)
  }
  inference <- metadata$inference
  if (length(inference) != 1L || is.na(inference) ||
      !inference %in% c("hmc", "pathfinder")) {
    stop("Fit provenance does not identify a supported inference method.",
         call. = FALSE)
  }
  if (isTRUE(require_hmc) && !identical(inference, "hmc")) {
    stop("LOO and bridge evidence require HMC posterior draws.",
         call. = FALSE)
  }
  kernel <- metadata$stan_kernel
  if (length(kernel) != 1L || is.na(kernel) ||
      !kernel %in% .bayes_bridge_kernel_names) {
    stop("Fit provenance does not identify a supported Stan kernel.",
         call. = FALSE)
  }
  saved_hash <- metadata$stan_source_hash
  if (!is.character(saved_hash) || length(saved_hash) != 1L ||
      !grepl("^[0-9a-f]{64}$", saved_hash)) {
    stop("Fit provenance has no valid Stan source SHA-256 hash.",
         call. = FALSE)
  }
  source_file <- tryCatch(.cmdstan_kernel_file(kernel), error = function(e) {
    stop("Cannot locate the current Stan source for ", kernel, ": ",
         conditionMessage(e), call. = FALSE)
  })
  current_hash <- digest::digest(file = source_file, algo = "sha256")
  if (!identical(saved_hash, current_hash)) {
    stop("Saved and installed Stan source hashes differ for ", kernel,
         "; evidence is unavailable for this fit.", call. = FALSE)
  }
  current_code <- paste0(paste(readLines(source_file, warn = FALSE),
                               collapse = "\n"), "\n")
  fit_code <- .bayes_stan_source_text(source_code)
  current_code <- .bayes_stan_source_text(current_code)
  code_hash <- digest::digest(fit_code, algo = "sha256", serialize = FALSE)
  current_code_hash <- digest::digest(current_code, algo = "sha256",
                                      serialize = FALSE)
  if (!identical(code_hash, current_code_hash)) {
    stop("Saved CmdStan fit source does not match the current kernel source; ",
         "evidence is unavailable for this fit.", call. = FALSE)
  }
  list(
    kernel = kernel,
    inference = inference,
    source_file = source_file,
    saved_hash = saved_hash,
    current_hash = current_hash,
    code_hash = code_hash
  )
}

.bayes_cmdstan_evidence_contract <- function(fit, require_hmc = TRUE) {
  stanfit <- fit@stanfit
  if (!inherits(stanfit, "CmdStanMCMC")) {
    stop("Evidence requires a CmdStanR CmdStanMCMC fit.", call. = FALSE)
  }
  source_code <- tryCatch(stanfit$code(), error = function(e) NULL)
  .bayes_cmdstan_source_contract(
    fit@metadata, source_code, require_hmc = require_hmc
  )
}

.bayes_bridge_restore_context <- function(fit, contract = NULL) {
  if (is.null(contract)) contract <- .bayes_cmdstan_evidence_contract(fit)
  fields <- fit@metadata$stan_data_fields
  if (is.null(fields)) {
    fields <- names(.cmdstan_model_for_kernel(contract$kernel)$variables()$data)
  }
  if (!is.character(fields) || !length(fields) || anyNA(fields) ||
      anyDuplicated(fields) || !all(fields %in% names(fit@stan_data))) {
    stop("Saved fit has no valid Stan data field contract for bridge sampling.",
         call. = FALSE)
  }
  stan_data <- fit@stan_data[fields]
  saved_data_hash <- fit@metadata$stan_data_hash
  if (!is.null(saved_data_hash) &&
      !identical(saved_data_hash,
                 digest::digest(stan_data, algo = "sha256"))) {
    stop("Saved Stan data differ from the data used for fitting; bridge evidence is unavailable.",
         call. = FALSE)
  }
  data_file <- tempfile("antitrustBayes-bridge-data-", fileext = ".json")
  tryCatch(cmdstanr::write_stan_json(stan_data, file = data_file),
           error = function(e) {
             unlink(data_file)
             stop("Cannot reconstruct CmdStan bridge data: ",
                  conditionMessage(e), call. = FALSE)
           })
  ## CmdStanR's saved fit still points at its original temporary JSON. Use an
  ## independent copy so bridge evaluation cannot mutate the user's fit.
  samples <- unserialize(serialize(fit@stanfit, NULL))
  samples$runset$args$data_file <- data_file
  list(
    samples = samples,
    source = "cmdstanr_saved_fit_reconstructed_data",
    kernel = contract$kernel,
    saved_hash = contract$saved_hash,
    current_hash = contract$current_hash,
    migration = FALSE,
    cleanup = function() unlink(data_file)
  )
}

.bayes_cmdstan_log_lik_array <- function(stanfit, parameter_name,
                                         expected_units = NULL) {
  if (!inherits(stanfit, "CmdStanMCMC")) {
    stop("Predictive scoring requires a CmdStanR CmdStanMCMC fit.",
         call. = FALSE)
  }
  if (!is.character(parameter_name) || length(parameter_name) != 1L ||
      is.na(parameter_name) || !nzchar(parameter_name)) {
    stop("A single log-likelihood variable name is required.",
         call. = FALSE)
  }
  draws <- tryCatch(
    stanfit$draws(variables = parameter_name, format = "draws_array"),
    error = function(e) stop("Cannot read generated quantity '",
                             parameter_name, "': ", conditionMessage(e),
                             call. = FALSE)
  )
  log_lik <- as.array(draws)
  if (length(dim(log_lik)) != 3L || any(dim(log_lik) < 1L) ||
      !is.numeric(log_lik)) {
    stop("Generated quantity '", parameter_name,
         "' did not produce an iteration-by-chain-by-unit array.",
         call. = FALSE)
  }
  variables <- dimnames(log_lik)[[3L]]
  if (is.null(variables) || length(variables) != dim(log_lik)[[3L]]) {
    stop("Generated quantity '", parameter_name,
         "' has no complete unit names.", call. = FALSE)
  }
  if (length(variables) == 1L && identical(variables, parameter_name)) {
    unit_index <- 1L
  } else {
    prefix <- paste0(parameter_name, "[")
    indexed <- startsWith(variables, prefix) & endsWith(variables, "]")
    indexes <- rep(NA_integer_, length(variables))
    indexes[indexed] <- suppressWarnings(as.integer(sub(
      "\\]$", "", substring(variables[indexed], nchar(prefix) + 1L),
      perl = TRUE
    )))
    unit_index <- indexes
    if (anyNA(unit_index) || anyDuplicated(unit_index) ||
        !identical(sort(unit_index), seq_along(unit_index))) {
      stop("Generated quantity '", parameter_name,
           "' has missing, duplicated, or noncontiguous unit indices.",
           call. = FALSE)
    }
    order_units <- order(unit_index)
    log_lik <- log_lik[, , order_units, drop = FALSE]
    unit_index <- unit_index[order_units]
  }
  if (!is.null(expected_units) &&
      (length(expected_units) != 1L || !is.finite(expected_units) ||
       expected_units < 1L || expected_units != as.integer(expected_units) ||
       dim(log_lik)[[3L]] != as.integer(expected_units))) {
    stop("Generated quantity '", parameter_name,
         "' does not match the expected predictive-unit count.",
         call. = FALSE)
  }
  if (any(!is.finite(log_lik))) {
    stop("Generated quantity '", parameter_name,
         "' contains non-finite values.", call. = FALSE)
  }
  array(as.numeric(log_lik), dim = dim(log_lik),
        dimnames = dimnames(log_lik))
}

.bayes_cmdstan_log_lik_units <- function(fit, parameter_name) {
  field <- switch(parameter_name,
    log_lik = "N",
    log_lik_period = "N_market_year",
    log_lik_market = "N_market",
    log_lik_cluster = "N_cluster",
    NULL
  )
  if (is.null(field)) return(NULL)
  value <- fit@stan_data[[field]]
  if (is.numeric(value) && length(value) == 1L && is.finite(value) &&
      value >= 1L && value == as.integer(value)) as.integer(value) else NULL
}

#' Compare Bayesian Models or Exact Ordered-k Restrictions
#'
#' Computes LOO-IC, WAIC, and optionally log marginal likelihoods via bridge
#' sampling for one or more fitted \code{\linkS4class{BayesLogit}} objects.
#' For one exact ordered-k fit, \code{restrictions} instead compares restrictions
#' on k with the exact encompassing-prior identity. Returns a comparison table
#' appropriate to the selected route.
#'
#' @param ... One or more \code{BayesLogit} objects, optionally named. With
#'   \code{universe}, a single named list of source fits is also accepted.
#' @param restrictions Optional named restrictions on ordered-k states. When
#'   supplied, \code{...} must contain exactly one exact ordered-k
#'   \code{\linkS4class{BayesFit}}. See
#'   \code{\link{bayes_restriction_compare}} for the one- and multimarket forms.
#' @param reference Optional restriction name or row index used as the Bayes
#'   factor reference. Only available when \code{restrictions} is supplied.
#' @param restriction_draws Optional generated-quantity draws passed to
#'   \code{\link{bayes_restriction_compare}}. This is mainly useful for testing
#'   or recovering summaries from serialized draws.
#' @param bridge Logical. If \code{TRUE} (default), compute log marginal
#'   likelihoods via \code{\link[bridgesampling]{bridge_sampler}} for Bayes
#'   factor model comparison. If \code{FALSE}, only LOO/WAIC are computed.
#' @param bridge_maxiter Maximum iterations for bridge sampling.
#'   Default 5000.
#' @param bridge_repetitions Number of independently initialized bridge runs
#'   per model. The reported log marginal likelihood is their median.
#' @param bridge_method Bridge-sampling estimator passed to
#'   \code{\link[bridgesampling]{bridge_sampler}}: \code{"normal"} (default)
#'   or \code{"warp3"}. \code{"warp3"} symmetrizes the posterior before
#'   fitting the proposal and can reduce across-repetition log marginal
#'   likelihood variance for skewed posteriors, at roughly double the
#'   per-repetition cost.
#' @param loo_unit Natural predictive unit: observation, period, market, or
#'   cluster. Cluster uses the identifiers supplied to \code{bayes_est}.
#' @param prior_weights Optional numeric vector of prior model probabilities.
#'   Must be the same length as the number of models and sum to 1. If
#'   \code{NULL} (default), recorded \code{model_prior} values are used when
#'   all are available; otherwise equal probabilities are used.
#' @param universe Optional \code{bayes_model_universe} descriptor. When
#'   supplied, the named fits in \code{...} are its source fits and the
#'   descriptor's explicit scientific prior is used. This route requires
#'   bridge-sampling evidence.
#'
#' @return With multiple fitted models, a data frame containing model name,
#'   supply model, LOOIC, WAIC, stacking weights, and (if
#'   \code{bridge = TRUE}) log marginal likelihoods, Bayes factors, and
#'   posterior model probabilities. With \code{restrictions}, a
#'   \code{bayes_restriction_comparison} table containing prior and posterior
#'   restriction masses and encompassing-prior Bayes factors.
#'
#' @details
#' Four model comparison approaches are provided:
#' \describe{
#'   \item{LOO-IC}{Leave-one-out cross-validation via \pkg{loo}. Stacking
#'     weights reflect predictive performance.}
#'   \item{WAIC}{Widely Applicable Information Criterion.}
#'   \item{Bridge Sampling}{Computes the log marginal likelihood
#'     \eqn{p(y | M_k)} via \code{\link[bridgesampling]{bridge_sampler}}.
#'     Posterior model probabilities use \code{prior_weights}, or the
#'     model-specific probabilities recorded by \code{bayes_est}.
#'     Bayes factors are computed relative to the best model.}
#'   \item{Exact ordered-k restrictions}{Computes
#'     \eqn{BF_{A,F}=P(A\mid y)/P(A)} from the generated conditional
#'     probabilities in one marginalized-k fit. This route does not use bridge
#'     sampling or a Gelfand-Dey estimator.}
#' }
#'
#' Requires that the Stan model includes a \code{log_lik} generated quantity
#' (for LOO/WAIC). All fits must use the same ordered observations, market
#' side, outside prices, margin mask, and current quantity-space likelihood.
#' This permits Logit/CES comparison without comparing likelihoods defined on
#' different observed outcomes.
#' The current Stan specification retains prior density constants and support
#' normalizers directly. Earlier specification fits must be refitted. Evidence
#' refers to the transformed-outcome working likelihood. Bridge comparison
#' with an active, non-normalized HMT penalty is unsupported. Exact ordered-k
#' LOO clusters must contain whole partition units; observation LOO and period
#' LOO with a market-shared count are unsupported for fringe and Stackelberg.
#'
#' @examples
#' \donttest{
#' # fit1 <- bayes_est(..., supply_model = "bertrand")
#' # fit2 <- bayes_est(..., supply_model = "cournot")
#' # bayes_compare(Bertrand = fit1, Cournot = fit2)
#' # bayes_compare(Bertrand = fit1, Cournot = fit2, bridge = FALSE)
#' # bayes_compare(fit_stack, restrictions = list(k1 = 1L, k2 = 2L,
#' #                                                k3_to_5 = 3:5),
#' #               reference = "k2")
#' }
#'
#' @export
bayes_compare <- function(..., restrictions = NULL, reference = NULL,
                          restriction_draws = NULL,
                          bridge = TRUE, bridge_maxiter = 5000L,
                          bridge_repetitions = 3L,
                          bridge_method = c("normal", "warp3"),
                          loo_unit = c("cluster", "period", "market", "observation"),
                          prior_weights = NULL, universe = NULL) {
  bridge_method <- match.arg(bridge_method)
  args <- list(...)
  if (!is.null(universe)) {
    # Accept a named collection as well as named fits in ... for this route.
    if (length(args) == 1L && is.list(args[[1L]]) &&
        !methods::is(args[[1L]], "BayesFit")) args <- args[[1L]]
    if (!is.null(restrictions) || !is.null(reference) ||
        !is.null(restriction_draws)) {
      stop("'universe' cannot be combined with ordered-k 'restrictions', ",
           "'reference', or 'restriction_draws'.", call. = FALSE)
    }
    if (!is.null(prior_weights)) {
      stop("'prior_weights' belongs in the explicit model-universe descriptor ",
           "when 'universe' is supplied.", call. = FALSE)
    }
    loo_unit <- match.arg(loo_unit)
    return(.bayes_compare_universe(
      args, universe = universe, bridge = bridge,
      bridge_maxiter = bridge_maxiter,
      bridge_repetitions = bridge_repetitions, bridge_method = bridge_method,
      loo_unit = loo_unit
    ))
  }
  if (!is.null(restrictions)) {
    is_fit <- vapply(args, methods::is, logical(1), "BayesFit")
    if (length(args) != 1L || sum(is_fit) != 1L) {
      stop("Restriction comparison requires exactly one BayesFit object in '...'.",
           call. = FALSE)
    }
    if (!is.null(prior_weights)) {
      stop("'prior_weights' compares separate fitted models and is unavailable ",
           "for ordered-k restrictions; the fit's recorded k prior is used.",
           call. = FALSE)
    }
    return(bayes_restriction_compare(
      args[[which(is_fit)]], restrictions = restrictions,
      reference = reference, draws = restriction_draws
    ))
  }
  if (!is.null(reference) || !is.null(restriction_draws)) {
    stop("'reference' and 'restriction_draws' require 'restrictions'.",
         call. = FALSE)
  }
  loo_unit <- match.arg(loo_unit)
  if (length(bridge_repetitions) != 1L || !is.finite(bridge_repetitions) ||
      bridge_repetitions < 1 || bridge_repetitions != as.integer(bridge_repetitions)) {
    stop("'bridge_repetitions' must be a positive integer.")
  }
  bridge_repetitions <- as.integer(bridge_repetitions)

  # Separate BayesLogit objects from named arguments
  fits <- args[vapply(args, is, logical(1), "BayesLogit")]
  if (length(fits) == 0L) {
    stop("At least one BayesLogit fit is required.")
  }
  invisible(lapply(fits, validate_fit_spec))
  validate_comparison_contract(fits)
  collapsed_quality <- any(vapply(fits, .bayes_is_collapsed_quality, logical(1)))
  if (collapsed_quality && !isTRUE(bridge)) {
    stop("LOO/WAIC comparison is unavailable for collapsed quality fits; ",
         "use bridge = TRUE for normalized evidence.", call. = FALSE)
  }
  if (!collapsed_quality) .bayes_comparison_units(fits, loo_unit)
  if (bridge) .bayes_validate_bridge(fits)

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
    recorded <- vapply(fits, function(x) {
      value <- x@metadata$model_prior
      if (length(value) == 1L && is.numeric(value)) as.numeric(value) else NA_real_
    }, numeric(1))
    if (all(is.finite(recorded)) && sum(recorded) > 0) {
      prior_weights <- recorded / sum(recorded)
    } else {
      prior_weights <- rep(1 / length(fits), length(fits))
    }
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

  supply_names <- c(
    "Bertrand", "Monopolistic", "Cournot", "Auction", "Bargaining"
  )

  loo_list  <- vector("list", length(fits))
  waic_list <- vector("list", length(fits))
  bridge_list <- list()
  bridge_context_info <- vector("list", length(fits))
  parameter_name <- switch(loo_unit,
    observation = "log_lik",
    period = "log_lik_period",
    market = "log_lik_market",
    cluster = "log_lik_cluster"
  )
  loo_available <- !collapsed_quality
  loo_failure <- NULL

  for (i in seq_along(fits)) {
    fit <- fits[[i]]
    evidence_contract <- NULL
    contract_error <- NULL
    if (loo_available || bridge) {
      evidence_contract <- tryCatch(
        .bayes_cmdstan_evidence_contract(fit),
        error = function(e) {
          contract_error <<- conditionMessage(e)
          NULL
        }
      )
    }

    if (loo_available) {
      loo_error <- contract_error
      loo_result <- if (!is.null(loo_error)) NULL else tryCatch({
        log_lik <- .bayes_cmdstan_log_lik_array(
          fit@stanfit, parameter_name,
          expected_units = .bayes_cmdstan_log_lik_units(fit, parameter_name)
        )
        n_units <- dim(log_lik)[[3L]]
        if (n_units < 2L) {
          stop("at least two ", loo_unit,
               " units are required for LOO/WAIC", call. = FALSE)
        }
        r_eff <- loo::relative_eff(exp(sweep(
          log_lik, 3L, apply(log_lik, 3L, max), FUN = "-"
        )))
        list(
          loo = loo::loo(log_lik, r_eff = r_eff),
          waic = loo::waic(log_lik)
        )
      }, error = function(e) {
        loo_error <<- conditionMessage(e)
        NULL
      })
      if (is.null(loo_result)) {
        loo_available <- FALSE
        loo_failure <- paste0("model ", nms[i], ": ", loo_error)
        warning("Predictive comparison failed for ", loo_failure,
                "; all LOO/WAIC summaries are marked NA.", call. = FALSE)
      } else {
        loo_list[[i]] <- loo_result$loo
        waic_list[[i]] <- loo_result$waic
      }
    }

    # Bridge sampling for marginal likelihood
    if (bridge) {
      bridge_error <- contract_error
      bridge_context <- if (is.null(bridge_error)) tryCatch(
        .bayes_bridge_restore_context(fit, evidence_contract),
        error = function(e) {
          bridge_error <<- conditionMessage(e)
          warning("Bridge context restoration failed for model ", nms[i],
                  "; evidence marked incomplete: ", bridge_error,
                  call. = FALSE)
          NULL
        }
      ) else {
        warning("Bridge context restoration failed for model ", nms[i],
                "; evidence marked incomplete: ", bridge_error,
                call. = FALSE)
        NULL
      }
      if (is.null(bridge_context)) {
        bridge_list[[i]] <- rep(list(NULL), bridge_repetitions)
        bridge_context_info[[i]] <- list(status = "incomplete")
      } else {
        bridge_list[[i]] <- lapply(seq_len(bridge_repetitions), function(rep_id) {
          tryCatch(
            bridgesampling::bridge_sampler(
              bridge_context$samples, repetitions = 1L, method = bridge_method,
              silent = TRUE, maxiter = bridge_maxiter
            ),
            error = function(e) {
              warning("Bridge sampling failed for model ", nms[i],
                      ", repetition ", rep_id, ": ", conditionMessage(e))
              NULL
            }
          )
        })
        bridge_context$cleanup()
        bridge_valid <- vapply(bridge_list[[i]], function(x) {
          !is.null(x) && is.numeric(x$logml) &&
            length(x$logml) == 1L && is.finite(x$logml)
        }, logical(1))
        bridge_context_info[[i]] <- list(
          status = if (all(bridge_valid))
            "validated" else "bridge_failed",
          completed_repetitions = sum(bridge_valid),
          source = bridge_context$source,
          kernel = bridge_context$kernel,
          saved_hash = bridge_context$saved_hash,
          current_hash = bridge_context$current_hash,
          migration = bridge_context$migration
        )
      }
    }
  }

  names(loo_list) <- names(waic_list) <- nms
  if (bridge) names(bridge_list) <- nms
  if (loo_available) {
    comparison <- tryCatch({
      if (length(loo_list) == 1L) {
        list(compare = NULL, weights = stats::setNames(1, nms))
      } else {
        list(
          compare = loo::loo_compare(loo_list),
          weights = loo::loo_model_weights(loo_list, method = "stacking")
        )
      }
    }, error = function(e) {
      loo_failure <<- conditionMessage(e)
      NULL
    })
    if (is.null(comparison)) {
      loo_available <- FALSE
      warning("LOO/WAIC comparison failed; all predictive summaries are ",
              "marked NA: ", loo_failure, call. = FALSE)
    } else {
      comp <- comparison$compare
      weights_stacking <- comparison$weights
    }
  }
  if (!loo_available && !collapsed_quality) {
    if (is.null(loo_failure)) {
      warning("LOO/WAIC require at least two ", loo_unit,
              " units; returning NA predictive summaries.", call. = FALSE)
    }
    loo_list <- waic_list <- rep(list(NULL), length(fits))
    comp <- NULL
    weights_stacking <- rep(NA_real_, length(fits))
  } else if (!collapsed_quality) {
    ## All fits support the requested predictive comparison.
  } else {
    loo_list <- waic_list <- rep(list(NULL), length(fits))
    comp <- NULL
    weights_stacking <- rep(NA_real_, length(fits))
  }

  result <- data.frame(
    model         = nms,
    supply_model  = vapply(fits, function(f)
      supply_names[f@supply_model], character(1)),
    demand_model  = vapply(fits, slot, character(1), "demand_model"),
    market_side   = vapply(fits, slot, character(1), "market_side"),
    share_likelihood = vapply(fits, function(x)
      x@metadata$share_likelihood %||% "student_t", character(1)),
    margin_likelihood = vapply(fits, function(x)
      x@metadata$margin_likelihood %||% "student_t", character(1)),
    quality_parameterization = vapply(fits, function(x)
      x@metadata$quality_parameterization %||% "sampled", character(1)),
    quality_scope = vapply(fits, .quality_scope_from_fit, character(1), strict = FALSE),
    prior_manifest = vapply(fits, function(x) x@metadata$prior_manifest,
                            character(1)),
    model_prior = prior_weights,
    loo_unit = loo_unit,
    looic = if (loo_available) vapply(
      loo_list, function(x) x$estimates["looic", "Estimate"], numeric(1)
    ) else rep(NA_real_, length(fits)),
    looic_se = if (loo_available) vapply(
      loo_list, function(x) x$estimates["looic", "SE"], numeric(1)
    ) else rep(NA_real_, length(fits)),
    p_loo = if (loo_available) vapply(
      loo_list, function(x) x$estimates["p_loo", "Estimate"], numeric(1)
    ) else rep(NA_real_, length(fits)),
    waic = if (loo_available) vapply(
      waic_list, function(x) x$estimates["waic", "Estimate"], numeric(1)
    ) else rep(NA_real_, length(fits)),
    weight_stacking = as.numeric(weights_stacking),
    stringsAsFactors = FALSE
  )

  # Add bridge sampling results
  if (bridge && length(bridge_list) > 0L) {
    bridge_logml <- lapply(bridge_list, function(reps) {
      vapply(reps, function(b) if (!is.null(b)) b$logml else NA_real_, numeric(1))
    })
    logml <- vapply(bridge_logml, function(x) {
      if (any(is.finite(x))) stats::median(x[is.finite(x)]) else NA_real_
    }, numeric(1))
    logml_sd <- vapply(bridge_logml, function(x) {
      if (sum(is.finite(x)) > 1L) stats::sd(x[is.finite(x)]) else NA_real_
    }, numeric(1))

    finite <- is.finite(logml)
    diff_logml <- bf_vs_best <- post_prob <- rep(NA_real_, length(logml))
    if (all(finite)) {
      max_logml <- max(logml)
      diff_logml <- logml - max_logml
      bf_vs_best <- exp(diff_logml)
      unnorm <- exp(diff_logml) * prior_weights
      if (sum(unnorm) > 0) {
        post_prob <- unnorm / sum(unnorm)
      }
    }

    result$logml        <- logml
    result$logml_sd     <- logml_sd
    result$bf_vs_best   <- bf_vs_best
    result$post_prob    <- post_prob
  }

  rownames(result) <- NULL

  attr(result, "loo_compare") <- comp
  attr(result, "loo_list") <- loo_list
  if (bridge) attr(result, "bridge_list") <- bridge_list
  if (bridge) attr(result, "bridge_context_info") <- bridge_context_info
  attr(result, "prior_weights") <- prior_weights
  attr(result, "loo_unit") <- loo_unit
  result
}

#' Decompose Joint and Conditional-Margin Model Evidence
#'
#' Since the conditional log quantity-share likelihood and shared-parameter
#' priors are common across conduct models, subtracting a separately fitted
#' share-only log marginal likelihood isolates
#' \eqn{\log p(M\mid S,M_k)}.
#'
#' @param comparison A bridge-sampling table returned by \code{bayes_compare}.
#' @param share_logml Scalar or model-aligned vector of log marginal
#'   likelihoods from share-only fits. A scalar is recycled when the share
#'   evidence is common across conduct models.
#' @return The comparison table with common-share and conditional-margin
#'   evidence columns.
#' @export
bayes_evidence_decomposition <- function(comparison, share_logml) {
  if (!"logml" %in% names(comparison)) {
    stop("'comparison' must contain bridge-sampling log marginal likelihoods.")
  }
  if (!is.numeric(share_logml) || !length(share_logml) ||
      any(!is.finite(share_logml)) ||
      !(length(share_logml) %in% c(1L, nrow(comparison)))) {
    stop("'share_logml' must be one finite value or one value per model.")
  }
  if (length(share_logml) == 1L) share_logml <- rep(share_logml, nrow(comparison))
  out <- comparison
  out$logml_share <- share_logml
  out$logml_margin_given_share <- out$logml - share_logml
  out
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
  source_position <- stats::ave(source_order, source_order, FUN = seq_along)
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

  pooled_metadata <- list(
    source_models = names(sims),
    weights = weights,
    component_metadata = lapply(sims, slot, "metadata")
  )

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
    call           = match.call(),
    metadata       = pooled_metadata
  )
}
