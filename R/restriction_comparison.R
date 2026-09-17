# Comparison of restrictions on exact ordered-k partition states.

.restriction_logmeanexp <- function(x) {
  if (!is.numeric(x) || !length(x) || any(is.nan(x)) ||
      any(is.infinite(x) & x > 0)) {
    stop("Log probabilities must be finite or negative infinity.",
         call. = FALSE)
  }
  finite <- is.finite(x)
  if (!any(finite)) return(-Inf)
  m <- max(x[finite])
  ## Keep the denominator equal to the total number of draws.  Dropping
  ## -Inf draws would turn a zero-probability draw into a positive mass.
  m + log(sum(exp(x - m)) / length(x))
}

.restriction_log_ratio <- function(log_post, log_prior) {
  if (length(log_post) != 1L || length(log_prior) != 1L) {
    stop("Restriction log masses must be scalar values.", call. = FALSE)
  }
  if (is.na(log_post) || is.na(log_prior)) return(NA_real_)
  ## A zero-prior restriction is not a defined model comparison.  Keep this
  ## explicit instead of returning NaN from (-Inf) - (-Inf).
  if (is.infinite(log_prior) && log_prior < 0) {
    if (is.infinite(log_post) && log_post < 0) return(NA_real_)
    return(Inf)
  }
  if (is.infinite(log_post) && log_post < 0) return(-Inf)
  log_post - log_prior
}

.restriction_validate_allowed_k <- function(x, label) {
  if (!is.numeric(x) || !length(x) || anyNA(x) || any(!is.finite(x)) ||
      any(x != floor(x))) {
    stop("Restriction '", label,
         "' must contain one or more finite integer k values.",
         call. = FALSE)
  }
  ## Restrictions are sets.  Removing repeated entries also makes their
  ## result independent of the order supplied by the caller.
  as.integer(unique(x))
}

.restriction_validate_names <- function(x, what) {
  nms <- names(x)
  if (is.null(nms) || length(nms) != length(x) || anyNA(nms) ||
      any(!nzchar(nms)) ||
      anyDuplicated(nms)) {
    stop("'", what, "' must be a named collection with unique non-empty names.",
         call. = FALSE)
  }
  nms
}

.restriction_parse <- function(restrictions, markets, k_by_unit) {
  if (is.null(restrictions) || !length(restrictions)) {
    stop("'restrictions' must contain at least one named restriction.",
         call. = FALSE)
  }

  ## A named atomic vector is a convenient one-market shorthand for scalar
  ## restrictions.  The documented general form remains a named list whose
  ## values are integer vectors (or, for multiple markets, named lists).
  if (!is.list(restrictions)) {
    if (!is.numeric(restrictions)) {
      stop("'restrictions' must be a named list or numeric vector.",
           call. = FALSE)
    }
    if (is.null(names(restrictions))) {
      restrictions <- list(restriction_1 = restrictions)
    } else {
      restrictions <- as.list(restrictions)
    }
  }
  restriction_names <- .restriction_validate_names(restrictions, "restrictions")
  n_unit <- length(k_by_unit)
  n_market <- length(markets)
  if (n_unit < 1L || n_market != n_unit || anyDuplicated(markets) ||
      anyNA(markets) || any(!nzchar(markets))) {
    stop("The exact-k market map must contain one unique label per market.",
         call. = FALSE)
  }

  parsed <- lapply(seq_along(restrictions), function(r) {
    spec <- restrictions[[r]]
    rname <- restriction_names[[r]]
    if (n_unit == 1L && !is.list(spec)) {
      allowed <- list(.restriction_validate_allowed_k(spec, rname))
      names(allowed) <- markets
      return(allowed)
    }
    if (!is.list(spec)) {
      stop("Multimarket restriction '", rname,
           "' must be a named list mapping market labels to allowed k values.",
           call. = FALSE)
    }
    supplied_markets <- .restriction_validate_names(spec,
                                                    paste0("restriction '", rname, "'"))
    unknown <- setdiff(supplied_markets, markets)
    if (length(unknown)) {
      stop("Restriction '", rname, "' names unknown market label(s): ",
           paste(unknown, collapse = ", "), ".", call. = FALSE)
    }
    allowed <- lapply(seq_along(markets), function(u) {
      if (markets[[u]] %in% supplied_markets) {
        .restriction_validate_allowed_k(
          spec[[match(markets[[u]], supplied_markets)]],
          paste0(rname, "/", markets[[u]])
        )
      } else {
        ## Unspecified markets are unrestricted.
        as.integer(k_by_unit[[u]])
      }
    })
    names(allowed) <- markets
    allowed
  })
  names(parsed) <- restriction_names
  parsed
}

.restriction_unit_masses <- function(prob, k_by_unit, component_indices,
                                     allowed) {
  n_draw <- nrow(prob)
  n_unit <- length(component_indices)
  out <- matrix(0, nrow = n_draw, ncol = n_unit,
                dimnames = list(NULL, names(allowed)))
  for (u in seq_len(n_unit)) {
    ii <- component_indices[[u]]
    take <- k_by_unit[[u]] %in% allowed[[u]]
    if (any(take)) {
      out[, u] <- rowSums(prob[, ii[take], drop = FALSE])
    }
  }
  out
}

.restriction_prepare_draws <- function(fit, draws) {
  structure_name <- fit@metadata$structure %||%
    attr(fit@stan_data, "structure", exact = TRUE)
  if (length(structure_name) != 1L ||
      !structure_name %in% c("fringe", "stackelberg")) {
    stop("bayes_restriction_compare() requires structure = 'fringe' or 'stackelberg'.",
         call. = FALSE)
  }
  partition_prior <- fit@metadata$partition_prior %||%
    attr(fit@stan_data, "partition_prior", exact = TRUE)
  if (!identical(partition_prior, "ordered_k")) {
    stop("bayes_restriction_compare() requires an exact ordered-k fit.",
         call. = FALSE)
  }

  if (methods::is(fit, "BayesSpecifiedFit")) {
    ## As in posterior_k(), supplied fits have no likelihood-generated
    ## quantities.  Ignore an optional draw argument and use their coherent
    ## supplied prior for both sides of the identity.
    maps <- .posterior_k_metadata_maps(fit, structure_name)
    if (maps$n_unit != 1L || length(maps$market) != 1L ||
        anyDuplicated(maps$market)) {
      stop("BayesSpecifiedFit exact ordered-k metadata must describe one market.",
           call. = FALSE)
    }
    prior <- matrix(as.numeric(maps$prior[[1L]]), nrow = 1L)
    posterior <- prior
    return(list(
      structure = structure_name, partition_prior = partition_prior,
      markets = as.character(maps$market), k_by_unit = maps$k,
      component_indices = maps$component_indices, prior = prior,
      posterior = posterior, n_draws = 1L, specified = TRUE
    ))
  }

  if (!methods::is(fit, "BayesLogit")) {
    stop("'fit' must be a BayesLogit or BayesSpecifiedFit object.",
         call. = FALSE)
  }
  maps <- .posterior_k_stan_maps(fit@stan_data, structure_name, fit@data)
  n_component <- length(maps$k)
  n_draw <- fit@metadata$posterior_draws
  if (length(n_draw) != 1L || !is.numeric(n_draw) || !is.finite(n_draw) ||
      n_draw < 1 || n_draw != as.integer(n_draw)) n_draw <- NULL
  if (is.null(draws)) {
    draws <- .posterior_extract(fit, "posterior_k")
  }
  if (!is.list(draws)) {
    stop("'draws' must be a list containing generated ordered-k probabilities.",
         call. = FALSE)
  }
  posterior <- .posterior_k_draw_matrix(
    draws$posterior_component_prob, n_component,
    "posterior_component_prob", n_draw
  )
  if (is.null(n_draw)) n_draw <- nrow(posterior)
  prior <- .posterior_k_draw_matrix(
    draws$prior_component_prob, n_component,
    "prior_component_prob", n_draw
  )
  posterior <- .posterior_k_validate_prob_matrix(
    posterior, "posterior_component_prob", maps$component_indices
  )
  prior <- .posterior_k_validate_prob_matrix(
    prior, "prior_component_prob", maps$component_indices
  )
  if (nrow(posterior) != nrow(prior) || nrow(posterior) < 1L) {
    stop("Generated ordered-k posterior and prior draws have different counts.",
         call. = FALSE)
  }
  component_to_unit <- rep.int(NA_integer_, n_component)
  for (u in seq_len(maps$n_unit)) {
    component_to_unit[maps$component_indices[[u]]] <- u
  }
  if (anyNA(component_to_unit)) {
    stop("Stan component map does not assign every ordered-k state to a unit.",
         call. = FALSE)
  }
  list(
    structure = structure_name, partition_prior = partition_prior,
    markets = as.character(maps$market), k_by_unit = lapply(
      maps$component_indices, function(ii) as.integer(maps$k[ii])
    ), component_indices = maps$component_indices, prior = prior,
    posterior = posterior, n_draws = nrow(posterior), specified = FALSE
  )
}

#' Compare exact ordered-k model restrictions using encompassing-prior masses
#'
#' This computes Bayes factors for restrictions on the discrete ordered-k
#' partition state.  A restriction is represented by an integer vector of
#' allowed k values for a one-market fit.  For a multimarket fit, each
#' restriction is a named list mapping market labels to allowed values; any
#' omitted market is unrestricted.  Composite, unsorted, and noncontiguous
#' sets are supported.
#'
#' @param fit An exact ordered-k \code{\linkS4class{BayesFit}} with
#'   \code{structure = "fringe"} or \code{"stackelberg"}.
#' @param restrictions A named list of restrictions.  Each one-market value
#'   is an integer vector; multimarket values are named market-to-k lists.
#' @param reference Optional restriction name or one-based row index.  When
#'   supplied, the result includes Bayes factors relative to that restriction.
#' @param draws Optional generated-quantity draw list containing
#'   \code{posterior_component_prob} and \code{prior_component_prob}.  If
#'   \code{NULL}, generated quantities are extracted from the fitted Stan
#'   object.
#' @return A data frame with one row per restriction and columns
#'   \code{restriction}, \code{prior_mass}, \code{posterior_mass},
#'   \code{log_evidence_vs_full}, and \code{method}.  If \code{reference} is
#'   supplied, \code{log_bf_vs_reference} and \code{bf_vs_reference} are also
#'   included.  The parsed restrictions and fit mapping are stored as
#'   attributes.
#'
#' @details
#' The calculation uses the exact encompassing-prior identity
#' \eqn{BF_{A,F} = P(A \mid y) / P(A)}.  For multiple markets, each draw's
#' market-conditional masses are multiplied first and then averaged with a
#' log-mean-exp.  This is an exact discrete-state identity, not a
#' Gelfand-Dey estimate.  It applies when the continuous-parameter prior is
#' common across the k restrictions.  A \code{BayesSpecifiedFit} has no
#' likelihood draws and therefore returns evidence one for coherent
#' restrictions.
#'
#' @export
bayes_restriction_compare <- function(fit, restrictions, reference = NULL,
                                      draws = NULL) {
  if (!methods::is(fit, "BayesFit")) {
    stop("'fit' must be a BayesFit object.", call. = FALSE)
  }
  prepared <- .restriction_prepare_draws(fit, draws)
  parsed <- .restriction_parse(
    restrictions, prepared$markets, prepared$k_by_unit
  )
  n_restriction <- length(parsed)
  log_prior <- log_posterior <- numeric(n_restriction)
  for (r in seq_len(n_restriction)) {
    allowed <- parsed[[r]]
    prior_unit <- .restriction_unit_masses(
      prepared$prior, prepared$k_by_unit, prepared$component_indices, allowed
    )
    posterior_unit <- .restriction_unit_masses(
      prepared$posterior, prepared$k_by_unit, prepared$component_indices, allowed
    )
    ## Product is evaluated on each draw before log-mean-exp.  This retains
    ## posterior dependence across markets in the generated quantities.
    log_prior_draw <- rowSums(suppressWarnings(log(prior_unit)))
    log_posterior_draw <- rowSums(suppressWarnings(log(posterior_unit)))
    log_prior[[r]] <- .restriction_logmeanexp(log_prior_draw)
    log_posterior[[r]] <- .restriction_logmeanexp(log_posterior_draw)
  }
  zero_prior <- is.infinite(log_prior) & log_prior < 0
  if (any(zero_prior)) {
    stop("Restriction(s) with zero prior mass cannot define a restricted ",
         "model or Bayes factor: ",
         paste(names(parsed)[zero_prior], collapse = ", "), ".",
         call. = FALSE)
  }
  log_evidence <- vapply(seq_len(n_restriction), function(r) {
    .restriction_log_ratio(log_posterior[[r]], log_prior[[r]])
  }, numeric(1))

  out <- data.frame(
    restriction = names(parsed),
    prior_mass = exp(log_prior),
    posterior_mass = exp(log_posterior),
    log_evidence_vs_full = log_evidence,
    method = rep("encompassing_prior", n_restriction),
    stringsAsFactors = FALSE
  )
  if (!is.null(reference)) {
    if (length(reference) != 1L || is.na(reference)) {
      stop("'reference' must be one restriction name or row index.",
           call. = FALSE)
    }
    ref <- if (is.character(reference)) {
      match(reference, names(parsed))
    } else if (is.numeric(reference) && is.finite(reference) &&
               reference == as.integer(reference)) {
      as.integer(reference)
    } else NA_integer_
    if (is.na(ref) || ref < 1L || ref > n_restriction) {
      stop("'reference' does not identify a supplied restriction.",
           call. = FALSE)
    }
    out$log_bf_vs_reference <- vapply(seq_len(n_restriction), function(r) {
      if (is.na(log_evidence[[r]]) || is.na(log_evidence[[ref]])) {
        NA_real_
      } else if (is.infinite(log_evidence[[r]]) &&
                 is.infinite(log_evidence[[ref]]) &&
                 log_evidence[[r]] == log_evidence[[ref]]) {
        NA_real_
      } else {
        log_evidence[[r]] - log_evidence[[ref]]
      }
    }, numeric(1))
    out$bf_vs_reference <- exp(out$log_bf_vs_reference)
    attr(out, "reference") <- names(parsed)[[ref]]
  }
  attr(out, "method") <- "encompassing_prior"
  attr(out, "restrictions") <- parsed
  attr(out, "metadata") <- list(
    structure = prepared$structure,
    partition_prior = prepared$partition_prior,
    markets = prepared$markets,
    k_support = prepared$k_by_unit,
    n_draws = prepared$n_draws,
    specified = prepared$specified
  )
  class(out) <- c("bayes_restriction_comparison", class(out))
  out
}
