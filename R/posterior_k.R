# Exact ordered-k posterior summaries.

## The exact partition is a market-level latent state.  This file deliberately
## does not use the cutoff summary's smooth role calculations: an ordered-k
## fit has one common k for the complete market and its generated quantity is
## already the conditional posterior over those states.

.posterior_k_scalar <- function(x, name, integer = FALSE) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x)) {
    stop("'", name, "' must be a single finite numeric value.", call. = FALSE)
  }
  if (integer && x != as.integer(x)) {
    stop("'", name, "' must be an integer.", call. = FALSE)
  }
  as.numeric(x)
}

.posterior_k_discrete_quantile <- function(k, p, prob) {
  ord <- order(k)
  k <- as.numeric(k[ord])
  p <- as.numeric(p[ord])
  k[[which(cumsum(p) >= prob - 1e-12)[[1L]]]]
}

.posterior_k_draw_matrix <- function(x, n_component, name, n_draw = NULL) {
  if (is.null(x)) {
    stop("Exact ordered-k fit is missing generated quantity '", name, "'.",
         call. = FALSE)
  }
  if (!is.numeric(x)) {
    stop("Generated quantity '", name, "' must be numeric.", call. = FALSE)
  }
  d <- dim(x)
  if (!is.null(d)) {
    if (length(d) != 2L) {
      stop("Generated quantity '", name,
           "' must have draws by components dimensions.", call. = FALSE)
    }
    if (d[[2L]] != n_component) {
      stop("Generated quantity '", name, "' has ", d[[2L]],
           " components; Stan data map ", n_component, ".", call. = FALSE)
    }
    if (!is.null(n_draw) && d[[1L]] != n_draw) {
      stop("Generated quantity '", name, "' has ", d[[1L]],
           " draws; the fitted object records ", n_draw, ".", call. = FALSE)
    }
    return(matrix(as.numeric(x), nrow = d[[1L]], ncol = d[[2L]]))
  }
  if (is.null(n_draw)) {
    if (length(x) %% n_component != 0L) {
      stop("Generated quantity '", name,
           "' has a length incompatible with the component map.", call. = FALSE)
    }
    n_draw <- length(x) / n_component
  }
  if (length(x) != n_draw * n_component) {
    stop("Generated quantity '", name, "' has length ", length(x),
         "; expected ", n_draw * n_component, ".", call. = FALSE)
  }
  matrix(as.numeric(x), nrow = n_draw, ncol = n_component, byrow = TRUE)
}

.posterior_k_validate_prob_matrix <- function(x, name, component_indices,
                                               tol = 1e-8) {
  if (any(!is.finite(x)) || any(x < -tol) || any(x > 1 + tol)) {
    stop("Generated quantity '", name,
         "' contains invalid probabilities.", call. = FALSE)
  }
  ## Each market has a separate categorical distribution. Validate and
  ## normalize each unit separately; summing the flattened row would conceal
  ## a malformed all-zero market whenever another market has positive mass.
  x[x < 0] <- 0
  for (u in seq_along(component_indices)) {
    ii <- component_indices[[u]]
    rs <- rowSums(x[, ii, drop = FALSE])
    if (any(!is.finite(rs)) || any(rs <= tol)) {
      stop("Generated quantity '", name,
           "' has a draw with no positive probability mass in partition unit ",
           u, ".", call. = FALSE)
    }
    x[, ii] <- x[, ii, drop = FALSE] / rs
  }
  x
}

.posterior_k_stan_maps <- function(stan_data, structure_name, fit_data = NULL) {
  required <- c("N_partition_unit", "partition_component_k",
                "partition_unit_market")
  missing <- setdiff(required, names(stan_data))
  if (length(missing)) {
    stop("Exact ordered-k fit is missing Stan partition map field(s): ",
         paste(missing, collapse = ", "), ".", call. = FALSE)
  }
  n_unit <- .posterior_k_scalar(stan_data$N_partition_unit,
                                "N_partition_unit", integer = TRUE)
  if (n_unit < 1L) {
    stop("Stan field 'N_partition_unit' must be positive.", call. = FALSE)
  }
  n_unit <- as.integer(n_unit)
  k <- as.numeric(stan_data$partition_component_k)
  if (!length(k) || any(!is.finite(k)) || any(k != as.integer(k))) {
    stop("Stan field 'partition_component_k' is malformed.", call. = FALSE)
  }
  k <- as.integer(k)
  if (structure_name == "stackelberg" && any(k < 1L)) {
    stop("Stackelberg ordered-k maps must contain only positive leader counts.",
         call. = FALSE)
  }
  if (structure_name == "fringe" && any(k < 0L)) {
    stop("Fringe ordered-k maps cannot contain negative core counts.", call. = FALSE)
  }
  unit_market <- as.integer(stan_data$partition_unit_market)
  if (length(unit_market) != n_unit || anyNA(unit_market) ||
      any(unit_market < 1L)) {
    stop("Stan field 'partition_unit_market' is malformed.", call. = FALSE)
  }
  component_unit <- if ("partition_component_unit" %in% names(stan_data)) {
    as.integer(stan_data$partition_component_unit)
  } else NULL
  if (!is.null(component_unit) &&
      (length(component_unit) != length(k) || anyNA(component_unit) ||
       any(component_unit < 1L | component_unit > n_unit))) {
    stop("Stan field 'partition_component_unit' is malformed.", call. = FALSE)
  }
  starts <- counts <- NULL
  if (all(c("partition_unit_start", "partition_unit_count") %in%
          names(stan_data))) {
    starts <- as.integer(stan_data$partition_unit_start)
    counts <- as.integer(stan_data$partition_unit_count)
    if (length(starts) != n_unit || length(counts) != n_unit ||
        anyNA(starts) || anyNA(counts) || any(starts < 1L) ||
        any(counts < 1L) || any(starts + counts - 1L > length(k))) {
      stop("Stan partition unit start/count maps are malformed.", call. = FALSE)
    }
  }
  if (is.null(component_unit) && is.null(starts)) {
    stop("Exact ordered-k fit has no usable Stan component-to-unit map.",
         call. = FALSE)
  }
  component_indices <- lapply(seq_len(n_unit), function(u) {
    by_unit <- if (!is.null(component_unit)) which(component_unit == u) else integer()
    by_range <- if (!is.null(starts)) {
      seq.int(starts[[u]], length.out = counts[[u]])
    } else integer()
    if (length(by_unit) && length(by_range) &&
        !identical(sort(by_unit), sort(by_range))) {
      stop("Stan partition component maps disagree for partition unit ", u,
           ".", call. = FALSE)
    }
    ans <- if (length(by_unit)) by_unit else by_range
    if (!length(ans) || anyDuplicated(ans)) {
      stop("Stan partition unit ", u, " has no unique component states.",
           call. = FALSE)
    }
    if (anyDuplicated(k[ans])) {
      stop("Stan partition unit ", u, " repeats a k state.", call. = FALSE)
    }
    ans
  })
  if (length(unique(unit_market)) != n_unit) {
    stop("Exact ordered-k requires one partition unit per market; the Stan map has ",
         "multiple units for at least one market.", call. = FALSE)
  }
  market_levels <- attr(stan_data, "market_levels", exact = TRUE)
  if (is.null(market_levels) && !is.null(fit_data) && "market" %in% names(fit_data)) {
    market_levels <- levels(as.factor(fit_data$market))
  }
  if (!is.null(market_levels) &&
      any(unit_market > length(market_levels))) {
    stop("Stan partition unit market indices exceed the fitted market map.",
         call. = FALSE)
  }
  market <- if (is.null(market_levels)) {
    as.character(unit_market)
  } else {
    as.character(market_levels[unit_market])
  }
  list(
    n_unit = n_unit, k = k, component_unit = component_unit,
    component_indices = component_indices, unit_market = unit_market,
    market = market
  )
}

.posterior_k_metadata_maps <- function(fit, structure_name) {
  md <- fit@metadata
  supports <- md$k_support_resolved %||% md$prior_normalization$k_support_resolved %||%
    md$prior_normalization$partition_support
  priors <- md$k_prior_resolved %||% md$prior_normalization$k_prior_resolved
  if (is.null(supports) || is.null(priors) || !is.list(supports) ||
      !is.list(priors) || length(supports) != length(priors) ||
      !length(supports)) {
    stop("Exact ordered-k fit is missing market-specific k support/prior metadata.",
         call. = FALSE)
  }
  n_unit <- length(supports)
  component_indices <- lapply(seq_len(n_unit), seq_along)
  k <- lapply(supports, function(x) {
    if (!is.numeric(x) || !length(x) || any(!is.finite(x)) ||
        any(x != as.integer(x)) || anyDuplicated(x) ||
        (structure_name == "stackelberg" && any(x < 1L)) ||
        (structure_name == "fringe" && any(x < 0L))) {
      stop("Exact ordered-k metadata contains a malformed k support.", call. = FALSE)
    }
    as.integer(x)
  })
  p <- lapply(priors, function(x) {
    if (!is.numeric(x) || !length(x) || any(!is.finite(x)) || any(x < 0)) {
      stop("Exact ordered-k metadata contains a malformed k prior.", call. = FALSE)
    }
    if (sum(x) <= 0 || abs(sum(x) - 1) > 1e-8) {
      stop("Exact ordered-k metadata contains a k prior that does not sum to one.",
           call. = FALSE)
    }
    as.numeric(x / sum(x))
  })
  if (any(lengths(k) != lengths(p))) {
    stop("Exact ordered-k support and prior metadata are not conformable.",
         call. = FALSE)
  }
  ## A BayesSpecifiedFit represents one market by construction.  Keep the
  ## label explicit instead of inferring it from factor levels or firm IDs.
  market <- "unit_1"
  list(n_unit = n_unit, k = k, prior = p, component_indices = component_indices,
       market = rep(market, n_unit))
}

#' Summarize posterior probabilities over exact ordered strategic-core sizes.
#'
#' For an exact ordered-k fit, one state determines the complete equilibrium
#' for a market.  This helper averages Stan's generated conditional state
#' probabilities over posterior draws and returns one row per market and k.
#' For a fringe fit, k is the number of strategic core firms; for a
#' Stackelberg fit, k is the number of simultaneous leaders.
#'
#' @param fit A \code{BayesFit} produced with \code{partition_prior =
#'   "ordered_k"}. A \code{BayesSpecifiedFit} reports its supplied prior as
#'   the posterior because it has no likelihood draws.
#' @param level Credible interval level for the discrete posterior summary.
#' @param draws Optional generated-quantity draw list. This is primarily for
#'   testing and serialized-fit recovery; when \code{NULL}, draws are obtained
#'   from \code{rstan::extract}.
#' @return A data frame with one row per market and supported k. Columns are
#'   \code{market}, \code{k}, \code{prior}, \code{posterior}, \code{role},
#'   \code{posterior_mean}, \code{posterior_median},
#'   \code{posterior_lower}, and \code{posterior_upper}. The \code{summary}
#'   attribute contains one row per market and \code{metadata} records the
#'   structure and component mapping used.
#' @export
posterior_k <- function(fit, level = 0.95, draws = NULL) {
  if (!methods::is(fit, "BayesFit")) {
    stop("'fit' must be a BayesFit object.", call. = FALSE)
  }
  structure_name <- fit@metadata$structure %||% attr(fit@stan_data,
                                                       "structure", exact = TRUE)
  if (length(structure_name) != 1L ||
      !structure_name %in% c("fringe", "stackelberg")) {
    stop("posterior_k() requires structure = 'fringe' or 'stackelberg'.",
         call. = FALSE)
  }
  partition_prior <- fit@metadata$partition_prior %||%
    attr(fit@stan_data, "partition_prior", exact = TRUE)
  if (!identical(partition_prior, "ordered_k")) {
    stop("posterior_k() requires an exact ordered-k fit.", call. = FALSE)
  }
  if (length(level) != 1L || !is.numeric(level) || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("'level' must lie in (0, 1).", call. = FALSE)
  }
  role <- if (identical(structure_name, "stackelberg")) "leaders" else "strategic_core"
  estimated <- !methods::is(fit, "BayesSpecifiedFit")

  if (estimated) {
    stan_data <- fit@stan_data
    maps <- .posterior_k_stan_maps(stan_data, structure_name, fit@data)
    n_component <- length(maps$k)
    n_draw <- fit@metadata$posterior_draws
    if (!is.numeric(n_draw) || length(n_draw) != 1L || !is.finite(n_draw) || n_draw < 1L ||
        n_draw != as.integer(n_draw)) n_draw <- NULL
    if (is.null(draws)) {
      draws <- rstan::extract(fit@stanfit, pars = c(
        "posterior_component_prob", "prior_component_prob"
      ))
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
    if (nrow(posterior) != nrow(prior)) {
      stop("Generated ordered-k posterior and prior draws have different counts.",
           call. = FALSE)
    }
    component_to_unit <- rep.int(NA_integer_, n_component)
    for (u in seq_len(maps$n_unit)) component_to_unit[maps$component_indices[[u]]] <- u
    if (anyNA(component_to_unit)) {
      stop("Stan component map does not assign every ordered-k state to a unit.",
           call. = FALSE)
    }
    supports <- lapply(maps$component_indices, function(ii) maps$k[ii])
    units <- lapply(seq_len(maps$n_unit), function(u) {
      ii <- maps$component_indices[[u]]
      list(k = maps$k[ii], prior = colMeans(prior[, ii, drop = FALSE]),
           posterior = colMeans(posterior[, ii, drop = FALSE]))
    })
    market <- maps$market
    mapping <- list(stan = TRUE, component_unit = component_to_unit,
                    component_k = maps$k, component_indices = maps$component_indices)
  } else {
    maps <- .posterior_k_metadata_maps(fit, structure_name)
    units <- lapply(seq_len(maps$n_unit), function(u) {
      list(k = maps$k[[u]], prior = maps$prior[[u]],
           posterior = maps$prior[[u]])
    })
    market <- maps$market
    mapping <- list(stan = FALSE, component_indices = maps$component_indices)
  }

  summary_rows <- lapply(seq_along(units), function(u) {
    x <- units[[u]]
    ord <- order(x$k)
    k <- as.numeric(x$k[ord]); p <- as.numeric(x$posterior[ord]);
    prior <- as.numeric(x$prior[ord])
    if (length(k) != length(p) || length(k) != length(prior) ||
        any(!is.finite(p)) || sum(p) <= 0) {
      stop("Ordered-k market component probabilities are malformed for market '",
           market[[u]], "'.", call. = FALSE)
    }
    p <- p / sum(p); prior <- prior / sum(prior)
    data.frame(
      market = as.character(market[[u]]), k = as.integer(k),
      prior = prior, posterior = p, role = role,
      posterior_mean = sum(k * p),
      posterior_median = .posterior_k_discrete_quantile(k, p, .5),
      posterior_lower = .posterior_k_discrete_quantile(k, p, (1 - level) / 2),
      posterior_upper = .posterior_k_discrete_quantile(k, p, 1 - (1 - level) / 2),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, summary_rows)
  rownames(out) <- NULL
  smry <- do.call(rbind, lapply(summary_rows, function(x) x[1L, c(
    "market", "role", "posterior_mean", "posterior_median",
    "posterior_lower", "posterior_upper"
  ), drop = FALSE]))
  rownames(smry) <- NULL
  attr(out, "summary") <- smry
  attr(out, "metadata") <- list(
    structure = structure_name, role = role, partition_prior = partition_prior,
    level = level, n_draws = if (estimated) nrow(posterior) else fit@params$n_draws,
    market = market, mapping = mapping,
    k_support = lapply(units, function(x) as.integer(x$k)),
    k_prior = lapply(units, function(x) as.numeric(x$prior))
  )
  class(out) <- c("posterior_k_summary", class(out))
  out
}
