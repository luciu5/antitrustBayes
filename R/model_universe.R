# Formal model-universe descriptors and evidence aggregation.

.mu_or <- function(x, y) if (is.null(x) || !length(x)) y else x

.mu_scalar_character <- function(x, name) {
  if (length(x) != 1L || !is.character(x) || is.na(x) || !nzchar(x)) {
    stop("'", name, "' must be one non-empty character value.", call. = FALSE)
  }
  x
}

.mu_validate_prior_weights <- function(prior_weights, hypothesis_names) {
  if (missing(prior_weights) || is.null(prior_weights)) {
    stop("'prior_weights' is required for a model universe; specify the full ",
         "scientific prior over the named hypotheses.", call. = FALSE)
  }
  if (!is.numeric(prior_weights) || length(prior_weights) != length(hypothesis_names) ||
      anyNA(prior_weights) || any(!is.finite(prior_weights)) ||
      any(prior_weights < 0) || abs(sum(prior_weights) - 1) > 1e-8) {
    stop("'prior_weights' must be finite, non-negative, and sum to 1 with one ",
         "entry per named hypothesis.", call. = FALSE)
  }
  if (!is.null(names(prior_weights))) {
    if (!identical(as.character(names(prior_weights)), hypothesis_names)) {
      stop("Named 'prior_weights' must use exactly the hypothesis names in the ",
           "same order.", call. = FALSE)
    }
  }
  stats::setNames(as.numeric(prior_weights), hypothesis_names)
}

.mu_parse_hypothesis <- function(spec, label) {
  if (is.character(spec)) {
    source <- .mu_scalar_character(spec, paste0("hypotheses[['", label, "']]") )
    return(list(
      hypothesis = label, source = source, kind = "atomic",
      fit = source, k = NULL, k_requested = NULL
    ))
  }
  if (!is.list(spec) || is.null(spec$fit)) {
    stop("Hypothesis '", label, "' must be a source fit name or a list with ",
         "a scalar 'fit' name.", call. = FALSE)
  }
  source <- .mu_scalar_character(spec$fit,
                                 paste0("hypotheses[['", label, "']]$fit"))
  unknown <- setdiff(names(spec), c("fit", "k"))
  if (length(unknown)) {
    stop("Hypothesis '", label, "' has unsupported field(s): ",
         paste(unknown, collapse = ", "), ".", call. = FALSE)
  }
  k <- spec$k
  if (!is.null(k) && !identical(k, "interior") &&
      !is.numeric(k) && !is.list(k)) {
    stop("Hypothesis '", label, "' has invalid 'k'; use 'interior', integer ",
         "values, or a named market-to-values list.", call. = FALSE)
  }
  if (is.character(k) && !identical(k, "interior")) {
    stop("Hypothesis '", label, "' has invalid character 'k'.", call. = FALSE)
  }
  list(
    hypothesis = label, source = source, kind = if (is.null(k)) "atomic" else "family",
    fit = source, k = k, k_requested = k
  )
}

#' Declare a non-overlapping Bayesian model universe
#'
#' The descriptor is consumed by \code{\link{bayes_compare}}.  Entries are
#' named atomic hypotheses: a character value names a supplied full fit, and
#' an entry such as \code{list(fit = "BF", k = "interior")} retains a strict
#' interior ordered-k family.  The prior weights are deliberately explicit;
#' this constructor never reallocates a family prior after endpoints are
#' removed.
#'
#' @param hypotheses Named list of atomic hypothesis specifications.
#' @param prior_weights Numeric vector of scientific prior probabilities, in
#'   the same order as \code{hypotheses}.
#' @return An object of class \code{bayes_model_universe}.
#' @export
bayes_model_universe <- function(hypotheses, prior_weights) {
  if (!is.list(hypotheses) || !length(hypotheses) ||
      is.null(names(hypotheses)) || anyNA(names(hypotheses)) ||
      any(!nzchar(names(hypotheses))) || anyDuplicated(names(hypotheses))) {
    stop("'hypotheses' must be a non-empty named list with unique names.",
         call. = FALSE)
  }
  hypothesis_names <- as.character(names(hypotheses))
  parsed <- lapply(seq_along(hypotheses), function(i) {
    .mu_parse_hypothesis(hypotheses[[i]], hypothesis_names[[i]])
  })
  names(parsed) <- hypothesis_names
  weights <- .mu_validate_prior_weights(prior_weights, hypothesis_names)
  structure(
    list(
      hypotheses = parsed,
      prior_weights = weights,
      prior_weight_provenance = list(
        supplied = TRUE,
        names = hypothesis_names,
        normalized = FALSE,
        rule = "explicit_named_scientific_prior"
      )
    ),
    class = "bayes_model_universe"
  )
}

.mu_validate_descriptor <- function(universe) {
  if (!inherits(universe, "bayes_model_universe") ||
      !is.list(universe$hypotheses) || !length(universe$hypotheses)) {
    stop("'universe' must be a bayes_model_universe descriptor.", call. = FALSE)
  }
  nms <- names(universe$hypotheses)
  if (is.null(nms) || anyNA(nms) || any(!nzchar(nms)) || anyDuplicated(nms)) {
    stop("The model-universe hypotheses must have unique non-empty names.",
         call. = FALSE)
  }
  weights <- .mu_validate_prior_weights(universe$prior_weights, nms)
  universe$prior_weights <- weights
  invisible(universe)
}

.mu_structure <- function(fit) {
  md <- fit@metadata
  stan <- fit@stan_data
  value <- .mu_or(md$structure, attr(stan, "structure", exact = TRUE))
  if (is.null(value) && isTRUE(stan$structure == 0L)) value <- "simultaneous"
  if (is.null(value) && isTRUE(stan$structure == 1L)) value <- "fringe"
  if (is.null(value) && isTRUE(stan$structure == 2L)) value <- "stackelberg"
  if (length(value) != 1L || !is.character(value) ||
      !value %in% c("simultaneous", "fringe", "stackelberg")) {
    stop("Each universe source fit must identify its structural family.",
         call. = FALSE)
  }
  value
}

.mu_partition_prior <- function(fit) {
  md <- fit@metadata
  value <- .mu_or(md$partition_prior,
                 attr(fit@stan_data, "partition_prior", exact = TRUE))
  if (is.null(value) && isTRUE(fit@stan_data$partition_prior_code == 1L)) {
    value <- "ordered_k"
  }
  .mu_or(value, "cutoff")
}

.mu_supply <- function(fit) {
  value <- .mu_or(fit@metadata$supply_model, fit@supply_model)
  if (is.numeric(value)) {
    names <- c("bertrand", "monopolistic", "cournot", "auction", "bargaining")
    if (length(value) != 1L || !is.finite(value) || value != as.integer(value) ||
        value < 1 || value > length(names)) {
      stop("A universe source fit has an invalid supply model.", call. = FALSE)
    }
    return(names[[as.integer(value)]])
  }
  canonical_supply_name(as.character(value))
}

.mu_scalar_state <- function(x, fallback = NULL) {
  if (is.null(x) || !length(x)) return(fallback)
  if (length(x) == 1L) return(x)
  paste(as.character(x), collapse = "|")
}

.mu_active_priors <- function(fit, gaussian_only) {
  metadata_priors <- fit@metadata$resolved_priors
  resolved <- if (is.list(metadata_priors)) metadata_priors else list()
  stan <- fit@stan_data
  n_market <- .mu_or(stan$N_market,
                     length(unique(as.character(fit@data$market))))
  alpha_pool <- isTRUE(stan$use_alpha_market_pool == 1L) && n_market > 1L
  s0_scale <- isTRUE(stan$use_s0_market_year_pool == 1L) ||
    (isTRUE(stan$use_s0_market_pool == 1L) && n_market > 1L)
  k_s0 <- .mu_or(stan$K_s0, if (is.matrix(stan$X_s0)) ncol(stan$X_s0) else 0L)
  inactive <- character()
  if (!.bayes_uses_cutoff_partition(fit))
    inactive <- c(inactive, "prior_cutoff_alpha", "prior_cutoff_beta")
  if (!.bayes_is_smooth_fringe(fit))
    inactive <- c(inactive, "prior_sigma_meanval_fringe")
  if (!alpha_pool) inactive <- c(inactive, "prior_sigma_alpha")
  if (!s0_scale) inactive <- c(inactive, "prior_sigma_s0_pool")
  if (k_s0 == 0L) inactive <- c(inactive, "prior_sigma_beta_s0")
  if (!isTRUE(stan$use_market_resid_scale == 1L))
    inactive <- c(inactive, "prior_sigma_resid_market")
  if (!isTRUE(stan$estimate_df == 1L))
    inactive <- c(inactive, "prior_nu_shape", "prior_nu_rate")
  if (!identical(.mu_supply(fit), "bargaining"))
    inactive <- c(inactive, grep("^bargaining_",
                                union(names(resolved), names(stan)), value = TRUE))
  if (!isTRUE(stan$use_rho == 1L))
    inactive <- c(inactive, "prior_lkj")
  if (gaussian_only) inactive <- c(inactive, "resid_df_fixed")
  authoritative <- stan[grep("^(prior_|s0_prior_|bargaining_prior_)", names(stan))]
  authoritative <- authoritative[setdiff(names(authoritative), inactive)]
  metadata_sources <- list(metadata_priors,
                          fit@metadata$prior_normalization$resolved_priors)
  for (source in metadata_sources) {
    if (!is.list(source)) next
    for (field in intersect(names(source), names(authoritative))) {
      if (!isTRUE(all.equal(as.numeric(source[[field]]),
                            as.numeric(authoritative[[field]]),
                            tolerance = 1e-12, check.attributes = FALSE))) {
        stop("Active prior metadata disagrees with authoritative Stan field '",
             field, "'.", call. = FALSE)
      }
    }
  }
  resolved <- c(authoritative,
                resolved[setdiff(names(resolved), names(authoritative))])
  resolved <- resolved[setdiff(names(resolved), inactive)]
  resolved[order(names(resolved))]
}

.mu_exact_prior_state <- function(fit, restriction = NULL) {
  if (!.mu_structure(fit) %in% c("fringe", "stackelberg") ||
      !identical(.mu_partition_prior(fit), "ordered_k")) return(NULL)
  maps <- .mu_exact_maps(fit)
  stan <- fit@stan_data
  probability <- as.numeric(stan$partition_prior_probability)
  log_probability <- as.numeric(stan$partition_prior_logprob)
  if (length(probability) != length(maps$k) ||
      length(log_probability) != length(probability) ||
      any(!is.finite(probability)) || any(probability < 0) ||
      anyNA(log_probability)) {
    stop("Ordered-k fingerprint requires valid authoritative Stan prior probabilities.",
         call. = FALSE)
  }
  positive <- probability > 0
  if (any(!is.finite(log_probability[positive])) ||
      any(abs(log_probability[positive] - log(probability[positive])) > 1e-10) ||
      any(log_probability[!positive] != -Inf)) {
    stop("Ordered-k Stan prior probabilities and log probabilities disagree.",
         call. = FALSE)
  }
  md <- fit@metadata
  normalization <- if (is.list(md$prior_normalization)) md$prior_normalization else list()
  align_metadata <- function(x, field) {
    if (is.null(x)) return(NULL)
    if (!is.list(x) || length(x) != maps$n_unit)
      stop("Resolved ordered-k ", field, " metadata is malformed.", call. = FALSE)
    if (!is.null(names(x))) {
      if (anyDuplicated(names(x)) || !setequal(names(x), as.character(maps$market)))
        stop("Resolved ordered-k ", field, " metadata market labels disagree.", call. = FALSE)
      x <- x[as.character(maps$market)]
    }
    x
  }
  metadata_candidates <- list(
    list(support = align_metadata(md$k_support_resolved, "support"),
         prior = align_metadata(md$k_prior_resolved, "prior")),
    list(support = align_metadata(normalization$k_support_resolved, "support"),
         prior = align_metadata(normalization$k_prior_resolved, "prior"))
  )
  state <- lapply(seq_len(maps$n_unit), function(u) {
    ii <- maps$component_indices[[u]]
    k <- as.integer(maps$k[ii]); p <- probability[ii]
    if (abs(sum(p) - 1) > 1e-8)
      stop("Ordered-k Stan prior must sum to one within each market.", call. = FALSE)
    p <- p / sum(p)
    for (candidate in metadata_candidates) {
      mk <- if (is.null(candidate$support)) k else candidate$support[[u]]
      if (!is.numeric(mk) || anyNA(mk) || any(!is.finite(mk)) ||
          any(mk != floor(mk)) || anyDuplicated(mk) || !setequal(mk, k))
        stop("Resolved ordered-k support metadata disagrees with Stan.", call. = FALSE)
      if (!is.null(candidate$prior)) {
        mp <- as.numeric(candidate$prior[[u]])
        if (length(mp) != length(mk) || any(!is.finite(mp)) || any(mp < 0) ||
            abs(sum(mp) - 1) > 1e-8 ||
            any(abs(mp[match(k, mk)] / sum(mp) - p) > 1e-12))
          stop("Resolved ordered-k prior metadata disagrees with Stan.", call. = FALSE)
      }
    }
    allowed <- if (is.null(restriction)) k else restriction$allowed[[u]]
    take <- k %in% allowed & p > 0
    k <- k[take]; p <- p[take]
    if (length(p)) p <- p / sum(p)
    oo <- order(k)
    list(k = unname(k[oo]), probability = unname(p[oo]),
         zero_prior_event = !length(p))
  })
  names(state) <- as.character(maps$market)
  state[order(names(state))]
}

.mu_fingerprint <- function(fit, restriction = NULL) {
  md <- fit@metadata
  dat <- fit@data
  stan <- fit@stan_data
  data_fields <- intersect(c(
    "prices", "quantities", "margins", "market", "year", "group",
    "margin_observed", "outside_price", "loo_cluster", "cluster"
  ), names(dat))
  data_state <- lapply(dat[data_fields], function(x) unname(as.vector(x)))
  names(data_state) <- data_fields
  stan_fields <- intersect(c(
    "shareIn", "clr_shareIn", "marginNormInv", "log_marginNormInv",
    "margin_obs", "market_year_idx", "cluster", "likelihood_mode",
    "outside_price", "X_s0", "s0_pooling",
    "s0_pooling_resolved", "use_rho", "use_alpha_market_pool",
    "use_s0_market_pool", "use_s0_market_year_pool", "use_market_resid_scale",
    "use_hmt", "estimate_df", "resid_df_fixed", "quality_scope",
    "N_quality", "obs_to_quality", "N_quality_block",
    "quality_block_start", "quality_block_n", "N_quality_coord",
    "quality_coord_start", "market_obs_order", "market_obs_start",
    "market_obs_n", "market_margin_order", "market_margin_start",
    "market_margin_n"
  ), names(stan))
  stan_state <- lapply(stan[stan_fields], function(x) {
    if (is.matrix(x) || is.data.frame(x)) unname(as.matrix(x)) else unname(x)
  })
  names(stan_state) <- stan_fields
  quality_scope <- .quality_scope_from_fit(fit, strict = FALSE)
  quality_incidence <- .quality_incidence_from_fit(fit, strict = FALSE)
  quality_lookup <- if (is.null(quality_incidence)) NULL else
    quality_incidence$quality_lookup
  n_market <- stan$N_market
  if (length(n_market) != 1L || is.na(n_market) ||
      !is.finite(n_market) || n_market < 1L) {
    n_market <- if ("market" %in% names(dat)) {
      length(unique(as.character(dat$market)))
    } else if (is.list(quality_incidence)) {
      length(quality_incidence$market_obs_n)
    } else {
      NA_integer_
    }
  }
  n_market <- suppressWarnings(as.integer(n_market[[1L]]))
  ## A market-scoped one-market block is the historical global block.  Keep
  ## scope-specific metadata and caches intact, while making the statistical
  ## model-universe identity use the global code and lookup convention so the
  ## same hypothesis is not counted twice.
  canonical_single_market <- isTRUE(n_market == 1L)
  if (canonical_single_market && "quality_scope" %in% names(stan_state))
    stan_state$quality_scope <- 1L
  if (canonical_single_market && is.data.frame(quality_lookup) &&
      nrow(quality_lookup)) {
    quality_lookup$market_index <- NA_integer_
    if ("market" %in% names(quality_lookup))
      quality_lookup$market <- NA_character_
    if ("market_label" %in% names(quality_lookup))
      quality_lookup$market_label <- NA_character_
    if ("quality_unit" %in% names(quality_lookup) &&
        "group_label" %in% names(quality_lookup))
      quality_lookup$quality_unit <- as.character(quality_lookup$group_label)
  }
  share_mode <- .mu_or(md$share_likelihood,
                       attr(stan, "share_likelihood", exact = TRUE))
  margin_mode <- .mu_or(md$margin_likelihood,
                        attr(stan, "margin_likelihood", exact = TRUE))
  gaussian_only <- identical(share_mode, "gaussian") &&
    identical(margin_mode, "gaussian")
  if (gaussian_only) stan_state$resid_df_fixed <- NULL
  rho_mode <- .mu_or(md$use_rho, stan$use_rho)
  rho_off <- length(rho_mode) == 1L &&
    (isFALSE(rho_mode) || (is.numeric(rho_mode) && isTRUE(rho_mode == 0)))
  quality_mode <- if (identical(share_mode, "gaussian") &&
                      (identical(margin_mode, "gaussian") || rho_off)) NULL else
    md$quality_parameterization
  likelihood_mode_metadata <- list(
    share_likelihood = share_mode,
    margin_likelihood = margin_mode,
    quality_parameterization = quality_mode,
    mode = paste(share_mode, if (is.null(quality_mode)) "gaussian_equivalent" else
                   quality_mode, margin_mode, sep = "/")
  )
  attr_state <- list(
    likelihood_mode_metadata = likelihood_mode_metadata,
    s0_pooling = attr(stan, "s0_pooling", exact = TRUE),
    market_levels = attr(stan, "market_levels", exact = TRUE),
    quality_scope = if (canonical_single_market) "global" else quality_scope,
    quality_lookup = quality_lookup
  )
  resolved <- .mu_active_priors(fit, gaussian_only)
  normalization <- md$prior_normalization
  model_spec_version <- .mu_or(
    md$model_spec_version,
    .mu_or(attr(stan, "model_spec_version", exact = TRUE),
           if (is.list(normalization)) normalization$model_spec_version else NULL)
  )
  state <- list(
    model_spec_version = model_spec_version,
    demand_model = .mu_or(md$demand_model, fit@demand_model),
    market_side = fit@market_side,
    likelihood = .mu_or(md$likelihood, NULL),
    share_likelihood = share_mode,
    margin_likelihood = margin_mode,
    ## Gaussian quality collapse preserves the statistical model. Independent
    ## non-Gaussian margins also preserve this equivalence when rho is off;
    ## other joint specifications retain their quality-mode fingerprint.
    quality_parameterization = {
      quality_mode
    },
    quality_scope = attr_state$quality_scope,
    likelihood_mode = .mu_or(md$likelihood_mode, stan$likelihood_mode),
    use_rho = .mu_or(md$use_rho, stan$use_rho),
    prior_manifest = if (length(resolved)) NULL else md$prior_manifest,
    resolved_priors = resolved,
    ordered_k_prior = .mu_exact_prior_state(fit, restriction),
    df_mode = md$df_mode,
    resid_df_fixed = if (gaussian_only) NULL else md$resid_df_fixed,
    data = data_state,
    stan = stan_state,
    attributes = attr_state
  )
  ## Version 2 materializes ALTREP vectors while preserving the complete state,
  ## names, classes, dimensions, and attributes.  Version 3 can retain the
  ## original ALTREP/reference graph, making equal M=1 global/market states
  ## hash differently solely because their incidence objects were constructed
  ## through different paths.
  digest::digest(serialize(state, NULL, version = 2), algo = "sha256")
}

.mu_source_fits <- function(args, descriptor) {
  if (!length(args) || is.null(names(args)) || any(!nzchar(names(args))) ||
      anyDuplicated(names(args))) {
    stop("Universe comparison requires every supplied source fit to have a ",
         "unique name matching a descriptor source.", call. = FALSE)
  }
  is_fit <- vapply(args, methods::is, logical(1), "BayesLogit")
  if (!all(is_fit)) {
    stop("Universe comparison currently requires named BayesLogit source fits.",
         call. = FALSE)
  }
  source_names <- unique(vapply(descriptor$hypotheses, `[[`, character(1), "source"))
  missing <- setdiff(source_names, names(args))
  if (length(missing)) {
    stop("Universe source fit(s) not supplied in '...': ",
         paste(missing, collapse = ", "), ".", call. = FALSE)
  }
  if (length(setdiff(names(args), source_names))) {
    stop("Supplied fit name(s) are not referenced by the model universe: ",
         paste(setdiff(names(args), source_names), collapse = ", "), ".",
         call. = FALSE)
  }
  stats::setNames(args[source_names], source_names)
}

.mu_exact_maps <- function(fit) {
  structure_name <- .mu_structure(fit)
  if (!structure_name %in% c("fringe", "stackelberg") ||
      !identical(.mu_partition_prior(fit), "ordered_k")) {
    stop("Ordered-k universe families require an exact fringe or Stackelberg ",
         "source fit with partition_prior = 'ordered_k'.", call. = FALSE)
  }
  scope <- .mu_or(fit@metadata$k_scope,
                  attr(fit@stan_data, "k_scope", exact = TRUE))
  if (!is.null(scope) && !identical(scope, "market")) {
    stop("Model-universe k restrictions are market-level; period-level k_scope ",
         "cannot be reinterpreted as a market restriction.", call. = FALSE)
  }
  maps <- if (methods::is(fit, "BayesSpecifiedFit")) {
    .posterior_k_metadata_maps(fit, structure_name)
  } else {
    .posterior_k_stan_maps(fit@stan_data, structure_name, fit@data)
  }
  if (length(maps$market) != maps$n_unit || anyDuplicated(as.character(maps$market))) {
    stop("Ordered-k universe source must have one unique market partition unit ",
         "per market.", call. = FALSE)
  }
  maps
}

.mu_ranking_counts <- function(fit, maps) {
  ranking <- fit@metadata$firm_market_ranking
  required <- c("firm_market_market", "firm_market_group_label", "firm_market_rank")
  if (!is.list(ranking) || !all(required %in% names(ranking))) {
    stop("Strict interior k requires the source fit's stored firm-market ranking.",
         call. = FALSE)
  }
  market_index <- as.integer(ranking$firm_market_market)
  labels <- as.character(ranking$firm_market_group_label)
  ranks <- as.integer(ranking$firm_market_rank)
  if (length(market_index) != length(labels) || length(labels) != length(ranks) ||
      !length(labels) || anyNA(market_index) || anyNA(labels) || any(!nzchar(labels)) ||
      anyNA(ranks)) {
    stop("The source fit's firm-market ranking metadata is malformed.", call. = FALSE)
  }
  if (!all(c("market", "group") %in% names(fit@data))) {
    stop("Strict interior k requires observed market and firm identifiers so ",
         "the ranked firm union can be checked against the fitted data.",
         call. = FALSE)
  }
  observed_market <- as.character(fit@data$market)
  observed_group <- as.character(fit@data$group)
  if (length(observed_market) != length(observed_group) ||
      anyNA(observed_market) || anyNA(observed_group)) {
    stop("The source fit's observed market and firm identifiers are malformed.",
         call. = FALSE)
  }
  out <- integer(maps$n_unit)
  for (u in seq_len(maps$n_unit)) {
    take <- market_index == as.integer(maps$unit_market[[u]])
    if (!any(take) || anyDuplicated(labels[take]) ||
        !identical(sort(ranks[take]), seq_len(sum(take)))) {
      stop("The source fit's firm-market ranking is malformed for market '",
           maps$market[[u]], "'.", call. = FALSE)
    }
    observed <- unique(observed_group[observed_market == as.character(maps$market[[u]])])
    if (!length(observed) || !setequal(observed, labels[take])) {
      stop("The source fit's ranked firm union disagrees with observed firms in ",
           "market '", maps$market[[u]], "'.", call. = FALSE)
    }
    out[[u]] <- length(unique(labels[take]))
  }
  names(out) <- as.character(maps$market)
  out
}

.mu_validate_k_values <- function(x, label, n_union) {
  if (!is.numeric(x) || !length(x) || anyNA(x) || any(!is.finite(x)) ||
      any(x != floor(x)) || anyDuplicated(x)) {
    stop("Hypothesis '", label, "' has malformed integer k values.", call. = FALSE)
  }
  x <- as.integer(x)
  if (any(x <= 0L) || any(x >= n_union)) {
    stop("Hypothesis '", label, "' must retain strictly interior k in every ",
         "market (0 < k < the ranked firm union); canonical endpoints are ",
         "represented by MC, B, or C.", call. = FALSE)
  }
  sort(x)
}

.mu_resolve_restriction <- function(spec, fit, label) {
  maps <- .mu_exact_maps(fit)
  n_union <- .mu_ranking_counts(fit, maps)
  support <- lapply(maps$component_indices, function(ii) as.integer(maps$k[ii]))
  interior <- lapply(seq_along(n_union), function(u) {
    support[[u]][support[[u]] > 0L & support[[u]] < n_union[[u]]]
  })
  names(interior) <- as.character(maps$market)
  if (any(lengths(interior) == 0L)) {
    stop("Hypothesis '", label, "' has no strict interior k in at least one ",
         "market.", call. = FALSE)
  }
  requested <- spec$k
  if (is.null(requested) || identical(requested, "interior")) {
    allowed <- interior
  } else if (is.numeric(requested)) {
    if (maps$n_unit != 1L) {
      stop("Multimarket hypothesis '", label,
           "' must provide a named market-to-k list.", call. = FALSE)
    }
    allowed <- list(.mu_validate_k_values(requested, label, n_union[[1L]]))
    names(allowed) <- as.character(maps$market)
  } else {
    if (!is.list(requested) || is.null(names(requested)) ||
        any(!nzchar(names(requested))) || anyDuplicated(names(requested))) {
      stop("Hypothesis '", label,
           "' must use a named market-to-k list for multiple markets.",
           call. = FALSE)
    }
    unknown <- setdiff(names(requested), names(interior))
    if (length(unknown)) {
      stop("Hypothesis '", label, "' names unknown market(s): ",
           paste(unknown, collapse = ", "), ".", call. = FALSE)
    }
    allowed <- interior
    for (u in seq_along(allowed)) {
      nm <- names(interior)[[u]]
      if (nm %in% names(requested)) {
        allowed[[u]] <- .mu_validate_k_values(requested[[nm]], label,
                                              n_union[[u]])
      }
    }
  }
  for (u in seq_along(allowed)) {
    unsupported <- setdiff(allowed[[u]], support[[u]])
    if (length(unsupported)) {
      stop("Hypothesis '", label, "' requests k unsupported by source fit in ",
           "market '", names(allowed)[[u]], "': ",
           paste(unsupported, collapse = ", "), ".", call. = FALSE)
    }
    if (!length(intersect(allowed[[u]], support[[u]]))) {
      stop("Hypothesis '", label, "' retains no supported strict interior k in ",
           "market '", names(allowed)[[u]], "'.", call. = FALSE)
    }
  }
  list(
    maps = maps, n_union = n_union, allowed = allowed,
    support = support, requested = requested,
    restriction = if (maps$n_unit == 1L) allowed[[1L]] else allowed
  )
}

.mu_overlap <- function(x, y) {
  length(x) == length(y) && all(vapply(seq_along(x), function(u) {
    length(intersect(x[[u]], y[[u]])) > 0L
  }, logical(1)))
}

.mu_restriction_label <- function(resolved) {
  if (is.null(resolved)) return(NA_character_)
  paste(vapply(seq_along(resolved$allowed), function(u) {
    paste0(names(resolved$allowed)[[u]], "={",
           paste(resolved$allowed[[u]], collapse = ","), "}")
  }, character(1)), collapse = ";")
}

.mu_fit_descriptor_rows <- function(fits, descriptor) {
  rows <- descriptor$hypotheses
  out <- vector("list", length(rows))
  for (i in seq_along(rows)) {
    spec <- rows[[i]]
    fit <- fits[[spec$source]]
    structure_name <- .mu_structure(fit)
    supply <- .mu_supply(fit)
    fp <- .mu_fingerprint(fit)
    if (identical(spec$kind, "atomic")) {
      if (structure_name %in% c("fringe", "stackelberg")) {
        resolved <- .mu_resolve_restriction(
          list(k = NULL), fit, spec$hypothesis
        )
        support_is_interior <- all(vapply(seq_along(resolved$support), function(u) {
          all(resolved$support[[u]] %in% resolved$allowed[[u]]) &&
            all(resolved$support[[u]] > 0L &
                  resolved$support[[u]] < resolved$n_union[[u]])
        }, logical(1)))
        if (!support_is_interior) {
          stop("Scalar full-family hypothesis '", spec$hypothesis,
               "' is only accepted when its exact source support is already ",
               "strictly interior; use canonical MC/B/C for endpoints and ",
               "k = 'interior' for a restricted family.", call. = FALSE)
        }
        spec$kind <- "family"
        spec$k <- "interior"
        resolved$requested <- NULL
      } else if (!identical(structure_name, "simultaneous")) {
        stop("Hypothesis '", spec$hypothesis,
             "' uses an unsupported structural family.", call. = FALSE)
      } else {
        resolved <- NULL
      }
    } else {
      resolved <- .mu_resolve_restriction(spec, fit, spec$hypothesis)
    }
    if (identical(spec$kind, "family") &&
        !structure_name %in% c("fringe", "stackelberg")) {
      stop("Ordered-k family '", spec$hypothesis,
           "' must use fringe or Stackelberg structure.", call. = FALSE)
    }
    if (identical(spec$kind, "family") &&
        !supply %in% c("bertrand", "cournot")) {
      stop("Exact ordered-k family '", spec$hypothesis,
           "' supports Bertrand or Cournot conduct only.", call. = FALSE)
    }
    out[[i]] <- list(
      hypothesis = spec$hypothesis, source = spec$source, kind = spec$kind,
      fit = fit, structure = structure_name, supply = supply,
      fingerprint = fp, resolved = resolved,
      restricted_fingerprint = if (is.null(resolved)) fp else
        .mu_fingerprint(fit, restriction = resolved),
      restriction_requested = !is.null(resolved) &&
        !is.null(spec$k) && !is.null(resolved$requested)
    )
  }
  names(out) <- names(rows)

  ## The same exact atomic state or overlapping rectangle may not be encoded
  ## twice.  These checks are conditional on a matching statistical state, so
  ## different likelihoods, priors, or demand systems remain separate models.
  for (i in seq_along(out)) for (j in seq_len(i - 1L)) {
    x <- out[[i]]; y <- out[[j]]
    if (identical(x$kind, "family") && identical(y$kind, "family") &&
        identical(x$structure, y$structure) && identical(x$supply, y$supply) &&
        identical(x$restricted_fingerprint, y$restricted_fingerprint)) {
      stop("Hypotheses '", x$hypothesis, "' and '", y$hypothesis,
           "' duplicate the same conditional ordered-k prior state.", call. = FALSE)
    }
    if (!identical(x$fingerprint, y$fingerprint)) next
    if (identical(x$kind, "atomic") && identical(y$kind, "atomic") &&
        identical(x$structure, "simultaneous") &&
        identical(y$structure, "simultaneous") &&
        identical(x$supply, y$supply)) {
      stop("Hypotheses '", x$hypothesis, "' and '", y$hypothesis,
           "' duplicate the same atomic conduct state under a matching ",
           "statistical contract.", call. = FALSE)
    }
    if (identical(x$kind, "family") && identical(y$kind, "family") &&
        identical(x$structure, y$structure) && identical(x$supply, y$supply) &&
        .mu_overlap(x$resolved$allowed, y$resolved$allowed)) {
      stop("Hypotheses '", x$hypothesis, "' and '", y$hypothesis,
           "' overlap on an exact ordered-k family rectangle. Restrict each ",
           "family to a disjoint state set.", call. = FALSE)
    }
  }
  out
}

.mu_restriction_result <- function(row, draws = NULL) {
  if (is.null(row$resolved)) {
    return(list(
      prior_mass = 1, posterior_mass = 1,
      log_evidence_vs_full = 0, method = "full_atomic"
    ))
  }
  fit <- row$fit
  restrictions <- list(universe_state = row$resolved$restriction)
  result <- bayes_restriction_compare(fit, restrictions = restrictions,
                                      draws = draws)
  result[1L, c("prior_mass", "posterior_mass", "log_evidence_vs_full", "method"),
         drop = FALSE]
}

.mu_logsumexp <- function(x) {
  if (!length(x)) return(-Inf)
  finite <- is.finite(x)
  if (!any(finite)) return(-Inf)
  m <- max(x[finite])
  m + log(sum(exp(x[finite] - m)))
}

.mu_result_from_evidence <- function(rows, source_evidence, prior_weights,
                                     status_override = NULL,
                                     restriction_results = NULL) {
  n <- length(rows)
  if (!is.data.frame(source_evidence) || !all(c("model", "logml") %in% names(source_evidence))) {
    stop("Source evidence must contain model and logml columns.", call. = FALSE)
  }
  records <- vector("list", n)
  evidence_bad <- FALSE
  for (i in seq_len(n)) {
    row <- rows[[i]]
    source_row <- source_evidence[match(row$source, source_evidence$model), , drop = FALSE]
    if (!nrow(source_row)) {
      source_logml <- NA_real_; source_sd <- NA_real_
    } else {
      source_logml <- source_row$logml[[1L]]
      source_sd <- if ("logml_sd" %in% names(source_row)) source_row$logml_sd[[1L]] else NA_real_
    }
    supplied_restriction <- if (!is.null(restriction_results)) {
      restriction_results[[i]]
    } else NULL
    if (!is.null(supplied_restriction)) {
      restriction <- supplied_restriction
    } else if (is.null(row$evidence_override)) {
      restriction <- if (is.null(row$resolved)) {
        list(prior_mass = 1, posterior_mass = 1,
             log_evidence_vs_full = 0, method = "full_atomic")
      } else {
        rr <- .mu_restriction_result(row)
        as.list(rr[1L, , drop = FALSE])
      }
    } else {
      restriction <- row$evidence_override
    }
    prior_mass <- as.numeric(restriction$prior_mass)
    posterior_mass <- as.numeric(restriction$posterior_mass)
    log_ratio <- as.numeric(restriction$log_evidence_vs_full)
    if (length(prior_mass) != 1L || !is.finite(prior_mass) || prior_mass <= 0 ||
        prior_mass > 1 + 1e-8) {
      stop("Hypothesis '", row$hypothesis,
           "' has invalid or zero prior restriction mass.", call. = FALSE)
    }
    if (length(posterior_mass) != 1L || is.na(posterior_mass) ||
        !is.finite(posterior_mass) || posterior_mass < 0 || posterior_mass > 1 + 1e-8) {
      stop("Hypothesis '", row$hypothesis,
           "' has invalid posterior restriction mass.", call. = FALSE)
    }
    linear_underflow <- !is.null(row$resolved) && identical(posterior_mass, 0)
    if (is.na(source_logml) || is.nan(source_logml) ||
        !is.finite(source_logml) || is.na(log_ratio) || is.nan(log_ratio) ||
        (!is.finite(log_ratio) && !linear_underflow)) {
      evidence_bad <- evidence_bad || prior_weights[[i]] > 0
      restricted <- NA_real_
    } else if (linear_underflow) {
      ## A sampled/generated probability of exactly zero is usually a linear
      ## underflow for a positive-prior continuous region. Preserve a finite
      ## log-space ratio when available. The log-space evidence remains
      ## usable; retain the flag so downstream policy code can avoid treating
      ## the displayed linear mass as an exact mathematical zero.
      restricted <- if (is.finite(log_ratio)) source_logml + log_ratio else NA_real_
    } else {
      restricted <- source_logml + log_ratio
    }
    if (is.nan(restricted) || is.infinite(restricted)) {
      evidence_bad <- evidence_bad || prior_weights[[i]] > 0
      if (is.finite(source_logml) && is.finite(log_ratio) &&
          !linear_underflow) restricted <- NA_real_
    }
    records[[i]] <- data.frame(
      hypothesis = row$hypothesis, model = row$hypothesis,
      source = row$source, kind = row$kind,
      supply_model = row$supply, structure = row$structure,
      statistical_fingerprint = row$fingerprint,
      k_restriction = .mu_restriction_label(row$resolved),
      prior_weight = prior_weights[[i]], prior_mass = prior_mass,
      posterior_mass = posterior_mass, logml_source = source_logml,
      logml_restricted = restricted, log_evidence_vs_source = log_ratio,
      logml_sd = source_sd,
      restriction_method = as.character(.mu_or(restriction$method, "encompassing_prior")),
      evidence_status = if (is.na(restricted)) "bridge_failed" else
        if (linear_underflow) "restriction_mass_underflow" else
        if (is.finite(restricted)) "complete" else "zero_restricted_mass",
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, records)
  out$direction <- vapply(out$supply_model, function(supply) switch(supply,
    bertrand = "price", cournot = "quantity", auction = "auction",
    monopolistic = "moncom", bargaining = "bargaining", "other"
  ), character(1))
  if (!is.null(status_override)) evidence_bad <- evidence_bad || status_override
  if (evidence_bad) {
    out$posterior_probability <- NA_real_
    status <- "incomplete_evidence"
  } else {
    log_weight <- log(out$prior_weight) + out$logml_restricted
    log_weight[out$prior_weight == 0] <- -Inf
    normalizer <- .mu_logsumexp(log_weight)
    if (!is.finite(normalizer)) {
      out$posterior_probability <- NA_real_
      status <- "incomplete_evidence"
    } else {
      out$posterior_probability <- exp(log_weight - normalizer)
      out$posterior_probability[out$prior_weight == 0] <- 0
      status <- "complete"
    }
  }
  rownames(out) <- NULL
  list(table = out, status = status)
}

.bayes_compare_universe <- function(args, universe, bridge, bridge_maxiter,
                                    bridge_repetitions, bridge_method, loo_unit) {
  descriptor <- .mu_validate_descriptor(universe)
  if (!isTRUE(bridge)) {
    stop("Model-universe comparison requires bridge = TRUE so every retained ",
         "hypothesis has normalized evidence.", call. = FALSE)
  }
  fits <- .mu_source_fits(args, descriptor)
  invisible(lapply(fits, validate_fit_spec))
  rows <- .mu_fit_descriptor_rows(fits, descriptor)
  ## Bridge each source fit exactly once.  A source may back several disjoint
  ## atomic restrictions, but unreferenced supplied fits are not part of this
  ## universe and must not affect its evidence contract.
  ordinary_args <- c(fits, list(
    restrictions = NULL, reference = NULL, restriction_draws = NULL,
    bridge = bridge, bridge_maxiter = bridge_maxiter,
    bridge_repetitions = bridge_repetitions, bridge_method = bridge_method,
    loo_unit = loo_unit,
    prior_weights = NULL, universe = NULL
  ))
  source_evidence <- do.call(bayes_compare, ordinary_args)
  result <- .mu_result_from_evidence(rows, source_evidence,
                                     descriptor$prior_weights)
  table <- result$table
  ## Keep the ordinary bayes_compare table contract available to downstream
  ## tools such as bayes_weights(), while retaining universe-specific fields.
  table$model_prior <- table$prior_weight
  table$logml <- table$logml_restricted
  table$post_prob <- table$posterior_probability
  table$bf_vs_best <- NA_real_
  finite <- is.finite(table$logml)
  if (any(finite)) table$bf_vs_best[finite] <- exp(
    table$logml[finite] - max(table$logml[finite])
  )
  attr(table, "status") <- result$status
  attr(table, "descriptor") <- descriptor
  attr(table, "source_comparison") <- source_evidence
  attr(table, "diagnostics") <- list(
    evidence_identity = "Z_restricted = Z_source * posterior_mass / prior_mass",
    multimarket_mass = "within_draw_market_product_then_draw_average",
    bridge = TRUE
  )
  attr(table, "restriction_details") <- lapply(rows, function(row) list(
    hypothesis = row$hypothesis, source = row$source,
    allowed_k = if (is.null(row$resolved)) NULL else row$resolved$allowed,
    n_ranked_firm_union = if (is.null(row$resolved)) NULL else row$resolved$n_union,
    requested_k = if (is.null(row$resolved)) NULL else row$resolved$requested
  ))
  attr(table, "prior_weights") <- descriptor$prior_weights
  class(table) <- c("bayes_universe_comparison", class(table))
  table
}

.mu_dimension_for_row <- function(row) {
  conduct <- switch(row$supply_model,
    bertrand = "price", cournot = "quantity", auction = "auction",
    monopolistic = "monopolistic", bargaining = "bargaining", "other"
  )
  core <- if (row$supply_model == "monopolistic") "none" else if (
    row$supply_model %in% c("bertrand", "cournot") && row$kind == "family"
  ) "interior" else if (row$supply_model %in% c("bertrand", "cournot") &&
                         row$structure == "simultaneous") "all" else "not_applicable"
  list(conduct = conduct, structure = row$structure, strategic_core = core)
}

#' Aggregate a model-universe posterior by scientific dimensions
#'
#' @param comparison A \code{bayes_universe_comparison} returned by
#'   \code{\link{bayes_compare}} with \code{universe}.
#' @param dimension One or more of \code{"conduct"}, \code{"structure"}, or
#'   \code{"strategic_core"}.
#' @return A data frame of exact sums of atomic posterior probabilities.
#'   Incomplete evidence returns the requested levels with missing probabilities.
#' @export
bayes_dimensions <- function(comparison,
                             dimension = c("conduct", "structure", "strategic_core")) {
  if (!inherits(comparison, "bayes_universe_comparison")) {
    stop("'comparison' must be a bayes_universe_comparison.", call. = FALSE)
  }
  dimension <- match.arg(dimension, several.ok = TRUE)
  tab <- comparison
  maps <- lapply(seq_len(nrow(tab)), function(i) .mu_dimension_for_row(tab[i, , drop = FALSE]))
  levels_by_dimension <- list(
    conduct = c("price", "quantity", "auction", "monopolistic", "bargaining", "other"),
    structure = c("simultaneous", "fringe", "stackelberg"),
    strategic_core = c("none", "interior", "all", "not_applicable")
  )
  rows <- lapply(dimension, function(d) {
    assigned <- vapply(maps, `[[`, character(1), d)
    p <- tab$posterior_probability
    prior <- tab$prior_weight
    complete <- identical(attr(comparison, "status"), "complete") &&
      all(is.finite(p))
    post <- prior_sum <- numeric(length(levels_by_dimension[[d]]))
    for (j in seq_along(post)) {
      take <- assigned == levels_by_dimension[[d]][[j]]
      prior_sum[[j]] <- sum(prior[take])
      post[[j]] <- if (complete) sum(p[take]) else NA_real_
    }
    data.frame(
      dimension = d, level = levels_by_dimension[[d]],
      direction = if (d == "conduct") ifelse(
        levels_by_dimension[[d]] == "monopolistic", "moncom",
        levels_by_dimension[[d]]
      ) else levels_by_dimension[[d]],
      posterior_probability = post, prior_probability = prior_sum,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  attr(out, "status") <- attr(comparison, "status")
  attr(out, "strategic_core_event") <-
    "interior means a strict k<n state in every market of that exact family; it is a global all-market event."
  out
}
