## Fit-level draw plans and demand replacement for latent-quality simulation.
## The observed-share route remains in simulate.R; these helpers are called
## only when demand_anchoring = "latent_quality" (or an explicit draw plan is
## supplied).

.lq_or <- function(x, y) if (is.null(x) || !length(x)) y else x[[1L]]

.lq_with_seed <- function(seed, fun) {
  if (is.null(seed)) return(fun())
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(as.integer(seed))
  fun()
}

.lq_validate_seed <- function(seed) {
  if (!is.null(seed) && (length(seed) != 1L || !is.numeric(seed) ||
      !is.finite(seed) || seed != as.integer(seed) || seed < 0 ||
      seed > .Machine$integer.max)) {
    stop("'seed' must be a nonnegative integer.", call. = FALSE)
  }
  if (is.null(seed)) NULL else as.integer(seed)
}

.lq_draw_matrix <- function(x, n_draw, n_col, name) {
  if (is.null(x)) return(NULL)
  d <- dim(x)
  if (is.null(d)) {
    if (n_col == 1L && length(x) == n_draw) {
      return(matrix(as.numeric(x), nrow = n_draw, ncol = 1L))
    }
    if (n_draw == 1L && length(x) == n_col) {
      return(matrix(as.numeric(x), nrow = 1L, ncol = n_col))
    }
    if (length(x) == n_draw * n_col) {
      return(matrix(as.numeric(x), nrow = n_draw, ncol = n_col, byrow = TRUE))
    }
    stop("Generated quantity '", name, "' has malformed dimensions.",
         call. = FALSE)
  }
  if (length(d) != 2L || d[[1L]] != n_draw || d[[2L]] != n_col) {
    stop("Generated quantity '", name, "' must have draws by columns dimensions.",
         call. = FALSE)
  }
  matrix(as.numeric(x), nrow = n_draw, ncol = n_col)
}

.lq_group_levels <- function(fit) {
  if (methods::is(fit, "BayesSpecifiedFit")) {
    return(sort(unique(as.character(fit@market$group))))
  }
  labels <- attr(fit@stan_data, "group_levels", exact = TRUE)
  if (is.null(labels)) labels <- sort(unique(as.character(fit@data$group)))
  labels <- as.character(labels)
  if (!length(labels) || anyNA(labels) || any(!nzchar(labels)) ||
      anyDuplicated(labels)) {
    stop("The fitted group-level map is malformed.", call. = FALSE)
  }
  labels
}

.lq_quality_scope <- function(fit) {
  if (methods::is(fit, "BayesSpecifiedFit")) return("global")
  .quality_scope_from_fit(fit, strict = FALSE)
}

.lq_quality_incidence <- function(fit) {
  if (methods::is(fit, "BayesSpecifiedFit")) return(NULL)
  .quality_incidence_from_fit(fit, strict = FALSE)
}

.lq_quality_mode <- function(fit) {
  if (methods::is(fit, "BayesSpecifiedFit")) return("none")
  stan_data <- fit@stan_data
  value <- .lq_or(fit@metadata$quality_parameterization,
                  attr(stan_data, "quality_parameterization", exact = TRUE))
  if (is.null(value)) value <- stan_data$quality_parameterization
  if (is.numeric(value) && length(value) == 1L) {
    value <- if (value == 1) "sampled" else if (value == 2) "collapsed" else NULL
  }
  if (is.character(value) && length(value) == 1L && value %in% c("sampled", "collapsed")) return(value)
  if (isTRUE(fit@metadata$bayes_is_collapsed_quality)) return("collapsed")
  stop("Unable to resolve the fitted quality parameterization.", call. = FALSE)
}

.lq_market_labels <- function(fit) {
  if (methods::is(fit, "BayesSpecifiedFit")) return(character())
  labels <- attr(fit@stan_data, "market_levels", exact = TRUE) %||%
    fit@metadata$market_levels %||%
    if ("market" %in% names(fit@data))
      sort(unique(as.character(fit@data$market))) else character()
  labels <- as.character(labels)
  n_market <- as.integer(fit@stan_data$N_market %||% length(labels))
  if (length(labels) != n_market || !length(labels) || anyNA(labels) ||
      any(!nzchar(labels)) || anyDuplicated(labels)) {
    stop("The fitted market-level map is malformed.", call. = FALSE)
  }
  labels
}

.lq_draw_count <- function(fit, draws = NULL) {
  if (methods::is(fit, "BayesSpecifiedFit")) return(as.integer(fit@params$n_draws))
  if (is.null(draws)) draws <- .posterior_extract(fit, "simulation")
  .posterior_draw_count_from_draws(draws, fit)
}

.lq_subset_draws <- function(draws, indices, n_total) {
  ## The package draw contract uses vectors, matrices, and higher-dimensional arrays with the
  ## draw dimension first. Subsetting only that dimension avoids constructing
  ## conditional quality draws for every available posterior row when a plan
  ## requests a small reusable subset.
  lapply(draws, function(x) {
    d <- dim(x)
    if (is.null(d)) {
      if (length(x) == n_total) x[indices] else x
    } else if (length(d) >= 1L && d[[1L]] == n_total) {
      idx <- c(list(indices), rep(list(TRUE), length(d) - 1L))
      do.call("[", c(list(x), idx, list(drop = FALSE)))
    } else {
      x
    }
  })
}

## The policy simulator stores one role state for every selected simulation
## unit (market-year in a panel), while exact ordered-k quality conditioning
## needs one state for every partition market.  Convert the former to the
## latter without changing the simulator-facing cache.  Shared market-scope
## role draws must agree across all years in a market; a period-scoped state
## therefore fails explicitly instead of silently selecting one year.
.lq_market_quality_k_cache <- function(fit, incidence, k_cache) {
  if (is.null(k_cache))
    return(k_cache)
  scope <- as.character(incidence$quality_scope)
  if (!scope %in% c("market", "global")) return(k_cache)
  setup <- .lq_all_units(fit)
  info_market <- vapply(setup$info, function(z) as.integer(z$market_index), integer(1))
  n_market <- as.integer(fit@stan_data$N_market %||% max(info_market))
  ## Global collapsed quality can use the exact-k/rho joint conditional only
  ## for one market. It still needs the same market-level role state when that
  ## market contains several years.
  if (scope == "global" && n_market != 1L) return(k_cache)
  if (length(info_market) != length(setup$info) || anyNA(info_market) ||
      any(info_market < 1L | info_market > n_market)) {
    stop("The exact-k quality cache cannot map simulation units to markets.", call. = FALSE)
  }
  out <- lapply(seq_along(k_cache), function(draw_pos) {
    vals <- numeric(n_market)
    for (m in seq_len(n_market)) {
      unit_pos <- which(info_market == m)
      if (!length(unit_pos))
        stop("The exact-k quality cache is missing market ", m, ".", call. = FALSE)
      vv <- vapply(unit_pos, function(pos) {
        .lq_plan_k_for_unit(list(k_cache = k_cache), draw_pos, setup$info[[pos]])
      }, integer(1))
      if (any(vv != vv[[1L]]))
        stop("Exact-k quality conditioning requires one exact-k state across all years in market ", m, ".", call. = FALSE)
      vals[[m]] <- vv[[1L]]
    }
    vals
  })
  names(out) <- names(k_cache)
  out
}

.lq_quality_draws <- function(fit, draws, draw_indices, seed = NULL, k_cache = NULL) {
  groups <- .lq_group_levels(fit)
  n_draw <- .lq_draw_count(fit, draws)
  n_group <- length(groups)
  mode <- .lq_quality_mode(fit)
  incidence <- .lq_quality_incidence(fit)
  scope <- .lq_quality_scope(fit)
  units <- if (is.null(incidence)) groups else as.character(incidence$quality_lookup$quality_unit)
  n_quality <- length(units)
  if (mode == "none") {
    out <- matrix(0, nrow = length(draw_indices), ncol = n_quality,
                  dimnames = list(as.character(draw_indices), units))
    return(list(draws = out, source = "none", conditioning = "none",
                mode = mode, group_labels = groups, quality_units = units,
                quality_scope = scope, quality_lookup = incidence$quality_lookup,
                incidence = incidence, market_intercept_draws = NULL,
                market_labels = NULL))
  }
  if (mode == "sampled") {
    q <- .lq_draw_matrix(draws$b_quality_all, n_draw, n_quality,
                         "b_quality_all")
    if (is.null(q)) {
      if (n_quality == 1L) q <- matrix(0, nrow = n_draw, ncol = 1L)
      else stop("Sampled-quality fit is missing generated quantity 'b_quality_all'.",
                call. = FALSE)
    }
    q <- q[draw_indices, , drop = FALSE]
    block_ok <- if (is.null(incidence)) abs(rowSums(q)) else {
      vapply(seq_len(incidence$N_quality_block), function(b) {
        abs(rowSums(q[, incidence$quality_block_start[[b]]:
          (incidence$quality_block_start[[b]] + incidence$quality_block_n[[b]] - 1L), drop = FALSE]))
      }, numeric(nrow(q)))
    }
    if (any(!is.finite(q)) || any(block_ok > 1e-8)) {
      stop("Sampled quality draws are non-finite or violate per-block centering.",
           call. = FALSE)
    }
    colnames(q) <- units
    rownames(q) <- as.character(draw_indices)
    return(list(draws = q, source = "sampled_stan", conditioning = "joint_draw",
                mode = mode, group_labels = groups, quality_units = units,
                quality_scope = scope, quality_lookup = incidence$quality_lookup,
                incidence = incidence, market_intercept_draws = NULL,
                market_labels = NULL))
  }
  ## posterior_quality() supplies the exact Gaussian conditional draw. Restrict
  ## every draw-indexed field first, then request all rows of that subset. This
  ## retains the selected original theta rows without allocating an all-draw
  ## quality matrix or performing a second independent draw selection.
  selected_draws <- .lq_subset_draws(draws, draw_indices, n_draw)
  selected_k_cache <- if (is.null(k_cache)) NULL else {
    z <- k_cache[as.character(draw_indices)]
    if (length(z) != length(draw_indices) || any(vapply(z, is.null, logical(1))))
      stop("The quality draw cache is missing an exact-k state for a selected draw.", call. = FALSE)
    names(z) <- as.character(seq_along(draw_indices)); z
  }
  selected_k_cache <- .lq_market_quality_k_cache(fit, incidence, selected_k_cache)
  quality_fit <- fit
  ## posterior_quality() deliberately checks the fit's recorded draw count.
  ## The shallow clone changes only that validation field; provenance in the
  ## returned plan continues to carry the original posterior indices and fit
  ## fingerprint.
  if (!is.null(quality_fit@metadata$posterior_draws)) {
    quality_fit@metadata$posterior_draws <- length(draw_indices)
  }
  pq <- posterior_quality(quality_fit, n_draws = length(draw_indices), seed = seed,
                          draws = selected_draws, conditioning = "auto",
                          k_cache = selected_k_cache)
  if (!identical(as.character(pq$group_labels), groups) ||
      !identical(as.character(pq$quality_scope), scope)) {
    stop("Conditional quality labels do not match the fitted group map.", call. = FALSE)
  }
  q <- pq$draws
  if (nrow(q) != length(draw_indices)) {
    stop("Conditional quality recovery returned an unexpected draw count.", call. = FALSE)
  }
  if (any(!is.finite(q))) {
    stop("Conditional quality draws are non-finite.", call. = FALSE)
  }
  dimnames(q) <- list(as.character(draw_indices), units)
  market_intercept_draws <- pq$market_intercept_draws
  market_labels <- pq$market_intercept_labels
  if (identical(.market_intercept_parameterization_from_fit(fit), "collapsed")) {
    if (!is.matrix(market_intercept_draws) ||
        nrow(market_intercept_draws) != length(draw_indices) ||
        ncol(market_intercept_draws) != length(market_labels) ||
        any(!is.finite(market_intercept_draws))) {
      stop("Joint market-intercept recovery returned malformed draws.",
           call. = FALSE)
    }
    rownames(market_intercept_draws) <- as.character(draw_indices)
    colnames(market_intercept_draws) <- as.character(market_labels)
  } else {
    market_intercept_draws <- NULL
    market_labels <- NULL
  }
  list(draws = q, source = "posterior_quality", conditioning =
         pq$metadata$conditioning, mode = mode, group_labels = groups,
       quality_units = units, quality_scope = scope,
       quality_lookup = incidence$quality_lookup, incidence = incidence,
       realized_k = pq$realized_k,
       market_intercept_draws = market_intercept_draws,
       market_labels = market_labels)
}

.lq_all_units <- function(fit) {
  if (methods::is(fit, "BayesSpecifiedFit")) {
    return(.specified_simulation_units(fit, NULL, NULL))
  }
  years <- unique(as.character(fit@data$year))
  if (length(years) > 1L) .estimated_simulation_units(fit, NULL, years)
  else .estimated_simulation_units(fit, NULL, NULL)
}

.lq_role_cache <- function(fit, draws, draw_indices,
                           role_mapping = "mixture", fixed_k = NULL,
                           role_offset = 0L) {
  if (is.null(.native_ordered_structure(fit))) return(NULL)
  setup <- .lq_all_units(fit)
  .stack_k_cache(fit, draws, setup, draw_indices,
                 role_mapping = role_mapping, fixed_k = fixed_k,
                 role_offset = role_offset)
}

## A single-year plan uses market IDs as unit keys; a caller that explicitly
## selects that year uses market_year IDs. Preserve the same realized k.
.lq_plan_k_for_unit <- function(plan, draw_pos, unit) {
  state <- plan$k_cache[[draw_pos]]
  if (is.null(state) || is.null(names(state))) {
    stop("The exact-k draw plan has no named role state.", call. = FALSE)
  }
  key <- if (unit$unit_id %in% names(state)) unit$unit_id else
    if (unit$market_id %in% names(state)) unit$market_id else NA_character_
  if (is.na(key)) {
    stop("The exact-k draw plan has no role state for simulation unit ",
         unit$unit_id, ".", call. = FALSE)
  }
  value <- state[[key]]
  if (length(value) != 1L || is.na(value) || !is.finite(value) ||
      value < 0L || value != as.integer(value)) {
    stop("The exact-k draw plan has an invalid role state for simulation unit ",
         unit$unit_id, ".", call. = FALSE)
  }
  as.integer(value)
}

.lq_validate_role_plan_args <- function(fit, role_mapping, fixed_k,
                                        role_offset) {
  role_mapping <- match.arg(role_mapping, c("mixture", "map", "threshold"))
  if (length(role_offset) != 1L || !is.numeric(role_offset) ||
      !is.finite(role_offset) || abs(role_offset) > .Machine$integer.max ||
      role_offset != as.integer(role_offset)) {
    stop("'role_offset' must be an integer scalar.", call. = FALSE)
  }
  role_offset <- as.integer(role_offset)
  native_structure <- .native_ordered_structure(fit)
  if (is.null(native_structure) &&
      (!identical(role_mapping, "mixture") || !is.null(fixed_k) ||
       !identical(role_offset, 0L))) {
    stop("Exact role controls require a Stackelberg or exact ordered-k fringe fit.",
         call. = FALSE)
  }
  if (!is.null(fixed_k)) {
    if (!is.numeric(fixed_k) || !length(fixed_k) ||
        any(!is.finite(fixed_k)) ||
        any(abs(fixed_k) > .Machine$integer.max) ||
        any(fixed_k != as.integer(fixed_k)) ||
        (length(fixed_k) > 1L &&
         (is.null(names(fixed_k)) || anyNA(names(fixed_k)) ||
          any(!nzchar(names(fixed_k))) || anyDuplicated(names(fixed_k))))) {
      stop("'fixed_k' must be an integer scalar or uniquely named integer vector.",
           call. = FALSE)
    }
    fixed_k_names <- names(fixed_k)
    fixed_k <- as.integer(fixed_k)
    names(fixed_k) <- fixed_k_names
  }
  list(role_mapping = role_mapping, fixed_k = fixed_k,
       role_offset = role_offset)
}

.lq_fixed_k_equal <- function(x, y) {
  if (is.null(x) || is.null(y)) return(is.null(x) && is.null(y))
  if (!is.numeric(x) || !is.numeric(y) || !length(x) || !length(y) ||
      any(!is.finite(x)) || any(!is.finite(y)) ||
      any(abs(x) > .Machine$integer.max) ||
      any(abs(y) > .Machine$integer.max) || any(x != as.integer(x)) ||
      any(y != as.integer(y))) return(FALSE)
  identical(as.integer(x), as.integer(y)) &&
    identical(names(x), names(y))
}

.lq_fit_fingerprint <- function(fit, draws = NULL, version = 2L,
                                omit_redundant_quality_alias = FALSE) {
  if (length(version) != 1L || !is.numeric(version) ||
      !is.finite(version) || version != as.integer(version) ||
      !version %in% c(1L, 2L, 3L, 4L)) {
    stop("Unsupported latent-quality fit fingerprint version.", call. = FALSE)
  }
  version <- as.integer(version)
  if (version == 1L && !methods::is(fit, "BayesSpecifiedFit") && is.null(draws)) {
    if (is.null(fit@stanfit)) stop("Cannot validate a legacy plan without its Stan fit.", call. = FALSE)
    draws <- .bayes_backend_extract(fit@stanfit)
  }
  if (version %in% c(2L, 3L, 4L) && !methods::is(fit, "BayesSpecifiedFit") && is.null(draws)) {
    draws <- .posterior_extract(fit, "policy_plan")
  }
  canonical_draws <- if (version == 1L || methods::is(fit, "BayesSpecifiedFit")) {
    draws
  } else {
    .posterior_canonical_draws(draws, fit, "policy_plan")
  }
  if (isTRUE(omit_redundant_quality_alias)) {
    if (!version %in% c(2L, 3L))
      stop("Historical quality-alias fingerprints require version 2 or 3.", call. = FALSE)
    if (!is.null(canonical_draws$sigma_b_strat_economic) &&
        identical(canonical_draws$sigma_b_strat_economic, canonical_draws$sigma_b_strat))
      canonical_draws$sigma_b_strat_economic <- NULL
  }
  payload <- list(
    data = fit@data,
    demand_model = as.character(fit@demand_model),
    market_side = as.character(fit@market_side),
    supply_model = as.integer(fit@supply_model),
    metadata = fit@metadata,
    stan_data = if (methods::is(fit, "BayesSpecifiedFit")) NULL else fit@stan_data,
    params = if (methods::is(fit, "BayesSpecifiedFit")) fit@params else NULL,
    market = if (methods::is(fit, "BayesSpecifiedFit")) fit@market else NULL,
    draws = if (methods::is(fit, "BayesSpecifiedFit")) NULL else canonical_draws
  )
  if (version %in% c(2L, 3L, 4L)) payload$fit_signature_version <- as.integer(version)
  digest::digest(payload, algo = "xxhash64", serialize = TRUE)
}

## Early v2/v3 canonical payloads did not retain the redundant economic alias.
## Accept that exact serialization variant, while keeping every economic draw,
## fitted source/data field and the separately checked latent state hash fixed.
.lq_fit_fingerprint_matches <- function(fit, draws, fingerprint, version) {
  current <- .lq_fit_fingerprint(fit, draws, version = version)
  if (identical(as.character(fingerprint), current)) return(TRUE)
  version %in% c(2L, 3L) && identical(as.character(fingerprint),
    .lq_fit_fingerprint(fit, draws, version = version,
                        omit_redundant_quality_alias = TRUE))
}

.lq_plan_signature <- function(fit, groups, n_total, mode, draws = NULL) {
  incidence <- .lq_quality_incidence(fit)
  scope <- .lq_quality_scope(fit)
  list(
    fit_signature_version = if (.posterior_demand_shock_model(fit)) 4L else 3L,
    demand_model = as.character(fit@demand_model),
    market_side = as.character(fit@market_side),
    posterior_draws = as.integer(n_total),
    group_labels = as.character(groups),
    quality_scope = scope,
    quality_units = if (is.null(incidence)) as.character(groups) else
      as.character(incidence$quality_lookup$quality_unit),
    quality_lookup = if (is.null(incidence)) NULL else incidence$quality_lookup,
    quality_mode = mode,
    structure = as.character(.lq_or(fit@metadata$structure, "simultaneous")),
    fit_fingerprint = .lq_fit_fingerprint(fit, draws, version =
      if (.posterior_demand_shock_model(fit)) 4L else 3L)
  )
}

.lq_plan_state_fingerprint <- function(plan) {
  ## Unversioned/v1 plans predate quality scope. Preserve their historical
  ## state hash rather than adding absent scoped fields to the payload.
  version <- if (is.list(plan$fit_signature) &&
                 !is.null(plan$fit_signature$fit_signature_version))
    as.integer(plan$fit_signature$fit_signature_version) else 1L
  keys <- c("draw_indices", "quality_draws", "group_labels",
            "quality_source", "quality_conditioning", "quality_mode",
            "k_cache", "role_mapping", "fixed_k", "role_offset")
  if (identical(version, 2L)) {
    keys <- append(keys, c("quality_scope", "quality_units", "quality_lookup"),
                   after = 3L)
  } else if (version %in% c(3L, 4L)) {
    keys <- append(keys, c("quality_scope", "quality_units", "quality_lookup"),
                   after = 3L)
    keys <- c(keys, "market_intercept_draws", "market_labels")
  }
  if (identical(version, 4L)) keys <- c(keys, "demand_shock_draws",
    "measurement_sd_draws", "demand_shock_sd_draws", "demand_shock_correlation_draws",
    "outside_share_draws", "share_residual_model")
  if (identical(plan$share_basis, "relative"))
    keys <- c(keys, "share_basis", "derived_market_year_intercept_draws", "market_year_labels")
  digest::digest(plan[keys], algo = "xxhash64", serialize = TRUE)
}

## Internal draw-plan constructor. `draws` is the policy-plan contract union;
## the public wrapper below performs one selected extraction and delegates here.
.lq_make_policy_plan <- function(fit, draws, n_draws = NULL, seed = NULL,
                                role_mapping = c("mixture", "map", "threshold"),
                                fixed_k = NULL, role_offset = 0L) {
  if (!methods::is(fit, "BayesFit")) stop("'fit' must be a BayesFit object.", call. = FALSE)
  role_args <- .lq_validate_role_plan_args(fit, role_mapping, fixed_k,
                                          role_offset)
  seed <- .lq_validate_seed(seed)
  n_total <- .lq_draw_count(fit, draws)
  if (is.null(n_draws)) n_draws <- min(200L, n_total)
  if (length(n_draws) != 1L || !is.finite(n_draws) || n_draws < 1L ||
      n_draws != as.integer(n_draws)) stop("'n_draws' must be a positive integer.", call. = FALSE)
  n_draws <- min(as.integer(n_draws), n_total)
  draw_indices <- .lq_with_seed(seed, function() {
    if (n_total == 1L) 1L else sample.int(n_total, n_draws)
  })
  role_seed <- if (is.null(seed)) NULL else as.integer((as.double(seed) + 2) %% .Machine$integer.max)
  k_cache <- .lq_with_seed(role_seed, function() {
    .lq_role_cache(fit, draws, draw_indices,
                   role_mapping = role_args$role_mapping,
                   fixed_k = role_args$fixed_k,
                   role_offset = role_args$role_offset)
  })
  quality_seed <- if (is.null(seed)) NULL else as.integer((as.double(seed) + 1) %% .Machine$integer.max)
  quality <- .lq_quality_draws(fit, draws, draw_indices, seed = quality_seed,
                               k_cache = k_cache)
  plan <- structure(list(
    draw_indices = as.integer(draw_indices), quality_draws = quality$draws,
    market_intercept_draws = quality$market_intercept_draws,
    market_labels = quality$market_labels,
    quality_source = quality$source, quality_conditioning = quality$conditioning,
    quality_mode = quality$mode, group_labels = quality$group_labels,
    quality_scope = quality$quality_scope, quality_units = quality$quality_units,
    quality_lookup = quality$quality_lookup,
    k_cache = k_cache, role_mapping = role_args$role_mapping,
    fixed_k = role_args$fixed_k, role_offset = role_args$role_offset,
    fit_signature = .lq_plan_signature(fit, quality$group_labels, n_total,
                                       quality$mode, draws = draws),
    metadata = list(
      posterior_draws_available = as.integer(n_total),
      posterior_draws_requested = as.integer(length(draw_indices)),
      draw_indices = as.integer(draw_indices), quality_source = quality$source,
      quality_conditioning = quality$conditioning,
      quality_parameterization = quality$mode, group_labels = quality$group_labels,
      quality_scope = quality$quality_scope, quality_units = quality$quality_units,
      quality_lookup = quality$quality_lookup,
      market_intercept_parameterization =
        .market_intercept_parameterization_from_fit(fit),
      market_labels = quality$market_labels,
      role_mapping = role_args$role_mapping, fixed_k = role_args$fixed_k,
      role_offset = role_args$role_offset, exact_k = !is.null(k_cache)
    )
  ), class = c("bayes_policy_draw_plan", "list"))
  if (.posterior_demand_shock_model(fit)) {
    shock_seed <- if (is.null(seed)) NULL else
      as.integer((as.double(seed) + 4) %% .Machine$integer.max)
    plan <- .lq_add_demand_shocks(fit, draws, plan, shock_seed)
  }
  if (identical(.share_basis_from_fit(fit), "relative")) {
    derived <- .lq_relative_derived_intercepts(fit, draws, draw_indices, quality$draws,
      demand_shock_draws = plan$demand_shock_draws)
    plan$share_basis <- "relative"
    plan$derived_market_year_intercept_draws <- derived$draws
    plan$market_year_labels <- derived$labels
    plan$metadata$share_basis <- "relative"
    plan$metadata$observation_basis <- "ilr_quantity_shares_v1"
    plan$metadata$market_year_labels <- derived$labels
  }
  plan$fit_signature$latent_state_fingerprint <- .lq_plan_state_fingerprint(plan)
  plan
}

#' Build a reusable posterior draw plan for policy simulation
#'
#' The plan samples posterior indices once and, for collapsed fits, draws one
#' conditional quality vector per selected parameter draw, centered globally
#' or within each market according to the fit's quality scope. Exact native
#' ordered-k fits also receive one role realization per market and draw. Reuse
#' the returned object in observed-anchor and latent-quality simulations.
#' Demand-shock fits additionally recover observation-level period demand
#' conditional on the same quality and role realization. Policy demand retains
#' quality plus period demand; independent measurement error is excluded.
#'
#' @param fit A \code{BayesFit}.
#' @param n_draws Number of posterior draws, defaulting to at most 200.
#' @param seed Optional nonnegative integer random seed.
#' @param role_mapping Exact role realization: posterior mixture, posterior
#'   MAP, or legacy cutoff threshold where supported by the fit. Exact
#'   ordered-k fits store the selected role state in the plan.
#' @param fixed_k Optional exact leader/core count, scalar or named by
#'   simulation unit. When supplied, every selected draw uses this count.
#' @param role_offset Integer adjacent-partition sensitivity offset.
#' @param allow_experimental Explicitly allow exploratory Pathfinder draws
#'   for policy work. Their approximation diagnostics remain attached to the
#'   originating fit and do not certify posterior accuracy.
#' @return A reusable object of class \code{bayes_policy_draw_plan}.
#' @details Demand-shock plans bind a version-4 fingerprint that includes
#' recovered period demand, measurement/demand scales and outside shares.
#' Combined-model plans retain version 3. Both bind a canonical, name-ordered
#' economic draw fingerprint, including reconstructed market-intercept states when those
#' effects were integrated during fitting. Unversioned historical plans are validated through their
#' version-1 unrestricted-draw fingerprint as an explicit compatibility path;
#' they are not reinterpreted as current plans. Reuse a plan with the same
#' \code{role_mapping}, \code{fixed_k}, and \code{role_offset} arguments used to
#' create it; simulation rejects conflicting controls.
#' @export
posterior_policy_draws <- function(fit, n_draws = NULL, seed = NULL,
                                   role_mapping = c("mixture", "map", "threshold"),
                                   fixed_k = NULL, role_offset = 0L,
                                   allow_experimental = FALSE) {
  if (!methods::is(fit, "BayesFit")) stop("'fit' must be a BayesFit object.", call. = FALSE)
  .cmdstan_require_policy_inference(fit, allow_experimental)
  draws <- if (methods::is(fit, "BayesSpecifiedFit")) NULL else
    .posterior_extract(fit, "policy_plan")
  .lq_make_policy_plan(fit, draws, n_draws = n_draws, seed = seed,
                       role_mapping = role_mapping, fixed_k = fixed_k,
                       role_offset = role_offset)
}

.lq_validate_plan <- function(fit, plan, n_draws = NULL,
                              role_mapping = "mixture", fixed_k = NULL,
                              role_offset = 0L, draws = NULL) {
  if (!inherits(plan, "bayes_policy_draw_plan")) {
    stop("'draw_plan' must be created by posterior_policy_draws().", call. = FALSE)
  }
  groups <- .lq_group_levels(fit)
  incidence <- .lq_quality_incidence(fit)
  scope <- .lq_quality_scope(fit)
  units <- if (is.null(incidence)) groups else as.character(incidence$quality_lookup$quality_unit)
  signature_version <- if (is.list(plan$fit_signature) &&
                           !is.null(plan$fit_signature$fit_signature_version))
    as.integer(plan$fit_signature$fit_signature_version) else 1L
  if (signature_version == 1L && !methods::is(fit, "BayesSpecifiedFit")) {
    ## Explicit compatibility path: v1 hashes the old unrestricted extract.
    if (is.null(fit@stanfit)) stop("Cannot validate a legacy plan without its Stan fit.", call. = FALSE)
    draws <- .bayes_backend_extract(fit@stanfit)
  } else if (!methods::is(fit, "BayesSpecifiedFit") && is.null(draws)) {
    draws <- .posterior_extract(fit, "policy_plan")
  }
  n_total <- .lq_draw_count(fit, draws)
  sig <- plan$fit_signature
  legacy_scope_missing <- is.list(sig) && is.null(sig$quality_scope) && signature_version == 1L
  legacy_units_missing <- is.list(sig) && is.null(sig$quality_units) && signature_version == 1L
  if (!is.list(sig) || !identical(as.integer(sig$posterior_draws), as.integer(n_total)) ||
      !identical(as.character(sig$group_labels), groups) ||
      (!legacy_scope_missing && !identical(as.character(sig$quality_scope), scope)) ||
      (!legacy_units_missing && !identical(as.character(sig$quality_units), units)) ||
      (!is.null(sig$quality_lookup) && !identical(sig$quality_lookup,
        if (is.null(incidence)) NULL else incidence$quality_lookup)) ||
      !identical(as.character(sig$demand_model), as.character(fit@demand_model)) ||
      !identical(as.character(sig$market_side), as.character(fit@market_side)) ||
      !identical(as.character(sig$quality_mode), as.character(.lq_quality_mode(fit))) ||
      !identical(as.character(sig$structure), as.character(.lq_or(fit@metadata$structure, "simultaneous"))) ||
      !.lq_fit_fingerprint_matches(fit, draws, sig$fit_fingerprint, signature_version) ||
      !identical(as.integer(signature_version), as.integer(
        if (is.null(sig$fit_signature_version)) 1L else sig$fit_signature_version))) {
    stop("'draw_plan' was created for a different fit or group map.", call. = FALSE)
  }
  if (.posterior_demand_shock_model(fit) && signature_version != 4L)
    stop("Demand-shock fits require a version-4 posterior state plan.", call. = FALSE)
  di <- as.integer(plan$draw_indices)
  if (!length(di) || anyNA(di) || any(di < 1L | di > n_total) || anyDuplicated(di)) {
    stop("'draw_plan' has malformed posterior draw indices.", call. = FALSE)
  }
  if (!identical(sig$latent_state_fingerprint, .lq_plan_state_fingerprint(plan)) ||
      !identical(as.integer(plan$metadata$draw_indices), di) ||
      !identical(rownames(plan$quality_draws), as.character(di))) {
    stop("'draw_plan' latent states or selected posterior rows were changed; ",
         "rebuild it with posterior_policy_draws().", call. = FALSE)
  }
  if (signature_version == 4L) {
    z <- plan$demand_shock_draws
    if (!is.matrix(z) || !identical(dim(z), c(length(di), nrow(fit@data))) ||
        any(!is.finite(z)) || !identical(rownames(z), as.character(di)) ||
        !identical(plan$share_residual_model, "demand_shock"))
      stop("'draw_plan' has malformed observation demand-shock states.", call. = FALSE)
  }
  metadata_role_offset <- plan$metadata$role_offset
  valid_metadata_role_offset <- length(metadata_role_offset) == 1L &&
    is.numeric(metadata_role_offset) && is.finite(metadata_role_offset) &&
    metadata_role_offset == as.integer(metadata_role_offset)
  if (!identical(plan$metadata$role_mapping, plan$role_mapping) ||
      !.lq_fixed_k_equal(plan$metadata$fixed_k, plan$fixed_k) ||
      !valid_metadata_role_offset ||
      !identical(as.integer(metadata_role_offset),
                 as.integer(plan$role_offset)) ||
      !identical(isTRUE(plan$metadata$exact_k), !is.null(plan$k_cache))) {
    stop("'draw_plan' role metadata does not match its role realization.",
         call. = FALSE)
  }
  if (!is.null(n_draws) &&
      (length(n_draws) != 1L || !is.numeric(n_draws) ||
       !is.finite(n_draws) || n_draws != as.integer(n_draws) ||
       n_draws != length(di))) {
    stop("'n_draws' conflicts with the selected draw_plan; use its draw count.", call. = FALSE)
  }
  if (!identical(match.arg(role_mapping, c("mixture", "map", "threshold")), plan$role_mapping)) {
    stop("'role_mapping' conflicts with the draw_plan role realization.", call. = FALSE)
  }
  if (!is.null(fixed_k) && !.lq_fixed_k_equal(fixed_k, plan$fixed_k)) {
    stop("'fixed_k' conflicts with the draw_plan role realization.", call. = FALSE)
  }
  if (!identical(as.integer(role_offset), as.integer(plan$role_offset))) {
    stop("'role_offset' conflicts with the draw_plan role realization.", call. = FALSE)
  }
  if (!is.matrix(plan$quality_draws) || nrow(plan$quality_draws) != length(di) ||
      ncol(plan$quality_draws) != length(units) ||
      !identical(colnames(plan$quality_draws), units)) {
    stop("'draw_plan' has malformed scoped quality draws.", call. = FALSE)
  }
  centered <- if (is.null(incidence)) abs(rowSums(plan$quality_draws)) else
    unlist(lapply(seq_len(incidence$N_quality_block), function(b) {
      first <- incidence$quality_block_start[[b]]
      last <- first + incidence$quality_block_n[[b]] - 1L
      abs(rowSums(plan$quality_draws[, first:last, drop = FALSE]))
    }), use.names = FALSE)
  if (any(!is.finite(plan$quality_draws)) || any(centered > 1e-8)) {
    stop("'draw_plan' quality draws are non-finite or violate centering.",
         call. = FALSE)
  }
  intercept_mode <- .market_intercept_parameterization_from_fit(fit)
  if (identical(intercept_mode, "collapsed")) {
    market_labels <- .lq_market_labels(fit)
    if (!is.matrix(plan$market_intercept_draws) ||
        nrow(plan$market_intercept_draws) != length(di) ||
        ncol(plan$market_intercept_draws) != length(market_labels) ||
        !identical(rownames(plan$market_intercept_draws), as.character(di)) ||
        !identical(colnames(plan$market_intercept_draws), market_labels) ||
        !identical(as.character(plan$market_labels), market_labels) ||
        !identical(as.character(plan$metadata$market_labels), market_labels) ||
        any(!is.finite(plan$market_intercept_draws))) {
      stop("'draw_plan' has malformed reconstructed market-intercept draws or labels.",
           call. = FALSE)
    }
  } else if (identical(signature_version, 3L) &&
             (!is.null(plan$market_intercept_draws) ||
              !is.null(plan$market_labels))) {
    stop("'draw_plan' contains market-intercept draws for a sampled-intercept fit.",
         call. = FALSE)
  }
  if (identical(.share_basis_from_fit(fit), "relative")) {
    if (!identical(plan$share_basis, "relative"))
      stop("Relative-share fit requires a matching relative-share draw plan.", call. = FALSE)
    expected <- .lq_relative_derived_intercepts(fit, draws, di, plan$quality_draws,
      demand_shock_draws = plan$demand_shock_draws)
    if (!identical(as.character(plan$market_year_labels), expected$labels) ||
        !isTRUE(all.equal(plan$derived_market_year_intercept_draws, expected$draws,
                         tolerance = 1e-12)))
      stop("'draw_plan' has inconsistent derived market-year intercepts.", call. = FALSE)
  }
  plan
}

.lq_market_year_labels <- function(fit) {
  d <- fit@stan_data
  n <- as.integer(d$N_market_year %||% max(d$market_year_idx %||% 1L))
  labels <- fit@metadata$market_year_labels %||%
    attr(d, "market_year_size", exact = TRUE)$market_year_labels %||%
    as.character(seq_len(n))
  labels <- as.character(labels)
  if (length(labels) != n || anyNA(labels) || any(!nzchar(labels)) || anyDuplicated(labels))
    stop("Relative-share fit has malformed market-year labels.", call. = FALSE)
  labels
}

## Retained predictor includes q only when quality was sampled in Stan.
## Recover q before deriving the level calibration; never use a zero-q GQ.
.lq_relative_full_predictor <- function(fit, draws, draw_indices, quality_draws) {
  n <- as.integer(fit@stan_data$N %||% nrow(fit@data))
  value <- .posterior_quality_matrix(draws$relative_share_predictor, n,
                                      "relative_share_predictor")
  if (nrow(value) != .lq_draw_count(fit, draws) || any(!is.finite(value)))
    stop("Relative-share predictor draws are malformed.", call. = FALSE)
  value <- value[draw_indices, , drop = FALSE]
  if (identical(.lq_quality_mode(fit), "collapsed")) {
    incidence <- .lq_quality_incidence(fit)
    index <- incidence$obs_to_quality
    if (length(index) != n || anyNA(index) || !is.matrix(quality_draws) ||
        nrow(quality_draws) != length(draw_indices) ||
        any(index < 1L | index > ncol(quality_draws)))
      stop("Relative-share predictor has malformed quality incidence.", call. = FALSE)
    value <- value + quality_draws[, index, drop = FALSE]
  }
  if (any(!is.finite(value))) stop("Relative-share predictor is nonfinite.", call. = FALSE)
  value
}

.lq_relative_derived_intercepts <- function(fit, draws, draw_indices, quality_draws,
                                           demand_shock_draws = NULL) {
  labels <- .lq_market_year_labels(fit)
  predictor <- .lq_relative_full_predictor(fit, draws, draw_indices, quality_draws)
  if (!is.null(demand_shock_draws)) predictor <- predictor + demand_shock_draws
  cells <- as.integer(fit@stan_data$market_year_idx %||% rep.int(1L, ncol(predictor)))
  if (length(cells) != ncol(predictor) || anyNA(cells) ||
      !setequal(cells, seq_along(labels)))
    stop("Relative-share fit has malformed market-year incidence.", call. = FALSE)
  outside <- .posterior_quality_matrix(draws$s0_prob, length(labels), "s0_prob")
  if (nrow(outside) != .lq_draw_count(fit, draws) ||
      any(!is.finite(outside)) || any(outside <= 0 | outside >= 1))
    stop("Relative-share outside-share draws are malformed.", call. = FALSE)
  out <- matrix(NA_real_, nrow = length(draw_indices), ncol = length(labels),
    dimnames = list(as.character(draw_indices), labels))
  for (cell in seq_along(labels)) {
    rows <- which(cells == cell)
    v <- predictor[, rows, drop = FALSE]
    top <- apply(v, 1L, max)
    out[, cell] <- -stats::qlogis(outside[draw_indices, cell]) - top -
      log(rowSums(exp(sweep(v, 1L, top, "-"))))
  }
  if (any(!is.finite(out))) stop("Derived market-year intercepts are nonfinite.", call. = FALSE)
  list(draws = out, labels = labels)
}

.lq_quality_for_draw <- function(plan, draw_index) {
  row <- match(as.integer(draw_index), plan$draw_indices)
  if (is.na(row)) stop("The draw_plan has no quality row for posterior draw ", draw_index, ".", call. = FALSE)
  as.numeric(plan$quality_draws[row, ])
}

.lq_reference_score <- function(object, draws, draw, unit, plan, market_object,
                                params) {
  p <- as.numeric(market_object$prices)
  t_obs <- as.numeric(market_object$quantities) / sum(market_object$quantities)
  ## .estimated_draw_params() returns a one-draw parameter object, whereas
  ## BayesSpecifiedFit stores its scalar/vector parameters directly.  Do not
  ## index the former by the original posterior row: for draw > 1 that would
  ## silently produce NA/length-zero outside shares.
  s0_raw <- params$s0
  s0 <- if (length(s0_raw) == 1L) as.numeric(s0_raw) else as.numeric(s0_raw[draw])
  if (!is.finite(s0) || s0 <= 0 || s0 >= 1) stop("Latent-quality simulation has an invalid outside share.", call. = FALSE)
  if (methods::is(object, "BayesSpecifiedFit")) {
    score <- log((1 - s0) * t_obs / s0)
    return(list(score = score, s0 = s0,
                quality = rep(0, length(plan$group_labels)),
                draw_index = as.integer(draw)))
  }
  if (identical(.share_basis_from_fit(object), "relative")) {
    position <- match(as.integer(draw), plan$draw_indices)
    cell <- as.integer(unit$market_year_index)
    if (is.na(position) || length(cell) != 1L || is.na(cell) ||
        cell < 1L || cell > ncol(plan$derived_market_year_intercept_draws))
      stop("Relative-share policy state has an invalid market-year map.", call. = FALSE)
    value <- .lq_relative_full_predictor(object, draws, draw,
      plan$quality_draws[position, , drop = FALSE])
    score <- value[1L, unit$idx] + plan$derived_market_year_intercept_draws[position, cell]
    if (!is.null(plan$demand_shock_draws))
      score <- score + plan$demand_shock_draws[position, unit$idx]
    return(list(score = as.numeric(score), s0 = s0,
      quality = .lq_quality_for_draw(plan, draw),
      demand_shock = if (is.null(plan$demand_shock_draws)) NULL else
        as.numeric(plan$demand_shock_draws[position, unit$idx]),
      draw_index = as.integer(draw)))
  }
  intercept_mode <- .market_intercept_parameterization_from_fit(object)
  mean_name <- if (identical(intercept_mode, "collapsed"))
    "share_mean_integrated" else "share_mean"
  mean_draws <- if (identical(intercept_mode, "collapsed"))
    draws$share_mean_integrated else draws$share_mean
  sm <- .lq_draw_matrix(mean_draws, .lq_draw_count(object, draws),
                        nrow(object@data), mean_name)
  score <- as.numeric(sm[draw, unit$idx])
  if (identical(intercept_mode, "collapsed")) {
    market_index <- as.integer(unit$market_index)
    if (length(market_index) != 1L || is.na(market_index) ||
        market_index < 1L || market_index > ncol(plan$market_intercept_draws)) {
      stop("Latent-quality market-intercept map is invalid.", call. = FALSE)
    }
    intercept_row <- match(as.integer(draw), plan$draw_indices)
    if (is.na(intercept_row))
      stop("The draw plan has no market-intercept state for posterior draw ",
           draw, ".", call. = FALSE)
    mu_b <- .lq_draw_matrix(draws$mu_b_market_vec,
                            .lq_draw_count(object, draws), 1L,
                            "mu_b_market_vec")[draw, 1L]
    intercept <- plan$market_intercept_draws[intercept_row, market_index]
    if (!is.finite(mu_b) || !is.finite(intercept))
      stop("Latent-quality market-intercept state is non-finite.", call. = FALSE)
    ## The integrated mean already contains mu_b. Restore only the realized
    ## centered offset represented by this plan's actual intercept draw.
    score <- score + intercept - mu_b
  }
  if (identical(plan$quality_mode, "collapsed")) {
    incidence <- .lq_quality_incidence(object)
    quality_index <- if (is.null(incidence)) match(
      as.character(object@data$group[unit$idx]), plan$group_labels
    ) else incidence$obs_to_quality[unit$idx]
    q <- .lq_quality_for_draw(plan, draw)
    if (anyNA(quality_index) || any(quality_index < 1L | quality_index > length(q)))
      stop("Latent-quality incidence labels are not in the fit map.", call. = FALSE)
    score <- score + q[quality_index]
  }
  list(score = score, s0 = s0, quality = .lq_quality_for_draw(plan, draw),
       draw_index = as.integer(draw))
}

.lq_replace_demand <- function(model, score, s0, prices_reference,
                               price_outside, quality, draw_index,
                               normalization = "preserve_s0_reanchor",
                               price_domain = "positive") {
  n <- length(prices_reference)
  if (length(score) != n || any(!is.finite(score))) stop("Latent demand score has invalid dimensions.", call. = FALSE)
  t <- exp(score - max(score)); t <- t / sum(t)
  shares <- (1 - s0) * t
  delta <- log(shares / s0)
  real_rates <- .model_uses_real_rate_domain(model) ||
    identical(price_domain, "real")
  old_price <- as.numeric(model@pricePre)
  if (length(old_price) != n || any(!is.finite(old_price)) ||
      (!real_rates && any(old_price <= 0))) {
    stop("Latent-quality preservation requires initialized prices valid for its domain.",
         call. = FALSE)
  }
  old_mkt_size <- if ("mktSize" %in% methods::slotNames(model)) {
    as.numeric(model@mktSize)
  } else {
    NA_real_
  }
  if (length(old_mkt_size) != 1L || !is.finite(old_mkt_size) || old_mkt_size <= 0) {
    old_mkt_size <- NA_real_
  }
  old_mc <- as.numeric(model@mcPre)
  if (length(old_mc) != n || any(!is.finite(old_mc))) {
    stop("Latent-quality preservation requires an initialized numeric mcPre.", call. = FALSE)
  }
  is_ces <- methods::is(model, "CES")
  is_auction <- methods::is(model, "Auction2ndLogit") && !is_ces
  if (methods::is(model, "BayesConductLogit") &&
      identical(as.integer(model@supply_model), 4L) && !is_auction) {
    stop("Latent-quality smooth-cutoff Logit auctions have no implemented cost-based value-index reconstruction; posted-price Logit reconstruction is invalid for this auction family.",
         call. = FALSE)
  }
  if (!is_ces && !methods::is(model, "Logit")) {
    stop("Latent-quality demand replacement currently supports Logit input/output and output CES only.",
         call. = FALSE)
  }
  if (is_ces) {
    if (!isTRUE(model@output)) {
      stop("Latent-quality demand replacement currently supports Logit input/output and output CES only.",
           call. = FALSE)
    }
    if (!is.finite(old_mkt_size) || old_mkt_size <= 0) {
      stop("Latent CES simulation requires the original positive mktSize budget.",
           call. = FALSE)
    }
    gamma <- as.numeric(model@slopes$gamma)
    if (length(gamma) != 1L || !is.finite(gamma) || gamma <= 1 ||
        any(prices_reference <= 0) || !is.finite(price_outside) || price_outside <= 0) {
      stop("Latent CES demand has invalid gamma or prices.", call. = FALSE)
    }
    model@slopes$meanval <- (shares / s0) * (prices_reference / price_outside)^gamma
  } else {
    alpha <- as.numeric(model@slopes$alpha)
    if (length(alpha) != 1L || !is.finite(alpha) || alpha == 0) stop("Latent Logit demand has an invalid signed price coefficient.", call. = FALSE)
    if (is_auction) {
      if ((isTRUE(model@output) && alpha >= 0) ||
          (!isTRUE(model@output) && alpha <= 0) || !is.na(model@normIndex) ||
          !isTRUE(all.equal(as.numeric(price_outside), 0)))
        stop("Latent auction quality requires a correctly signed Logit coefficient, a positive outside share and zero outside cost shock.", call. = FALSE)
      # Auction2ndLogit stores eta = v + alpha * mcPre, not the ordinary
      # price-based intercept. Its public share method adds alpha times the
      # additive cost change after policy. Keep mcPre and v = eta - alpha*mcPre
      # fixed while realizing probabilities from this conditional quality draw.
      # The input-market version uses positive alpha and mcPre represents the
      # value primitive; its native additive value shock follows the same index.
      model@slopes$meanval <- delta
      model@mcDelta <- rep(0, n)
    } else {
      model@slopes$meanval <- delta - alpha * (prices_reference - price_outside)
    }
  }
  model@priceStart <- prices_reference
  model@pricePre <- prices_reference
  ## calcShares() and calcMargins() read @pricePre. Capture the explicit
  ## reference state before the root changes that slot.
  reference_quantity <- as.numeric(calcShares(model, preMerger = TRUE,
                                               revenue = FALSE))
  reference_revenue <- as.numeric(calcShares(model, preMerger = TRUE,
                                              revenue = TRUE))
  reference_markup <- tryCatch(
    as.numeric(calcMargins(model, preMerger = TRUE, level = TRUE)),
    error = function(e) rep(NA_real_, n)
  )
  foc <- if (isTRUE(model@output)) {
    ## Output demand reports a positive dollar markup p - mc.
    function(price, markup) price - old_mc - markup
  } else {
    ## Input demand reports a positive dollar markdown mc - p (the value or
    ## procurement-price primitive is p + markdown).
    function(price, markup) old_mc - price - markup
  }
  foc_reference <- foc(prices_reference, reference_markup)
  ## This is deliberately the native/antitrust price root, not calcMC().
  solved <- calcPrices(model, preMerger = TRUE)
  if (length(solved) != n || any(!is.finite(solved)) ||
      (!real_rates && any(solved <= 0))) stop("Latent-quality pre-policy equilibrium returned invalid prices.", call. = FALSE)
  if (!identical(as.numeric(model@mcPre), old_mc)) {
    stop("Latent-quality pre-policy solve modified the preserved mcPre cost primitive.",
         call. = FALSE)
  }
  model@pricePre <- as.numeric(solved)
  latent_quantity <- as.numeric(calcShares(model, preMerger = TRUE, revenue = FALSE))
  latent_revenue <- as.numeric(calcShares(model, preMerger = TRUE, revenue = TRUE))
  if (methods::is(model, "CES")) {
    ## mktSize is the original total expenditure budget. The inside slot is
    ## derived from it after demand replacement; it must not redefine budget.
    if (is.finite(old_mkt_size) && old_mkt_size > 0) {
      model@mktSize <- old_mkt_size
      model@insideSize <- old_mkt_size * sum(latent_revenue)
    }
    model@shares <- latent_revenue
    model@shareInside <- sum(latent_revenue)
  } else {
    model@shares <- latent_quantity
    model@shareInside <- sum(latent_quantity)
  }
  markup_pre <- tryCatch(
    as.numeric(calcMargins(model, preMerger = TRUE, level = TRUE)),
    error = function(e) rep(NA_real_, n)
  )
  foc_solved <- foc(model@pricePre, markup_pre)
  if (any(!is.finite(markup_pre)) || any(!is.finite(foc_solved)) ||
      max(abs(foc_solved)) > 1e-6) {
    stop("Latent-quality pre-policy root failed the marginal-cost FOC check.",
         call. = FALSE)
  }
  list(
    model = model,
    quality = quality,
    draw_index = as.integer(draw_index),
    normalization = normalization,
    reference_price = prices_reference,
    solved_price = as.numeric(model@pricePre),
    reference_shares = if (methods::is(model, "CES")) reference_revenue else reference_quantity,
    solved_shares = if (methods::is(model, "CES")) latent_revenue else latent_quantity,
    reference_quantity_shares = reference_quantity,
    reference_revenue_shares = reference_revenue,
    solved_quantity_shares = latent_quantity,
    solved_revenue_shares = latent_revenue,
    mc_pre = old_mc,
    mkt_size = old_mkt_size,
    demand_representation = if (is_auction) "auction_cost_based_value_index" else
      if (is_ces) "ces_mean_value" else "logit_price_based_mean_value",
    auction_value = if (is_auction) as.numeric(model@slopes$meanval) -
      as.numeric(model@slopes$alpha) * old_mc else NULL,
    foc_reference = foc_reference,
    foc_solved = foc_solved,
    price_discrepancy = as.numeric(model@pricePre) - prices_reference
  )
}

.lq_metadata <- function(info, quality_source, quality_conditioning) {
  list(
    demand_anchoring = "latent_quality",
    demand_representation = info$demand_representation,
    auction_value = info$auction_value,
    normalization = info$normalization,
    quality_source = quality_source,
    quality_conditioning = quality_conditioning,
    posterior_draw_index = info$draw_index,
    quality = info$quality,
    demand_shock = info$demand_shock,
    reference_shares = info$reference_shares,
    latent_pre_shares = info$solved_shares,
    reference_quantity_shares = info$reference_quantity_shares,
    reference_revenue_shares = info$reference_revenue_shares,
    latent_pre_quantity_shares = info$solved_quantity_shares,
    latent_pre_revenue_shares = info$solved_revenue_shares,
    reference_inside_share = sum(info$reference_shares),
    latent_pre_inside_share = sum(info$solved_shares),
    mc_primitive = "observed-share FOC-inverted mcPre held fixed",
    foc_definition = if (isTRUE(info$model@output)) {
      "price - mc - positive dollar markup"
    } else {
      "mc - price - positive dollar markdown"
    },
    mkt_size_budget = info$mkt_size,
    reference_price = info$reference_price,
    latent_pre_price = info$solved_price,
    price_domain = if (identical(info$price_domain, "real") ||
                       .model_uses_real_rate_domain(info$model)) "real" else "positive",
    nonpositive_pre_rate = info$solved_price <= 0,
    pre_price_discrepancy = info$price_discrepancy,
    mc_pre = info$mc_pre,
    foc_at_reference = info$foc_reference,
    foc_at_latent_pre = info$foc_solved,
    max_abs_foc_at_latent_pre = max(abs(info$foc_solved), na.rm = TRUE),
    max_abs_pre_price_discrepancy = max(abs(info$price_discrepancy), na.rm = TRUE)
  )
}

.lq_apply_to_fit <- function(fit, score_info, plan,
                             price_domain = "positive") {
  info <- .lq_replace_demand(
    fit@model, score_info$score, score_info$s0,
    as.numeric(fit@model@pricePre), as.numeric(fit@model@priceOutside),
    score_info$quality, score_info$draw_index,
    price_domain = price_domain
  )
  info$demand_shock <- score_info$demand_shock
  info$price_domain <- price_domain
  fit@model <- info$model
  ## Public structural consumers also read @parameters (for example demand
  ## translations). Keep that reusable parameter state at the same realized
  ## quality as the native model rather than the observed-share inversion.
  fit@parameters <- info$model@slopes
  list(fit = fit, info = info,
       metadata = .lq_metadata(info, plan$quality_source,
                               plan$quality_conditioning))
}

.lq_apply_to_native <- function(native, score_info, plan) {
  info <- .lq_replace_demand(
    native, score_info$score, score_info$s0,
    as.numeric(native@pricePre), as.numeric(native@priceOutside),
    score_info$quality, score_info$draw_index
  )
  info$demand_shock <- score_info$demand_shock
  list(native = info$model, info = info,
       metadata = .lq_metadata(info, plan$quality_source,
                               plan$quality_conditioning))
}

## Recover the observation shock after the SAME sampled/conditional q and role
## state. Matheron's update conditions only the fitted ILR contrasts. Independent
## measurement error participates in conditioning but never enters demand.
.lq_conditional_demand_shock <- function(residual, shift, variance,
                                         measurement_sd, cell,
                                         condition_shares = TRUE) {
  n <- length(residual)
  if (length(shift) != n || length(variance) != n || length(cell) != n ||
      any(!is.finite(c(residual, shift, variance, measurement_sd))) ||
      any(variance < 0) || length(measurement_sd) != 1L || measurement_sd <= 0)
    stop("Malformed conditional demand-shock inputs.", call. = FALSE)
  out <- numeric(n)
  for (idx in split(seq_len(n), cell)) {
    basis <- .relative_share_basis(rep.int(1L, length(idx)))
    if (!ncol(basis)) next
    x <- shift[idx] + sqrt(variance[idx]) * stats::rnorm(length(idx))
    if (condition_shares) {
      error <- measurement_sd * stats::rnorm(length(idx))
      covariance <- crossprod(basis, (variance[idx] + measurement_sd^2) * basis)
      discrepancy <- crossprod(basis, residual[idx] - x - error)
      x <- x + as.numeric(variance[idx] * basis %*% solve(covariance, discrepancy))
    }
    out[idx] <- as.numeric(basis %*% crossprod(basis, x))
  }
  out
}

.lq_add_demand_shocks <- function(fit, draws, plan, seed = NULL) {
  info <- .posterior_predict_mask_and_mode(fit)
  means <- .posterior_predict_mean_state(fit, draws, plan, info)
  indices <- plan$draw_indices
  n_total <- .lq_draw_count(fit, draws)
  measurement <- .posterior_quality_vector(draws$sigma_measurement, n_total,
                                           "sigma_measurement")[indices]
  omega <- .posterior_quality_vector(draws$sigma_demand_shock, n_total,
                                     "sigma_demand_shock")[indices]
  rho <- .posterior_quality_vector(draws$rho_demand_shock, n_total,
                                   "rho_demand_shock")[indices]
  if (any(!is.finite(c(measurement, omega, rho))) || any(measurement <= 0) ||
      any(omega <= 0) || any(abs(rho) >= 1))
    stop("Invalid posterior measurement or demand-shock scales/correlations.", call. = FALSE)
  omega_by_obs <- if (isTRUE(fit@stan_data$use_market_resid_scale == 1L)) {
    .posterior_quality_matrix(draws$sigma_demand_shock_market, info$n_market,
      "sigma_demand_shock_market")[indices, info$market, drop = FALSE]
  } else matrix(omega, nrow = length(indices), ncol = info$n)
  if (!is.matrix(omega_by_obs) ||
      !identical(dim(omega_by_obs), c(length(indices), info$n)) ||
      any(!is.finite(omega_by_obs)) || any(omega_by_obs <= 0))
    stop("Invalid market-specific posterior demand-shock scales.", call. = FALSE)
  scales <- if (.posterior_contract_rho_active(fit))
    .posterior_predict_scales(fit, draws, indices, info) else NULL
  calculate <- function() {
    out <- matrix(0, length(indices), info$n,
      dimnames = list(as.character(indices), as.character(seq_len(info$n))))
    for (r in seq_along(indices)) {
      shift <- numeric(info$n)
      variance <- omega_by_obs[r, ]^2
      if (!is.null(scales)) {
        idx <- info$margin_idx
        shift[idx] <- rho[[r]] * omega_by_obs[r, idx] / scales$margin[r, idx] *
          (info$observed_margin[idx] - means$margin_mean[r, idx])
        variance[idx] <- variance[idx] * (1 - rho[[r]]^2)
      }
      out[r, ] <- .lq_conditional_demand_shock(
        info$observed_share - means$share_mean[r, ], shift, variance,
        measurement[[r]], fit@stan_data$market_year_idx,
        condition_shares = info$likelihood_mode != 3L)
    }
    out
  }
  plan$demand_shock_draws <- .lq_with_seed(seed, calculate)
  plan$measurement_sd_draws <- measurement
  plan$demand_shock_sd_draws <- omega_by_obs
  plan$demand_shock_correlation_draws <- rho
  plan$outside_share_draws <- .posterior_quality_matrix(draws$s0_prob,
    fit@stan_data$N_market_year, "s0_prob")[indices, , drop = FALSE]
  plan$share_residual_model <- "demand_shock"
  plan$metadata$share_residual_model <- "demand_shock"
  plan$metadata$demand_state <- "persistent quality plus observation demand shock"
  plan
}
