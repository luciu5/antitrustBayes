#' Simulate counterfactuals over Bayesian structural draws
#'
#' Applies an \pkg{antitrust} counterfactual to draws stored in a
#' \code{BayesFit}. Legacy structures delegate to \code{antitrust::simulate()}.
#' Exact Stackelberg and ordered-k core-fringe structures reconstruct native
#' coordination models and route supported one-step ownership, cost, role, and
#' exit changes through the corresponding deterministic simulator.
#'
#' \code{bayes_sim()} is retained as a convenience wrapper around this S4
#' method. New code may call \code{simulate()} directly.
#'
#' @param object,fit A \code{BayesFit}: either a Stan-estimated
#'   \code{BayesLogit} from \code{bayes_est()} or a fit created by
#'   \code{bayes_specify()}.
#' @param ownerPost A post-counterfactual ownership vector or matrix, an
#'   \code{antitrust::Counterfactual}, or a named list of either, one per
#'   simulation unit. It may be omitted for native exact hard-role cost-only
#'   counterfactuals supplied through \code{mcDelta}.
#' @param n_draws Number of posterior draws to simulate. By default, uses at
#'   most 200 draws.
#' @param mcDelta Optional vector of proportional marginal-cost changes for
#'   the legacy ownership interface. Put cost changes in \code{ownerPost} when
#'   it is a \code{Counterfactual}.
#' @param subset Optional logical vector selecting products for the legacy
#'   ownership interface. Put exits in \code{ownerPost} when it is a
#'   \code{Counterfactual}.
#' @param insideSize Total inside quantity. For an estimated fit the default is
#'   one, preserving the historical \code{bayes_sim()} behavior. A specified
#'   fit defaults to the value stored by \code{bayes_market()}.
#' @param priceOutside Optional positive outside-good price. The fitted or
#'   specified market value is used by default.
#' @param labels Optional product labels.
#' @param market Which markets to simulate. The specified-fit unit is named
#'   \code{"unit_1"}.
#' @param year Which years to simulate. Panel fits require an explicit year;
#'   the specified-fit year is \code{"1"}.
#' @param role_mapping Exact role realization: posterior mixture draw,
#'   posterior MAP, or legacy cutoff threshold. The default is
#'   \code{"mixture"}; threshold is unavailable for ordered-k fringe.
#' @param fixed_k Optional explicit leader or strategic-core count (scalar or
#'   named by simulation unit). Fringe permits zero.
#' @param role_offset Integer adjacent-partition sensitivity offset.
#' @param leadersPost Optional Stackelberg post-counterfactual leader override.
#' @param corePost Optional exact-fringe post-counterfactual strategic-core
#'   override. An empty vector gives a wholly MonCom post state.
#' @param demand_anchoring Demand baseline for policy simulation. The default
#'   \code{"observed"} retains the historical observed-share reconstruction.
#'   \code{"latent_quality"} draws the fitted global quality state, normalizes
#'   latent relative shares while retaining the posterior outside share, and
#'   solves a new pre-policy equilibrium using the observed-share FOC cost
#'   primitive. It supports output Logit/CES and input Logit. For input Logit,
#'   the native positive dollar markdown is retained and the observed
#'   value/procurement-price primitive is held fixed.
#' @param draw_plan Optional object returned by
#'   \code{\link{posterior_policy_draws}}. A plan fixes posterior draw,
#'   quality, and exact-k role realizations so multiple scenarios can reuse one
#'   fit-level draw plan.
#' @param mc.cores Number of parallel fork workers. The default is the number
#'   of detected cores capped at four; Windows always uses one.
#' @param ... Additional arguments passed to \code{antitrust::simulate()}.
#'
#' @return A \code{BayesSim} containing one \pkg{antitrust} result (or
#'   \code{CounterfactualPath}) per selected draw and simulation unit.
#'
#' @seealso \code{bayes_est()}, \code{bayes_specify()},
#'   \code{\link[antitrust]{counterfactual}},
#'   \code{\link[antitrust]{simulate}}
#'
#' @examples
#' data(pnb)
#' params <- bayes_params(alpha = 0.65, s0 = 0.015, supply_model = 1L)
#' market <- bayes_market(
#'   prices = pnb$prices, quantities = pnb$quantities, margins = pnb$margins,
#'   group = pnb$group, outside_price = 1,
#'   insideSize = attr(pnb, "insideSize"), labels = pnb$bank
#' )
#' fit <- bayes_specify(params, market)
#' owner_post <- pnb$group
#' owner_post[owner_post == attr(pnb, "target")] <- attr(pnb, "acquirer")
#' result <- simulate(fit, antitrust::counterfactual(ownership = owner_post))
#'
#' @name simulate-BayesFit
NULL

## Keep parallel defaults conservative on shared machines. mclapply() cannot
## fork on Windows, so explicitly normalize that platform to one worker.
.default_mc_cores <- function() {
  if (identical(.Platform$OS.type, "windows")) return(1L)
  cores <- tryCatch(
    parallel::detectCores(logical = TRUE),
    error = function(e) NA_integer_
  )
  if (is.na(cores) || cores < 1L) cores <- 1L
  min(4L, cores)
}

.validate_mc_cores <- function(mc.cores) {
  if (length(mc.cores) != 1L || !is.finite(mc.cores) || mc.cores < 1 ||
      mc.cores != as.integer(mc.cores)) {
    stop("'mc.cores' must be a positive integer.")
  }
  if (identical(.Platform$OS.type, "windows")) return(1L)
  as.integer(mc.cores)
}

.select_unit_arg <- function(x, idx, unit_id, market_id, n_rows,
                             scalar = FALSE, arg = "argument") {
  if (is.list(x) && !methods::is(x, "Counterfactual")) {
    keys <- c(unit_id, market_id)
    keys <- keys[nzchar(keys)]
    named_key <- intersect(keys, names(x))
    if (length(named_key)) return(x[[named_key[[1L]]]])
    if (length(x) == 1L && is.null(names(x))) return(x[[1L]])
    stop("No value for simulation unit '", unit_id, "' in '", arg, "'.")
  }
  if (is.matrix(x)) return(x)
  if (scalar && length(x) == 1L) return(x)
  if (length(x) == n_rows) return(x[idx])
  x
}

.estimated_simulation_units <- function(object, market, year) {
  validate_fit_spec(object)
  df <- object@data
  required <- c("prices", "quantities", "margins", "market", "year",
                "group", "outside_price")
  missing_fields <- setdiff(required, names(df))
  if (length(missing_fields)) {
    stop("The fitted data are missing required field(s): ",
         paste(missing_fields, collapse = ", "), ".")
  }

  market_factor <- as.factor(df$market)
  selected_markets <- levels(market_factor)
  if (!is.null(market)) {
    selected_markets <- intersect(as.character(market), selected_markets)
    if (!length(selected_markets)) {
      stop("None of the specified markets found in the fitted data.")
    }
  }

  year_factor <- as.factor(df$year)
  selected_years <- levels(year_factor)
  if (!is.null(year)) {
    selected_years <- intersect(as.character(year), selected_years)
    if (!length(selected_years)) {
      stop("None of the specified years found in the fitted data.")
    }
  }

  if (is.null(year)) {
    year_counts <- tapply(as.character(df$year), as.character(df$market),
                          function(x) length(unique(x)))
    if (any(year_counts[selected_markets] > 1L, na.rm = TRUE)) {
      stop("year must be supplied when simulating a panel fit; otherwise ",
           "different years would be treated as simultaneous products.",
           call. = FALSE)
    }
    units <- data.frame(
      unit_id = selected_markets,
      market = selected_markets,
      year = NA_character_,
      stringsAsFactors = FALSE
    )
  } else {
    units <- unique(data.frame(
      market = as.character(df$market),
      year = as.character(df$year),
      stringsAsFactors = FALSE
    ))
    units <- units[
      units$market %in% selected_markets & units$year %in% selected_years,
      , drop = FALSE
    ]
    if (!nrow(units)) {
      stop("No market-year cells match the requested market/year filters.")
    }
    units$unit_id <- paste(units$market, units$year, sep = "_")
    units <- units[, c("unit_id", "market", "year"), drop = FALSE]
    rownames(units) <- NULL
  }

  market_year_labels <- paste(
    as.integer(market_factor), as.integer(year_factor), sep = "_"
  )
  market_year_factor <- as.factor(market_year_labels)
  unit_info <- lapply(seq_len(nrow(units)), function(i) {
    idx <- if (is.na(units$year[[i]])) {
      which(as.character(df$market) == units$market[[i]])
    } else {
      which(as.character(df$market) == units$market[[i]] &
              as.character(df$year) == units$year[[i]])
    }
    list(
      unit_id = units$unit_id[[i]],
      market_id = units$market[[i]],
      year_id = units$year[[i]],
      idx = idx,
      market_index = as.integer(market_factor)[idx[[1L]]],
      market_year_index = as.integer(market_year_factor)[idx[[1L]]]
    )
  })
  list(data = df, units = units, info = unit_info)
}

.specified_simulation_units <- function(object, market, year) {
  if (!is.null(market) && !"unit_1" %in% as.character(market)) {
    stop("None of the specified markets found in the specified fit.")
  }
  if (!is.null(year) && !"1" %in% as.character(year)) {
    stop("None of the specified years found in the specified fit.")
  }
  list(
    data = object@data,
    units = data.frame(
      unit_id = "unit_1", market = "unit_1", year = "1",
      stringsAsFactors = FALSE
    ),
    info = list(list(
      unit_id = "unit_1", market_id = "unit_1", year_id = "1",
      idx = seq_len(object@market$n_prod),
      market_index = NA_integer_, market_year_index = NA_integer_
    ))
  )
}

.posterior_draw_count <- function(object, draws = NULL) {
  if (methods::is(object, "BayesSpecifiedFit")) {
    return(as.integer(object@params$n_draws))
  }
  a_market <- draws$a_market
  n <- if (is.null(dim(a_market))) length(a_market) else dim(a_market)[[1L]]
  if (length(n) != 1L || !is.finite(n) || n < 1L) {
    stop("The fitted object contains no posterior draws.", call. = FALSE)
  }
  as.integer(n)
}

.indexed_draw <- function(x, draw, index) {
  if (is.null(dim(x))) x[[draw]] else x[draw, index]
}

.estimated_draw_params <- function(object, draws, draw, unit) {
  is_ces <- identical(object@demand_model, "ces")
  a_raw <- .indexed_draw(draws$a_market, draw, unit$market_index)
  mky <- unit$market_year_index
  s0 <- if ("s0_prob" %in% names(draws)) {
    .indexed_draw(draws$s0_prob, draw, mky)
  } else if ("s0" %in% names(draws)) {
    stats::plogis(.indexed_draw(draws$s0, draw, mky))
  } else {
    stats::plogis(.indexed_draw(draws$s0_logit, draw, mky))
  }
  structure <- object@metadata$structure %||%
    attr(object@stan_data, "structure", exact = TRUE) %||%
    if (isTRUE(object@stan_data$use_cutoff == 1L)) "fringe" else "simultaneous"
  use_cutoff <- identical(structure, "fringe") || identical(structure, "stackelberg")
  cutoff_mode <- identical(object@metadata$partition_prior %||%
    attr(object@stan_data, "partition_prior", exact = TRUE) %||% "cutoff", "cutoff")
  cutoff <- if (use_cutoff && cutoff_mode) {
    .indexed_draw(draws$cutoff_share, draw, unit$market_index)
  } else {
    0
  }
  bargaining <- if (object@supply_model == 5L) {
    as.numeric(draws$bargaining_power[[draw]])
  } else {
    0
  }
  slope <- if (is_ces) {
    1 + a_raw / object@price_sd
  } else {
    a_raw / object@price_sd
  }
  .bayes_params_raw(
    alpha = if (is_ces) NULL else slope,
    gamma = if (is_ces) slope else NULL,
    s0 = as.numeric(s0),
    bargaining_power = bargaining,
    cutoff_share = as.numeric(cutoff),
    cutoff_scale = object@metadata$cutoff_scale %||%
      attr(object@stan_data, "cutoff_scale", exact = TRUE) %||%
      if (identical(structure, "stackelberg")) 0.01 else 0.2,
    use_cutoff = use_cutoff,
    supply_model = object@supply_model,
    structure = structure,
    partition_prior = object@metadata$partition_prior %||%
      attr(object@stan_data, "partition_prior", exact = TRUE) %||% "cutoff",
    k_support = .stack_whole_or(object@metadata$k_support, NULL),
    k_prior = .stack_whole_or(object@metadata$k_prior, NULL),
    k_scope = object@metadata$k_scope %||%
      attr(object@stan_data, "k_scope", exact = TRUE) %||% "market"
  )
}

## Translate the one-step antitrust counterfactual contract into the native
## coordination Stackelberg simulator.  The native API works with firm IDs,
## product-level proportional cost changes, and an active-product mask.
.stack_native_owner <- function(value, default, n) {
  if (is.null(value)) return(as.character(default))
  if (is.matrix(value)) {
    stop("Stackelberg counterfactual ownership must be a firm-ID vector.",
         call. = FALSE)
  }
  value <- as.character(value)
  if (length(value) != n || anyNA(value) || any(!nzchar(value))) {
    stop("Stackelberg ownership must contain one non-empty firm ID per product.",
         call. = FALSE)
  }
  value
}

.stack_native_subset <- function(value, labels, n) {
  if (is.null(value)) return(rep(TRUE, n))
  if (is.logical(value)) {
    if (length(value) != n || anyNA(value)) {
      stop("'exit' must be a logical vector matching the product count.",
           call. = FALSE)
    }
    return(as.logical(value))
  }
  if (is.character(value)) {
    if (anyNA(value) || any(!value %in% labels)) {
      stop("Character 'exit' values must match product labels.", call. = FALSE)
    }
    value <- match(value, labels)
  }
  if (!is.numeric(value) || anyNA(value) || any(value != as.integer(value)) ||
      any(value < 1L) || any(value > n)) {
    stop("'exit' must contain valid product indices or labels.", call. = FALSE)
  }
  out <- rep(TRUE, n)
  out[as.integer(value)] <- FALSE
  out
}

.stack_native_cost_delta <- function(value, labels, n) {
  if (is.null(value)) return(rep(0, n))
  if (!is.numeric(value) || anyNA(value) || any(!is.finite(value))) {
    stop("Stackelberg costs must be finite numeric proportional changes.",
         call. = FALSE)
  }
  if (!is.null(names(value))) {
    if (any(!nzchar(names(value))) || anyDuplicated(names(value)) ||
        any(!names(value) %in% labels)) {
      stop("Named Stackelberg costs must use unique product labels.",
           call. = FALSE)
    }
    out <- stats::setNames(rep(0, n), labels)
    out[names(value)] <- as.numeric(value)
    return(unname(out))
  }
  if (length(value) != n) {
    stop("Un-named Stackelberg costs must match the product count.",
         call. = FALSE)
  }
  as.numeric(value)
}

.stack_native_leaders <- function(value, owner_post, labels, default, active = NULL) {
  n <- length(owner_post)
  if (is.null(value)) return(unique(as.character(default)))
  if (is.logical(value)) {
    if (!is.null(names(value))) {
      if (any(!names(value) %in% labels) || anyDuplicated(names(value))) {
        stop("Named post-leader indicators must use product labels.", call. = FALSE)
      }
      mark <- rep(FALSE, n); names(mark) <- labels
      mark[names(value)] <- value
      value <- mark
    }
    if (length(value) != n || anyNA(value)) {
      stop("Logical post-leader indicators must match the product count.",
           call. = FALSE)
    }
    return(unique(owner_post[which(value)]))
  }
  if (is.numeric(value)) {
    if (anyNA(value) || any(value != as.integer(value)) ||
        any(value < 1L) || any(value > n)) {
      stop("Numeric post-leader roles must be valid product indices.",
           call. = FALSE)
    }
    return(unique(owner_post[as.integer(value)]))
  }
  value <- as.character(value)
  if (!length(value) || anyNA(value) || any(!nzchar(value))) {
    stop("Post-leader roles must identify at least one firm.", call. = FALSE)
  }
  if (all(value %in% unique(owner_post))) return(unique(value))
  if (all(value %in% labels)) return(unique(owner_post[match(value, labels)]))
  stop("Post-leader roles must be firm IDs or product labels.", call. = FALSE)
}

.stack_native_scenario <- function(native, scenario = NULL, owner_post = NULL,
                                   mc_delta = NULL, subset = NULL,
                                   leaders_post = NULL) {
  pre_owner <- as.character(native@firmOwnerPre)
  labels <- as.character(native@labels)
  n <- length(pre_owner)
  if (length(labels) != n || anyNA(labels) || any(!nzchar(labels)) ||
      anyDuplicated(labels)) {
    stop("Native Stackelberg counterfactuals require unique product labels.",
         call. = FALSE)
  }
  if (methods::is(scenario, "Counterfactual")) {
    if (length(scenario@steps) != 1L) {
      stop("Bayesian Stackelberg simulation supports exactly one counterfactual step.",
           call. = FALSE)
    }
    changes <- scenario@steps[[1L]]@changes
    supported <- c("ownership", "costs", "exit", "leader")
    unsupported <- setdiff(names(changes), supported)
    if (length(unsupported)) {
      stop("Native Stackelberg simulation does not support counterfactual field(s): ",
           paste(unsupported, collapse = ", "), call. = FALSE)
    }
    if (!is.null(mc_delta) || !is.null(subset) || !is.null(leaders_post)) {
      stop("A Counterfactual cannot be combined with legacy scenario arguments.",
           call. = FALSE)
    }
    owner_value <- changes$ownership
    cost_value <- changes$costs
    exit_value <- changes$exit
    leader_value <- changes$leader
  } else {
    owner_value <- scenario
    cost_value <- mc_delta
    exit_value <- NULL
    leader_value <- leaders_post
  }
  owner_changed <- !is.null(owner_value)
  owner <- .stack_native_owner(owner_value, pre_owner, n)
  active <- if (is.null(exit_value)) {
    if (is.null(subset)) rep(TRUE, n) else .stack_native_subset(subset, labels, n)
  } else .stack_native_subset(exit_value, labels, n)
  costs <- .stack_native_cost_delta(cost_value, labels, n)
  if (is.null(leader_value)) {
    if (owner_changed && !identical(owner, pre_owner)) {
      pre_leaders <- as.character(native@leadersPre)
      ## A merged post-merger owner leads whenever it contains a pre-merger
      ## leader product.  This preserves all constituent leader roles.
      leaders <- unique(owner[pre_owner %in% pre_leaders & active])
    } else {
      leaders <- as.character(native@leadersPre)
    }
  } else {
    leaders <- .stack_native_leaders(leader_value, owner, labels,
                                     native@leadersPre, active)
  }
  if (any(!leaders %in% unique(owner[active]))) {
    stop("Post-leader roles must identify active post-counterfactual firms.",
         call. = FALSE)
  }
  if (!any(active)) stop("The counterfactual must retain at least one product.", call. = FALSE)
  list(ownerPost = owner, leadersPost = leaders, mcDelta = costs, subset = active)
}

.core_native_roles <- function(value, owner_post, labels) {
  n <- length(owner_post)
  if (is.logical(value)) {
    if (!is.null(names(value))) {
      if (any(!names(value) %in% labels) || anyDuplicated(names(value))) {
        stop("Named post-core indicators must use product labels.", call. = FALSE)
      }
      mark <- rep(FALSE, n); names(mark) <- labels
      mark[names(value)] <- value
      value <- mark
    }
    if (length(value) != n || anyNA(value)) {
      stop("Logical post-core indicators must match the product count.",
           call. = FALSE)
    }
    return(unique(owner_post[which(value)]))
  }
  if (is.numeric(value)) {
    if (!length(value)) return(character())
    if (anyNA(value) || any(value != as.integer(value)) ||
        any(value < 1L) || any(value > n)) {
      stop("Numeric post-core roles must be valid product indices.",
           call. = FALSE)
    }
    return(unique(owner_post[as.integer(value)]))
  }
  value <- as.character(value)
  if (!length(value)) return(character())
  if (anyNA(value) || any(!nzchar(value))) {
    stop("Post-core roles must contain non-missing firm IDs.", call. = FALSE)
  }
  if (all(value %in% unique(owner_post))) return(unique(value))
  if (all(value %in% labels)) return(unique(owner_post[match(value, labels)]))
  stop("Post-core roles must be firm IDs or product labels.", call. = FALSE)
}

.core_native_scenario <- function(native, scenario = NULL, mc_delta = NULL,
                                  subset = NULL, core_post = NULL) {
  pre_owner <- as.character(native@firmOwnerPre)
  labels <- as.character(native@labels)
  n <- length(pre_owner)
  if (length(labels) != n || anyNA(labels) || any(!nzchar(labels)) ||
      anyDuplicated(labels)) {
    stop("Native core-fringe counterfactuals require unique product labels.",
         call. = FALSE)
  }
  if (methods::is(scenario, "Counterfactual")) {
    if (length(scenario@steps) != 1L) {
      stop("Bayesian core-fringe simulation supports exactly one counterfactual step.",
           call. = FALSE)
    }
    changes <- scenario@steps[[1L]]@changes
    supported <- c("ownership", "costs", "exit", "core")
    unsupported <- setdiff(names(changes), supported)
    if (length(unsupported)) {
      stop("Native core-fringe simulation does not support counterfactual field(s): ",
           paste(unsupported, collapse = ", "), call. = FALSE)
    }
    if (!is.null(mc_delta) || !is.null(subset) || !is.null(core_post)) {
      stop("A Counterfactual cannot be combined with legacy scenario arguments.",
           call. = FALSE)
    }
    owner_value <- changes$ownership
    cost_value <- changes$costs
    exit_value <- changes$exit
    core_value <- changes$core
  } else {
    owner_value <- scenario
    cost_value <- mc_delta
    exit_value <- NULL
    core_value <- core_post
  }
  owner_changed <- !is.null(owner_value)
  owner <- .stack_native_owner(owner_value, pre_owner, n)
  active <- if (is.null(exit_value)) {
    if (is.null(subset)) rep(TRUE, n) else .stack_native_subset(subset, labels, n)
  } else .stack_native_subset(exit_value, labels, n)
  if (!any(active)) stop("The counterfactual must retain at least one product.",
                         call. = FALSE)
  costs <- .stack_native_cost_delta(cost_value, labels, n)
  if (is.null(core_value)) {
    if (owner_changed && !identical(owner, pre_owner)) {
      core <- unique(owner[pre_owner %in% as.character(native@corePre) & active])
    } else {
      core <- as.character(native@corePre)
    }
  } else {
    core <- .core_native_roles(core_value, owner, labels)
  }
  if (any(!core %in% unique(owner[active]))) {
    stop("Post-core roles must identify active post-counterfactual firms.",
         call. = FALSE)
  }
  list(ownerPost = owner, corePost = core, mcDelta = costs, subset = active)
}

.native_ordered_structure <- function(object) {
  structure_name <- .stack_bridge_or(object@metadata$structure, "simultaneous")
  partition_prior <- .stack_bridge_or(object@metadata$partition_prior, "cutoff")
  if (identical(structure_name, "stackelberg")) return("stackelberg")
  if (identical(structure_name, "fringe") &&
      identical(partition_prior, "ordered_k")) return("fringe")
  NULL
}

.stack_fixed_k_value <- function(fixed_k, unit, n_rows) {
  if (is.null(fixed_k)) return(NULL)
  if (length(fixed_k) == 1L) return(fixed_k)
  if (!is.null(names(fixed_k))) {
    keys <- c(unit$unit_id, unit$market_id)
    hit <- keys[keys %in% names(fixed_k)]
    if (length(hit)) return(fixed_k[[hit[[1L]]]])
  }
  stop("'fixed_k' must be scalar or named by simulation market/unit.", call. = FALSE)
}

.stack_k_cache <- function(object, draws, setup, draw_indices,
                           role_mapping, fixed_k, role_offset) {
  info <- setup$info
  structure_name <- .native_ordered_structure(object)
  if (is.null(structure_name)) return(NULL)
  cache <- vector("list", length(draw_indices)); names(cache) <- as.character(draw_indices)
  for (di in seq_along(draw_indices)) {
    draw <- draw_indices[[di]]
    by_market <- split(seq_along(info), vapply(info, `[[`, character(1), "market_id"))
    one_draw <- integer(length(info)); names(one_draw) <- vapply(info, `[[`, character(1), "unit_id")
    for (market_units in by_market) {
      pscope <- .stack_bridge_or(object@metadata$partition_prior, "cutoff")
      kscope <- .stack_bridge_or(object@metadata$k_scope, "market")
      shared <- identical(pscope, "ordered_k") && identical(kscope, "market")
      groups <- if (shared) list(market_units) else lapply(market_units, identity)
      for (which_units in groups) {
        rep_unit <- info[[which_units[[1L]]]]
        bridge_unit <- list(
          data = setup$data[rep_unit$idx, , drop = FALSE],
          idx = rep_unit$idx, market = rep_unit$market_id,
          year = rep_unit$year_id, market_index = rep_unit$market_index,
          market_year_index = rep_unit$market_year_index
        )
        params <- if (methods::is(object, "BayesSpecifiedFit")) {
          .stack_bridge_params(object, bridge_unit, draw)
        } else .estimated_draw_params(object, draws, draw, rep_unit)
        sizes <- .stack_bridge_ranked_sizes(
          object, bridge_unit, bridge_unit$data, object@demand_model
        )
        comp <- .stack_bridge_component_prob(object, bridge_unit, params, sizes, draw,
                                             draws = draws)
        fk <- .stack_fixed_k_value(fixed_k, rep_unit, nrow(setup$data))
        ## The flattened component support is authoritative.  For a
        ## shared-market ordered-k unit it already reflects the minimum firm
        ## count across every market-year, including years omitted by the
        ## simulation selector.
        min_k <- min(comp$k)
        max_k <- max(comp$k)
        if (!is.null(fk) && (length(fk) != 1L || !is.finite(fk) ||
                             fk != as.integer(fk) ||
                             !as.integer(fk) %in% comp$k)) {
          stop("'fixed_k' must lie in the shared market support {",
               paste(comp$k, collapse = ", "), "}",
               " for simulation unit '", rep_unit$unit_id, "'.")
        }
        selected <- .stack_bridge_select_k(comp, params, sizes, role_mapping,
                                            fk, role_offset, support_min = min_k,
                                            support_max = max_k)
        if (selected > max_k) {
          stop("Selected role count ", selected,
               " exceeds the shared market support ", min_k, ":", max_k,
               " for simulation unit '", rep_unit$unit_id, "'.")
        }
        one_draw[which_units] <- selected
      }
    }
    cache[[di]] <- one_draw
  }
  cache
}

.simulate_native_role_fit <- function(object, ownerPost, n_draws, mcDelta, subset,
                                      insideSize, priceOutside, labels, market,
                                      year, mc.cores, role_mapping, fixed_k,
                                      role_offset, leadersPost, corePost, dots,
                                      demand_anchoring, draw_plan, .call, setup,
                                      draws, draw_indices) {
  if (!requireNamespace("coordination", quietly = TRUE) ||
      utils::packageVersion("coordination") < "0.0.0.9003") {
    stop("native Bayesian reconstruction requires coordination >= 0.0.0.9003.")
  }
  if (length(dots)) {
    stop("Unsupported native hard-role simulation argument(s): ",
         paste(names(dots), collapse = ", "), call. = FALSE)
  }
  role_mapping <- match.arg(role_mapping, c("mixture", "map", "threshold"))
  structure_name <- .native_ordered_structure(object)
  if (is.null(structure_name)) {
    stop("native hard-role simulation requires Stackelberg or exact ordered-k fringe.",
         call. = FALSE)
  }
  if (identical(structure_name, "fringe") && role_mapping == "threshold") {
    stop("role_mapping = 'threshold' is unavailable for exact ordered-k fringe.",
         call. = FALSE)
  }
  if (identical(structure_name, "fringe") && !is.null(leadersPost)) {
    stop("'leadersPost' is only available for Stackelberg simulation.",
         call. = FALSE)
  }
  if (identical(structure_name, "stackelberg") && !is.null(corePost)) {
    stop("'corePost' is only available for exact core-fringe simulation.",
         call. = FALSE)
  }
  if (length(role_offset) != 1L || !is.finite(role_offset) ||
      role_offset != as.integer(role_offset)) {
    stop("'role_offset' must be an integer scalar.", call. = FALSE)
  }
  cache <- if (!is.null(draw_plan)) draw_plan$k_cache else {
    .stack_k_cache(object, draws, setup, draw_indices, role_mapping,
                   fixed_k, role_offset)
  }
  if (is.null(cache)) {
    stop("The native role draw plan did not contain exact-k realizations.",
         call. = FALSE)
  }
  all_sims <- price_pre <- price_post <- price_change <- list()
  shares_pre <- shares_post <- product_labels <- failures <- list()
  latent_diagnostics <- list()
  n_valid <- integer(length(setup$info)); names(n_valid) <- setup$units$unit_id
  for (u in seq_along(setup$info)) {
    unit <- setup$info[[u]]; idx <- unit$idx; unit_id <- unit$unit_id
    unit_data <- setup$data[idx, , drop = FALSE]
    owner_pre <- as.character(unit_data$group)
    labels_unit <- if (is.null(labels)) {
      if (methods::is(object, "BayesSpecifiedFit")) {
        object@market$labels
      } else {
        label_col <- intersect(c("labels", "product", "product_id"), names(unit_data))
        candidate <- if (length(label_col)) as.character(unit_data[[label_col[[1L]]]]) else {
          paste0("product_", seq_along(unit_data$group))
        }
        if (anyNA(candidate) || any(!nzchar(candidate)) || anyDuplicated(candidate)) {
          paste0("product_", seq_along(candidate))
        } else candidate
      }
    } else {
      .select_unit_arg(labels, idx, unit_id, unit$market_id, nrow(setup$data), arg = "labels")
    }
    if (length(labels_unit) != length(idx)) stop("labels must match the simulation unit.", call. = FALSE)
    p0 <- if (is.null(priceOutside)) unique(unit_data$outside_price) else {
      .select_unit_arg(priceOutside, idx, unit_id, unit$market_id, nrow(setup$data),
                       scalar = TRUE, arg = "priceOutside")
    }
    inside <- if (is.null(insideSize)) {
      if (methods::is(object, "BayesSpecifiedFit")) object@market$insideSize else 1
    } else .select_unit_arg(insideSize, idx, unit_id, unit$market_id, nrow(setup$data),
                            scalar = TRUE, arg = "insideSize")
    market_object <- bayes_market(
      prices = as.numeric(unit_data$prices), quantities = as.numeric(unit_data$quantities),
      margins = as.numeric(unit_data$margins), group = owner_pre,
      outside_price = p0, market_side = object@market_side,
      insideSize = inside, labels = labels_unit
    )
    bridge_unit <- list(
      data = unit_data, idx = idx, market = unit$market_id,
      year = unit$year_id, market_index = unit$market_index,
      market_year_index = unit$market_year_index
    )
    scenario <- .select_unit_arg(ownerPost, idx, unit_id, unit$market_id,
                                 nrow(setup$data), arg = "ownerPost")
    mc_delta_unit <- if (is.null(mcDelta)) NULL else .select_unit_arg(
      mcDelta, idx, unit_id, unit$market_id, nrow(setup$data), arg = "mcDelta")
    subset_unit <- if (is.null(subset)) NULL else .select_unit_arg(
      subset, idx, unit_id, unit$market_id, nrow(setup$data), arg = "subset")
    leaders_unit <- if (is.null(leadersPost)) NULL else .select_unit_arg(
      leadersPost, idx, unit_id, unit$market_id, nrow(setup$data), arg = "leadersPost")
    core_unit <- if (is.null(corePost)) NULL else .select_unit_arg(
      corePost, idx, unit_id, unit$market_id, nrow(setup$data), arg = "corePost")
    run_one <- function(draw_index, draw_pos) {
      k <- cache[[draw_pos]][[unit_id]]
      if (identical(structure_name, "stackelberg")) {
        native <- .stack_bridge_model_draw(
          object, unit = bridge_unit, draw = draw_index,
          fixed_k = k, role_mapping = "map", role_offset = 0L,
          insideSize = inside, priceOutside = p0, labels = labels_unit,
          draws = draws, dots = list()
        )
        if (identical(demand_anchoring, "latent_quality")) {
          params <- if (methods::is(object, "BayesSpecifiedFit")) {
            object@params
          } else {
            .estimated_draw_params(object, draws, draw_index, unit)
          }
          score_info <- .lq_reference_score(
            object, draws, draw_index, bridge_unit, draw_plan,
            market_object, params
          )
          latent <- .lq_apply_to_native(native, score_info, draw_plan)
          native <- latent$native
          latent_metadata <- c(latent$metadata, list(k = as.integer(k)))
        } else {
          latent_metadata <- NULL
        }
        sc <- .stack_native_scenario(
          native, scenario = scenario,
          owner_post = if (methods::is(scenario, "Counterfactual")) NULL else scenario,
          mc_delta = mc_delta_unit, subset = subset_unit,
          leaders_post = leaders_unit
        )
        result <- coordination::stackelberg_simulate(
          native, ownerPost = sc$ownerPost, leadersPost = sc$leadersPost,
          mcDelta = sc$mcDelta, subset = sc$subset
        )
      } else {
        native <- .core_fringe_bridge_model_draw(
          object, unit = bridge_unit, draw = draw_index,
          fixed_k = k, role_mapping = "map", role_offset = 0L,
          insideSize = inside, priceOutside = p0, labels = labels_unit,
          draws = draws, dots = list()
        )
        if (identical(demand_anchoring, "latent_quality")) {
          params <- if (methods::is(object, "BayesSpecifiedFit")) {
            object@params
          } else {
            .estimated_draw_params(object, draws, draw_index, unit)
          }
          score_info <- .lq_reference_score(
            object, draws, draw_index, bridge_unit, draw_plan,
            market_object, params
          )
          latent <- .lq_apply_to_native(native, score_info, draw_plan)
          native <- latent$native
          latent_metadata <- c(latent$metadata, list(k = as.integer(k)))
        } else {
          latent_metadata <- NULL
        }
        sc <- .core_native_scenario(
          native, scenario = scenario, mc_delta = mc_delta_unit,
          subset = subset_unit, core_post = core_unit
        )
        result <- coordination::core_fringe_simulate(
          native, ownerPost = sc$ownerPost, corePost = sc$corePost,
          mcDelta = sc$mcDelta, subset = sc$subset
        )
      }
      if (!is.null(latent_metadata)) attr(result, "bayes_latent_quality") <- latent_metadata
      result
    }
    run_safe <- function(draw_pos) {
      draw_index <- draw_indices[[draw_pos]]
      tryCatch(list(result = run_one(draw_index, draw_pos), error = NA_character_),
               error = function(e) list(result = NULL, error = conditionMessage(e)))
    }
    raw <- if (length(draw_indices) == 1L || mc.cores == 1L) {
      lapply(seq_along(draw_indices), run_safe)
    } else {
      parallel::mclapply(seq_along(draw_indices), run_safe, mc.cores = mc.cores)
    }
    assembled <- .assemble_unit(raw, market_object)
    all_sims[[unit_id]] <- assembled$results
    price_pre[[unit_id]] <- assembled$price_pre
    price_post[[unit_id]] <- assembled$price_post
    price_change[[unit_id]] <- assembled$price_change
    shares_pre[[unit_id]] <- assembled$shares_pre
    shares_post[[unit_id]] <- assembled$shares_post
    product_labels[[unit_id]] <- assembled$labels
    failures[[unit_id]] <- assembled$failures
    latent_diagnostics[[unit_id]] <- assembled$latent_diagnostics
    n_valid[[unit_id]] <- assembled$n_valid
  }
  n_failed <- sum(vapply(failures, function(x) sum(!is.na(x)), integer(1)))
  if (n_failed) warning(n_failed, " of ", length(draw_indices) * length(setup$info),
                        " native hard-role simulations failed. Inspect result@failures.",
                        call. = FALSE)
  metadata <- object@metadata
  metadata$simulation_route <- if (identical(structure_name, "stackelberg")) {
    "coordination::stackelberg_simulate"
  } else {
    "coordination::core_fringe_simulate"
  }
  metadata$role_mapping <- role_mapping
  metadata$role_offset <- as.integer(role_offset)
  metadata$k_cache <- cache
  metadata$posterior_draws_available <- .posterior_draw_count(object, draws)
  metadata$posterior_draws_requested <- length(draw_indices)
  metadata$draw_indices <- as.integer(draw_indices)
  metadata$demand_anchoring <- demand_anchoring
  metadata$draw_plan <- if (is.null(draw_plan)) NULL else draw_plan$metadata
  metadata$latent_diagnostics <- latent_diagnostics
  methods::new("BayesSim", sims = all_sims, price_pre = price_pre,
               price_post = price_post, price_change = price_change,
               shares_pre = shares_pre, shares_post = shares_post,
               market_ids = as.character(setup$units$unit_id),
               product_labels = product_labels, n_draws = as.integer(length(draw_indices)),
               n_valid = n_valid, failures = failures, call = .call,
               metadata = metadata)
}

.terminal_model <- function(result) {
  if (methods::is(result, "CounterfactualPath")) {
    antitrust::final_result(result)
  } else {
    result
  }
}

.result_values <- function(result, fallback_labels) {
  model <- .terminal_model(result)
  price_pre <- as.numeric(model@pricePre)
  price_post <- as.numeric(model@pricePost)
  shares_pre <- as.numeric(antitrust::calcShares(model, preMerger = TRUE))
  shares_post <- as.numeric(antitrust::calcShares(model, preMerger = FALSE))
  k <- length(price_post)
  lengths <- c(length(price_pre), length(shares_pre), length(shares_post))
  if (!k || any(lengths != k)) {
    stop("The simulated result has inconsistent product dimensions.")
  }
  labels <- as.character(model@labels)
  if (length(labels) != k) labels <- names(model@pricePost)
  if (is.null(labels) || length(labels) != k) {
    labels <- if (length(fallback_labels) == k) {
      fallback_labels
    } else {
      paste0("product_", seq_len(k))
    }
  }
  list(
    price_pre = price_pre, price_post = price_post,
    shares_pre = shares_pre, shares_post = shares_post,
    labels = labels,
    latent_quality = attr(result, "bayes_latent_quality", exact = TRUE)
  )
}

.assemble_unit <- function(raw, market_object) {
  results <- lapply(raw, `[[`, "result")
  failures <- vapply(raw, `[[`, character(1), "error")
  values <- vector("list", length(results))
  reference_labels <- NULL

  for (i in seq_along(results)) {
    if (is.null(results[[i]])) next
    values[[i]] <- tryCatch(
      .result_values(results[[i]], market_object$labels),
      error = function(e) {
        failures[[i]] <<- conditionMessage(e)
        results[[i]] <<- NULL
        NULL
      }
    )
    if (is.null(values[[i]])) next
    if (is.null(reference_labels)) {
      reference_labels <- values[[i]]$labels
    } else if (!identical(values[[i]]$labels, reference_labels)) {
      failures[[i]] <- "The simulated result has inconsistent product labels or dimensions."
      results[[i]] <- NULL
      values[[i]] <- NULL
    }
  }

  if (is.null(reference_labels)) reference_labels <- market_object$labels
  n <- length(results)
  k <- length(reference_labels)
  empty <- matrix(NA_real_, nrow = n, ncol = k,
                  dimnames = list(NULL, reference_labels))
  price_pre <- price_post <- shares_pre <- shares_post <- empty
  for (i in seq_along(values)) {
    if (is.null(values[[i]])) next
    price_pre[i, ] <- values[[i]]$price_pre
    price_post[i, ] <- values[[i]]$price_post
    shares_pre[i, ] <- values[[i]]$shares_pre
    shares_post[i, ] <- values[[i]]$shares_post
  }
  list(
    results = results,
    failures = failures,
    labels = reference_labels,
    price_pre = price_pre,
    price_post = price_post,
    price_change = price_post / price_pre - 1,
    shares_pre = shares_pre,
    shares_post = shares_post,
    latent_diagnostics = lapply(values, function(x) {
      if (is.null(x)) NULL else x$latent_quality
    }),
    n_valid = as.integer(sum(vapply(results, Negate(is.null), logical(1))))
  )
}

.simulate_bayes_fit <- function(object, ownerPost = NULL, n_draws = NULL,
                                mcDelta = NULL, subset = NULL,
                                insideSize = NULL, priceOutside = NULL,
                                labels = NULL, market = NULL, year = NULL,
                                mc.cores = .default_mc_cores(),
                                role_mapping = c("mixture", "map", "threshold"),
                                fixed_k = NULL, role_offset = 0L,
                                leadersPost = NULL, corePost = NULL,
                                demand_anchoring = c("observed", "latent_quality"),
                                draw_plan = NULL, ..., .call) {
  demand_anchoring <- match.arg(demand_anchoring)
  if (identical(demand_anchoring, "latent_quality")) {
    if (as.integer(object@supply_model) %in% c(4L, 5L)) {
      stop("demand_anchoring = 'latent_quality' is unavailable for auction or bargaining conduct.",
           call. = FALSE)
    }
    if (!identical(as.character(object@market_side), "output") &&
        !identical(as.character(object@demand_model), "logit")) {
      stop("demand_anchoring = 'latent_quality' currently supports input Logit and output CES/Logit markets only.",
           call. = FALSE)
    }
  }
  mc.cores <- .validate_mc_cores(mc.cores)
  specified <- methods::is(object, "BayesSpecifiedFit")
  legacy_plan <- !is.null(draw_plan) &&
    (!is.list(draw_plan) || !is.list(draw_plan$fit_signature) ||
     is.null(draw_plan$fit_signature$fit_signature_version) ||
     as.integer(draw_plan$fit_signature$fit_signature_version) == 1L)
  draws <- if (specified) NULL else if (legacy_plan) {
    ## Version-1 plans hash the old unrestricted extract.  This is the one
    ## deliberate compatibility exception to the selective draw contract.
    rstan::extract(object@stanfit)
  } else .posterior_extract(
    object,
    if (!is.null(draw_plan) || identical(demand_anchoring, "latent_quality"))
      "policy_plan" else "simulation"
  )
  setup <- if (specified) {
    .specified_simulation_units(object, market, year)
  } else {
    .estimated_simulation_units(object, market, year)
  }
  n_total <- .posterior_draw_count(object, draws)
  plan <- draw_plan
  role_mapping_selected <- match.arg(role_mapping,
                                      c("mixture", "map", "threshold"))
  if (identical(demand_anchoring, "latent_quality") && is.null(plan) &&
      (!is.null(fixed_k) || !identical(role_mapping_selected, "mixture") ||
       !identical(as.integer(role_offset), 0L))) {
    stop("latent-quality simulation currently supports only the automatic mixture role plan; fixed_k, non-mixture role_mapping, and role_offset are unsupported.",
         call. = FALSE)
  }
  if (!is.null(plan)) {
    plan <- .lq_validate_plan(object, plan, n_draws = n_draws,
                              role_mapping = role_mapping,
                              fixed_k = fixed_k, role_offset = role_offset,
                              draws = draws)
    draw_indices <- as.integer(plan$draw_indices)
    n_draws <- length(draw_indices)
  } else if (identical(demand_anchoring, "latent_quality")) {
    plan <- .lq_make_policy_plan(object, draws, n_draws = n_draws)
    draw_indices <- as.integer(plan$draw_indices)
    n_draws <- length(draw_indices)
  } else {
    if (is.null(n_draws)) n_draws <- min(200L, n_total)
    if (length(n_draws) != 1L || !is.finite(n_draws) || n_draws < 1 ||
        n_draws != as.integer(n_draws)) {
      stop("'n_draws' must be a positive integer.")
    }
    n_draws <- min(as.integer(n_draws), n_total)
    draw_indices <- if (n_total == 1L) {
      1L
    } else {
      sample.int(n_total, n_draws)
    }
  }

  native_structure <- .native_ordered_structure(object)
  if (!is.null(native_structure)) {
    return(.simulate_native_role_fit(
      object = object, ownerPost = ownerPost, n_draws = n_draws,
      mcDelta = mcDelta, subset = subset, insideSize = insideSize,
      priceOutside = priceOutside, labels = labels, market = market,
      year = year, mc.cores = mc.cores, role_mapping = role_mapping,
      fixed_k = fixed_k, role_offset = role_offset, leadersPost = leadersPost,
      corePost = corePost, dots = list(...), .call = .call,
      demand_anchoring = demand_anchoring, draw_plan = plan,
      setup = setup, draws = draws,
      draw_indices = draw_indices
    ))
  }

  all_sims <- price_pre <- price_post <- price_change <- list()
  shares_pre <- shares_post <- product_labels <- failures <- list()
  latent_diagnostics <- list()
  n_valid <- integer(length(setup$info))
  names(n_valid) <- setup$units$unit_id

  for (u in seq_along(setup$info)) {
    unit <- setup$info[[u]]
    idx <- unit$idx
    unit_id <- unit$unit_id
    market_id <- unit$market_id
    df <- setup$data

    if (specified) {
      base_market <- object@market
      owner_pre <- base_market$group
      labels_unit <- if (is.null(labels)) {
        base_market$labels
      } else {
        .select_unit_arg(labels, idx, unit_id, market_id, nrow(df),
                         arg = "labels")
      }
      outside_price <- if (is.null(priceOutside)) {
        base_market$outside_price
      } else {
        .select_unit_arg(priceOutside, idx, unit_id, market_id, nrow(df),
                         scalar = TRUE, arg = "priceOutside")
      }
      inside_size <- if (is.null(insideSize)) {
        base_market$insideSize
      } else {
        .select_unit_arg(insideSize, idx, unit_id, market_id, nrow(df),
                         scalar = TRUE, arg = "insideSize")
      }
      market_object <- bayes_market(
        prices = base_market$prices,
        quantities = base_market$quantities,
        margins = base_market$margins,
        group = owner_pre,
        outside_price = outside_price,
        market_side = object@market_side,
        insideSize = inside_size,
        labels = labels_unit
      )
    } else {
      unit_data <- df[idx, , drop = FALSE]
      owner_pre <- as.character(unit_data$group)
      labels_unit <- if (is.null(labels)) {
        as.character(unit_data$group)
      } else {
        .select_unit_arg(labels, idx, unit_id, market_id, nrow(df),
                         arg = "labels")
      }
      outside_price <- if (is.null(priceOutside)) {
        unique(unit_data$outside_price)
      } else {
        .select_unit_arg(priceOutside, idx, unit_id, market_id, nrow(df),
                         scalar = TRUE, arg = "priceOutside")
      }
      inside_size <- if (is.null(insideSize)) {
        1
      } else {
        .select_unit_arg(insideSize, idx, unit_id, market_id, nrow(df),
                         scalar = TRUE, arg = "insideSize")
      }
      market_object <- bayes_market(
        prices = as.numeric(unit_data$prices),
        quantities = as.numeric(unit_data$quantities),
        margins = as.numeric(unit_data$margins),
        group = owner_pre,
        outside_price = outside_price,
        market_side = object@market_side,
        insideSize = inside_size,
        labels = labels_unit
      )
    }

    scenario <- .select_unit_arg(
      ownerPost, idx, unit_id, market_id, nrow(df), arg = "ownerPost"
    )
    is_counterfactual <- methods::is(scenario, "Counterfactual")
    if (is_counterfactual && (!is.null(mcDelta) || !is.null(subset))) {
      stop("cannot combine a Counterfactual with legacy scenario argument(s): ",
           paste(c(if (!is.null(mcDelta)) "mcDelta",
                   if (!is.null(subset)) "subset"), collapse = ", "))
    }

    if (!is_counterfactual) {
      if (is.null(scenario)) {
        stop("No post-counterfactual ownership supplied for simulation unit '",
             unit_id, "'.")
      }
      if (!is.matrix(scenario)) scenario <- as.character(scenario)
      if (is.matrix(scenario)) {
        expected <- c(market_object$n_prod, market_object$n_prod)
        if (!identical(dim(scenario), expected)) {
          stop("An ownership matrix for '", unit_id, "' must be ",
               expected[[1L]], " by ", expected[[2L]], ".")
        }
      } else if (length(scenario) != market_object$n_prod) {
        stop("'ownerPost' must have one owner per product in '", unit_id, "'.")
      }
    }

    mc_delta_unit <- if (is.null(mcDelta)) {
      rep(0, market_object$n_prod)
    } else {
      .select_unit_arg(mcDelta, idx, unit_id, market_id, nrow(df),
                       arg = "mcDelta")
    }
    subset_unit <- if (is.null(subset)) {
      rep(TRUE, market_object$n_prod)
    } else {
      as.logical(.select_unit_arg(
        subset, idx, unit_id, market_id, nrow(df), arg = "subset"
      ))
    }
    if (!is_counterfactual) {
      if (length(mc_delta_unit) != market_object$n_prod ||
          any(!is.finite(mc_delta_unit)) || any(mc_delta_unit <= -1)) {
        stop("'mcDelta' must contain finite proportional changes greater than -1 ",
             "for every product in '", unit_id, "'.")
      }
      if (length(subset_unit) != market_object$n_prod || anyNA(subset_unit) ||
          !any(subset_unit)) {
        stop("'subset' must select at least one product in '", unit_id, "'.")
      }
    }

    run_one <- function(draw_index) {
      params <- if (specified) {
        object@params
      } else {
        .estimated_draw_params(object, draws, draw_index, unit)
      }
      parameter_index <- if (specified) draw_index else 1L
      baseline <- build_conduct_object(
        params, parameter_index, market_object,
        ownerPre = owner_pre, ownerPost = owner_pre,
        mcDelta = rep(0, market_object$n_prod),
        subset = rep(TRUE, market_object$n_prod)
      )
      fit <- as_antitrust_fit(baseline)
      latent_metadata <- NULL
      if (identical(demand_anchoring, "latent_quality")) {
        score_info <- .lq_reference_score(
          object, draws, draw_index, unit, plan, market_object, params
        )
        latent <- .lq_apply_to_fit(fit, score_info, plan)
        fit <- latent$fit
        latent_metadata <- latent$metadata
      }
      result <- if (is_counterfactual) {
        antitrust::simulate(fit, scenario, ...)
      } else {
        antitrust::simulate(
          fit, scenario, mcDelta = mc_delta_unit, subset = subset_unit, ...
        )
      }
      if (!is.null(latent_metadata)) attr(result, "bayes_latent_quality") <- latent_metadata
      result
    }
    run_one_safe <- function(draw_index) {
      tryCatch(
        list(result = run_one(draw_index), error = NA_character_),
        error = function(e) list(result = NULL, error = conditionMessage(e))
      )
    }
    raw <- if (length(draw_indices) == 1L || mc.cores == 1L) {
      lapply(draw_indices, run_one_safe)
    } else {
      parallel::mclapply(draw_indices, run_one_safe, mc.cores = mc.cores)
    }
    malformed <- !vapply(raw, function(x) {
      is.list(x) && identical(names(x), c("result", "error"))
    }, logical(1))
    raw[malformed] <- list(list(
      result = NULL, error = "worker process did not return a result"
    ))

    assembled <- .assemble_unit(raw, market_object)
    all_sims[[unit_id]] <- assembled$results
    price_pre[[unit_id]] <- assembled$price_pre
    price_post[[unit_id]] <- assembled$price_post
    price_change[[unit_id]] <- assembled$price_change
    shares_pre[[unit_id]] <- assembled$shares_pre
    shares_post[[unit_id]] <- assembled$shares_post
    product_labels[[unit_id]] <- assembled$labels
    failures[[unit_id]] <- assembled$failures
    latent_diagnostics[[unit_id]] <- assembled$latent_diagnostics
    n_valid[[unit_id]] <- assembled$n_valid
  }

  n_failed <- sum(vapply(failures, function(x) sum(!is.na(x)), integer(1)))
  if (n_failed) {
    warning(
      n_failed, " of ", length(draw_indices) * length(setup$info),
      " unit-draw simulations failed. Inspect result@failures for details.",
      call. = FALSE
    )
  }
  metadata <- object@metadata
  metadata$simulation_route <- "antitrust::simulate"
  metadata$fit_class <- class(object)[[1L]]
  metadata$source_supply_model <- object@supply_model
  metadata$posterior_draws_available <- n_total
  metadata$posterior_draws_requested <- length(draw_indices)
  metadata$draw_indices <- as.integer(draw_indices)
  metadata$demand_anchoring <- demand_anchoring
  metadata$draw_plan <- if (is.null(draw_plan)) NULL else draw_plan$metadata
  metadata$latent_diagnostics <- latent_diagnostics

  methods::new(
    "BayesSim",
    sims = all_sims,
    price_pre = price_pre,
    price_post = price_post,
    price_change = price_change,
    shares_pre = shares_pre,
    shares_post = shares_post,
    market_ids = as.character(setup$units$unit_id),
    product_labels = product_labels,
    n_draws = as.integer(length(draw_indices)),
    n_valid = n_valid,
    failures = failures,
    call = .call,
    metadata = metadata
  )
}

#' @rdname simulate-BayesFit
#' @exportMethod simulate
setMethod(
  "simulate", "BayesFit",
  function(object, ownerPost = NULL, n_draws = NULL, mcDelta = NULL, subset = NULL,
           insideSize = NULL, priceOutside = NULL, labels = NULL,
           market = NULL, year = NULL,
           mc.cores = .default_mc_cores(),
           role_mapping = c("mixture", "map", "threshold"), fixed_k = NULL,
           role_offset = 0L, leadersPost = NULL, corePost = NULL,
           demand_anchoring = c("observed", "latent_quality"),
           draw_plan = NULL, ...) {
    native_structure <- .native_ordered_structure(object)
    if (missing(ownerPost) && is.null(native_structure)) {
      stop("'ownerPost' must be supplied for simulate().")
    }
    .simulate_bayes_fit(
      object = object, ownerPost = ownerPost, n_draws = n_draws,
      mcDelta = mcDelta, subset = subset, insideSize = insideSize,
      priceOutside = priceOutside, labels = labels, market = market,
      year = year, mc.cores = mc.cores, role_mapping = role_mapping,
      fixed_k = fixed_k, role_offset = role_offset, leadersPost = leadersPost,
      corePost = corePost, demand_anchoring = demand_anchoring,
      draw_plan = draw_plan,
      ..., .call = match.call()
    )
  }
)

#' @rdname simulate-BayesFit
#' @export
bayes_sim <- function(fit, ownerPost = NULL, n_draws = NULL, mcDelta = NULL,
                      subset = NULL, insideSize = NULL,
                      priceOutside = NULL, labels = NULL,
                      market = NULL, year = NULL,
                      mc.cores = .default_mc_cores(),
                      role_mapping = c("mixture", "map", "threshold"),
                      fixed_k = NULL, role_offset = 0L, leadersPost = NULL,
                      corePost = NULL,
                      demand_anchoring = c("observed", "latent_quality"),
                      draw_plan = NULL, ...) {
  native_structure <- .native_ordered_structure(fit)
  if (missing(ownerPost) && is.null(native_structure)) {
    stop("'ownerPost' must be supplied for bayes_sim().")
  }
  result <- simulate(
    fit, ownerPost = ownerPost, n_draws = n_draws,
    mcDelta = mcDelta, subset = subset, insideSize = insideSize,
    priceOutside = priceOutside, labels = labels,
    market = market, year = year, mc.cores = mc.cores,
    role_mapping = role_mapping, fixed_k = fixed_k,
    role_offset = role_offset, leadersPost = leadersPost,
    corePost = corePost, demand_anchoring = demand_anchoring,
    draw_plan = draw_plan, ...
  )
  result@call <- match.call()
  result
}
