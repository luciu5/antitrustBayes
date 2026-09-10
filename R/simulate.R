#' Simulate counterfactuals over Bayesian structural draws
#'
#' Applies an \pkg{antitrust} counterfactual to draws stored in a
#' \code{BayesFit}. For every draw and simulation unit, the method constructs
#' the pre-counterfactual structural realization, wraps it as an
#' \code{antitrust::AntitrustFit}, and delegates the scenario solve to
#' \code{antitrust::simulate()}. The result therefore uses the same ownership,
#' cost, entry, exit, quality, and multi-step machinery as \pkg{antitrust}.
#'
#' \code{bayes_sim()} is retained as a convenience wrapper around this S4
#' method. New code may call \code{simulate()} directly.
#'
#' @param object,fit A \code{BayesFit}: either a Stan-estimated
#'   \code{BayesLogit} from \code{bayes_est()} or a fit created by
#'   \code{bayes_specify()}.
#' @param ownerPost A post-counterfactual ownership vector or matrix, an
#'   \code{antitrust::Counterfactual}, or a named list of either, one per
#'   simulation unit.
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
  use_cutoff <- isTRUE(object@stan_data$use_cutoff == 1L)
  cutoff <- if (use_cutoff) {
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
    cutoff_scale = 0.2,
    use_cutoff = use_cutoff,
    supply_model = object@supply_model
  )
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
    labels = labels
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
    n_valid = as.integer(sum(vapply(results, Negate(is.null), logical(1))))
  )
}

.simulate_bayes_fit <- function(object, ownerPost, n_draws = NULL,
                                mcDelta = NULL, subset = NULL,
                                insideSize = NULL, priceOutside = NULL,
                                labels = NULL, market = NULL, year = NULL,
                                mc.cores = .default_mc_cores(), ..., .call) {
  mc.cores <- .validate_mc_cores(mc.cores)
  specified <- methods::is(object, "BayesSpecifiedFit")
  draws <- if (specified) NULL else rstan::extract(object@stanfit)
  setup <- if (specified) {
    .specified_simulation_units(object, market, year)
  } else {
    .estimated_simulation_units(object, market, year)
  }
  n_total <- .posterior_draw_count(object, draws)
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

  all_sims <- price_pre <- price_post <- price_change <- list()
  shares_pre <- shares_post <- product_labels <- failures <- list()
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
      if (is_counterfactual) {
        antitrust::simulate(fit, scenario, ...)
      } else {
        antitrust::simulate(
          fit, scenario, mcDelta = mc_delta_unit, subset = subset_unit, ...
        )
      }
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
  function(object, ownerPost, n_draws = NULL, mcDelta = NULL, subset = NULL,
           insideSize = NULL, priceOutside = NULL, labels = NULL,
           market = NULL, year = NULL,
           mc.cores = .default_mc_cores(), ...) {
    if (missing(ownerPost)) stop("'ownerPost' must be supplied for simulate().")
    .simulate_bayes_fit(
      object = object, ownerPost = ownerPost, n_draws = n_draws,
      mcDelta = mcDelta, subset = subset, insideSize = insideSize,
      priceOutside = priceOutside, labels = labels, market = market,
      year = year, mc.cores = mc.cores, ..., .call = match.call()
    )
  }
)

#' @rdname simulate-BayesFit
#' @export
bayes_sim <- function(fit, ownerPost, n_draws = NULL, mcDelta = NULL,
                      subset = NULL, insideSize = NULL,
                      priceOutside = NULL, labels = NULL,
                      market = NULL, year = NULL,
                      mc.cores = .default_mc_cores(), ...) {
  if (missing(ownerPost)) stop("'ownerPost' must be supplied for bayes_sim().")
  result <- simulate(
    fit, ownerPost = ownerPost, n_draws = n_draws,
    mcDelta = mcDelta, subset = subset, insideSize = insideSize,
    priceOutside = priceOutside, labels = labels,
    market = market, year = year, mc.cores = mc.cores, ...
  )
  result@call <- match.call()
  result
}
