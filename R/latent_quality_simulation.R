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

.lq_quality_mode <- function(fit) {
  if (methods::is(fit, "BayesSpecifiedFit")) return("none")
  stan_data <- fit@stan_data
  value <- .lq_or(fit@metadata$quality_parameterization,
                  attr(stan_data, "quality_parameterization", exact = TRUE))
  if (is.numeric(value) && length(value) == 1L) {
    value <- if (value == 1) "sampled" else if (value == 2) "collapsed" else NULL
  }
  if (is.character(value) && length(value) == 1L && value %in% c("sampled", "collapsed")) return(value)
  if (isTRUE(fit@metadata$bayes_is_collapsed_quality)) return("collapsed")
  if (identical(stan_data$quality_parameterization, 1L)) return("sampled")
  if (identical(stan_data$quality_parameterization, 2L)) return("collapsed")
  stop("Unable to resolve the fitted quality parameterization.", call. = FALSE)
}

.lq_draw_count <- function(fit, draws = NULL) {
  if (methods::is(fit, "BayesSpecifiedFit")) return(as.integer(fit@params$n_draws))
  if (is.null(draws)) draws <- rstan::extract(fit@stanfit)
  x <- draws$a_market
  n <- if (is.null(dim(x))) length(x) else dim(x)[[1L]]
  if (length(n) != 1L || !is.finite(n) || n < 1L || n != as.integer(n)) {
    stop("The fitted object contains no posterior draws.", call. = FALSE)
  }
  as.integer(n)
}

.lq_subset_draws <- function(draws, indices, n_total) {
  ## rstan returns vectors, matrices, and higher-dimensional arrays with the
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

.lq_quality_draws <- function(fit, draws, draw_indices, seed = NULL) {
  groups <- .lq_group_levels(fit)
  n_draw <- .lq_draw_count(fit, draws)
  n_group <- length(groups)
  mode <- .lq_quality_mode(fit)
  if (mode == "none") {
    out <- matrix(0, nrow = length(draw_indices), ncol = n_group,
                  dimnames = list(as.character(draw_indices), groups))
    return(list(draws = out, source = "none", conditioning = "none",
                mode = mode, group_labels = groups))
  }
  if (mode == "sampled") {
    q <- .lq_draw_matrix(draws$b_quality_all, n_draw, n_group,
                         "b_quality_all")
    if (is.null(q)) {
      if (n_group == 1L) q <- matrix(0, nrow = n_draw, ncol = 1L)
      else stop("Sampled-quality fit is missing generated quantity 'b_quality_all'.",
                call. = FALSE)
    }
    q <- q[draw_indices, , drop = FALSE]
    if (any(!is.finite(q)) || any(abs(rowSums(q)) > 1e-8)) {
      stop("Sampled global quality draws are non-finite or violate centering.",
           call. = FALSE)
    }
    colnames(q) <- groups
    rownames(q) <- as.character(draw_indices)
    return(list(draws = q, source = "sampled_stan", conditioning = "joint_draw",
                mode = mode, group_labels = groups))
  }
  ## posterior_quality() supplies the exact Gaussian conditional draw. Restrict
  ## every draw-indexed field first, then request all rows of that subset. This
  ## retains the selected original theta rows without allocating an all-draw
  ## quality matrix or performing a second independent draw selection.
  selected_draws <- .lq_subset_draws(draws, draw_indices, n_draw)
  quality_fit <- fit
  ## posterior_quality() deliberately checks the fit's recorded draw count.
  ## The shallow clone changes only that validation field; provenance in the
  ## returned plan continues to carry the original posterior indices and fit
  ## fingerprint.
  if (!is.null(quality_fit@metadata$posterior_draws)) {
    quality_fit@metadata$posterior_draws <- length(draw_indices)
  }
  pq <- posterior_quality(quality_fit, n_draws = length(draw_indices), seed = seed,
                          draws = selected_draws, conditioning = "auto")
  if (!identical(as.character(pq$group_labels), groups)) {
    stop("Conditional quality labels do not match the fitted group map.", call. = FALSE)
  }
  q <- pq$draws
  if (nrow(q) != length(draw_indices)) {
    stop("Conditional quality recovery returned an unexpected draw count.", call. = FALSE)
  }
  if (any(!is.finite(q)) || any(abs(rowSums(q)) > 1e-8)) {
    stop("Conditional quality draws are non-finite or violate centering.", call. = FALSE)
  }
  dimnames(q) <- list(as.character(draw_indices), groups)
  list(draws = q, source = "posterior_quality", conditioning =
         pq$metadata$conditioning, mode = mode, group_labels = groups)
}

.lq_all_units <- function(fit) {
  if (methods::is(fit, "BayesSpecifiedFit")) {
    return(.specified_simulation_units(fit, NULL, NULL))
  }
  years <- unique(as.character(fit@data$year))
  if (length(years) > 1L) .estimated_simulation_units(fit, NULL, years)
  else .estimated_simulation_units(fit, NULL, NULL)
}

.lq_role_cache <- function(fit, draws, draw_indices) {
  if (is.null(.native_ordered_structure(fit))) return(NULL)
  setup <- .lq_all_units(fit)
  .stack_k_cache(fit, draws, setup, draw_indices,
                 role_mapping = "mixture", fixed_k = NULL, role_offset = 0L)
}

.lq_fit_fingerprint <- function(fit, draws = NULL) {
  payload <- list(
    data = fit@data,
    demand_model = as.character(fit@demand_model),
    market_side = as.character(fit@market_side),
    supply_model = as.integer(fit@supply_model),
    metadata = fit@metadata,
    stan_data = if (methods::is(fit, "BayesSpecifiedFit")) NULL else fit@stan_data,
    params = if (methods::is(fit, "BayesSpecifiedFit")) fit@params else NULL,
    market = if (methods::is(fit, "BayesSpecifiedFit")) fit@market else NULL,
    draws = if (methods::is(fit, "BayesSpecifiedFit")) NULL else draws
  )
  digest::digest(payload, algo = "xxhash64", serialize = TRUE)
}

.lq_plan_signature <- function(fit, groups, n_total, mode, draws = NULL) {
  list(
    demand_model = as.character(fit@demand_model),
    market_side = as.character(fit@market_side),
    posterior_draws = as.integer(n_total),
    group_labels = as.character(groups),
    quality_mode = mode,
    structure = as.character(.lq_or(fit@metadata$structure, "simultaneous")),
    fit_fingerprint = .lq_fit_fingerprint(fit, draws)
  )
}

.lq_plan_state_fingerprint <- function(plan) {
  digest::digest(plan[c("draw_indices", "quality_draws", "group_labels",
                       "quality_source", "quality_conditioning", "quality_mode",
                       "k_cache", "role_mapping", "fixed_k", "role_offset")],
                 algo = "xxhash64", serialize = TRUE)
}

#' Build a reusable posterior draw plan for policy simulation
#'
#' The plan samples posterior indices once and, for collapsed fits, draws one
#' conditional global quality vector per selected parameter draw. Exact native
#' ordered-k fits also receive one role realization per market and draw. Reuse
#' the returned object in observed-anchor and latent-quality simulations.
#'
#' @param fit A \code{BayesFit}.
#' @param n_draws Number of posterior draws, defaulting to at most 200.
#' @param seed Optional nonnegative integer random seed.
#' @return A reusable object of class \code{bayes_policy_draw_plan}.
#' @export
posterior_policy_draws <- function(fit, n_draws = NULL, seed = NULL) {
  if (!methods::is(fit, "BayesFit")) stop("'fit' must be a BayesFit object.", call. = FALSE)
  seed <- .lq_validate_seed(seed)
  draws <- if (methods::is(fit, "BayesSpecifiedFit")) NULL else rstan::extract(fit@stanfit)
  n_total <- .lq_draw_count(fit, draws)
  if (is.null(n_draws)) n_draws <- min(200L, n_total)
  if (length(n_draws) != 1L || !is.finite(n_draws) || n_draws < 1L ||
      n_draws != as.integer(n_draws)) stop("'n_draws' must be a positive integer.", call. = FALSE)
  n_draws <- min(as.integer(n_draws), n_total)
  draw_indices <- .lq_with_seed(seed, function() {
    if (n_total == 1L) 1L else sample.int(n_total, n_draws)
  })
  quality_seed <- if (is.null(seed)) NULL else as.integer((as.double(seed) + 1) %% .Machine$integer.max)
  quality <- .lq_quality_draws(fit, draws, draw_indices, seed = quality_seed)
  role_seed <- if (is.null(seed)) NULL else as.integer((as.double(seed) + 2) %% .Machine$integer.max)
  k_cache <- .lq_with_seed(role_seed, function() .lq_role_cache(fit, draws, draw_indices))
  plan <- structure(list(
    draw_indices = as.integer(draw_indices),
    quality_draws = quality$draws,
    quality_source = quality$source,
    quality_conditioning = quality$conditioning,
    quality_mode = quality$mode,
    group_labels = quality$group_labels,
    k_cache = k_cache,
    role_mapping = "mixture",
    fixed_k = NULL,
    role_offset = 0L,
    fit_signature = .lq_plan_signature(fit, quality$group_labels, n_total,
                                       quality$mode, draws = draws),
    metadata = list(
      posterior_draws_available = as.integer(n_total),
      posterior_draws_requested = as.integer(length(draw_indices)),
      draw_indices = as.integer(draw_indices),
      quality_source = quality$source,
      quality_conditioning = quality$conditioning,
      quality_parameterization = quality$mode,
      group_labels = quality$group_labels,
      role_mapping = "mixture", role_offset = 0L,
      exact_k = !is.null(k_cache)
    )
  ), class = c("bayes_policy_draw_plan", "list"))
  plan$fit_signature$latent_state_fingerprint <- .lq_plan_state_fingerprint(plan)
  plan
}

.lq_validate_plan <- function(fit, plan, n_draws = NULL,
                              role_mapping = "mixture", fixed_k = NULL,
                              role_offset = 0L, draws = NULL) {
  if (!inherits(plan, "bayes_policy_draw_plan")) {
    stop("'draw_plan' must be created by posterior_policy_draws().", call. = FALSE)
  }
  groups <- .lq_group_levels(fit)
  draws <- if (methods::is(fit, "BayesSpecifiedFit")) NULL else {
    if (is.null(draws)) rstan::extract(fit@stanfit) else draws
  }
  n_total <- .lq_draw_count(fit, draws)
  sig <- plan$fit_signature
  if (!is.list(sig) || !identical(as.integer(sig$posterior_draws), as.integer(n_total)) ||
      !identical(as.character(sig$group_labels), groups) ||
      !identical(as.character(sig$demand_model), as.character(fit@demand_model)) ||
      !identical(as.character(sig$market_side), as.character(fit@market_side)) ||
      !identical(as.character(sig$quality_mode), as.character(.lq_quality_mode(fit))) ||
      !identical(as.character(sig$structure), as.character(.lq_or(fit@metadata$structure, "simultaneous"))) ||
      !identical(as.character(sig$fit_fingerprint),
                 as.character(.lq_fit_fingerprint(fit, draws)))) {
    stop("'draw_plan' was created for a different fit or group map.", call. = FALSE)
  }
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
  if (!is.null(n_draws) &&
      (length(n_draws) != 1L || !is.numeric(n_draws) ||
       !is.finite(n_draws) || n_draws != as.integer(n_draws) ||
       n_draws != length(di))) {
    stop("'n_draws' conflicts with the selected draw_plan; use its draw count.", call. = FALSE)
  }
  if (!identical(match.arg(role_mapping, c("mixture", "map", "threshold")), plan$role_mapping)) {
    stop("'role_mapping' conflicts with the draw_plan role realization.", call. = FALSE)
  }
  if (!is.null(fixed_k) || !identical(as.integer(role_offset), as.integer(plan$role_offset))) {
    stop("'fixed_k' or 'role_offset' conflicts with the draw_plan role realization.", call. = FALSE)
  }
  if (!is.matrix(plan$quality_draws) || nrow(plan$quality_draws) != length(di) ||
      ncol(plan$quality_draws) != length(groups) ||
      !identical(colnames(plan$quality_draws), groups)) {
    stop("'draw_plan' has malformed global quality draws.", call. = FALSE)
  }
  if (any(!is.finite(plan$quality_draws)) ||
      any(abs(rowSums(plan$quality_draws)) > 1e-8)) {
    stop("'draw_plan' quality draws are non-finite or violate centering.",
         call. = FALSE)
  }
  plan
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
  sm <- .lq_draw_matrix(draws$share_mean, .lq_draw_count(object, draws),
                        nrow(object@data), "share_mean")
  score <- as.numeric(sm[draw, unit$idx])
  if (identical(plan$quality_mode, "collapsed")) {
    group <- as.character(object@data$group[unit$idx])
    q <- .lq_quality_for_draw(plan, draw)
    map <- match(group, plan$group_labels)
    if (anyNA(map)) stop("Latent-quality group labels are not in the fit map.", call. = FALSE)
    score <- score + q[map]
  }
  list(score = score, s0 = s0, quality = .lq_quality_for_draw(plan, draw),
       draw_index = as.integer(draw))
}

.lq_replace_demand <- function(model, score, s0, prices_reference,
                               price_outside, quality, draw_index,
                               normalization = "preserve_s0_reanchor") {
  n <- length(prices_reference)
  if (length(score) != n || any(!is.finite(score))) stop("Latent demand score has invalid dimensions.", call. = FALSE)
  t <- exp(score - max(score)); t <- t / sum(t)
  shares <- (1 - s0) * t
  delta <- log(shares / s0)
  old_price <- as.numeric(model@pricePre)
  if (length(old_price) != n || any(!is.finite(old_price)) || any(old_price <= 0)) {
    stop("Latent-quality preservation requires initialized positive baseline prices.",
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
    model@slopes$meanval <- delta - alpha * (prices_reference - price_outside)
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
  if (length(solved) != n || any(!is.finite(solved)) || any(solved <= 0)) stop("Latent-quality pre-policy equilibrium returned invalid prices.", call. = FALSE)
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
    foc_reference = foc_reference,
    foc_solved = foc_solved,
    price_discrepancy = as.numeric(model@pricePre) - prices_reference
  )
}

.lq_metadata <- function(info, quality_source, quality_conditioning) {
  list(
    demand_anchoring = "latent_quality",
    normalization = info$normalization,
    quality_source = quality_source,
    quality_conditioning = quality_conditioning,
    posterior_draw_index = info$draw_index,
    quality = info$quality,
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
    pre_price_discrepancy = info$price_discrepancy,
    mc_pre = info$mc_pre,
    foc_at_reference = info$foc_reference,
    foc_at_latent_pre = info$foc_solved,
    max_abs_foc_at_latent_pre = max(abs(info$foc_solved), na.rm = TRUE),
    max_abs_pre_price_discrepancy = max(abs(info$price_discrepancy), na.rm = TRUE)
  )
}

.lq_apply_to_fit <- function(fit, score_info, plan) {
  info <- .lq_replace_demand(
    fit@model, score_info$score, score_info$s0,
    as.numeric(fit@model@pricePre), as.numeric(fit@model@priceOutside),
    score_info$quality, score_info$draw_index
  )
  fit@model <- info$model
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
  list(native = info$model, info = info,
       metadata = .lq_metadata(info, plan$quality_source,
                               plan$quality_conditioning))
}
