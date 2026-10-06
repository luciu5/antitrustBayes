# Calibrated one-market output-Logit designs for end-to-end estimation checks.

.bayes_calibrated_conducts <- function() {
  list(
    MC = list(label = "MC", supply_model = "monopolistic",
              antitrust_conduct = "moncom", structure = "simultaneous",
              route = "antitrust"),
    B = list(label = "B", supply_model = "bertrand",
             antitrust_conduct = "bertrand", structure = "simultaneous",
             route = "antitrust"),
    C = list(label = "C", supply_model = "cournot",
             antitrust_conduct = "cournot", structure = "simultaneous",
             route = "antitrust"),
    A = list(label = "A", supply_model = "auction",
             antitrust_conduct = "auction2nd", structure = "simultaneous",
             route = "antitrust"),
    BF = list(label = "BF", supply_model = "bertrand",
              antitrust_conduct = "bertrand", structure = "fringe",
              route = "core_fringe"),
    CF = list(label = "CF", supply_model = "cournot",
              antitrust_conduct = "cournot", structure = "fringe",
              route = "core_fringe"),
    BS = list(label = "BS", supply_model = "bertrand",
              antitrust_conduct = "bertrand", structure = "stackelberg",
              route = "stackelberg"),
    CS = list(label = "CS", supply_model = "cournot",
              antitrust_conduct = "cournot", structure = "stackelberg",
              route = "stackelberg")
  )
}

.bayes_calibrated_geometry <- function(J, concentration, inside_size,
                                       conditional_shares = NULL) {
  shape <- if (length(concentration) == 1L) {
    rep(concentration, J)
  } else {
    concentration
  }
  if (is.null(conditional_shares)) {
    weights <- stats::rgamma(J, shape = shape, rate = 1)
    if (any(!is.finite(weights)) || any(weights <= 0) ||
        !is.finite(sum(weights)) || sum(weights) <= 0) {
      stop("Dirichlet draw was numerically degenerate; increase the concentration.",
           call. = FALSE)
    }
    conditional_shares <- unname(weights / sum(weights))
  } else {
    conditional_shares <- unname(as.numeric(conditional_shares))
  }
  mean_log_share <- mean(log(conditional_shares))
  s0_formula <- stats::plogis(mean_log_share)
  shares <- unname((1 - s0_formula) * conditional_shares)
  s0_identity <- 1 - sum(shares)
  if (any(!is.finite(conditional_shares)) ||
      any(conditional_shares <= 0) ||
      any(!is.finite(shares)) || any(shares <= 0) ||
      !is.finite(s0_formula) || s0_formula <= 0 || s0_formula >= 1) {
    stop("Dirichlet draw did not retain positive finite shares; increase the concentration.",
         call. = FALSE)
  }
  prices <- rep(1, J)
  products <- paste0("product_", seq_len(J))
  names(prices) <- names(shares) <- names(conditional_shares) <- products
  b <- log(conditional_shares) - mean_log_share
  names(b) <- products
  ranking <- order(-conditional_shares, seq_len(J))
  list(
    prices = prices,
    shares = shares,
    conditional_shares = conditional_shares,
    s0 = s0_formula,
    s0_identity = s0_identity,
    mean_log_share = mean_log_share,
    b = b,
    quantities = unname(inside_size * conditional_shares),
    group = products,
    labels = products,
    concentration = shape,
    ranking = ranking,
    anchor_index = ranking[[1L]]
  )
}

.bayes_calibrated_specify <- function(conduct, geometry, alpha, inside_size,
                                      price_start) {
  meanval <- if (identical(conduct$antitrust_conduct, "auction2nd")) {
    geometry$b
  } else {
    geometry$b + alpha * geometry$prices
  }
  antitrust::specify(
    demand = "logit", conduct = conduct$antitrust_conduct,
    prices = unname(geometry$prices),
    parameters = list(alpha = -alpha, meanval = unname(meanval)),
    ownerPre = geometry$group,
    insideSize = inside_size, labels = geometry$labels,
    priceOutside = 0, priceStart = price_start, output = TRUE
  )
}

.bayes_calibrated_public_model <- function(conduct, geometry, alpha, k,
                                           inside_size, price_start) {
  J <- length(geometry$prices)
  margins_observed <- rep(NA_real_, J)
  if (identical(conduct$route, "core_fringe")) {
    model <- coordination::core_fringe(
      prices = unname(geometry$prices), shares = unname(geometry$shares),
      margins = margins_observed, ownerPre = geometry$group,
      corePre = if (k == 0L) character() else
        geometry$group[geometry$ranking[seq_len(k)]],
      demand = "logit", conduct = conduct$antitrust_conduct,
      output = TRUE, insideSize = inside_size, normIndex = NA_real_,
      priceOutside = 0, priceStart = price_start, labels = geometry$labels,
      alpha = alpha, control.equ = list(implicitCheck = FALSE)
    )
    fitted_alpha <- as.numeric(model@diagnostics$alpha)
  } else if (identical(conduct$route, "stackelberg")) {
    model <- coordination::stackelberg(
      prices = unname(geometry$prices), shares = unname(geometry$shares),
      margins = margins_observed, ownerPre = geometry$group,
      leadersPre = geometry$group[geometry$ranking[seq_len(k)]],
      demand = "logit", conduct = conduct$antitrust_conduct,
      output = TRUE, insideSize = inside_size, normIndex = NA_real_,
      priceOutside = 0, priceStart = price_start, labels = geometry$labels,
      alpha = alpha, control.equ = list(implicitCheck = FALSE)
    )
    fitted_alpha <- as.numeric(model@diagnostics$alpha)
  } else {
    source_fit <- .bayes_calibrated_specify(
      conduct, geometry, alpha, inside_size, price_start
    )
    model <- source_fit@model
    fitted_alpha <- alpha
  }
  if (length(fitted_alpha) != 1L || !is.finite(fitted_alpha) ||
      abs(fitted_alpha - alpha) > 1e-8 * max(1, alpha)) {
    stop("Public model did not retain the supplied Logit alpha.", call. = FALSE)
  }
  list(model = model, source_fit = if (identical(conduct$route, "antitrust"))
    source_fit else NULL)
}

.bayes_calibrated_unit_markup <- function(conduct, geometry, k, inside_size) {
  ## Logit dollar markups depend on shares and alpha, not on the price level.
  ## Use public FOCs at alpha = 1 to obtain the conduct-specific coefficient.
  ## The probe's implied costs are never used as market primitives.
  probe <- geometry
  probe$prices <- stats::setNames(rep(100, length(geometry$prices)),
                                  geometry$labels)
  constructed <- .bayes_calibrated_public_model(
    conduct, probe, alpha = 1, k = k, inside_size = inside_size,
    price_start = unname(probe$prices * 1.13 + 0.037)
  )
  markup <- as.numeric(antitrust::calcMargins(
    constructed$model, preMerger = TRUE, level = TRUE
  ))
  if (length(markup) != length(geometry$prices) ||
      any(!is.finite(markup)) || any(markup <= 0)) {
    stop("Public conduct model returned invalid structural markups.",
         call. = FALSE)
  }
  stats::setNames(markup, geometry$labels)
}

.bayes_calibrated_model <- function(conduct, geometry, alpha, costs, k,
                                    inside_size, price_tolerance,
                                    anchor_margin = NULL) {
  J <- length(geometry$prices)
  price_start <- unname(geometry$prices * 1.13 + 0.037)
  names(price_start) <- geometry$labels
  constructed <- .bayes_calibrated_public_model(
    conduct, geometry, alpha, k, inside_size, price_start
  )
  model <- constructed$model
  fitted_margins <- tryCatch(
    as.numeric(antitrust::calcMargins(model, preMerger = TRUE, level = FALSE)),
    error = function(e) stop("Could not recover public-model margins: ",
                             conditionMessage(e), call. = FALSE)
  )
  recovered_costs <- tryCatch(
    as.numeric(antitrust::calcMC(model, preMerger = TRUE)),
    error = function(e) stop("Could not recover public-model costs: ",
                             conditionMessage(e), call. = FALSE)
  )
  if (length(fitted_margins) != J || any(!is.finite(fitted_margins)) ||
      any(fitted_margins <= 0 | fitted_margins >= 1)) {
    stop("Public model returned invalid full-market margins.", call. = FALSE)
  }
  if (length(recovered_costs) != J || any(!is.finite(recovered_costs))) {
    stop("Public model returned invalid marginal costs.", call. = FALSE)
  }
  cost_error <- max(abs(recovered_costs - costs) / pmax(1, abs(costs)))
  if (!is.finite(cost_error) || cost_error > price_tolerance) {
    stop("Public model did not recover the supplied physical costs (error = ",
         format(cost_error, digits = 6), ").", call. = FALSE)
  }
  if (!is.null(anchor_margin) &&
      !isTRUE(all.equal(fitted_margins[[geometry$anchor_index]], anchor_margin,
                        tolerance = 1e-7))) {
    stop("Public model did not reproduce the supplied anchor margin.",
         call. = FALSE)
  }

  check_model <- model
  if (!("priceStart" %in% methods::slotNames(check_model))) {
    stop("Public model does not expose a perturbed priceStart slot.",
         call. = FALSE)
  }
  check_model@priceStart <- price_start
  if (!isTRUE(all.equal(as.numeric(check_model@priceStart),
                        unname(price_start), tolerance = 0))) {
    stop("Public model did not retain the requested perturbed price start.",
         call. = FALSE)
  }
  equilibrium_prices <- tryCatch(
    as.numeric(antitrust::calcPrices(check_model, preMerger = TRUE)),
    error = function(e) stop("Fresh equilibrium solve failed: ",
                             conditionMessage(e), call. = FALSE)
  )
  if (length(equilibrium_prices) != J ||
      any(!is.finite(equilibrium_prices) | equilibrium_prices <= 0)) {
    stop("Fresh equilibrium solve returned invalid prices.", call. = FALSE)
  }
  max_relative_error <- max(abs(equilibrium_prices - geometry$prices) /
                              pmax(abs(geometry$prices), 1e-12))
  if (!is.finite(max_relative_error) || max_relative_error > price_tolerance) {
    stop("Fresh equilibrium solve does not reproduce the market prices (relative error = ",
         format(max_relative_error, digits = 6), ").", call. = FALSE)
  }
  check_model@pricePre <- equilibrium_prices
  public_shares <- tryCatch(
    as.numeric(antitrust::calcShares(check_model, preMerger = TRUE)),
    error = function(e) stop("Could not recompute public-model shares at solved prices: ",
                             conditionMessage(e), call. = FALSE)
  )
  if (length(public_shares) != J || any(!is.finite(public_shares))) {
    stop("Public model returned invalid baseline shares.", call. = FALSE)
  }
  share_error <- max(abs(public_shares - geometry$shares))
  if (!is.finite(share_error) || share_error > price_tolerance) {
    stop("Public model did not reproduce the generated unconditional shares (error = ",
         format(share_error, digits = 6), ").", call. = FALSE)
  }
  role_ids <- if (is.na(k)) character() else
    geometry$group[geometry$ranking[seq_len(k)]]
  if (identical(conduct$route, "core_fringe") &&
      !identical(as.character(model@corePre), as.character(role_ids))) {
    stop("Public model core roles disagree with the fixed share ranking.",
         call. = FALSE)
  }
  if (identical(conduct$route, "stackelberg") &&
      !identical(as.character(model@leadersPre), as.character(role_ids))) {
    stop("Public model leader roles disagree with the fixed share ranking.",
         call. = FALSE)
  }
  foc_residual <- NA_real_
  if (identical(conduct$route, "core_fringe")) {
    foc_residual <- as.numeric(coordination::core_fringe_residuals(
      model, preMerger = TRUE
    )$max)
  } else if (identical(conduct$route, "stackelberg")) {
    foc_residual <- as.numeric(coordination::stackelberg_residuals(
      model, preMerger = TRUE
    )$max)
  }
  foc_available <- conduct$route %in% c("core_fringe", "stackelberg")
  if (foc_available && (!is.finite(foc_residual) || foc_residual > 1e-6)) {
    stop("Public model baseline FOC residual exceeds tolerance.",
         call. = FALSE)
  }
  gates <- list(
    positive_costs = all(is.finite(costs) & costs > 0),
    primitive_costs = cost_error <= price_tolerance,
    anchor_margin = if (is.null(anchor_margin)) TRUE else
      isTRUE(all.equal(fitted_margins[[geometry$anchor_index]], anchor_margin,
                       tolerance = 1e-7)),
    perturbed_price_solve = max_relative_error <= price_tolerance,
    solved_share_reproduction = share_error <= price_tolerance,
    baseline_foc = if (foc_available) foc_residual <= 1e-6 else NA
  )
  known_gates <- !vapply(gates, is.na, logical(1))
  overall_passed <- all(unlist(gates[known_gates], use.names = FALSE))
  list(
    alpha = alpha,
    margins = stats::setNames(fitted_margins, geometry$labels),
    costs = stats::setNames(as.numeric(costs), geometry$labels),
    model = model,
    source_fit = constructed$source_fit,
    calibration_source = "cost primitives and public conduct FOCs",
    equilibrium = list(
      prices = stats::setNames(equilibrium_prices, geometry$labels),
      shares = stats::setNames(public_shares, geometry$labels),
      max_relative_error = max_relative_error,
      max_share_error = share_error,
      max_cost_error = cost_error,
      max_foc_residual = foc_residual,
      gates = gates,
      foc_status = if (foc_available) "checked" else
        "not exposed by public antitrust API",
      used_perturbed_start = TRUE,
      tolerance = price_tolerance,
      passed = overall_passed,
      start = price_start
    )
  )
}

.bayes_calibrated_data <- function(geometry, margins, conduct, k,
                                   inside_size) {
  J <- length(geometry$prices)
  observed <- seq_len(J) == geometry$anchor_index
  observed_margins <- rep(NA_real_, J)
  observed_margins[observed] <- margins[observed]
  market <- rep("market_1", J)
  year <- rep("1", J)
  args <- list(
    prices = unname(geometry$prices),
    quantities = unname(geometry$quantities),
    margins = observed_margins,
    market = market,
    year = year,
    group = geometry$group,
    demand_model = "logit",
    market_side = "output",
    supply_model = conduct$supply_model,
    margin_observed = observed
  )
  if (identical(conduct$structure, "fringe")) {
    args$structure <- "fringe"
    args$partition_prior <- "ordered_k"
    args$k_support <- 0:J
    args$k_scope <- "market"
  } else if (identical(conduct$structure, "stackelberg")) {
    args$structure <- "stackelberg"
    args$partition_prior <- "ordered_k"
    args$k_support <- seq_len(J)
    args$k_scope <- "market"
  }
  args
}

#' Generate a calibrated one-market output-Logit design
#'
#' Draws one conduct-independent Logit share geometry and physical cost vector,
#' then sets each price to its physical cost plus the structural markup implied
#' by the supplied conduct. It either calibrates the structural price
#' coefficient from one observed margin or uses a supplied coefficient. The
#' conditional inside shares are drawn from a `J`-dimensional Dirichlet with
#' the required concentration parameter. Their geometric mean sets the
#' outside share through `plogis(mean(log(r)))`; centered log shares are kept
#' as the baseline quality coordinates.
#' Calibration and equilibrium checks use the public `antitrust` or
#' `coordination` model constructors. The returned estimation inputs keep
#' only the highest-share anchor margin; every other margin is `NA`, while
#' `truth$data` retains the complete model-implied margin vector. The
#' highest conditional-share product is the anchor, and the same ranking
#' assigns core firms or leaders before any conduct-specific calculations.
#'
#' The supported conduct labels are `MC` (monopolistic competition), `B`
#' (Bertrand), `C` (Cournot), `A` (second-price auction), `BF` and `CF`
#' (Bertrand and Cournot core-fringe), and `BS` and `CS` (Bertrand and
#' Cournot Stackelberg). For role models, `k` fixes the number of core firms or
#' leaders in the truth. The estimation inputs retain the full supported
#' ordered-k range so the estimator does not learn the true role count.
#'
#' @param conduct One of `MC`, `B`, `C`, `A`, `BF`, `CF`, `BS`, or `CS`.
#' @param J Number of inside products; each product belongs to its own firm.
#' @param dirichlet_concentration Required positive Dirichlet concentration,
#'   either a scalar recycled across the `J` inside products or a vector of
#'   length `J`, unless `conditional_shares` is supplied.
#' @param conditional_shares Optional explicit positive inside-share vector
#'   of length `J` summing to one. Supplying it skips the Dirichlet draw.
#' @param anchor_margin Optional observed Lerner margin for the prespecified
#'   highest-share product in calibrated mode. If omitted, `margin_range` is
#'   required and supplies an open interval for a seeded draw.
#' @param margin_range Optional open interval used to draw the anchor margin.
#'   It is required when `anchor_margin` is omitted; both endpoints must lie
#'   strictly between zero and one.
#' @param seed Optional integer seed. The generator restores the caller's
#'   random-number state after drawing the design.
#' @param k Optional number of core firms or Stackelberg leaders. The default
#'   is `max(1, floor(J / 2))`. It is only used by `BF`, `CF`, `BS`, and `CS`.
#' @param mode Either `"calibrated"`, which recovers alpha from the anchor
#'   margin, or `"parameter"`, which uses the supplied `alpha`.
#' @param alpha Positive absolute Logit price coefficient required in
#'   parameter mode.
#' @param costs Optional explicit positive physical marginal costs, one per
#'   product. If omitted, costs are drawn independently of conduct.
#' @param cost_rule Draw `"uniform"` costs from `cost_range`, or use the common
#'   `cost_level`. Ignored when `costs` is supplied.
#' @param cost_level Positive common cost for `cost_rule = "common"`.
#' @param cost_range Positive increasing endpoints for uniform cost draws.
#' @param dirichlet_alpha Deprecated alias for `dirichlet_concentration`.
#' @param inside_size Positive total inside quantity used to form estimation
#'   quantities.
#' @param equilibrium_tolerance Maximum relative error allowed when checking
#'   the public-model equilibrium from a perturbed price start.
#' @return A list of class `bayes_calibrated_market` with primitive costs,
#'   full structural truth, estimation inputs and margin mask, and a public-model
#'   equilibrium solve from a perturbed starting price vector. Auction truth
#'   includes `auction_value`, its cost-based latent valuation.
#' @export
bayes_calibrated_market <- function(
    conduct = c("MC", "B", "C", "A", "BF", "CF", "BS", "CS"),
    J = 5L, dirichlet_concentration = NULL, anchor_margin = NULL,
    margin_range = NULL, seed = NULL, k = NULL,
    mode = c("calibrated", "parameter"), alpha = NULL,
    dirichlet_alpha = NULL,
    inside_size = 10000, equilibrium_tolerance = 2e-6,
    costs = NULL, cost_rule = c("uniform", "common"),
    cost_level = 1, cost_range = c(0.7, 1.3),
    conditional_shares = NULL) {
  mode <- match.arg(mode)
  cost_rule <- match.arg(cost_rule)
  conduct <- match.arg(conduct)
  conducts <- .bayes_calibrated_conducts()
  if (!is.character(conduct) || length(conduct) != 1L || is.na(conduct)) {
    stop("'conduct' must be one supported conduct label.", call. = FALSE)
  }
  conduct <- toupper(trimws(conduct))
  if (!conduct %in% names(conducts)) {
    stop("'conduct' must be one of: ", paste(names(conducts), collapse = ", "),
         ".", call. = FALSE)
  }
  spec <- conducts[[conduct]]
  if (length(J) != 1L || !is.numeric(J) || !is.finite(J) ||
      J < 2 || J != floor(J) || J > .Machine$integer.max - 1L) {
    stop("'J' must be an integer of at least two.", call. = FALSE)
  }
  J <- as.integer(J)
  if (!is.null(dirichlet_alpha)) {
    if (!is.null(dirichlet_concentration)) {
      stop("Supply only one of 'dirichlet_concentration' and deprecated 'dirichlet_alpha'.",
           call. = FALSE)
    }
    dirichlet_concentration <- dirichlet_alpha
  }
  if (is.null(dirichlet_concentration) && is.null(conditional_shares)) {
    stop("'dirichlet_concentration' is required.", call. = FALSE)
  }
  if (!is.null(conditional_shares) &&
      (!is.numeric(conditional_shares) || length(conditional_shares) != J ||
       any(!is.finite(conditional_shares)) || any(conditional_shares <= 0) ||
       abs(sum(conditional_shares) - 1) > 1e-10)) {
    stop("'conditional_shares' must be a positive finite vector of length J summing to one.",
         call. = FALSE)
  }
  if (!is.null(dirichlet_concentration) && (!is.numeric(dirichlet_concentration) ||
      !length(dirichlet_concentration) %in% c(1L, J) ||
      any(!is.finite(dirichlet_concentration)) ||
      any(dirichlet_concentration <= 0))) {
    stop("'dirichlet_concentration' must be positive and have length one or J.",
         call. = FALSE)
  }
  if (!is.numeric(inside_size) || length(inside_size) != 1L ||
      !is.finite(inside_size) || inside_size <= 0) {
    stop("'inside_size' must be a positive finite scalar.", call. = FALSE)
  }
  if (!is.null(costs) && (!is.numeric(costs) || length(costs) != J ||
                          any(!is.finite(costs)) || any(costs <= 0))) {
    stop("'costs' must be a positive finite vector of length J.", call. = FALSE)
  }
  if (is.null(costs) && identical(cost_rule, "common") &&
      (!is.numeric(cost_level) || length(cost_level) != 1L ||
       !is.finite(cost_level) || cost_level <= 0)) {
    stop("'cost_level' must be a positive finite scalar.", call. = FALSE)
  }
  if (is.null(costs) && identical(cost_rule, "uniform") &&
      (!is.numeric(cost_range) || length(cost_range) != 2L ||
       any(!is.finite(cost_range)) || any(cost_range <= 0) ||
       cost_range[[1L]] >= cost_range[[2L]])) {
    stop("'cost_range' must be two increasing positive finite values.",
         call. = FALSE)
  }
  if (!is.numeric(equilibrium_tolerance) ||
      length(equilibrium_tolerance) != 1L ||
      !is.finite(equilibrium_tolerance) || equilibrium_tolerance <= 0) {
    stop("'equilibrium_tolerance' must be a positive finite scalar.",
         call. = FALSE)
  }
  if (identical(mode, "calibrated")) {
    if (!is.null(alpha)) {
      stop("Supply 'alpha' only when mode = 'parameter'.", call. = FALSE)
    }
    if (is.null(anchor_margin)) {
      if (!is.numeric(margin_range) || length(margin_range) != 2L ||
          any(!is.finite(margin_range)) || margin_range[[1L]] <= 0 ||
          margin_range[[2L]] >= 1 || margin_range[[1L]] >= margin_range[[2L]]) {
        stop("'margin_range' must be an increasing pair strictly inside (0, 1).",
             call. = FALSE)
      }
    } else if (!is.numeric(anchor_margin) || length(anchor_margin) != 1L ||
               !is.finite(anchor_margin) || anchor_margin <= 0 ||
               anchor_margin >= 1) {
      stop("'anchor_margin' must lie strictly inside (0, 1).", call. = FALSE)
    }
    if (!is.null(anchor_margin) && !is.null(margin_range) &&
        (!is.numeric(margin_range) || length(margin_range) != 2L ||
         any(!is.finite(margin_range)) || margin_range[[1L]] <= 0 ||
         margin_range[[2L]] >= 1 || margin_range[[1L]] >= margin_range[[2L]])) {
      stop("When supplied, 'margin_range' must be an increasing pair strictly inside (0, 1).",
           call. = FALSE)
    }
  } else {
    if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
        alpha <= 0) {
      stop("'alpha' must be a positive finite scalar when mode = 'parameter'.",
           call. = FALSE)
    }
    if (!is.null(anchor_margin)) {
      stop("'anchor_margin' is available only when mode = 'calibrated'.",
           call. = FALSE)
    }
  }
  uses_roles <- spec$route %in% c("core_fringe", "stackelberg")
  if (!uses_roles && !is.null(k)) {
    stop("'k' is only used for BF, CF, BS, and CS conduct.", call. = FALSE)
  }
  if (uses_roles) {
    if (is.null(k)) k <- max(1L, floor(J / 2))
    if (!is.numeric(k) || length(k) != 1L || !is.finite(k) ||
        k != floor(k)) {
      stop("'k' must be a finite integer.", call. = FALSE)
    }
    if (k < 0 || k > J ||
        (identical(spec$route, "stackelberg") && (k < 1 || k >= J))) {
      stop("'k' must lie in 0:J for fringe conduct and 1:(J - 1) for Stackelberg conduct.",
           call. = FALSE)
    }
    k <- as.integer(k)
  } else {
    k <- NA_integer_
  }

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  max_seed <- .Machine$integer.max - 1L
  if (is.null(seed)) {
    seed <- sample.int(max_seed, 1L)
  } else if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
             seed != floor(seed) || seed < 1 || seed > max_seed) {
    stop("'seed' must be an integer between 1 and ", max_seed, ".",
         call. = FALSE)
  }
  seed <- as.integer(seed)
  set.seed(seed)
  concentration <- if (is.null(dirichlet_concentration)) rep(NA_real_, J) else
    if (length(dirichlet_concentration) == 1L) {
    rep(dirichlet_concentration, J)
  } else as.numeric(dirichlet_concentration)
  geometry <- .bayes_calibrated_geometry(
    J, concentration, inside_size, conditional_shares
  )
  if (identical(mode, "calibrated")) {
    if (is.null(anchor_margin)) {
      anchor_margin <- stats::runif(
        1L, min = margin_range[[1L]], max = margin_range[[2L]]
      )
    }
  }
  primitive_costs <- if (!is.null(costs)) as.numeric(costs) else
    if (identical(cost_rule, "common")) rep(cost_level, J) else
      stats::runif(J, cost_range[[1L]], cost_range[[2L]])
  names(primitive_costs) <- geometry$labels
  unit_markup <- .bayes_calibrated_unit_markup(spec, geometry, k, inside_size)
  fitted_alpha <- if (identical(mode, "parameter")) alpha else {
    anchor <- geometry$anchor_index
    unit_markup[[anchor]] * (1 - anchor_margin) /
      (anchor_margin * primitive_costs[[anchor]])
  }
  if (!is.finite(fitted_alpha) || fitted_alpha <= 0) {
    stop("Structural alpha implied by primitive costs is invalid.", call. = FALSE)
  }
  geometry$prices <- stats::setNames(
    as.numeric(primitive_costs + unit_markup / fitted_alpha),
    geometry$labels
  )
  primary <- .bayes_calibrated_model(
    spec, geometry, alpha = fitted_alpha, costs = primitive_costs,
    k = k, inside_size = inside_size,
    price_tolerance = equilibrium_tolerance,
    anchor_margin = if (identical(mode, "calibrated")) anchor_margin else NULL
  )
  observed <- seq_len(J) == geometry$anchor_index
  data_frame <- function(g, margins, mask) {
    data.frame(
      product = g$labels,
      prices = unname(g$prices),
      quantities = unname(g$quantities),
      shares = unname(g$shares),
      conditional_shares = unname(g$conditional_shares),
      b = unname(g$b),
      margins = if (is.null(mask)) unname(margins) else {
        replace(unname(margins), !mask, NA_real_)
      },
      market = rep("market_1", J),
      year = rep("1", J),
      group = g$group,
      margin_observed = if (is.null(mask)) rep(TRUE, J) else mask,
      stringsAsFactors = FALSE
    )
  }
  truth_data <- data_frame(geometry, primary$margins, NULL)
  estimation_data <- data_frame(geometry, primary$margins, observed)
  estimation_args <- .bayes_calibrated_data(
    geometry, primary$margins, spec, k, inside_size
  )
  out <- list(
    design = list(
      conduct = conduct,
      demand_model = "logit",
      supply_model = spec$supply_model,
      structure = spec$structure,
      mode = mode,
      J = J,
      k = k,
      anchor_product = geometry$labels[[geometry$anchor_index]],
      anchor_index = geometry$anchor_index,
      anchor_price = geometry$prices[[geometry$anchor_index]],
      anchor_margin = if (identical(mode, "calibrated")) anchor_margin else NA_real_,
      realized_anchor_margin = primary$margins[[geometry$anchor_index]],
      dirichlet_concentration = concentration,
      cost_rule = if (is.null(costs)) cost_rule else "supplied",
      cost_level = if (is.null(costs) && identical(cost_rule, "common"))
        cost_level else NA_real_,
      cost_range = if (is.null(costs) && identical(cost_rule, "uniform"))
        cost_range else NULL,
      role_ranking = geometry$group[geometry$ranking],
      inside_size = inside_size,
      seed = seed,
      fresh_design = "same generated market; re-solved from a perturbed start"
    ),
    truth = list(
      data = truth_data,
      alpha = primary$alpha,
      s0 = geometry$s0,
      s0_identity = geometry$s0_identity,
      conditional_shares = geometry$conditional_shares,
      centered_log_shares = geometry$b,
      shares = geometry$shares,
      margins = primary$margins,
      costs = primary$costs,
      unit_markup = unit_markup,
      absolute_margins = stats::setNames(
        as.numeric(geometry$prices - primitive_costs), geometry$labels
      ),
      meanval = stats::setNames(
        as.numeric(if (identical(spec$antitrust_conduct, "auction2nd"))
          geometry$b else
          log(geometry$shares / geometry$s0) +
            fitted_alpha * geometry$prices), geometry$labels
      ),
      auction_value = if (identical(spec$antitrust_conduct, "auction2nd"))
        stats::setNames(as.numeric(geometry$b + fitted_alpha * primitive_costs),
                        geometry$labels) else NULL,
      ownership = stats::setNames(geometry$group, geometry$labels),
      core = if (identical(spec$route, "core_fringe"))
        geometry$group[geometry$ranking[seq_len(k)]] else character(),
      leaders = if (identical(spec$route, "stackelberg"))
        geometry$group[geometry$ranking[seq_len(k)]] else character(),
      calibration_source = primary$calibration_source,
      equilibrium = primary$equilibrium,
      model = primary$model,
      source_fit = primary$source_fit
    ),
    estimation = list(
      data = estimation_data,
      margin_observed = observed,
      s0_true = geometry$s0,
      conditional_shares = geometry$conditional_shares,
      inputs = estimation_args
    ),
    fresh_oracle = c(
      primary$equilibrium,
      list(source = "same-market public model solve from perturbed starting prices")
    )
  )
  class(out) <- c("bayes_calibrated_market", "list")
  out
}

#' Re-solve a calibrated market under new costs and product quality
#'
#' Keeps conduct, alpha, ownership, and any core or leader roles fixed. For
#' auction conduct, `meanval` denotes the cost-independent latent valuation
#' (`market$truth$auction_value`); for other conduct it denotes the Logit
#' quality intercept (`market$truth$meanval`). The public model is freshly
#' solved from perturbed baseline prices, then checked against its FOCs.
#'
#' @param market A result from [bayes_calibrated_market()].
#' @param costs Positive physical marginal costs, in product order.
#' @param meanval Product quality or auction valuation, in product order.
#' @param equilibrium_tolerance Maximum relative cost/FOC error.
#' @return A `bayes_calibrated_market` with updated standardized truth and
#'   fresh-equilibrium diagnostics.
#' @export
bayes_resolve_market <- function(market, costs = market$truth$costs,
                                 meanval = if (identical(market$design$conduct, "A"))
                                   market$truth$auction_value else market$truth$meanval,
                                 equilibrium_tolerance = 2e-6) {
  if (!inherits(market, "bayes_calibrated_market"))
    stop("'market' must be a calibrated market.", call. = FALSE)
  J <- market$design$J
  if (!is.numeric(costs) || length(costs) != J ||
      any(!is.finite(costs) | costs <= 0))
    stop("'costs' must be a positive finite vector of length J.", call. = FALSE)
  if (!is.numeric(meanval) || length(meanval) != J || any(!is.finite(meanval)))
    stop("'meanval' must be a finite vector of length J.", call. = FALSE)
  if (!is.numeric(equilibrium_tolerance) ||
      length(equilibrium_tolerance) != 1L ||
      !is.finite(equilibrium_tolerance) || equilibrium_tolerance <= 0)
    stop("'equilibrium_tolerance' must be positive and finite.", call. = FALSE)
  labels <- as.character(market$truth$data$product)
  costs <- stats::setNames(as.numeric(costs), labels)
  quality <- stats::setNames(as.numeric(meanval), labels)
  model <- market$truth$model
  alpha <- market$truth$alpha
  is_auction <- identical(market$design$conduct, "A")
  model@mcPre <- model@mcPost <- unname(costs)
  model <- antitrust::initialize_cost_state(model)
  model@slopes$meanval <- if (is_auction) unname(quality - alpha * costs) else
    unname(quality)
  start <- as.numeric(market$truth$data$prices) * 1.13 + 0.037
  model@priceStart <- start
  prices <- as.numeric(antitrust::calcPrices(model, preMerger = TRUE))
  if (length(prices) != J || any(!is.finite(prices) | prices <= 0))
    stop("Fresh equilibrium solve returned invalid prices.", call. = FALSE)
  model@pricePre <- prices
  model@prices <- prices
  shares <- as.numeric(antitrust::calcShares(model, preMerger = TRUE))
  if (length(shares) != J || any(!is.finite(shares) | shares <= 0) ||
      sum(shares) >= 1)
    stop("Fresh equilibrium solve returned invalid shares.", call. = FALSE)
  model@shares <- shares
  model@shareInside <- sum(shares)
  model@margins <- (prices - costs) / prices
  verify <- model
  verify@priceStart <- prices * 0.87 + 0.019
  check_prices <- as.numeric(antitrust::calcPrices(verify, preMerger = TRUE))
  verify@pricePre <- check_prices
  check_shares <- as.numeric(antitrust::calcShares(verify, preMerger = TRUE))
  price_error <- max(abs(check_prices - prices) / pmax(prices, 1e-12))
  share_error <- max(abs(check_shares - shares))
  margins <- as.numeric(antitrust::calcMargins(model, preMerger = TRUE,
                                               level = FALSE))
  markup <- as.numeric(antitrust::calcMargins(model, preMerger = TRUE,
                                              level = TRUE))
  recovered <- as.numeric(antitrust::calcMC(model, preMerger = TRUE))
  cost_error <- max(abs(recovered - costs) / pmax(1, abs(costs)))
  foc_error <- max(abs(prices - costs - markup))
  if (length(margins) != J || any(!is.finite(margins) | margins <= 0 | margins >= 1) ||
      !is.finite(cost_error) || !is.finite(foc_error) ||
      cost_error > equilibrium_tolerance || foc_error > equilibrium_tolerance ||
      !is.finite(price_error) || price_error > equilibrium_tolerance ||
      !is.finite(share_error) || share_error > equilibrium_tolerance)
    stop("Fresh equilibrium failed margin, cost, or FOC validation.", call. = FALSE)
  role_error <- NA_real_
  if (market$design$conduct %in% c("BF", "CF"))
    role_error <- as.numeric(coordination::core_fringe_residuals(model, preMerger = TRUE)$max)
  if (market$design$conduct %in% c("BS", "CS"))
    role_error <- as.numeric(coordination::stackelberg_residuals(model, preMerger = TRUE)$max)
  if (!is.na(role_error) && (!is.finite(role_error) || role_error > 1e-6))
    stop("Fresh equilibrium rich-model FOC residual exceeds tolerance.", call. = FALSE)
  out <- market
  s0 <- 1 - sum(shares)
  r <- shares / sum(shares)
  names(prices) <- names(shares) <- names(margins) <- names(r) <- labels
  out$truth$model <- model
  out$truth$source_fit <- NULL
  out$truth$costs <- costs
  out$truth$shares <- shares
  out$truth$conditional_shares <- r
  out$truth$s0 <- out$truth$s0_identity <- s0
  out$truth$margins <- margins
  out$truth$unit_markup <- stats::setNames(markup * alpha, labels)
  out$truth$absolute_margins <- stats::setNames(prices - costs, labels)
  out$truth$meanval <- stats::setNames(as.numeric(model@slopes$meanval), labels)
  out$truth$auction_value <- if (is_auction) quality else NULL
  out$truth$centered_log_shares <- stats::setNames(log(r) - mean(log(r)), labels)
  out$truth$data$prices <- unname(prices)
  out$truth$data$shares <- unname(shares)
  out$truth$data$conditional_shares <- unname(r)
  out$truth$data$quantities <- unname(market$design$inside_size * r)
  out$truth$data$margins <- unname(margins)
  out$truth$data$b <- unname(out$truth$centered_log_shares)
  out$truth$equilibrium <- list(
    prices = prices, shares = shares, max_cost_error = cost_error,
    max_relative_error = price_error, max_share_error = share_error,
    max_foc_residual = if (is.na(role_error)) foc_error else max(foc_error, role_error),
    used_perturbed_start = TRUE, start = start, passed = TRUE,
    gates = list(primitive_costs = TRUE, perturbed_price_solve = TRUE,
                 solved_share_reproduction = TRUE, baseline_foc = TRUE)
  )
  out$fresh_oracle <- out$truth$equilibrium
  out$design$realized_anchor_margin <- margins[[market$design$anchor_index]]
  out$design$anchor_price <- prices[[market$design$anchor_index]]
  out$estimation$data <- out$truth$data
  out$estimation$data$margin_observed <- market$estimation$margin_observed
  out$estimation$data$margins[!market$estimation$margin_observed] <- NA_real_
  out$estimation$conditional_shares <- r
  out$estimation$s0_true <- s0
  out$estimation$inputs$prices <- unname(prices)
  out$estimation$inputs$quantities <- unname(market$design$inside_size * r)
  out$estimation$inputs$margins <- out$estimation$data$margins
  out
}
