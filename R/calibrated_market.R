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

.bayes_calibrated_geometry <- function(J, concentration, inside_size) {
  shape <- if (length(concentration) == 1L) {
    rep(concentration, J)
  } else {
    concentration
  }
  weights <- stats::rgamma(J, shape = shape, rate = 1)
  if (any(!is.finite(weights)) || any(weights <= 0) ||
      !is.finite(sum(weights)) || sum(weights) <= 0) {
    stop("Dirichlet draw was numerically degenerate; increase the concentration.",
         call. = FALSE)
  }
  conditional_shares <- unname(weights / sum(weights))
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
    geometry$b + alpha
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

.bayes_calibrated_public_margin <- function(conduct, geometry, alpha,
                                            inside_size, price_start) {
  fit <- .bayes_calibrated_specify(
    conduct, geometry, alpha, inside_size, price_start
  )
  margins <- as.numeric(antitrust::calcMargins(
    fit@model, preMerger = TRUE, level = FALSE
  ))
  margins[[geometry$anchor_index]]
}

.bayes_calibrated_root_alpha <- function(conduct, geometry, target_margin,
                                         initial_alpha, inside_size,
                                         price_start) {
  if (!is.finite(initial_alpha) || initial_alpha <= 0) {
    initial_alpha <- 1 / target_margin
  }
  objective <- function(log_alpha) {
    alpha <- exp(log_alpha)
    .bayes_calibrated_public_margin(
      conduct, geometry, alpha, inside_size, price_start
    ) - target_margin
  }
  center <- log(initial_alpha)
  f_center <- objective(center)
  if (is.finite(f_center) && abs(f_center) <= 1e-11) {
    return(exp(center))
  }
  lower <- upper <- center
  f_lower <- f_upper <- f_center
  for (i in seq_len(80L)) {
    if (is.finite(f_lower) && is.finite(f_upper) &&
        sign(f_lower) != sign(f_upper)) break
    if (!is.finite(f_lower) || f_lower < 0) {
      lower <- lower - log(2)
      f_lower <- tryCatch(objective(lower), error = function(e) NA_real_)
    }
    if (!is.finite(f_upper) || f_upper > 0) {
      upper <- upper + log(2)
      f_upper <- tryCatch(objective(upper), error = function(e) NA_real_)
    }
  }
  if (!is.finite(f_lower) || !is.finite(f_upper) ||
      sign(f_lower) == sign(f_upper)) {
    stop("Could not bracket positive alpha for the requested public-model margin.",
         call. = FALSE)
  }
  root <- stats::uniroot(
    objective, interval = sort(c(lower, upper)), tol = 1e-11
  )$root
  alpha <- exp(root)
  residual <- objective(root)
  if (!is.finite(alpha) || alpha <= 0 || !is.finite(residual) ||
      abs(residual) > 1e-8) {
    stop("Public-model alpha solve did not reproduce the anchor margin.",
         call. = FALSE)
  }
  alpha
}

.bayes_calibrated_model <- function(conduct, geometry, mode, anchor_margin,
                                    alpha, k, inside_size, price_tolerance) {
  J <- length(geometry$prices)
  margins_observed <- rep(NA_real_, J)
  anchor_index <- geometry$anchor_index
  if (identical(mode, "calibrated")) {
    margins_observed[[anchor_index]] <- anchor_margin
  }
  price_start <- unname(geometry$prices * 1.13 + 0.037)
  names(price_start) <- geometry$labels

  if (identical(conduct$route, "core_fringe")) {
    model <- coordination::core_fringe(
      prices = unname(geometry$prices),
      shares = unname(geometry$shares),
      margins = margins_observed,
      ownerPre = geometry$group,
      corePre = if (k == 0L) character() else
        geometry$group[geometry$ranking[seq_len(k)]],
      demand = "logit", conduct = conduct$antitrust_conduct,
      output = TRUE, insideSize = inside_size, normIndex = NA_real_,
      priceOutside = 0, priceStart = price_start, labels = geometry$labels,
      alpha = if (identical(mode, "parameter")) alpha else NULL,
      control.equ = list(implicitCheck = FALSE)
    )
    fitted_alpha <- as.numeric(model@diagnostics$alpha)
    calibration_source <- if (identical(mode, "calibrated")) {
      "coordination::core_fringe observed-margin calibration"
    } else "coordination::core_fringe fixed alpha"
  } else if (identical(conduct$route, "stackelberg")) {
    model <- coordination::stackelberg(
      prices = unname(geometry$prices),
      shares = unname(geometry$shares),
      margins = margins_observed,
      ownerPre = geometry$group,
      leadersPre = geometry$group[geometry$ranking[seq_len(k)]],
      demand = "logit", conduct = conduct$antitrust_conduct,
      output = TRUE, insideSize = inside_size, normIndex = NA_real_,
      priceOutside = 0, priceStart = price_start, labels = geometry$labels,
      alpha = if (identical(mode, "parameter")) alpha else NULL,
      control.equ = list(implicitCheck = FALSE)
    )
    fitted_alpha <- as.numeric(model@diagnostics$alpha)
    calibration_source <- if (identical(mode, "calibrated")) {
      "coordination::stackelberg observed-margin calibration"
    } else "coordination::stackelberg fixed alpha"
  } else if (identical(mode, "calibrated")) {
    fit <- antitrust::calibrate(
      demand = "logit", conduct = conduct$antitrust_conduct,
      prices = unname(geometry$prices), shares = unname(geometry$shares),
      margins = margins_observed, ownerPre = geometry$group,
      output = TRUE, insideSize = inside_size, labels = geometry$labels
    )
    model <- fit@model
    candidates <- list(fit@parameters$alpha, model@slopes$alpha)
    valid <- vapply(candidates, function(x) {
      is.numeric(x) && length(x) == 1L && is.finite(x) && x != 0
    }, logical(1))
    if (!any(valid)) {
      stop("antitrust::calibrate did not return a usable Logit alpha.",
           call. = FALSE)
    }
    initial_alpha <- abs(as.numeric(candidates[[which(valid)[[1L]]]]))
    ## Keep antitrust::calibrate as the public initializer, but solve the
    ## one-margin target against antitrust::specify/calcMargins so optimizer
    ## stopping conventions cannot leave the observed margin off-target.
    fitted_alpha <- .bayes_calibrated_root_alpha(
      conduct, geometry, anchor_margin, initial_alpha,
      inside_size, price_start
    )
    fit <- .bayes_calibrated_specify(
      conduct, geometry, fitted_alpha, inside_size, price_start
    )
    model <- fit@model
    calibration_source <- paste(
      "antitrust::calibrate initializer plus public antitrust::specify/calcMargins alpha root"
    )
  } else {
    ## antitrust Logit uses a negative price slope in output markets. Auction
    ## Logit keeps supplied shares fixed while solving, so its mean values are
    ## represented directly by log share odds.
    fit <- .bayes_calibrated_specify(
      conduct, geometry, alpha, inside_size, price_start
    )
    model <- fit@model
    fitted_alpha <- alpha
    calibration_source <- "antitrust::specify fixed alpha"
  }

  if (length(fitted_alpha) != 1L || !is.finite(fitted_alpha) ||
      fitted_alpha <= 0) {
    stop("Public model did not return a positive Logit price coefficient.",
         call. = FALSE)
  }
  fitted_margins <- tryCatch(
    as.numeric(antitrust::calcMargins(model, preMerger = TRUE, level = FALSE)),
    error = function(e) stop("Could not recover public-model margins: ",
                             conditionMessage(e), call. = FALSE)
  )
  costs <- tryCatch(
    as.numeric(antitrust::calcMC(model, preMerger = TRUE)),
    error = function(e) stop("Could not recover public-model costs: ",
                             conditionMessage(e), call. = FALSE)
  )
  if (length(fitted_margins) != J || any(!is.finite(fitted_margins)) ||
      any(fitted_margins <= 0)) {
    stop("Public model returned invalid full-market margins.", call. = FALSE)
  }
  if (length(costs) != J || any(!is.finite(costs))) {
    stop("Public model returned invalid marginal costs.", call. = FALSE)
  }
  if (any(!is.finite(costs)) || any(costs <= 0)) {
    stop("Public model implied non-positive physical marginal costs.",
         call. = FALSE)
  }
  if (identical(mode, "calibrated") &&
      !isTRUE(all.equal(fitted_margins[[anchor_index]], anchor_margin,
                        tolerance = 1e-7))) {
    stop("Public model did not reproduce the supplied anchor margin (got ",
         format(fitted_margins[[anchor_index]], digits = 10),
         ", expected ", format(anchor_margin, digits = 10), ").",
         call. = FALSE)
  }

  ## Re-solve from a perturbed start. This catches stale stored equilibria and
  ## confirms the generated costs support the generated price vector.
  check_model <- model
  if (!("priceStart" %in% methods::slotNames(check_model))) {
    stop("Public model does not expose a perturbed priceStart slot.",
         call. = FALSE)
  }
  check_model@priceStart <- price_start
  used_start <- as.numeric(check_model@priceStart)
  if (length(used_start) != J ||
      !isTRUE(all.equal(used_start, unname(price_start), tolerance = 0))) {
    stop("Public model did not retain the requested perturbed price start.",
         call. = FALSE)
  }
  equilibrium_prices <- tryCatch(
    as.numeric(antitrust::calcPrices(check_model, preMerger = TRUE)),
    error = function(e) stop("Fresh equilibrium solve failed: ",
                             conditionMessage(e), call. = FALSE)
  )
  if (length(equilibrium_prices) != J || any(!is.finite(equilibrium_prices))) {
    stop("Fresh equilibrium solve returned invalid prices.", call. = FALSE)
  }
  max_relative_error <- max(abs(equilibrium_prices - geometry$prices) /
                              pmax(abs(geometry$prices), 1e-12))
  if (!is.finite(max_relative_error) || max_relative_error > price_tolerance) {
    stop("Fresh equilibrium solve does not reproduce the market prices (relative error = ",
         format(max_relative_error, digits = 6), ").", call. = FALSE)
  }
  solved_geometry <- geometry
  solved_geometry$prices <- stats::setNames(equilibrium_prices, geometry$labels)
  demand_check <- .bayes_calibrated_specify(
    conduct, solved_geometry, fitted_alpha, inside_size,
    unname(equilibrium_prices * 1.13 + 0.037)
  )
  public_shares <- tryCatch(
    as.numeric(antitrust::calcShares(
      demand_check@model, preMerger = TRUE
    )),
    error = function(e) stop("Could not recompute public-model shares at the solved prices: ",
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
    foc <- coordination::core_fringe_residuals(model, preMerger = TRUE)
    foc_residual <- as.numeric(foc$max)
  } else if (identical(conduct$route, "stackelberg")) {
    foc <- coordination::stackelberg_residuals(model, preMerger = TRUE)
    foc_residual <- as.numeric(foc$max)
  }
  foc_available <- identical(conduct$route, "core_fringe") ||
    identical(conduct$route, "stackelberg")
  if (foc_available && !is.finite(foc_residual)) {
    stop("Coordination baseline FOC residual was not available from its public API.",
         call. = FALSE)
  }
  if (foc_available && foc_residual > 1e-6) {
    stop("Public model baseline FOC residual exceeds tolerance.",
         call. = FALSE)
  }

  gates <- list(
    positive_costs = all(is.finite(costs) & costs > 0),
    anchor_margin = if (identical(mode, "calibrated")) {
      isTRUE(all.equal(fitted_margins[[anchor_index]], anchor_margin,
                       tolerance = 1e-7))
    } else TRUE,
    perturbed_price_solve = max_relative_error <= price_tolerance,
    solved_share_reproduction = share_error <= price_tolerance,
    baseline_foc = if (foc_available) foc_residual <= 1e-6 else NA
  )
  known_gates <- !vapply(gates, is.na, logical(1))
  overall_passed <- if (all(known_gates)) {
    all(unlist(gates, use.names = FALSE))
  } else NA

  list(
    alpha = fitted_alpha,
    margins = stats::setNames(fitted_margins, geometry$labels),
    costs = stats::setNames(costs, geometry$labels),
    model = model,
    calibration_source = calibration_source,
    equilibrium = list(
      prices = stats::setNames(equilibrium_prices, geometry$labels),
      shares = stats::setNames(public_shares, geometry$labels),
      max_relative_error = max_relative_error,
      max_share_error = share_error,
      max_foc_residual = foc_residual,
      gates = gates,
      foc_status = if (foc_available) "checked" else "not exposed by public antitrust API",
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
#' Draws one conduct-independent Logit share geometry, fixes every inside
#' product's price at exactly one, and either calibrates the structural price
#' coefficient from one observed margin or supplies a known coefficient. The
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
#'   length `J`.
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
#' @param dirichlet_alpha Deprecated alias for `dirichlet_concentration`.
#' @param inside_size Positive total inside quantity used to form estimation
#'   quantities.
#' @param equilibrium_tolerance Maximum relative error allowed when checking
#'   the public-model equilibrium from a perturbed price start.
#' @return A list of class `bayes_calibrated_market` with the draw design,
#'   full truth data, estimation inputs and margin mask, and a public-model
#'   equilibrium solve from a perturbed starting price vector for the same
#'   generated baseline.
#' @export
bayes_calibrated_market <- function(
    conduct = c("MC", "B", "C", "A", "BF", "CF", "BS", "CS"),
    J = 5L, dirichlet_concentration = NULL, anchor_margin = NULL,
    margin_range = NULL, seed = NULL, k = NULL,
    mode = c("calibrated", "parameter"), alpha = NULL,
    dirichlet_alpha = NULL,
    inside_size = 10000, equilibrium_tolerance = 2e-6) {
  mode <- match.arg(mode)
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
  if (is.null(dirichlet_concentration)) {
    stop("'dirichlet_concentration' is required.", call. = FALSE)
  }
  if (!is.numeric(dirichlet_concentration) ||
      !length(dirichlet_concentration) %in% c(1L, J) ||
      any(!is.finite(dirichlet_concentration)) ||
      any(dirichlet_concentration <= 0)) {
    stop("'dirichlet_concentration' must be positive and have length one or J.",
         call. = FALSE)
  }
  if (!is.numeric(inside_size) || length(inside_size) != 1L ||
      !is.finite(inside_size) || inside_size <= 0) {
    stop("'inside_size' must be a positive finite scalar.", call. = FALSE)
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
  concentration <- if (length(dirichlet_concentration) == 1L) {
    rep(dirichlet_concentration, J)
  } else as.numeric(dirichlet_concentration)
  geometry <- .bayes_calibrated_geometry(
    J, concentration, inside_size
  )
  if (identical(mode, "calibrated")) {
    if (is.null(anchor_margin)) {
      anchor_margin <- stats::runif(
        1L, min = margin_range[[1L]], max = margin_range[[2L]]
      )
    }
  }
  primary <- .bayes_calibrated_model(
    spec, geometry, mode, anchor_margin,
    alpha = if (identical(mode, "parameter")) alpha else NULL,
    k = k, inside_size = inside_size,
    price_tolerance = equilibrium_tolerance
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
      calibration_source = primary$calibration_source,
      equilibrium = primary$equilibrium,
      model = primary$model
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
