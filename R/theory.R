# Shape-based conduct schedules and profiled distinguishability metrics.

## Absolute tolerance on the aggregate feasibility condition
## sum(q) <= 1 - s0. It absorbs floating-point error when q enumerates every
## inside firm (sum(q) == 1 - s0) without admitting materially infeasible data.
.conduct_share_tolerance <- 1e-8

#' Conduct Schedule
#'
#' Evaluates the inverse-margin schedule \eqn{g_k(q,s_0)} used by the package.
#' For Logit, inverse absolute margins are \eqn{alpha g_k}.  The CES schedules
#' are returned on their natural inverse-Lerner scale.
#'
#' @param q Numeric vector of unconditional firm shares. The vector may describe
#'   a subset of the market's firms, but the shares must be jointly feasible:
#'   \code{sum(q) <= 1 - s0} (up to a small numerical tolerance), with equality
#'   when \code{q} enumerates every inside firm.
#' @param model Conduct model name.
#' @param s0 Outside quantity share for Logit schedules and outside revenue
#'   share for CES schedules.
#' @param demand_model \code{"logit"} or \code{"ces"}.
#' @param gamma CES elasticity, required for CES schedules.
#' @param bargaining_power Buyer-side Nash bargaining weight.
#' @return A numeric vector with the same length as \code{q}.
#' @export
conduct_schedule <- function(q, model, s0 = 0.2,
                             demand_model = c("logit", "ces"),
                             gamma = 3, bargaining_power = 0.5) {
  demand_model <- match.arg(demand_model)
  model <- canonical_supply_name(model)
  if (!is.numeric(q) || any(!is.finite(q)) || any(q <= 0) || any(q >= 1)) {
    stop("'q' must contain finite shares strictly between zero and one.")
  }
  if (length(s0) != 1L || !is.finite(s0) || s0 <= 0 || s0 >= 1) {
    stop("'s0' must be a finite scalar strictly between zero and one.")
  }
  if (sum(q) > 1 - s0 + .conduct_share_tolerance) {
    stop("Firm shares must satisfy sum(q) <= 1 - s0; they cannot exceed ",
         "the total inside share.")
  }
  if (length(bargaining_power) != 1L || !is.finite(bargaining_power) ||
      bargaining_power < 0 || bargaining_power >= 1) {
    stop("'bargaining_power' must lie in [0, 1).")
  }

  if (demand_model == "logit") {
    return(switch(model,
      bertrand = 1 - q,
      monopolistic = rep(1, length(q)),
      cournot = s0 / (s0 + q),
      auction = q / -log1p(-q),
      bargaining = {
        kappa <- bargaining_power / (1 - bargaining_power)
        (kappa * q - (1 - q) * log1p(-q)) / -log1p(-q)
      }
    ))
  }

  if (length(gamma) != 1L || !is.finite(gamma) || gamma <= 1) {
    stop("'gamma' must exceed one for CES schedules.")
  }
  switch(model,
    bertrand = gamma - (gamma - 1) * q,
    monopolistic = rep(gamma, length(q)),
    cournot = 1 / (1 / gamma +
      (gamma - 1) * q / (gamma * (1 + (gamma - 1) * s0))),
    auction = 1 / (1 - (1 - q)^(1 / (gamma - 1))),
    bargaining = (gamma - (gamma - 1) * q) / (1 - bargaining_power)
  )
}

#' Local Conduct-Schedule Slopes
#'
#' Returns the small-share slopes that govern local conduct separation under
#' Logit demand.
#'
#' @param s0 Outside share.
#' @return A named numeric vector.
#' @export
conduct_local_slopes <- function(s0 = 0.2) {
  if (length(s0) != 1L || !is.finite(s0) || s0 <= 0 || s0 >= 1) {
    stop("'s0' must be a finite scalar strictly between zero and one.")
  }
  c(monopolistic = 0, auction = -0.5, bertrand = -1, cournot = -1 / s0)
}

#' Profiled Conduct Information Distance
#'
#' Computes the weighted Gaussian/KL approximation to pairwise conduct
#' separation after profiling the comparison model over its scale and,
#' optionally, its outside share.
#'
#' @param q Numeric vector of unconditional firm shares under \code{model_k}.
#'   It may describe a subset of the firms. A profiled outside share must keep
#'   \code{sum(q) <= 1 - s0}, so it is bounded above by \code{1 - sum(q)}.
#' @param model_k,model_l Conduct model names.
#' @param alpha_k Scale for the data-generating schedule.
#' @param s0_k Outside share for the data-generating schedule.
#' @param sigma_margin Margin-error standard deviation.
#' @param weights Optional non-negative observation weights.
#' @param demand_model \code{"logit"} or \code{"ces"}.
#' @param gamma CES elasticity when \code{demand_model = "ces"}.
#' @param profile_s0 Logical; profile the comparison model over \code{s0}.
#' @param s0_bounds Two-element admissible interval for profiled outside share,
#'   further truncated above at \code{1 - sum(q)}. Ignored when
#'   \code{profile_s0 = FALSE}.
#' @param bargaining_power_k,bargaining_power_l Bargaining weights when used.
#' @return A list containing the minimized distance and profiled parameters.
#' @export
conduct_distance <- function(q, model_k, model_l, alpha_k = 1, s0_k = 0.2,
                             sigma_margin = 1, weights = NULL,
                             profile_s0 = TRUE,
                             s0_bounds = c(0.001, 0.95),
                             demand_model = c("logit", "ces"), gamma = 3,
                             bargaining_power_k = 0.5,
                             bargaining_power_l = 0.5) {
  demand_model <- match.arg(demand_model)
  model_k <- canonical_supply_name(model_k)
  model_l <- canonical_supply_name(model_l)
  if (length(alpha_k) != 1L || !is.finite(alpha_k) || alpha_k <= 0) {
    stop("'alpha_k' must be a positive finite scalar.")
  }
  if (length(sigma_margin) != 1L || !is.finite(sigma_margin) ||
      sigma_margin <= 0) {
    stop("'sigma_margin' must be a positive finite scalar.")
  }
  if (is.null(weights)) weights <- rep(1, length(q))
  if (length(weights) != length(q) || any(!is.finite(weights)) ||
      any(weights < 0) || sum(weights) <= 0) {
    stop("'weights' must be finite, non-negative, and conformable with 'q'.")
  }
  if (length(s0_bounds) != 2L || any(!is.finite(s0_bounds)) ||
      s0_bounds[1] <= 0 || s0_bounds[2] >= 1 ||
      s0_bounds[1] >= s0_bounds[2]) {
    stop("'s0_bounds' must be an increasing interval inside (0, 1).")
  }

  target <- alpha_k * conduct_schedule(
    q, model_k, s0_k, demand_model = demand_model, gamma = gamma,
    bargaining_power = bargaining_power_k
  )
  ## Unconditional shares of any subset of firms must fit in the inside share:
  ## sum(q) <= 1 - s0, so the profiled s0 is bounded by 1 - sum(q), not by the
  ## largest single share. q need not enumerate every firm.
  upper_s0 <- min(s0_bounds[2], 1 - sum(q))
  if (profile_s0 && upper_s0 <= s0_bounds[1]) {
    stop("The profiled outside-share interval is incompatible with 'q'.")
  }

  objective <- function(par) {
    alpha_l <- exp(par[1])
    s0_l <- if (profile_s0) stats::plogis(par[2]) else s0_k
    fitted <- alpha_l * conduct_schedule(
      q, model_l, s0_l, demand_model = demand_model, gamma = gamma,
      bargaining_power = bargaining_power_l
    )
    0.5 * sum(weights * ((target - fitted) / sigma_margin)^2)
  }

  start_s0 <- if (profile_s0) {
    min(max(s0_k, s0_bounds[1] + 1e-8), upper_s0 - 1e-8)
  } else s0_k
  start_schedule <- conduct_schedule(
    q, model_l, start_s0, demand_model = demand_model, gamma = gamma,
    bargaining_power = bargaining_power_l
  )
  start_alpha <- sum(weights * target) / sum(weights * start_schedule)
  ## With s0 fixed the second coordinate is inert; give it harmless bounds.
  start <- c(log(max(start_alpha, 1e-8)),
             if (profile_s0) stats::qlogis(start_s0) else 0)
  lower <- c(log(1e-8), if (profile_s0) stats::qlogis(s0_bounds[1]) else -1)
  upper <- c(log(1e8), if (profile_s0) stats::qlogis(upper_s0) else 1)
  fit <- stats::optim(
    start, objective, method = "L-BFGS-B", lower = lower, upper = upper
  )

  structure(list(
    distance = unname(fit$value),
    alpha_l = exp(fit$par[1]),
    s0_l = if (profile_s0) stats::plogis(fit$par[2]) else s0_k,
    model_k = model_k,
    model_l = model_l,
    demand_model = demand_model,
    gamma = gamma,
    alpha_k = alpha_k,
    s0_k = s0_k,
    profile_s0 = profile_s0,
    convergence = fit$convergence
  ), class = "conduct_distance")
}

#' Pairwise Conduct Information Matrix
#'
#' @param q Firm shares.
#' @param models Conduct model names.
#' @param ... Arguments passed to \code{\link{conduct_distance}}.
#' @return A square matrix of directed profiled distances.
#' @export
conduct_information_matrix <- function(
    q, models = c("bertrand", "cournot", "auction", "monopolistic"), ...) {
  models <- canonical_supply_name(models)
  out <- matrix(0, length(models), length(models), dimnames = list(models, models))
  for (i in seq_along(models)) {
    for (j in seq_along(models)) {
      if (i != j) {
        out[i, j] <- conduct_distance(q, models[i], models[j], ...)$distance
      }
    }
  }
  out
}
