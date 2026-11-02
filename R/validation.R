# Fast schedule-level validation aligned with the distinguishability theory.

log_sum_exp_r <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

profile_schedule_fit <- function(y, q, model, sigma_margin, profile_s0,
                                 s0_start, bargaining_power = 0.5) {
  ## Feasibility requires sum(q) <= 1 - s0 for the enumerated firms.
  upper_s0 <- min(0.95, 1 - sum(q))
  objective <- function(par) {
    alpha <- exp(par[1])
    s0 <- if (profile_s0) min(stats::plogis(par[2]), upper_s0) else s0_start
    if (s0 <= 0.001) return(1e100)
    mu <- alpha * conduct_schedule(
      q, model, s0 = s0, bargaining_power = bargaining_power
    )
    sum((y - mu)^2) / sigma_margin^2
  }
  fit <- stats::optim(
    c(0, stats::qlogis(min(s0_start, upper_s0 - 1e-6))),
    objective,
    method = "L-BFGS-B",
    lower = c(log(1e-6), stats::qlogis(0.001)),
    upper = c(log(1e6), stats::qlogis(upper_s0))
  )
  list(
    rss_scaled = fit$value,
    alpha = exp(fit$par[1]),
    s0 = if (profile_s0) stats::plogis(fit$par[2]) else s0_start
  )
}

#' Fast Theory-Guided Conduct Validation
#'
#' Simulates inverse margins under one of the paper's targeted information
#' designs and compares profiled schedule fits. This is the inexpensive first
#' validation layer; representative cells can subsequently be run through
#' \code{bayes_est} for end-to-end Stan validation.
#'
#' @param design Weak, strong, Cournot-separating, outside-share,
#'   bargaining, or misspecified design.
#' @param true_model Data-generating conduct model.
#' @param models Candidate model set.
#' @param n_rep Number of replications.
#' @param seed Random seed.
#' @param model_priors Optional candidate-model probabilities.
#' @param bargaining_power Bargaining weight for bargaining schedules.
#' @return A data frame with posterior-score weights, correct-selection
#'   indicators, and the design's minimum theoretical information distance.
#' @export
conduct_validation <- function(
    design = c("weak", "strong", "cournot", "outside_share",
               "bargaining", "misspecified"),
    true_model = "bertrand",
    models = c("bertrand", "cournot", "auction", "monopolistic"),
    n_rep = 100L, seed = 1L, model_priors = NULL,
    bargaining_power = 0.5) {
  design <- match.arg(design)
  true_model <- canonical_supply_name(true_model)
  models <- canonical_supply_name(models)
  if (length(n_rep) != 1L || !is.finite(n_rep) || n_rep < 1 ||
      n_rep != as.integer(n_rep)) {
    stop("'n_rep' must be a positive integer.")
  }
  if (is.null(model_priors)) model_priors <- rep(1 / length(models), length(models))
  if (length(model_priors) != length(models) || any(!is.finite(model_priors)) ||
      any(model_priors <= 0)) {
    stop("'model_priors' must contain one positive finite value per model.")
  }
  model_priors <- model_priors / sum(model_priors)

  spec <- switch(design,
    weak = list(q = c(0.01, 0.02, 0.035, 0.05, 0.075), s0 = 0.20,
                sigma = 0.15, profile_s0 = TRUE),
    strong = list(q = c(0.03, 0.07, 0.12, 0.20, 0.28), s0 = 0.30,
                  sigma = 0.04, profile_s0 = FALSE),
    cournot = list(q = c(0.03, 0.05, 0.08, 0.12, 0.20), s0 = 0.05,
                   sigma = 0.06, profile_s0 = FALSE),
    outside_share = list(q = c(0.03, 0.07, 0.12, 0.20, 0.28), s0 = 0.30,
                         sigma = 0.08, profile_s0 = TRUE),
    bargaining = list(q = c(0.03, 0.08, 0.15, 0.25), s0 = 0.35,
                       sigma = 0.06, profile_s0 = TRUE),
    misspecified = list(q = c(0.03, 0.07, 0.12, 0.20, 0.28), s0 = 0.25,
                        sigma = 0.08, profile_s0 = TRUE)
  )
  if (sum(spec$q) > 1 - spec$s0 + .conduct_share_tolerance)
    stop("Internal design has infeasible shares.")

  true_mu <- conduct_schedule(
    spec$q, true_model, spec$s0, bargaining_power = bargaining_power
  )
  if (design == "misspecified") {
    true_mu <- 0.5 * conduct_schedule(spec$q, "bertrand", spec$s0) +
      0.5 * conduct_schedule(spec$q, "cournot", spec$s0)
  }
  distances <- vapply(models, function(model) {
    if (model == true_model && design != "misspecified") return(0)
    conduct_distance(
      spec$q, true_model, model, s0_k = spec$s0,
      sigma_margin = spec$sigma, profile_s0 = spec$profile_s0,
      bargaining_power_k = bargaining_power,
      bargaining_power_l = bargaining_power
    )$distance
  }, numeric(1))
  positive_distance <- distances[distances > 1e-10]
  min_information <- if (length(positive_distance)) min(positive_distance) else 0

  set.seed(seed)
  rows <- vector("list", as.integer(n_rep))
  for (r in seq_len(n_rep)) {
    y <- true_mu + stats::rnorm(length(true_mu), 0, spec$sigma)
    fits <- lapply(models, function(model) profile_schedule_fit(
      y, spec$q, model, spec$sigma, spec$profile_s0, spec$s0,
      bargaining_power = bargaining_power
    ))
    score <- vapply(fits, function(x) -0.5 * x$rss_scaled, numeric(1)) +
      log(model_priors)
    probability <- exp(score - log_sum_exp_r(score))
    rows[[r]] <- data.frame(
      replication = r,
      design = design,
      true_model = true_model,
      selected_model = models[which.max(probability)],
      true_model_probability = if (true_model %in% models) {
        probability[match(true_model, models)]
      } else NA_real_,
      correct_selection = models[which.max(probability)] == true_model,
      min_information = min_information,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}
