#' Draw posterior likelihood outcomes from a fitted Bayesian model
#'
#' @param fit An estimated \code{BayesFit} returned by \code{bayes_est()}.
#' @param draw_plan Optional reusable state plan returned by
#'   \code{posterior_policy_draws()}. Reusing a plan keeps collapsed quality,
#'   integrated market-intercept, and exact-k states aligned with other
#'   posterior policy work.
#' @param n_draws Number of posterior draws to use when creating a plan. The
#'   default is up to 200 draws. It must agree with a supplied \code{draw_plan}.
#' @param seed Optional nonnegative integer controlling replicated residuals,
#'   and, when creating a new plan, draw selection and latent-state
#'   reconstruction. The prior R random state is restored after a seeded call.
#' @param allow_experimental Explicitly allow exploratory Pathfinder draws.
#' @param share_conditioning \code{"replicate"} (default) draws new log-share
#'   outcomes. \code{"observed"} fixes log shares at their observed values and
#'   draws margins conditional on their residuals. The latter requires the
#'   joint likelihood (\code{likelihood_mode = 1}).
#'
#' @return A \code{posterior_predictive} list containing observed outcomes and
#'   draw-aligned \code{log_shareIn} and \code{log_marginNormInv} matrices,
#'   original posterior draw indices, the state plan, and provenance metadata.
#'   The observation-level \code{margin_observed} mask distinguishes observed
#'   margins from generated imputations. Replicated share entries are
#'   \code{NA} when the fitted likelihood omitted the share equation.
#'   With \code{share_conditioning = "observed"}, the share matrix repeats the
#'   observed log shares for each draw; \code{metadata$share_conditioning} is
#'   \code{"observed"} and \code{metadata$share_replicated} is \code{FALSE}.
#'
#' @details Each row of the outcome matrices corresponds to one posterior
#'   draw. Shared observed predictors remain fixed. Gaussian and Student-t
#'   residual laws follow the Stan likelihood exactly, including the joint
#'   Gaussian or Student-t rho conditional for margins. Collapsed quality and
#'   market-intercept states are reconstructed with the same helpers used by
#'   \code{posterior_policy_draws()}; exact role models reuse the plan's sampled
#'   component state and calculate its state-specific FOC margin mean. The
#'   returned object contains data and metadata only; it does not create plots.
#'   The raw \code{log_shareIn} values are likelihood outcomes. The
#'   \code{share_rep_normalized} matrix is a derived within-market-year softmax;
#'   it is not a likelihood variable. Prediction samples residual outcomes at
#'   the fitted latent states and does not solve a new structural equilibrium.
#'   Normalizing shares therefore does not force the raw log shares to match
#'   the observed values in \code{"replicate"} mode.
#'
#' @seealso \code{posterior_policy_draws()}, \code{posterior_quality()}
#' @export
posterior_predict <- function(fit, draw_plan = NULL, n_draws = NULL,
                              seed = NULL, allow_experimental = FALSE,
                              share_conditioning = c("replicate", "observed")) {
  share_conditioning <- match.arg(share_conditioning)
  if (!methods::is(fit, "BayesFit") || methods::is(fit, "BayesSpecifiedFit")) {
    stop("'fit' must be an estimated BayesFit object.", call. = FALSE)
  }
  if (length(allow_experimental) != 1L || is.na(allow_experimental) ||
      !is.logical(allow_experimental)) {
    stop("'allow_experimental' must be TRUE or FALSE.", call. = FALSE)
  }
  if (share_conditioning == "observed" &&
      .posterior_predict_mask_and_mode(fit)$likelihood_mode != 1L) {
    stop("'share_conditioning = \"observed\"' requires the joint likelihood (likelihood_mode = 1).",
         call. = FALSE)
  }
  .cmdstan_require_policy_inference(fit, allow_experimental)
  if (!is.null(seed)) seed <- .lq_validate_seed(seed)

  ## CmdStanR's draws() lazily caches a draws object in its R6 instance. Work
  ## on an isolated serialized copy so this read-only API leaves the user's
  ## fit (including its reference-class stanfit slot) unchanged.
  prediction_fit <- .posterior_predict_isolated_fit(fit)
  draws <- .posterior_extract(prediction_fit,
                              purpose = "posterior_predictive")
  plan <- if (is.null(draw_plan)) {
    .lq_make_policy_plan(prediction_fit, draws, n_draws = n_draws, seed = seed,
                         role_mapping = "mixture")
  } else {
    .lq_validate_plan(prediction_fit, draw_plan, n_draws = n_draws, draws = draws)
  }
  .posterior_predict_from_draws(prediction_fit, draws, plan, seed = seed,
                                share_conditioning = share_conditioning)
}

.posterior_predict_isolated_fit <- function(fit) {
  out <- fit
  if (!is.null(fit@stanfit)) {
    isolated <- tryCatch(
      unserialize(serialize(fit@stanfit, NULL, version = 3L)),
      error = function(e) {
        stop("Could not isolate the CmdStanR fit for posterior prediction: ",
             conditionMessage(e), call. = FALSE)
      }
    )
    if (!identical(class(isolated), class(fit@stanfit))) {
      stop("Isolating the CmdStanR fit changed its class; refusing to predict.",
           call. = FALSE)
    }
    out@stanfit <- isolated
  }
  out
}

.posterior_predict_code_name <- function(fit, field, choices, codes) {
  value <- fit@metadata[[field]]
  if (is.null(value)) {
    value <- attr(fit@stan_data, field, exact = TRUE)
  }
  if (is.null(value)) {
    code <- fit@stan_data[[field]]
    if (!is.null(code) && length(code) == 1L && is.numeric(code) &&
        is.finite(code) && code == floor(code)) {
      value <- unname(codes[as.character(as.integer(code))])
    }
  }
  if (length(value) != 1L || is.na(value) ||
      !as.character(value) %in% choices) {
    stop("Unable to resolve fitted ", field, " metadata.", call. = FALSE)
  }
  as.character(value)
}

.posterior_predict_mask_and_mode <- function(fit) {
  d <- fit@stan_data
  n <- as.integer(d$N %||% length(d$log_shareIn))
  if (length(n) != 1L || is.na(n) || n < 1L ||
      length(d$log_shareIn) != n || !is.numeric(d$log_shareIn) ||
      any(!is.finite(d$log_shareIn))) {
    stop("Stan data must contain N finite log_shareIn observations.",
         call. = FALSE)
  }
  mode <- as.integer(d$likelihood_mode)
  if (length(mode) != 1L || is.na(mode) || !mode %in% 1:3) {
    stop("Stan data contain a malformed likelihood_mode.", call. = FALSE)
  }
  margin_mask <- d$margin_obs
  if (is.null(margin_mask)) {
    stop("Stan data are missing the margin_obs mask.", call. = FALSE)
  }
  if (length(margin_mask) != n || anyNA(margin_mask) ||
      any(!as.numeric(margin_mask) %in% c(0, 1))) {
    stop("Stan data contain a malformed margin_obs mask.", call. = FALSE)
  }
  margin_mask <- as.logical(as.integer(margin_mask))
  margin_idx <- as.integer(d$margin_idx %||% integer())
  if (anyNA(margin_idx) || anyDuplicated(margin_idx) ||
      any(margin_idx < 1L | margin_idx > n) ||
      !identical(sort(margin_idx), which(margin_mask))) {
    stop("Stan data margin_idx and margin_obs do not identify the same rows.",
         call. = FALSE)
  }
  if (length(d$log_marginNormInv) != n ||
      !is.numeric(d$log_marginNormInv) ||
      any(!is.finite(d$log_marginNormInv[margin_mask]))) {
    stop("Stan data must contain finite observed log_marginNormInv values.",
         call. = FALSE)
  }
  market <- as.integer(d$market %||% rep.int(1L, n))
  n_market <- as.integer(d$N_market %||% max(market))
  if (length(market) != n || anyNA(market) || length(n_market) != 1L ||
      is.na(n_market) || n_market < 1L || any(market < 1L | market > n_market)) {
    stop("Stan data contain a malformed market map.", call. = FALSE)
  }
  list(n = n, likelihood_mode = mode, margin_mask = margin_mask,
       margin_idx = margin_idx, market = market, n_market = n_market,
       observed_share = as.numeric(d$log_shareIn),
       observed_margin = ifelse(margin_mask, as.numeric(d$log_marginNormInv), NA_real_))
}

.posterior_predict_observation_ids <- function(fit, data_info, market_year_idx) {
  data <- fit@data
  take <- function(candidates, fallback) {
    field <- candidates[candidates %in% names(data)]
    if (length(field)) as.character(data[[field[[1L]]]]) else fallback
  }
  if (nrow(data) && nrow(data) != data_info$n) {
    stop("Fitted observation IDs do not align with Stan observations.",
         call. = FALSE)
  }
  product <- take(c("product_id", "product", "group", "upc"),
                  as.character(seq_len(data_info$n)))
  market <- take("market", as.character(data_info$market))
  year <- take(c("year", "market_year"), as.character(market_year_idx))
  out <- data.frame(observation = seq_len(data_info$n), product_id = product,
                    market_id = market, year_id = year,
                    market_year_index = as.integer(market_year_idx),
                    stringsAsFactors = FALSE)
  out
}

.posterior_predict_exact_k_by_market_year <- function(fit, plan) {
  if (is.null(plan$k_cache)) return(NULL)
  setup <- .lq_all_units(fit)
  info <- setup$info
  if (!length(info)) stop("The exact-k draw plan has no simulation units.", call. = FALSE)
  unit_mky <- vapply(info, function(x) as.integer(x$market_year_index), integer(1))
  n_mky <- as.integer(fit@stan_data$N_market_year %||% max(unit_mky))
  if (length(unit_mky) != length(info) || anyNA(unit_mky) ||
      any(unit_mky < 1L | unit_mky > n_mky)) {
    stop("The exact-k draw plan cannot be mapped to market-year states.",
         call. = FALSE)
  }
  if (anyDuplicated(unit_mky)) {
    stop("The exact-k draw plan has multiple simulation units for a market-year.",
         call. = FALSE)
  }
  mky_unit <- match(seq_len(n_mky), unit_mky)
  if (anyNA(mky_unit)) {
    stop("The exact-k draw plan is missing a market-year state.", call. = FALSE)
  }
  out <- matrix(NA_integer_, nrow = length(plan$draw_indices), ncol = n_mky,
                dimnames = list(as.character(plan$draw_indices),
                                as.character(seq_len(n_mky))))
  for (rr in seq_along(plan$draw_indices)) {
    out[rr, ] <- vapply(mky_unit, function(unit_pos) {
      as.integer(.lq_plan_k_for_unit(plan, rr, info[[unit_pos]]))
    }, integer(1))
  }
  if (anyNA(out) || any(out < 0L)) {
    stop("The exact-k draw plan contains invalid states.", call. = FALSE)
  }
  out
}

## Exact posterior prediction uses the shared state-specific FOC helper.
.posterior_predict_exact_margin_mean <- function(fit, draws, draw_indices,
                                                k_by_market_year) {
  n_total <- .posterior_draw_count_from_draws(draws, fit)
  parsed <- list(n_draw = n_total, y = as.numeric(fit@stan_data$log_shareIn))
  .stack_pure_foc_margin_prediction(
    fit, parsed, draws, draw_indices, k_by_market_year,
    realized_k_scope = "market_year"
  )
}
.posterior_predict_mean_state <- function(fit, draws, plan, data_info) {
  d <- fit@stan_data
  n_total <- .posterior_draw_count_from_draws(draws, fit)
  draw_indices <- as.integer(plan$draw_indices)
  if (!length(draw_indices) || anyNA(draw_indices) || any(draw_indices < 1L) ||
      any(draw_indices > n_total) || anyDuplicated(draw_indices)) {
    stop("The draw plan contains invalid posterior draw indices.", call. = FALSE)
  }
  intercept_collapsed <- identical(
    .market_intercept_parameterization_from_fit(fit), "collapsed"
  )
  mean_name <- if (intercept_collapsed) "share_mean_integrated" else "share_mean"
  share_mean <- .posterior_quality_matrix(draws[[mean_name]], data_info$n,
                                           mean_name)
  if (nrow(share_mean) != n_total || any(!is.finite(share_mean))) {
    stop("Posterior share means are malformed.", call. = FALSE)
  }
  share_mean <- share_mean[draw_indices, , drop = FALSE]

  quality_mode <- .lq_quality_mode(fit)
  if (identical(quality_mode, "collapsed")) {
    incidence <- .lq_quality_incidence(fit)
    if (is.null(incidence)) {
      stop("Collapsed quality prediction has no quality incidence map.",
           call. = FALSE)
    }
    quality_index <- as.integer(incidence$obs_to_quality)
    quality <- plan$quality_draws
    if (!is.matrix(quality) || nrow(quality) != length(draw_indices) ||
        length(quality_index) != data_info$n || anyNA(quality_index) ||
        any(quality_index < 1L | quality_index > ncol(quality)) ||
        any(!is.finite(quality))) {
      stop("The draw plan has malformed collapsed-quality states.",
           call. = FALSE)
    }
    share_mean <- share_mean + quality[, quality_index, drop = FALSE]
  }
  if (intercept_collapsed) {
    mu_b <- .posterior_quality_vector(draws$mu_b_market_vec, n_total,
                                      "mu_b_market_vec")
    market_intercepts <- plan$market_intercept_draws
    if (!is.matrix(market_intercepts) ||
        nrow(market_intercepts) != length(draw_indices) ||
        ncol(market_intercepts) != data_info$n_market ||
        any(!is.finite(market_intercepts))) {
      stop("The draw plan has malformed collapsed market-intercept states.",
           call. = FALSE)
    }
    share_mean <- share_mean +
      market_intercepts[, data_info$market, drop = FALSE] -
      matrix(mu_b[draw_indices], nrow = length(draw_indices),
             ncol = data_info$n, byrow = FALSE)
  }

  exact <- isTRUE(.posterior_contract_native(fit))
  if (exact) {
    k_by_mky <- .posterior_predict_exact_k_by_market_year(fit, plan)
    if (is.null(k_by_mky)) {
      stop("Exact role prediction requires a draw plan with realized components.",
           call. = FALSE)
    }
  } else {
    k_by_mky <- NULL
  }
  if (exact && data_info$likelihood_mode != 2L) {
    margin_prediction <- .posterior_predict_exact_margin_mean(
      fit, draws, draw_indices, k_by_mky
    )
    margin_mean <- log(margin_prediction + 1e-9)
    margin_mean_source <- "selected exact-k component FOC"
  } else if (data_info$likelihood_mode == 2L) {
    margin_mean <- matrix(0, nrow = length(draw_indices), ncol = data_info$n)
    margin_mean_source <- "omitted by share-only likelihood"
  } else {
    pred <- .posterior_quality_matrix(draws$pred_marginNormInv, data_info$n,
                                      "pred_marginNormInv")
    if (nrow(pred) != n_total || any(!is.finite(pred)) || any(pred <= 0)) {
      stop("Posterior structural margin means are malformed.", call. = FALSE)
    }
    margin_mean <- log(pred[draw_indices, , drop = FALSE] + 1e-9)
    margin_mean_source <- "Stan pred_marginNormInv"
  }
  if (nrow(margin_mean) != length(draw_indices) ||
      ncol(margin_mean) != data_info$n || any(!is.finite(margin_mean))) {
    stop("Posterior margin means are malformed.", call. = FALSE)
  }
  rownames(share_mean) <- rownames(margin_mean) <- as.character(draw_indices)
  list(draw_indices = draw_indices, share_mean = share_mean,
       margin_mean = margin_mean, margin_mean_source = margin_mean_source,
       exact_k = k_by_mky)
}

.posterior_predict_scales <- function(fit, draws, draw_indices, data_info) {
  n_total <- .posterior_draw_count_from_draws(draws, fit)
  sigma_share <- .posterior_quality_vector(draws$sigma_share, n_total,
                                           "sigma_share")
  sigma_margin <- .posterior_quality_vector(draws$sigma_margin, n_total,
                                            "sigma_margin")
  if (any(!is.finite(sigma_share)) || any(sigma_share <= 0) ||
      any(!is.finite(sigma_margin)) || any(sigma_margin <= 0)) {
    stop("Posterior residual scales must be finite and positive.", call. = FALSE)
  }
  market_scales <- isTRUE(fit@stan_data$use_market_resid_scale == 1L)
  if (market_scales) {
    share_market <- .posterior_quality_matrix(
      draws$sigma_share_market, data_info$n_market, "sigma_share_market"
    )
    margin_market <- .posterior_quality_matrix(
      draws$sigma_margin_market, data_info$n_market, "sigma_margin_market"
    )
    if (nrow(share_market) != n_total || nrow(margin_market) != n_total ||
        any(!is.finite(share_market)) || any(share_market <= 0) ||
        any(!is.finite(margin_market)) || any(margin_market <= 0)) {
      stop("Posterior market-specific residual scales are malformed.",
           call. = FALSE)
    }
    share_scale <- share_market[draw_indices, data_info$market, drop = FALSE]
    margin_scale <- margin_market[draw_indices, data_info$market, drop = FALSE]
  } else {
    share_scale <- matrix(sigma_share[draw_indices], nrow = length(draw_indices),
                          ncol = data_info$n, byrow = FALSE)
    margin_scale <- matrix(sigma_margin[draw_indices], nrow = length(draw_indices),
                           ncol = data_info$n, byrow = FALSE)
  }
  list(share = share_scale, margin = margin_scale,
       sigma_share = sigma_share[draw_indices],
       sigma_margin = sigma_margin[draw_indices],
       market_specific = market_scales)
}

.posterior_predict_residual_draws <- function(share_mean, margin_mean,
                                              share_scale, margin_scale,
                                              share_family, margin_family,
                                              nu, rho, rho_active,
                                              likelihood_mode, margin_mask,
                                              seed = NULL,
                                              share_conditioning = c("replicate", "observed"),
                                              observed_share = NULL) {
  share_conditioning <- match.arg(share_conditioning)
  if (share_conditioning == "observed") {
    if (likelihood_mode != 1L) {
      stop("'share_conditioning = \"observed\"' requires the joint likelihood (likelihood_mode = 1).",
           call. = FALSE)
    }
    if (!is.numeric(observed_share) || length(observed_share) != ncol(share_mean) ||
        any(!is.finite(observed_share))) {
      stop("Conditional prediction requires finite observed log shares for every observation.",
           call. = FALSE)
    }
  }
  n_draw <- nrow(share_mean)
  n_obs <- ncol(share_mean)
  share_rep <- matrix(NA_real_, nrow = n_draw, ncol = n_obs)
  margin_rep <- matrix(NA_real_, nrow = n_draw, ncol = n_obs)
  calculate <- function() {
    for (r in seq_len(n_draw)) {
      nu_r <- if (length(nu) == 1L) nu[[1L]] else nu[[r]]
      rho_r <- if (length(rho) == 1L) rho[[1L]] else rho[[r]]
      if (likelihood_mode != 3L) {
        share_error <- if (share_conditioning == "observed") {
          observed_share - share_mean[r, ]
        } else if (identical(share_family, "gaussian")) {
          stats::rnorm(n_obs, sd = share_scale[r, ])
        } else {
          stats::rt(n_obs, df = nu_r) * share_scale[r, ] *
            sqrt((nu_r - 2) / nu_r)
        }
        share_rep[r, ] <- if (share_conditioning == "observed") {
          observed_share
        } else {
          share_mean[r, ] + share_error
        }
      }
      margin_predict_mask <- if (likelihood_mode != 2L) {
        rep.int(TRUE, n_obs)
      } else {
        rep.int(FALSE, n_obs)
      }
      if (any(margin_predict_mask)) {
        idx <- which(margin_predict_mask)
        if (isTRUE(rho_active)) {
          err_share <- share_rep[r, idx] - share_mean[r, idx]
          if (identical(margin_family, "gaussian")) {
            cond_mean <- margin_mean[r, idx] + rho_r *
              (margin_scale[r, idx] / share_scale[r, idx]) * err_share
            cond_scale <- margin_scale[r, idx] * sqrt(1 - rho_r^2)
            margin_rep[r, idx] <- stats::rnorm(
              length(idx), mean = cond_mean, sd = cond_scale
            )
          } else {
            z_share <- err_share / share_scale[r, idx]
            cond_mean <- margin_mean[r, idx] + rho_r *
              (margin_scale[r, idx] / share_scale[r, idx]) * err_share
            cond_scale <- margin_scale[r, idx] * sqrt(1 - rho_r^2) *
              sqrt((nu_r - 2 + z_share^2) / (nu_r + 1))
            margin_rep[r, idx] <- cond_mean + stats::rt(
              length(idx), df = nu_r + 1
            ) * cond_scale
          }
        } else if (identical(margin_family, "gaussian")) {
          margin_rep[r, idx] <- stats::rnorm(
            length(idx), mean = margin_mean[r, idx], sd = margin_scale[r, idx]
          )
        } else {
          margin_rep[r, idx] <- margin_mean[r, idx] + stats::rt(
            length(idx), df = nu_r
          ) * margin_scale[r, idx] * sqrt((nu_r - 2) / nu_r)
        }
      }
    }
    list(log_shareIn = share_rep, log_marginNormInv = margin_rep)
  }
  .posterior_quality_with_seed(seed, calculate)
}

.posterior_predict_from_draws <- function(fit, draws, plan, seed = NULL,
                                          share_conditioning = c("replicate", "observed")) {
  share_conditioning <- match.arg(share_conditioning)
  if (!is.list(draws)) stop("Posterior extraction did not return a draw list.",
                            call. = FALSE)
  data_info <- .posterior_predict_mask_and_mode(fit)
  if (share_conditioning == "observed" && data_info$likelihood_mode != 1L) {
    stop("'share_conditioning = \"observed\"' requires the joint likelihood (likelihood_mode = 1).",
         call. = FALSE)
  }
  share_family <- .posterior_predict_code_name(
    fit, "share_likelihood", c("gaussian", "student_t"),
    c(`1` = "student_t", `2` = "gaussian")
  )
  margin_family <- .posterior_predict_code_name(
    fit, "margin_likelihood", c("gaussian", "student_t"),
    c(`1` = "student_t", `2` = "gaussian")
  )
  means <- .posterior_predict_mean_state(fit, draws, plan, data_info)
  scales <- .posterior_predict_scales(
    fit, draws, means$draw_indices, data_info
  )
  n_total <- .posterior_draw_count_from_draws(draws, fit)
  need_t <- data_info$likelihood_mode != 3L && share_family == "student_t" ||
    data_info$likelihood_mode != 2L && margin_family == "student_t"
  nu <- if (need_t) {
    if (!is.null(draws$nu)) {
      .posterior_quality_vector(draws$nu, n_total, "nu")[means$draw_indices]
    } else if (length(fit@stan_data$resid_df_fixed) == 1L &&
               is.finite(fit@stan_data$resid_df_fixed)) {
      rep(as.numeric(fit@stan_data$resid_df_fixed), length(means$draw_indices))
    } else {
      stop("Student-t posterior prediction requires the fitted degrees of freedom.",
           call. = FALSE)
    }
  } else 2.5
  if (need_t && (any(!is.finite(nu)) || any(nu <= 2))) {
    stop("Student-t degrees-of-freedom draws must be finite and exceed 2.",
         call. = FALSE)
  }
  rho_active <- .posterior_contract_rho_active(fit)
  rho <- if (rho_active) {
    value <- .posterior_quality_vector(draws$rho, n_total, "rho")[means$draw_indices]
    if (any(!is.finite(value)) || any(abs(value) >= 1)) {
      stop("Posterior rho draws must be finite and lie in (-1, 1).",
           call. = FALSE)
    }
    value
  } else rep(0, length(means$draw_indices))
  replicate_seed <- if (is.null(seed)) NULL else
    as.integer((as.double(seed) + 3) %% .Machine$integer.max)
  replicated <- .posterior_predict_residual_draws(
    means$share_mean, means$margin_mean,
    scales$share, scales$margin,
    share_family, margin_family, nu, rho, rho_active,
    data_info$likelihood_mode, data_info$margin_mask, seed = replicate_seed,
    share_conditioning = share_conditioning,
    observed_share = if (share_conditioning == "observed") {
      data_info$observed_share
    } else NULL
  )
  rownames(replicated$log_shareIn) <-
    rownames(replicated$log_marginNormInv) <- as.character(means$draw_indices)
  colnames(replicated$log_shareIn) <- colnames(replicated$log_marginNormInv) <-
    as.character(seq_len(data_info$n))
  mky <- as.integer(fit@stan_data$market_year_idx %||%
    rep.int(1L, data_info$n))
  if (length(mky) != data_info$n || anyNA(mky) || any(mky < 1L)) {
    stop("Stan data contain a malformed market-year map.", call. = FALSE)
  }
  share_normalized <- matrix(NA_real_, nrow = nrow(replicated$log_shareIn),
                             ncol = ncol(replicated$log_shareIn),
                             dimnames = dimnames(replicated$log_shareIn))
  if (data_info$likelihood_mode != 3L) {
    for (rows in split(seq_len(data_info$n), mky)) {
      z <- replicated$log_shareIn[, rows, drop = FALSE]
      z <- sweep(z, 1L, apply(z, 1L, max), "-")
      p <- exp(z)
      share_normalized[, rows] <- p / rowSums(p)
    }
  }
  observation_ids <- .posterior_predict_observation_ids(fit, data_info, mky)
  long_index <- rep(seq_len(data_info$n), each = length(means$draw_indices))
  long_draw <- rep(means$draw_indices, times = data_info$n)
  long <- cbind(
    data.frame(draw_index = long_draw, observation_ids[long_index, , drop = FALSE],
               stringsAsFactors = FALSE),
    data.frame(
      log_shareIn_rep = as.vector(replicated$log_shareIn),
      share_rep_normalized = as.vector(share_normalized),
      log_marginNormInv_rep = as.vector(replicated$log_marginNormInv),
      margin_observed = data_info$margin_mask[long_index],
      exact_k = if (is.null(means$exact_k)) NA_integer_ else
        as.integer(means$exact_k[cbind(match(long_draw, means$draw_indices),
                                       mky[long_index])]),
      state_plan_fingerprint = rep(
        as.character(plan$fit_signature$latent_state_fingerprint %||% NA_character_),
        length(long_index)
      ),
      stringsAsFactors = FALSE
    )
  )
  metadata <- list(
      source = "posterior_predictive_from_fitted_Stan_likelihood",
      demand_model = as.character(fit@demand_model),
      supply_model = as.integer(fit@supply_model),
      share_likelihood = share_family,
      margin_likelihood = margin_family,
      likelihood_mode = data_info$likelihood_mode,
      rho_active = rho_active,
      market_specific_scales = scales$market_specific,
      quality_parameterization = .lq_quality_mode(fit),
      quality_scope = .lq_quality_scope(fit),
      market_intercept_parameterization =
        .market_intercept_parameterization_from_fit(fit),
      conditioning = if (.lq_quality_mode(fit) == "collapsed") {
        if (rho_active) "observed log shares and matched Gaussian margins" else
          "observed log shares"
      } else "posterior parameter draw",
      margin_mean_source = means$margin_mean_source,
      exact_k = !is.null(means$exact_k),
      exact_k_by_market_year = means$exact_k,
      state_plan_fingerprint = plan$fit_signature$latent_state_fingerprint,
      posterior_draws_available = as.integer(n_total),
      posterior_draws_used = as.integer(length(means$draw_indices)),
      margin_observed = data_info$margin_mask,
      margin_replicated = if (data_info$likelihood_mode != 2L) {
        rep.int(TRUE, data_info$n)
      } else {
        rep.int(FALSE, data_info$n)
      },
      share_rep_normalized_definition = "softmax of replicated log_shareIn within market-year; derived view, not a Stan likelihood variable",
      seed = if (is.null(seed)) NULL else as.integer(seed)
    )
  if (share_conditioning == "observed") {
    metadata$share_conditioning <- "observed"
    metadata$share_replicated <- FALSE
    metadata$share_rep_normalized_definition <-
      "softmax of fixed observed log_shareIn within market-year; derived view, not a Stan likelihood variable"
  }
  structure(list(
    observed = list(
      log_shareIn = data_info$observed_share,
      log_marginNormInv = data_info$observed_margin,
      margin_observed = data_info$margin_mask
    ),
    replicated = c(replicated, list(share_rep_normalized = share_normalized)),
    observation_ids = observation_ids,
    draws = long,
    draw_indices = means$draw_indices,
    draw_plan = plan,
    metadata = metadata
  ), class = c("posterior_predictive", "list"))
}
