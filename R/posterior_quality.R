## Conditional quality reconstruction for collapsed Gaussian fits.
##
## In collapsed mode Stan integrates the centered group-quality vector out of
## the share likelihood.  This helper reconstructs its conditional
## posterior after sampling the common parameters.  It intentionally works
## from the generated share_mean quantity rather than rebuilding the demand
## equation, so that all fitted intercepts, pooling, and conduct terms are
## inherited from the fit.

.posterior_quality_mode_value <- function(x, name, choices) {
  if (is.null(x)) return(NULL)
  if (length(x) != 1L || is.na(x)) {
    stop("Collapsed quality fit has malformed '", name, "' metadata.",
         call. = FALSE)
  }
  x <- as.character(x)
  if (!x %in% choices) {
    stop("posterior_quality() does not support '", x, "' for ", name,
         ".", call. = FALSE)
  }
  x
}

.posterior_quality_mode_code <- function(x, name, values) {
  if (is.null(x)) return(NULL)
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) ||
      x > .Machine$integer.max || x < -.Machine$integer.max ||
      x != floor(x) || !as.character(x) %in% names(values)) {
    stop("Collapsed quality fit has malformed Stan mode field '", name, "'.",
         call. = FALSE)
  }
  unname(values[[as.character(as.integer(x))]])
}

.posterior_quality_resolve_modes <- function(fit, stan_data) {
  metadata <- fit@metadata

  quality_attr <- attr(stan_data, "quality_parameterization", exact = TRUE)
  share_attr <- attr(stan_data, "share_likelihood", exact = TRUE)
  quality <- .posterior_quality_mode_value(
    metadata$quality_parameterization, "quality_parameterization",
    c("sampled", "collapsed")
  )
  if (is.null(quality)) {
    quality <- .posterior_quality_mode_value(
      quality_attr, "quality_parameterization", c("sampled", "collapsed")
    )
  }
  if (is.null(quality) && isTRUE(metadata$bayes_is_collapsed_quality)) {
    quality <- "collapsed"
  }
  quality_code <- .posterior_quality_mode_code(
    stan_data$quality_parameterization, "quality_parameterization",
    c(`1` = "sampled", `2` = "collapsed")
  )
  if (!is.null(quality) && !is.null(quality_code) &&
      !identical(quality, quality_code)) {
    stop("Collapsed quality fit has inconsistent quality parameterization metadata.",
         call. = FALSE)
  }
  quality <- quality %||% quality_code
  if (!identical(quality, "collapsed")) {
    stop("posterior_quality() requires quality_parameterization = 'collapsed'; ",
         "sampled quality fits are not supported.", call. = FALSE)
  }

  share <- .posterior_quality_mode_value(
    metadata$share_likelihood, "share_likelihood", c("student_t", "gaussian")
  )
  if (is.null(share)) {
    share <- .posterior_quality_mode_value(
      share_attr, "share_likelihood", c("student_t", "gaussian")
    )
  }
  share_code <- .posterior_quality_mode_code(
    stan_data$share_likelihood, "share_likelihood",
    c(`1` = "student_t", `2` = "gaussian")
  )
  if (!is.null(share) && !is.null(share_code) &&
      !identical(share, share_code)) {
    stop("Collapsed quality fit has inconsistent share likelihood metadata.",
         call. = FALSE)
  }
  ## The collapsed quality contract itself requires Gaussian shares.  Infer
  ## this only when no explicit field exists, which keeps small serialized
  ## fixtures useful while still rejecting an explicit Student-t fit.
  share <- share %||% share_code %||% "gaussian"
  if (!identical(share, "gaussian")) {
    stop("posterior_quality() requires share_likelihood = 'gaussian'; ",
         "non-Gaussian fits are not supported.", call. = FALSE)
  }

  likelihood_mode <- stan_data$likelihood_mode
  if (!is.null(likelihood_mode)) {
    if (length(likelihood_mode) != 1L || !is.numeric(likelihood_mode) ||
        !is.finite(likelihood_mode) || likelihood_mode > .Machine$integer.max ||
        likelihood_mode < -.Machine$integer.max || likelihood_mode != floor(likelihood_mode) ||
        !likelihood_mode %in% 1:3) {
      stop("Collapsed quality fit has malformed Stan field 'likelihood_mode'.",
           call. = FALSE)
    }
    if (as.integer(likelihood_mode) == 3L) {
      stop("posterior_quality() requires an active Gaussian share likelihood; ",
           "margin-only fits are not supported.", call. = FALSE)
    }
  }

  list(quality_parameterization = quality, share_likelihood = share,
       likelihood_mode = if (is.null(likelihood_mode)) NULL else as.integer(likelihood_mode))
}

.posterior_quality_matrix <- function(x, n_obs, name) {
  if (is.null(x) || !is.numeric(x)) {
    stop("Posterior draws field '", name, "' must be numeric.", call. = FALSE)
  }
  d <- dim(x)
  if (is.null(d)) {
    if (!length(x) || length(x) %% n_obs != 0L) {
      stop("Posterior draws field '", name,
           "' must contain complete draws by observations.", call. = FALSE)
    }
    ## A scalar generated vector is returned by rstan::extract as a vector;
    ## with one observation its length is therefore the number of draws.
    n_draw <- length(x) / n_obs
    return(matrix(as.numeric(x), nrow = n_draw, ncol = n_obs, byrow = TRUE))
  }
  if (length(d) != 2L || d[[2L]] != n_obs) {
    stop("Posterior draws field '", name,
         "' must have draws by observations dimensions.", call. = FALSE)
  }
  matrix(as.numeric(x), nrow = d[[1L]], ncol = d[[2L]])
}

.posterior_quality_vector <- function(x, n_draw, name) {
  if (is.null(x) || !is.numeric(x)) {
    stop("Posterior draws field '", name, "' must be numeric.", call. = FALSE)
  }
  d <- dim(x)
  if (!is.null(d)) {
    if (!(length(d) == 1L || (length(d) == 2L && any(d == 1L)))) {
      stop("Posterior draws field '", name,
           "' must be one value per draw.", call. = FALSE)
    }
    x <- as.numeric(x)
  }
  if (length(x) != n_draw) {
    stop("Posterior draws field '", name, "' has ", length(x),
         " values; expected ", n_draw, ".", call. = FALSE)
  }
  as.numeric(x)
}

.posterior_quality_draws <- function(fit, draws, conditioning = "auto") {
  stan_data <- fit@stan_data
  required <- c("group", "log_shareIn")
  missing <- setdiff(required, names(stan_data))
  if (length(missing)) {
    stop("Collapsed quality fit is missing Stan data field(s): ",
         paste(missing, collapse = ", "), ".", call. = FALSE)
  }
  y <- stan_data$log_shareIn
  group <- stan_data$group
  if (!is.numeric(y) || !length(y) || any(!is.finite(y))) {
    stop("Stan field 'log_shareIn' must be a nonempty finite numeric vector.",
         call. = FALSE)
  }
  if (!is.numeric(group) || length(group) != length(y) ||
      any(!is.finite(group)) || any(group > .Machine$integer.max) ||
      any(group < -.Machine$integer.max) || any(group != floor(group))) {
    stop("Stan field 'group' must be a finite integer vector matching log_shareIn.",
         call. = FALSE)
  }
  group <- as.integer(group)
  n_group <- stan_data$N_group %||% max(group)
  if (length(n_group) != 1L || !is.numeric(n_group) || !is.finite(n_group) ||
      n_group > .Machine$integer.max || n_group < 1 || n_group != floor(n_group)) {
    stop("Stan field 'N_group' must be a positive integer.", call. = FALSE)
  }
  n_group <- as.integer(n_group)
  if (any(group < 1L | group > n_group)) {
    stop("Stan field 'group' contains an index outside 1:N_group.", call. = FALSE)
  }
  n_obs <- length(y)
  if (!is.null(stan_data$N) && (length(stan_data$N) != 1L ||
      !is.numeric(stan_data$N) || !is.finite(stan_data$N) ||
      stan_data$N > .Machine$integer.max || stan_data$N < 1 ||
      stan_data$N != floor(stan_data$N) || as.integer(stan_data$N) != n_obs)) {
    stop("Stan field 'N' does not match log_shareIn.", call. = FALSE)
  }

  labels <- attr(stan_data, "group_levels", exact = TRUE)
  if (is.null(labels)) labels <- fit@metadata$group_levels
  if (is.null(labels) && length(fit@data$group) == n_obs) {
    ## prepare_stan_data uses factor()'s sorted levels for its integer map.
    original <- as.character(fit@data$group)
    labels <- vapply(seq_len(n_group), function(g) {
      vals <- original[group == g]
      if (length(vals)) vals[[1L]] else paste0("group_", g)
    }, character(1))
  }
  if (is.null(labels)) labels <- paste0("group_", seq_len(n_group))
  labels <- as.character(labels)
  if (length(labels) != n_group || anyNA(labels) || any(!nzchar(labels)) ||
      anyDuplicated(labels)) {
    stop("Collapsed quality fit has malformed or non-unique group labels.",
         call. = FALSE)
  }
  if (is.null(draws)) {
    draws <- .posterior_extract(fit, "posterior_quality", conditioning)
  }
  if (!is.list(draws)) {
    stop("'draws' must be a list containing share_mean, sigma_b_strat, and sigma_share.",
         call. = FALSE)
  }
  needed_draws <- c("share_mean", "sigma_b_strat", "sigma_share")
  missing_draws <- setdiff(needed_draws, names(draws))
  if (length(missing_draws)) {
    stop("Posterior draws are missing field(s): ",
         paste(missing_draws, collapse = ", "), ".", call. = FALSE)
  }
  share_mean <- .posterior_quality_matrix(draws$share_mean, n_obs, "share_mean")
  n_draw <- nrow(share_mean)
  if (!n_draw) stop("Posterior draws contain zero draws.", call. = FALSE)
  metadata_draws <- fit@metadata$posterior_draws
  if (!is.null(metadata_draws) && length(metadata_draws) == 1L &&
      is.numeric(metadata_draws) && is.finite(metadata_draws) &&
      metadata_draws <= .Machine$integer.max && metadata_draws >= 1 &&
      metadata_draws == floor(metadata_draws) && as.integer(metadata_draws) != n_draw) {
    stop("Posterior draws have ", n_draw,
         " draws; the fitted object records ", as.integer(metadata_draws), ".",
         call. = FALSE)
  }
  tau <- .posterior_quality_vector(draws$sigma_b_strat, n_draw, "sigma_b_strat")
  sigma <- .posterior_quality_vector(draws$sigma_share, n_draw, "sigma_share")
  if (any(!is.finite(share_mean)) || any(!is.finite(tau)) || any(!is.finite(sigma))) {
    stop("Posterior draws contain non-finite values.", call. = FALSE)
  }
  if (any(tau < 0) || any(sigma <= 0)) {
    stop("Posterior draws require sigma_b_strat >= 0 and sigma_share > 0.",
         call. = FALSE)
  }
  if (any(sigma > sqrt(.Machine$double.xmax)) ||
      any(tau > sqrt(.Machine$double.xmax))) {
    stop("Posterior scale draws are too large for stable Gaussian calculations.",
         call. = FALSE)
  }
  list(y = as.numeric(y), group = group, counts = tabulate(group, nbins = n_group),
       n_group = n_group, labels = labels,
       share_mean = share_mean, tau = tau, sigma = sigma, n_draw = n_draw,
       modes = .posterior_quality_resolve_modes(fit, stan_data))
}

.posterior_quality_with_seed <- function(seed, fun) {
  if (is.null(seed)) return(fun())
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  fun()
}

## The triangular Gaussian-rho transformation supplies extra observations
## of the same group quality. Only sufficient statistics are needed.
.posterior_quality_joint_draws <- function(fit, parsed, draws, conditioning) {
  data <- fit@stan_data
  active <- isTRUE(data$use_rho == 1L) &&
    isTRUE(data$likelihood_mode == 1L) && length(data$margin_idx) > 0L
  if (conditioning == "share" || !active) return(NULL)
  if (isTRUE(data$N_partition_component > 1L) ||
      identical(fit@metadata$partition_prior, "ordered_k") ||
      isTRUE(data$partition_prior_code == 1L)) {
    stop("Joint posterior_quality() recovery with rho and exact ordered-k ",
         "requires component-specific margin means; use conditioning = 'share' ",
         "for the explicitly share-only conditional distribution.", call. = FALSE)
  }
  needed <- c("rho", "sigma_margin", "pred_marginNormInv")
  ## The public caller extracts the complete quality/joint union once.  Keep
  ## this fallback for internal callers that explicitly invoke the helper.
  if (is.null(draws)) draws <- .posterior_extract(fit, "posterior_quality", conditioning)
  missing <- setdiff(needed, names(draws))
  if (length(missing)) stop("Joint quality recovery is missing draw field(s): ",
    paste(missing, collapse = ", "), ".", call. = FALSE)
  rho <- .posterior_quality_vector(draws$rho, parsed$n_draw, "rho")
  sigma <- .posterior_quality_vector(draws$sigma_margin, parsed$n_draw, "sigma_margin")
  prediction <- .posterior_quality_matrix(draws$pred_marginNormInv, length(parsed$y),
    "pred_marginNormInv")
  idx <- data$margin_idx
  if (!is.numeric(idx) || any(!is.finite(idx)) || any(idx != floor(idx)) ||
      any(idx < 1L | idx > length(parsed$y)) || anyDuplicated(idx))
    stop("Malformed share-margin matching in margin_idx.", call. = FALSE)
  idx <- as.integer(idx)
  y <- data$log_marginNormInv
  if (!is.numeric(y) || length(y) != length(parsed$y) || any(!is.finite(y[idx])))
    stop("Joint quality recovery requires matched finite log margins.", call. = FALSE)
  if (nrow(prediction) != parsed$n_draw || any(!is.finite(prediction[, idx])) ||
      any(prediction[, idx] <= 0) || any(!is.finite(rho)) || any(abs(rho) >= 1) ||
      any(!is.finite(sigma)) || any(sigma <= 0) ||
      any(sigma > sqrt(.Machine$double.xmax)))
    stop("Unusable Gaussian joint quality draw values.", call. = FALSE)
  list(rho = rho, sigma = sigma, prediction = prediction,
       idx = idx, y = as.numeric(y))
}

#' Reconstruct conditional group quality after a collapsed Gaussian fit
#'
#' @param fit A \code{\linkS4class{BayesLogit}} fit using Gaussian shares and
#'   \code{quality_parameterization = "collapsed"}.
#' @param n_draws Optional positive integer number of posterior common-
#'   parameter draws to use. If smaller than the available count, draws are
#'   selected without replacement. The default uses all available draws.
#' @param seed Optional integer seed controlling posterior draw selection and
#'   conditional quality simulation. The prior random-number state is restored
#'   after the function returns.
#' @param draws Optional list containing numeric \code{share_mean} (draws by
#'   observations), \code{sigma_b_strat}, and \code{sigma_share}. When omitted,
#'   these are extracted from the fit's Stan object. Joint recovery also needs
#'   \code{rho}, \code{sigma_margin}, and \code{pred_marginNormInv}.
#' @param conditioning \code{"auto"} (default) conditions on shares and
#'   matched margins when Gaussian rho is active, and on shares otherwise.
#'   \code{"share"} explicitly returns the share-only conditional distribution.
#'
#' @details
#' The quality vector is the group-level quality in the centered global prior
#' used by the collapsed Stan model. For each common-parameter draw this
#' function conditions that prior on the observed \code{log_shareIn} and,
#' with active Gaussian rho, the matched margins. Full joint recovery for
#' exact ordered-k rho fits is currently unsupported because the stored margin
#' prediction averages distinct component means; explicit share-only recovery
#' is available. The returned quality draw is generated
#' after the conditioning step, so HMC does not need to sample quality
#' coordinates. Computation is linear in the number of observations and groups
#' and never constructs a dense group covariance matrix.
#'
#' @return A \code{posterior_quality} list. \code{$draws} is a draws by groups
#' numeric matrix of conditional quality draws; \code{$conditional_means} and
#' \code{$conditional_sds} have the corresponding per-draw conditional moments.
#' Column names and \code{$group_labels} give stable group labels, while
#' \code{$draw_indices} records the original posterior draw number for each
#' row. \code{$draws_long} is a tidy version with columns \code{draw_index},
#' \code{group}, \code{quality}, \code{conditional_mean}, and
#' \code{conditional_sd}. \code{$summary} gives mixture posterior means and
#' SDs by group (including both common-parameter and conditional quality
#' uncertainty). \code{$metadata} records the selected draw count and modes.
#'
#' @export
posterior_quality <- function(fit, n_draws = NULL, seed = NULL, draws = NULL,
                              conditioning = c("auto", "share")) {
  conditioning <- match.arg(conditioning)
  if (!methods::is(fit, "BayesLogit")) {
    stop("'fit' must be a BayesLogit object.", call. = FALSE)
  }
  if (!is.null(n_draws) &&
      (length(n_draws) != 1L || !is.numeric(n_draws) || !is.finite(n_draws) ||
       n_draws > .Machine$integer.max || n_draws < 1 || n_draws != floor(n_draws))) {
    stop("'n_draws' must be a positive integer.", call. = FALSE)
  }
  if (!is.null(seed) &&
      (length(seed) != 1L || !is.numeric(seed) || !is.finite(seed) ||
       seed != floor(seed) || seed < 0 || seed > .Machine$integer.max)) {
    stop("'seed' must be a nonnegative integer.", call. = FALSE)
  }
  if (is.null(draws)) draws <- .posterior_extract(fit, "posterior_quality", conditioning)
  parsed <- .posterior_quality_draws(fit, draws, conditioning)
  joint <- .posterior_quality_joint_draws(fit, parsed, draws, conditioning)
  n_available <- parsed$n_draw
  n_use <- if (is.null(n_draws)) n_available else as.integer(n_draws)
  if (n_use > n_available) {
    stop("'n_draws' requests ", n_use, " draws but only ", n_available,
         " are available.", call. = FALSE)
  }
  seed <- if (is.null(seed)) NULL else as.integer(seed)

  calculate <- function() {
    draw_indices <- if (n_use == n_available) {
      seq_len(n_available)
    } else {
      sort(sample.int(n_available, n_use))
    }
    quality <- matrix(0, nrow = n_use, ncol = parsed$n_group,
                      dimnames = list(as.character(draw_indices), parsed$labels))
    means <- quality
    sds <- quality
    for (row in seq_len(n_use)) {
      draw <- draw_indices[[row]]
      sigma_sq <- parsed$sigma[[draw]]^2
      if (!is.finite(sigma_sq) || sigma_sq <= 0) {
        stop("Posterior draw ", draw,
             " has an unusable sigma_share squared.", call. = FALSE)
      }
      tau_sq <- parsed$tau[[draw]]^2
      if (!is.finite(tau_sq)) {
        stop("Posterior draw ", draw,
             " has an unusable sigma_b_strat squared.", call. = FALSE)
      }
      a <- parsed$counts / sigma_sq
      cstat <- numeric(parsed$n_group)
      csum <- rowsum(
        (parsed$y - parsed$share_mean[draw, ]) / sigma_sq,
        parsed$group, reorder = FALSE
      )
      cstat[as.integer(rownames(csum))] <- as.numeric(csum[, 1L])
      if (!is.null(joint)) {
        idx <- joint$idx
        slope <- joint$rho[[draw]] * joint$sigma[[draw]] / parsed$sigma[[draw]]
        variance <- joint$sigma[[draw]]^2 * (1 - joint$rho[[draw]]^2)
        if (!is.finite(variance) || variance <= 0 || !is.finite(slope))
          stop("Unusable Gaussian rho sufficient statistics.", call. = FALSE)
        residual <- joint$y[idx] - log(joint$prediction[draw, idx] + 1e-9) -
          slope * (parsed$y[idx] - parsed$share_mean[draw, idx])
        a <- a + tabulate(parsed$group[idx], nbins = parsed$n_group) * slope^2 / variance
        extra <- rowsum(-slope * residual / variance, parsed$group[idx], reorder = FALSE)
        which_group <- as.integer(rownames(extra))
        cstat[which_group] <- cstat[which_group] + as.numeric(extra[, 1L])
      }
      if (any(!is.finite(a)) || any(!is.finite(cstat))) {
        stop("Posterior draw ", draw,
             " produced non-finite Gaussian sufficient statistics.", call. = FALSE)
      }
      ## This form avoids tau^-2 overflow for tiny tau.  tau == 0 is the
      ## degenerate zero-quality limit and is handled without a 0/0 projection.
      d <- if (tau_sq == 0) rep(0, parsed$n_group) else {
        tau_sq / (1 + tau_sq * a)
      }
      if (any(!is.finite(d))) {
        stop("Posterior draw ", draw,
             " produced non-finite conditional quality variances.", call. = FALSE)
      }
      raw_mean <- d * cstat
      if (any(!is.finite(raw_mean))) {
        stop("Posterior draw ", draw,
             " produced non-finite conditional quality means.", call. = FALSE)
      }
      den <- sum(d)
      if (!is.finite(den)) {
        stop("Posterior draw ", draw,
             " produced a non-finite conditioning denominator.", call. = FALSE)
      }
      if (den > 0) {
        means[row, ] <- raw_mean - d * (sum(raw_mean) / den)
        variances <- pmax(0, d - d * d / den)
        sds[row, ] <- sqrt(variances)
        u <- stats::rnorm(parsed$n_group, mean = raw_mean, sd = sqrt(d))
        quality[row, ] <- u - d * (sum(u) / den)
        ## The projection is mathematically centered; remove the tiny final
        ## floating-point residual so the returned draw obeys the contract.
        quality[row, parsed$n_group] <- -sum(quality[row, -parsed$n_group, drop = TRUE])
      }
    }
    mixture_mean <- colMeans(means)
    mixture_variance <- colMeans(sds^2 + means^2) - mixture_mean^2
    mixture_variance <- pmax(0, mixture_variance)
    summary <- data.frame(
      group = parsed$labels,
      posterior_mean = as.numeric(mixture_mean),
      posterior_sd = sqrt(as.numeric(mixture_variance)),
      stringsAsFactors = FALSE
    )
    draws_long <- data.frame(
      draw_index = rep(as.integer(draw_indices), each = parsed$n_group),
      group = rep(parsed$labels, times = n_use),
      quality = as.vector(t(quality)),
      conditional_mean = as.vector(t(means)),
      conditional_sd = as.vector(t(sds)),
      stringsAsFactors = FALSE
    )
    structure(list(
      draws = quality,
      conditional_means = means,
      conditional_sds = sds,
      group_labels = parsed$labels,
      draw_indices = as.integer(draw_indices),
      draws_long = draws_long,
      summary = summary,
      metadata = list(
        available_draws = as.integer(n_available),
        selected_draws = as.integer(n_use),
        quality_parameterization = parsed$modes$quality_parameterization,
        share_likelihood = parsed$modes$share_likelihood,
        likelihood_mode = parsed$modes$likelihood_mode,
        conditioning = if (is.null(joint)) "observed log shares only" else
          "observed log shares and matched Gaussian margins",
        centered_constraint = "sum(group quality) = 0"
      )
    ), class = c("posterior_quality", "list"))
  }
  .posterior_quality_with_seed(seed, calculate)
}
