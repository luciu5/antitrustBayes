## Conditional quality reconstruction for collapsed Gaussian fits.
##
## In collapsed mode Stan integrates the centered group-quality vector out of
## the share likelihood.  This helper reconstructs its conditional
## posterior after sampling the common parameters.  It intentionally works
## from retained predictors rather than rebuilding the demand equation.
## Absolute-share fits retain their share_mean/intercept contract. Relative
## fits condition only on within-market-year contrasts of the retained
## relative_share_predictor and derive calibration levels after recovering q.

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
    if (as.integer(likelihood_mode) == 3L && !.posterior_demand_shock_model(fit)) {
      stop("posterior_quality() requires an active Gaussian share likelihood; ",
           "margin-only fits are not supported.", call. = FALSE)
    }
  }

  quality_scope <- .quality_scope_from_fit(fit, strict = FALSE)
  market_intercept <- .market_intercept_parameterization_from_fit(fit)
  if (identical(market_intercept, "collapsed") &&
      !identical(quality_scope, "market")) {
    stop("Collapsed market-intercept recovery requires market-scoped quality.",
         call. = FALSE)
  }
  list(quality_parameterization = quality, share_likelihood = share,
       likelihood_mode = if (is.null(likelihood_mode)) NULL else as.integer(likelihood_mode),
       quality_scope = quality_scope,
       market_intercept_parameterization = market_intercept)
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
    ## A scalar generated vector can be represented as a plain vector;
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

  incidence <- .quality_incidence_from_fit(fit, strict = FALSE)
  if (is.null(incidence)) stop("Collapsed quality fit has no quality incidence map.", call. = FALSE)
  scope <- as.character(incidence$quality_scope)
  n_quality <- as.integer(incidence$N_quality)
  if (length(n_quality) != 1L || is.na(n_quality) || n_quality < 1L)
    stop("Collapsed quality fit has malformed N_quality.", call. = FALSE)
  quality <- as.integer(incidence$obs_to_quality)
  if (length(quality) != n_obs || anyNA(quality) || any(quality < 1L | quality > n_quality))
    stop("Collapsed quality fit has malformed obs_to_quality.", call. = FALSE)
  blocks <- as.integer(incidence$quality_block_n)
  starts <- as.integer(incidence$quality_block_start)
  if (length(blocks) != incidence$N_quality_block || length(starts) != length(blocks) ||
      any(blocks < 1L) || any(starts < 1L) ||
      any(starts + blocks - 1L > n_quality) ||
      !identical(starts, cumsum(c(1L, utils::head(blocks, -1L)))))
    stop("Collapsed quality fit has malformed centering blocks.", call. = FALSE)
  lookup <- incidence$quality_lookup
  if (!is.data.frame(lookup) || nrow(lookup) != n_quality ||
      !all(c("quality_index", "group_label", "quality_unit") %in% names(lookup)))
    stop("Collapsed quality fit has malformed quality lookup metadata.", call. = FALSE)
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
  supplied_draws <- !is.null(draws)
  if (!supplied_draws) {
    draws <- .posterior_extract(fit, "posterior_quality", conditioning)
  }
  if (!is.list(draws)) {
    stop("'draws' must be a list containing a share-mean field, an economic quality scale, and sigma_share.",
         call. = FALSE)
  }
  modes <- .posterior_quality_resolve_modes(fit, stan_data)
  intercept_collapsed <- identical(
    modes$market_intercept_parameterization, "collapsed"
  )
  relative <- identical(.share_basis_from_fit(fit), "relative")
  share_mean_name <- if (relative) "relative_share_predictor" else if (intercept_collapsed)
    "share_mean_integrated" else "share_mean"
  share_mean_value <- draws[[share_mean_name]]
  needed_draws <- c("sigma_share")
  if (intercept_collapsed)
    needed_draws <- c(needed_draws, "mu_b_market_vec", "sigma_b_market_vec")
  missing_draws <- setdiff(needed_draws, names(draws))
  if (length(missing_draws)) {
    stop("Posterior draws are missing field(s): ",
         paste(missing_draws, collapse = ", "), ".", call. = FALSE)
  }
  if (is.null(share_mean_value)) missing_draws <- c(missing_draws, share_mean_name)
  if (length(missing_draws)) {
    stop("Posterior draws are missing field(s): ",
         paste(unique(missing_draws), collapse = ", "), ".")
  }
  share_mean <- .posterior_quality_matrix(share_mean_value, n_obs,
                                           share_mean_name)
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
  scale_mode <- fit@metadata$quality_scale_parameterization_effective %||%
    fit@metadata$quality_scale_parameterization %||%
    attr(stan_data, "quality_scale_parameterization", exact = TRUE) %||%
    if (isTRUE(stan_data$quality_scale_parameterization == 2L)) {
      "total_fraction"
    } else {
      "legacy"
    }
  if (!supplied_draws) {
    ## .posterior_extract() canonicalizes the preferred economic alias to this
    ## historical name, so it is already tau rather than a raw sampler field.
    tau_value <- draws$sigma_b_strat
    tau_name <- "sigma_b_strat"
  } else if ("sigma_b_strat_economic" %in% names(draws)) {
    tau_value <- draws$sigma_b_strat_economic
    tau_name <- "sigma_b_strat_economic"
  } else if (identical(scale_mode, "total_fraction")) {
    if (!"sigma_b_strat" %in% names(draws)) {
      stop("Posterior draws are missing the economic quality scale ",
           "'sigma_b_strat_economic'.", call. = FALSE)
    }
    total <- .posterior_quality_vector(draws$sigma_b_strat, n_draw,
                                       "sigma_b_strat")
    if ("log_sigma_share" %in% names(draws)) {
      fraction_logit <- .posterior_quality_vector(
        draws$log_sigma_share, n_draw, "log_sigma_share"
      )
      tau_value <- exp(log(total) +
        .5 * stats::plogis(-fraction_logit, log.p = TRUE))
    } else if ("share_variance_fraction" %in% names(draws)) {
      fraction <- .posterior_quality_vector(
        draws$share_variance_fraction, n_draw, "share_variance_fraction"
      )
      if (any(!is.finite(fraction)) || any(fraction < 0) || any(fraction > 1)) {
        stop("Posterior share_variance_fraction draws must lie in [0, 1].",
             call. = FALSE)
      }
      tau_value <- total * sqrt(1 - fraction)
    } else {
      stop("Under quality_scale_parameterization = 'total_fraction', supplied ",
           "draws must include 'sigma_b_strat_economic' or the raw pair ",
           "'sigma_b_strat' and 'log_sigma_share'.", call. = FALSE)
    }
    tau_name <- "reconstructed sigma_b_strat_economic"
  } else {
    if (!"sigma_b_strat" %in% names(draws)) {
      stop("Posterior draws are missing field(s): sigma_b_strat.",
           call. = FALSE)
    }
    tau_value <- draws$sigma_b_strat
    tau_name <- "sigma_b_strat"
  }
  tau <- .posterior_quality_vector(tau_value, n_draw, tau_name)
  sigma <- .posterior_quality_vector(draws$sigma_share, n_draw, "sigma_share")
  n_market <- as.integer(stan_data$N_market %||%
    if (!is.null(stan_data$market)) max(as.integer(stan_data$market)) else 1L)
  market_labels <- attr(stan_data, "market_levels", exact = TRUE) %||%
    fit@metadata$market_levels %||% if ("market" %in% names(fit@data)) {
      unique(as.character(fit@data$market))
    } else {
      as.character(seq_len(n_market))
    }
  market_labels <- as.character(market_labels)
  if (length(market_labels) != n_market || anyNA(market_labels) ||
      any(!nzchar(market_labels)) || anyDuplicated(market_labels)) {
    stop("Collapsed quality fit has malformed market labels.", call. = FALSE)
  }
  mu_b_market <- sigma_b_market <- sigma_share_market <- sigma_margin_market <- NULL
  if (intercept_collapsed) {
    mu_b_market <- .posterior_quality_vector(
      draws$mu_b_market_vec, n_draw, "mu_b_market_vec"
    )
    sigma_b_market <- .posterior_quality_vector(
      draws$sigma_b_market_vec, n_draw, "sigma_b_market_vec"
    )
    if (any(!is.finite(mu_b_market)) || any(!is.finite(sigma_b_market)) ||
        any(sigma_b_market <= 0)) {
      stop("Collapsed market-intercept hyperparameter draws are unusable.",
           call. = FALSE)
    }
    if (isTRUE(stan_data$use_market_resid_scale == 1L)) {
      sigma_share_market <- .posterior_quality_matrix(
        draws$sigma_share_market, n_market, "sigma_share_market"
      )
      if (nrow(sigma_share_market) != n_draw ||
          any(!is.finite(sigma_share_market)) || any(sigma_share_market <= 0)) {
        stop("Market-specific share residual scales are malformed.",
             call. = FALSE)
      }
      if (any(sigma_share_market > sqrt(.Machine$double.xmax)))
        stop("Market-specific share residual scales are too large for stable Gaussian calculations.",
             call. = FALSE)
    }
  }
  if (relative && isTRUE(stan_data$use_market_resid_scale == 1L)) {
    sigma_share_market <- .posterior_quality_matrix(
      draws$sigma_share_market, n_market, "sigma_share_market")
    if (nrow(sigma_share_market) != n_draw || any(!is.finite(sigma_share_market)) ||
        any(sigma_share_market <= 0))
      stop("Market-specific share residual scales are malformed.", call. = FALSE)
  }
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
  if (scope == "global") {
    if (n_quality != n_group || !identical(quality, group))
      stop("Global quality incidence disagrees with the economic group map.", call. = FALSE)
  }
  market <- if (!is.null(stan_data$market)) {
    as.integer(stan_data$market)
  } else {
    ## Legacy objects may carry only the market CSR ordering.  Rebuild the
    ## observation labels from each contiguous CSR segment; the inverse
    ## permutation alone is not a market identifier.
    out <- integer(n_obs)
    for (m in seq_along(incidence$market_obs_n)) {
      n_m <- as.integer(incidence$market_obs_n[[m]])
      if (!n_m) next
      first <- as.integer(incidence$market_obs_start[[m]])
      last <- first + n_m - 1L
      out[incidence$market_obs_order[first:last]] <- m
    }
    out
  }
  if (length(market) != n_obs || anyNA(market) ||
      any(market < 1L | market > n_market)) {
    stop("Collapsed quality fit has malformed market observation indices.",
         call. = FALSE)
  }
  margin_scale_active <- (intercept_collapsed || relative) &&
    isTRUE(stan_data$use_market_resid_scale == 1L) &&
    isTRUE(stan_data$use_rho == 1L) &&
    isTRUE(stan_data$likelihood_mode == 1L) &&
    length(stan_data$margin_idx) > 0L && identical(conditioning, "auto")
  if (margin_scale_active) {
    sigma_margin_market <- .posterior_quality_matrix(
      draws$sigma_margin_market, n_market, "sigma_margin_market"
    )
    if (nrow(sigma_margin_market) != n_draw ||
        any(!is.finite(sigma_margin_market)) || any(sigma_margin_market <= 0) ||
        any(sigma_margin_market > sqrt(.Machine$double.xmax))) {
      stop("Market-specific margin residual scales are malformed.",
           call. = FALSE)
    }
  }
  list(y = as.numeric(y), group = group, quality = quality,
       counts = tabulate(quality, nbins = n_quality),
       n_group = n_group, n_quality = n_quality, labels = labels,
       quality_labels = as.character(lookup$group_label),
       quality_units = as.character(lookup$quality_unit),
       quality_lookup = lookup, quality_scope = scope,
       blocks = blocks, block_starts = starts,
       incidence = incidence,
       share_mean = share_mean, tau = tau, sigma = sigma, n_draw = n_draw,
       n_market = n_market, market = market, market_labels = market_labels,
       mu_b_market = mu_b_market, sigma_b_market = sigma_b_market,
       sigma_share_market = sigma_share_market,
       sigma_margin_market = sigma_margin_market,
       relative = relative, modes = modes)
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

.posterior_quality_exact_k <- function(fit, parsed, draws, draw_indices,
                                       seed = NULL, k_cache = NULL) {
  d <- fit@stan_data
  required <- c("N_partition_unit", "partition_unit_market",
                "partition_component_k", "partition_component_unit",
                "partition_unit_count")
  if (!all(required %in% names(d)) || !is.null(k_cache) && !is.list(k_cache))
    return(NULL)
  n_unit <- as.integer(d$N_partition_unit)
  if (length(n_unit) != 1L || is.na(n_unit) || n_unit < 1L)
    stop("Exact ordered-k quality recovery has malformed N_partition_unit.", call. = FALSE)
  unit_market <- as.integer(d$partition_unit_market)
  component_k <- as.integer(d$partition_component_k)
  component_unit <- as.integer(d$partition_component_unit)
  unit_count <- as.integer(d$partition_unit_count)
  if (length(unit_market) != n_unit || anyNA(unit_market) || any(unit_market < 1L) ||
      length(unit_count) != n_unit || anyNA(unit_count) || any(unit_count < 1L) ||
      length(component_k) != sum(unit_count) ||
      length(component_k) != length(component_unit) ||
      anyNA(component_unit) || any(component_unit < 1L | component_unit > n_unit) ||
      anyNA(component_k) || any(component_k < 0L))
    stop("Exact ordered-k quality recovery has malformed partition maps.", call. = FALSE)
  n_market <- as.integer(d$N_market %||% max(unit_market))
  if (length(n_market) != 1L || is.na(n_market) || n_market < 1L ||
      any(unit_market > n_market) || n_unit != n_market ||
      !identical(sort(unique(unit_market)), seq_len(n_market)))
    stop("Exact ordered-k quality recovery requires one unit per market.", call. = FALSE)
  if (is.null(draws) || !is.list(draws) || is.null(draws$posterior_component_prob))
    stop("Exact ordered-k quality recovery requires posterior_component_prob draws.", call. = FALSE)
  prob <- .posterior_quality_matrix(draws$posterior_component_prob,
                                    length(component_k), "posterior_component_prob")
  if (nrow(prob) != parsed$n_draw || any(!is.finite(prob)) || any(prob < 0))
    stop("Posterior ordered-k probabilities are malformed.", call. = FALSE)
  unit_components <- lapply(seq_len(n_unit), function(u) which(component_unit == u))
  if (any(!lengths(unit_components))) stop("Exact ordered-k quality recovery has an empty partition unit.", call. = FALSE)
  if (any(vapply(seq_len(n_unit), function(u) {
    length(unit_components[[u]]) != unit_count[[u]]
  }, logical(1))))
    stop("Exact ordered-k quality recovery has inconsistent partition component counts.", call. = FALSE)
  market_unit <- match(seq_len(n_market), unit_market)
  if (anyNA(market_unit)) stop("Exact ordered-k quality recovery has an incomplete market map.", call. = FALSE)
  get_cached <- function(draw) {
    if (is.null(k_cache)) return(NULL)
    z <- k_cache[[as.character(draw)]]
    if (is.null(z) && length(k_cache) >= draw) z <- k_cache[[draw]]
    if (is.null(z)) return(NULL)
    z <- as.integer(z)
    if (length(z) != n_unit) stop("The ordered-k cache must contain one state per market partition unit.", call. = FALSE)
    z[market_unit]
  }
  calculate <- function() {
    out <- matrix(NA_integer_, nrow = length(draw_indices), ncol = n_market,
                  dimnames = list(as.character(draw_indices),
                                  as.character(attr(d, "market_levels", exact = TRUE) %||%
                                    seq_len(n_market))))
    for (rr in seq_along(draw_indices)) {
      draw <- as.integer(draw_indices[[rr]])
      cached <- get_cached(draw)
      if (!is.null(cached)) {
        if (length(cached) != n_market || anyNA(cached))
          stop("The ordered-k cache contains malformed market states.", call. = FALSE)
        if (any(vapply(seq_len(n_market), function(m) {
          cached[[m]] %in% component_k[unit_components[[market_unit[[m]]]]]
        }, logical(1)) == FALSE))
          stop("The ordered-k cache contains a state outside its partition support.", call. = FALSE)
        out[rr, ] <- cached
      } else {
        for (m in seq_len(n_market)) {
          cc <- unit_components[[market_unit[[m]]]]
          pp <- prob[draw, cc]
          if (!length(pp) || sum(pp) <= 0)
            stop("Ordered-k posterior probabilities have zero mass in a market.", call. = FALSE)
          pp <- pp / sum(pp)
          pick <- sample.int(length(cc), 1L, prob = pp)
          out[rr, m] <- component_k[cc[[pick]]]
        }
      }
    }
    out
  }
  .posterior_quality_with_seed(seed, calculate)
}

## Component-specific inverse-margin means used when a market-scoped quality
## fit has active rho and an exact ordered-k likelihood.  These formulas mirror
## the existing Stan Stackelberg kernels and use one firm statistic per firm-
## market-year, so multiproduct firms are never counted once per product.
.stack_pure_foc_margin_prediction <- function(fit, parsed, draws,
                                              draw_indices, realized_k,
                                              realized_k_scope = c("market",
                                                                   "market_year")) {
  realized_k_scope <- match.arg(realized_k_scope)
  d <- fit@stan_data
  needed <- c("firm_mky_start", "firm_mky_n",
              "firm_mky_rank", "obs_to_firm_mky", "market_year_idx",
              "mky_to_market")
  if (identical(as.character(fit@demand_model), "ces"))
    needed <- c(needed, "firm_mky_revShare", "outside_price", "inside_price_index")
  else needed <- c(needed, "firm_mky_shareIn", "price")
  missing <- setdiff(needed, names(d))
  if (length(missing)) stop("Exact ordered-k quality recovery cannot compute component margin means; missing field(s): ",
    paste(missing, collapse = ", "), ".", call. = FALSE)
  n_market <- max(as.integer(d$market))
  n_mky <- max(as.integer(d$market_year_idx))
  n <- length(parsed$y)
  mky_market <- as.integer(d$mky_to_market)
  if (length(mky_market) != n_mky || anyNA(mky_market) ||
      any(mky_market < 1L | mky_market > n_market))
    stop("Exact ordered-k market-year map is malformed.", call. = FALSE)
  expected_k <- if (realized_k_scope == "market") n_market else n_mky
  if (!is.matrix(realized_k) || nrow(realized_k) != length(draw_indices) ||
      ncol(realized_k) != expected_k || anyNA(realized_k) ||
      any(realized_k < 0L) || any(realized_k != floor(realized_k))) {
    stop("Exact ordered-k states do not match the selected scope.",
         call. = FALSE)
  }
  k_market_year <- if (realized_k_scope == "market") {
    realized_k[, mky_market, drop = FALSE]
  } else {
    realized_k
  }
  amat <- .posterior_quality_matrix(draws$a_market, n_market, "a_market")
  smat <- .posterior_quality_matrix(draws$s0_prob, n_mky, "s0_prob")
  if (nrow(amat) != parsed$n_draw || nrow(smat) != parsed$n_draw)
    stop("Exact ordered-k structural draws have inconsistent draw counts.", call. = FALSE)
  if (any(!is.finite(amat)) || any(amat <= 0) || any(!is.finite(smat)) ||
      any(smat <= 0) || any(smat >= 1))
    stop("Exact ordered-k structural draws contain unusable demand parameters.", call. = FALSE)
  structure_name <- as.character(fit@metadata$structure %||%
    attr(d, "structure", exact = TRUE) %||%
    if (isTRUE(d$structure == 2L)) "stackelberg" else "fringe")
  supply <- as.integer(fit@supply_model[[1L]])
  firm_start <- as.integer(d$firm_mky_start)
  firm_n <- as.integer(d$firm_mky_n)
  firm_rank <- as.integer(d$firm_mky_rank)
  obs_firm <- as.integer(d$obs_to_firm_mky)
  mky_market <- as.integer(d$mky_to_market)
  pred <- matrix(NA_real_, nrow = length(draw_indices), ncol = n)
  for (rr in seq_along(draw_indices)) {
    draw <- as.integer(draw_indices[[rr]])
    k_by_market_year <- k_market_year[rr, ]
    response <- numeric(n_mky)
    if (identical(structure_name, "stackelberg")) for (mky in seq_len(n_mky)) {
      m <- mky_market[[mky]]; k <- k_by_market_year[[mky]]
      st <- firm_start[[mky]]; nf <- firm_n[[mky]]
      ff <- seq.int(st, length.out = nf)
      if (identical(as.character(fit@demand_model), "ces")) {
        R0 <- d$outside_price[[mky]] * smat[draw, mky] /
          (d$outside_price[[mky]] * smat[draw, mky] +
             d$inside_price_index[[mky]] * (1 - smat[draw, mky]))
        gam <- 1 + amat[draw, m] / as.numeric(d$log_price_sd %||% fit@price_sd)
        h <- gam - 1; C <- 0
        for (fj in ff[firm_rank[ff] > k]) {
          R <- d$firm_mky_revShare[[fj]] * (1 - R0)
          if (supply == 1L) {
            t <- gam - h * R
            af <- h^2 * R / (t * (t - 1) + h^2 * R)
            C <- C + R * af
          } else {
            H <- 1 + h * R0
            af <- h * (H - 2 * R + gam * h * R0 * R / H) /
              (H * (H + (gam - 2) * R))
            C <- C + R * af
          }
        }
      } else {
        S <- d$firm_mky_shareIn[ff] * (1 - smat[draw, mky])
        if (supply == 1L) C <- sum(S[firm_rank[ff] > k]^2 /
          (1 - S[firm_rank[ff] > k] + S[firm_rank[ff] > k]^2)) else C <- sum(S[firm_rank[ff] > k])
      }
      response[[mky]] <- if (supply == 1L) 1 / (1 - C) else C
    }
    for (i in seq_len(n)) {
      m <- as.integer(d$market[[i]]); mky <- as.integer(d$market_year_idx[[i]])
      k <- as.integer(k_by_market_year[[mky]]); me <- obs_firm[[i]]
      rank <- firm_rank[[me]]
      if (identical(as.character(fit@demand_model), "ces")) {
        s0 <- smat[draw, mky]
        R0 <- d$outside_price[[mky]] * s0 /
          (d$outside_price[[mky]] * s0 + d$inside_price_index[[mky]] * (1 - s0))
        gam <- 1 + amat[draw, m] / as.numeric(d$log_price_sd %||% fit@price_sd)
        R <- d$firm_mky_revShare[[me]] * (1 - R0); h <- gam - 1
        if (identical(structure_name, "fringe")) {
          mm <- if (rank > k) 1 / gam else if (supply == 1L) 1 / (gam - h * R) else 1 / gam + h / (gam * (1 + h * R0)) * R
        } else if (supply == 1L) {
          mm <- if (rank > k) 1 / (gam - h * R) else 1 / (gam - h * response[[mky]] * R)
        } else {
          mm <- if (rank > k) 1 / gam + h / (gam * (1 + h * R0)) * R else
            1 / gam + h / (gam * (1 + h * R0)) * R / (1 + response[[mky]])
        }
        pred[rr, i] <- 1 / mm
      } else {
        alpha <- amat[draw, m] / as.numeric(d$price_sd %||% fit@price_sd)
        S <- d$firm_mky_shareIn[[me]] * (1 - smat[draw, mky])
        h <- if (identical(structure_name, "fringe")) {
          if (rank > k) 1 else if (supply == 1L) 1 / (1 - S) else 1 + S / smat[draw, mky]
        } else if (supply == 1L) {
          if (rank > k) 1 / (1 - S) else 1 / (1 - response[[mky]] * S)
        } else if (rank > k) 1 + S / smat[draw, mky] else
          1 + S / (smat[draw, mky] + response[[mky]])
        pred[rr, i] <- d$price[[i]] * alpha / h
      }
    }
  }
  if (any(!is.finite(pred)) || any(pred <= 0))
    stop("Exact ordered-k FOC margin prediction returned non-finite values.", call. = FALSE)
  pred
}

## The triangular Gaussian-rho transformation supplies extra observations
## of the same group quality. Only sufficient statistics are needed.
.posterior_quality_joint_draws <- function(fit, parsed, draws, conditioning,
                                           draw_indices = seq_len(parsed$n_draw),
                                           seed = NULL, k_cache = NULL) {
  data <- fit@stan_data
  active <- isTRUE(data$use_rho == 1L) &&
    isTRUE(data$likelihood_mode == 1L) && length(data$margin_idx) > 0L
  if (conditioning == "share" || !active) return(NULL)
  exact <- identical(.posterior_contract_partition(fit), "ordered_k") ||
    isTRUE(data$partition_prior_code == 1L)
  realized_k <- NULL
  n_market <- as.integer(data$N_market %||%
    if (!is.null(data$market)) max(as.integer(data$market)) else 1L)
  exact_quality_ok <- exact && (identical(parsed$quality_scope, "market") ||
    (identical(parsed$quality_scope, "global") && n_market == 1L))
  if (exact_quality_ok) {
    realized_k <- .posterior_quality_exact_k(fit, parsed, draws, draw_indices,
                                             seed = seed, k_cache = k_cache)
    if (is.null(realized_k))
      stop("Exact ordered-k joint quality recovery requires component-specific margin means plus posterior component probabilities and partition maps.", call. = FALSE)
  } else if (exact) {
    stop("Joint posterior_quality() recovery with rho and exact ordered-k ",
         "requires component-specific margin means; use conditioning = 'share' ",
         "for the explicitly share-only conditional distribution.", call. = FALSE)
  }
  needed <- c("rho", "sigma_margin")
  if (is.null(realized_k)) needed <- c(needed, "pred_marginNormInv")
  ## The public caller extracts the complete quality/joint union once.  Keep
  ## this fallback for internal callers that explicitly invoke the helper.
  if (is.null(draws)) draws <- .posterior_extract(fit, "posterior_quality", conditioning)
  missing <- setdiff(needed, names(draws))
  if (length(missing)) stop("Joint quality recovery is missing draw field(s): ",
    paste(missing, collapse = ", "), ".", call. = FALSE)
  rho <- .posterior_quality_vector(draws$rho, parsed$n_draw, "rho")
  sigma <- .posterior_quality_vector(draws$sigma_margin, parsed$n_draw, "sigma_margin")
  sigma_by_market <- if (!is.null(parsed$sigma_margin_market)) {
    parsed$sigma_margin_market
  } else {
    matrix(rep(sigma, parsed$n_market), nrow = parsed$n_draw,
           ncol = parsed$n_market)
  }
  prediction <- if (is.null(realized_k)) {
    .posterior_quality_matrix(draws$pred_marginNormInv, length(parsed$y),
      "pred_marginNormInv")
  } else {
    .stack_pure_foc_margin_prediction(fit, parsed, draws,
      draw_indices, realized_k)
  }
  idx <- data$margin_idx
  if (!is.numeric(idx) || any(!is.finite(idx)) || any(idx != floor(idx)) ||
      any(idx < 1L | idx > length(parsed$y)) || anyDuplicated(idx))
    stop("Malformed share-margin matching in margin_idx.", call. = FALSE)
  idx <- as.integer(idx)
  y <- data$log_marginNormInv
  if (!is.numeric(y) || length(y) != length(parsed$y) || any(!is.finite(y[idx])))
    stop("Joint quality recovery requires matched finite log margins.", call. = FALSE)
  expected_rows <- if (is.null(realized_k)) parsed$n_draw else length(draw_indices)
  if (nrow(prediction) != expected_rows || any(!is.finite(prediction[, idx])) ||
      any(prediction[, idx] <= 0) || any(!is.finite(rho)) || any(abs(rho) >= 1) ||
      any(!is.finite(sigma)) || any(sigma <= 0) ||
      any(sigma > sqrt(.Machine$double.xmax)))
    stop("Unusable Gaussian joint quality draw values.", call. = FALSE)
  if (any(!is.finite(sigma_by_market)) || any(sigma_by_market <= 0) ||
      any(sigma_by_market > sqrt(.Machine$double.xmax)))
    stop("Unusable market-specific Gaussian margin scale draw values.",
         call. = FALSE)
  if (!identical(as.integer(draw_indices), seq_len(parsed$n_draw))) {
    rho <- rho[draw_indices]
    sigma <- sigma[draw_indices]
    sigma_by_market <- sigma_by_market[draw_indices, , drop = FALSE]
    if (is.null(realized_k)) prediction <- prediction[draw_indices, , drop = FALSE]
  }
  rho_by_market <- if (.posterior_demand_shock_model(fit) &&
      isTRUE(data$use_market_resid_scale == 1L)) {
    z <- .posterior_quality_matrix(draws$rho_market, n_market, "rho_market")
    if (nrow(z) != parsed$n_draw || any(!is.finite(z)) || any(abs(z) >= 1))
      stop("Unusable market-specific Gaussian correlation draws.", call. = FALSE)
    z[draw_indices, , drop = FALSE]
  } else matrix(rho, nrow = length(draw_indices), ncol = n_market)
  list(rho = rho, rho_by_market = rho_by_market,
       sigma = sigma, sigma_by_market = sigma_by_market,
       prediction = prediction,
       idx = idx, y = as.numeric(y), realized_k = realized_k)
}

## Recover one market's joint Gaussian posterior using only group sufficient
## statistics. First draw the centered intercept offset from its marginal
## posterior, then draw centered quality conditional on that same offset.
.posterior_quality_market_intercept_conditional <- function(
    z, variance, group, tau, mu_b, sigma_b) {
  if (!length(z) || length(z) != length(variance) ||
      length(z) != length(group) || any(!is.finite(z)) ||
      any(!is.finite(variance)) || any(variance <= 0) ||
      anyNA(group) || any(group < 1L) || any(group != floor(group)) ||
      length(tau) != 1L || !is.finite(tau) || tau < 0 ||
      length(mu_b) != 1L || !is.finite(mu_b) ||
      length(sigma_b) != 1L || !is.finite(sigma_b) || sigma_b <= 0) {
    stop("Collapsed market-intercept recovery received malformed Gaussian inputs.",
         call. = FALSE)
  }
  group <- as.integer(group)
  n_group <- max(group)
  ## A single quality group is fixed at zero by the within-market centering
  ## constraint.  Match the Stan kernel by removing its otherwise spurious
  ## random-effect variance before recovering the market intercept.
  if (n_group == 1L) tau <- 0
  a <- cstat <- numeric(n_group)
  for (i in seq_along(z)) {
    g <- group[[i]]
    a[[g]] <- a[[g]] + 1 / variance[[i]]
    cstat[[g]] <- cstat[[g]] + z[[i]] / variance[[i]]
  }
  if (any(!is.finite(a)) || any(a <= 0) || any(!is.finite(cstat))) {
    stop("Collapsed market-intercept sufficient statistics are unusable.",
         call. = FALSE)
  }
  tau_sq <- tau^2
  sigma_b_sq <- sigma_b^2
  if (!is.finite(tau_sq) || !is.finite(sigma_b_sq)) {
    stop("Collapsed market-intercept scale draws are too large.", call. = FALSE)
  }
  h <- 1 / (1 + tau_sq * a)
  ah <- a * h
  A <- sum(ah)
  d_bar <- mean(h)
  den_b <- d_bar + sigma_b_sq * A
  if (!is.finite(A) || A <= 0 || !is.finite(den_b) || den_b <= 0) {
    stop("Collapsed market-intercept marginal variance is unusable.",
         call. = FALSE)
  }
  c_bar <- sum(h * cstat) / A
  b0_mean <- sigma_b_sq * A * c_bar / den_b
  b0_variance <- sigma_b_sq * d_bar / den_b
  if (!is.finite(b0_mean) || !is.finite(b0_variance) || b0_variance <= 0) {
    stop("Collapsed market-intercept posterior moments are unusable.",
         call. = FALSE)
  }

  d <- if (tau_sq == 0) rep(0, n_group) else tau_sq / (1 + tau_sq * a)
  d_sum <- sum(d)
  if (d_sum > 0) {
    slope_q <- -d * a + d * sum(d * a) / d_sum
    q_conditional_mean <- function(b0) {
      raw <- d * (cstat - b0 * a)
      raw - d * sum(raw) / d_sum
    }
    q_mean <- q_conditional_mean(b0_mean)
    q_variance_given_b <- pmax(0, d - d^2 / d_sum)
    q_variance <- q_variance_given_b + slope_q^2 * b0_variance
  } else {
    slope_q <- rep(0, n_group)
    q_conditional_mean <- function(b0) rep(0, n_group)
    q_mean <- rep(0, n_group)
    q_variance <- rep(0, n_group)
  }
  b0_draw <- stats::rnorm(1L, mean = b0_mean, sd = sqrt(b0_variance))
  q_draw <- if (d_sum > 0) {
    raw <- stats::rnorm(n_group, mean = d * (cstat - b0_draw * a),
                        sd = sqrt(d))
    centered <- raw - d * sum(raw) / d_sum
    centered[[n_group]] <- -sum(centered[-n_group])
    centered
  } else {
    rep(0, n_group)
  }
  if (any(!is.finite(q_mean)) || any(!is.finite(q_variance)) ||
      any(q_variance < 0) || any(!is.finite(q_draw))) {
    stop("Collapsed market-intercept quality recovery produced non-finite values.",
         call. = FALSE)
  }
  list(intercept = mu_b + b0_draw,
       intercept_mean = mu_b + b0_mean,
       intercept_sd = sqrt(b0_variance),
       quality_draw = q_draw,
       quality_mean = q_mean,
       quality_sd = sqrt(q_variance),
       intercept_offset_mean = b0_mean,
       intercept_offset_variance = b0_variance)
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
#' @param draws Optional list containing numeric \code{share_mean} for
#'   absolute-share fits or \code{relative_share_predictor} for relative-share
#'   fits (draws by observations), \code{sigma_share}, and the economic quality scale
#'   \code{sigma_b_strat_economic} (or historical \code{sigma_b_strat}). Under
#'   \code{quality_scale_parameterization = "total_fraction"}, unrestricted
#'   Stan extraction supplies the economic alias; a restricted raw extraction
#'   must instead include both \code{sigma_b_strat} (the total scale) and
#'   \code{log_sigma_share} (the fraction logit) for reconstruction. When
#'   omitted, these are extracted from the fit's Stan object. Joint recovery also needs
#'   \code{rho} and \code{sigma_margin}; non-exact fits use
#'   \code{pred_marginNormInv}, while exact market fits use the existing
#'   component maps and structural draws.
#' @param conditioning \code{"auto"} (default) conditions on shares and
#'   matched margins when Gaussian rho is active, and on shares otherwise.
#'   \code{"share"} explicitly returns the share-only conditional distribution.
#' @param k_cache Optional internal exact ordered-k cache from
#'   \code{\link{posterior_policy_draws}}. Reusing it makes quality recovery
#'   and policy simulation condition on the same realized market state.
#'
#' @details
#' The quality vector is the centered group-level quality in the collapsed Stan
#' model. Global fits retain one block across groups; market fits have one
#' persistent block for every market and use the same economic group labels in
#' the separate ownership/ranking map. For each common-parameter draw this
#' function conditions that prior on the observed \code{log_shareIn} under
#' \code{share_basis = "absolute"}, or its orthonormal within-market-year
#' contrasts under \code{share_basis = "relative"}. The relative conditional
#' projects out every market-year level, including when matched Gaussian
#' margins are observed and rho is active. Exact market ordered-k fits first
#' realize one existing component state per market and use its FOC margin mean.
#' The returned quality draw is generated
#' after the conditioning step, so HMC does not need to sample quality
#' coordinates. Reconstruction avoids a dense group covariance matrix.
#' Relative-share fits derive one calibration intercept per market-year only
#' after drawing this same conditional quality state. Their native collapsed
#' \code{pred_shareIn} is a zero-quality reference composition, not a recovered
#' posterior quality realization; use the returned quality state or a reusable
#' policy draw plan for posterior reconstruction and prediction.
#'
#' @return A \code{posterior_quality} list. \code{$draws} is a draws by scoped
#' units numeric matrix of conditional quality draws; \code{$conditional_means}
#' and \code{$conditional_sds} have the corresponding per-draw conditional
#' moments. Column names and \code{$quality_units} give stable scoped units, while
#' \code{$draw_indices} records the original posterior draw number for each
#' row. \code{$draws_long} is a tidy version with columns \code{draw_index},
#' \code{market}, \code{group}, \code{quality_unit}, \code{quality},
#' \code{conditional_mean}, and \code{conditional_sd}. \code{$summary} gives
#' mixture posterior means and SDs by scoped unit (including both common-
#' parameter and conditional quality uncertainty). \code{$realized_k} records
#' exact market states when applicable, and \code{$metadata} records whether
#' moments are conditional on a realized state.
#' Relative-share results additionally contain
#' \code{$derived_market_year_intercept_draws} and
#' \code{$derived_market_year_labels}; these are calibration levels computed
#' from the same recovered quality draw and outside-share state, with no
#' independent intercept prior.
#'
#' @export
posterior_quality <- function(fit, n_draws = NULL, seed = NULL, draws = NULL,
                              conditioning = c("auto", "share"),
                              k_cache = NULL) {
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
    joint <- .posterior_quality_joint_draws(
      fit, parsed, draws, conditioning, draw_indices = draw_indices,
      seed = if (is.null(seed)) NULL else as.integer((as.double(seed) + 1) %% .Machine$integer.max),
      k_cache = k_cache
    )
    quality <- matrix(0, nrow = n_use, ncol = parsed$n_quality,
                      dimnames = list(as.character(draw_indices), parsed$quality_units))
    means <- quality
    sds <- quality
    market_intercepts <- market_intercept_means <- market_intercept_sds <- NULL
    intercept_collapsed <- identical(
      parsed$modes$market_intercept_parameterization, "collapsed"
    )
    if (intercept_collapsed) {
      dn <- list(as.character(draw_indices), parsed$market_labels)
      market_intercepts <- matrix(NA_real_, nrow = n_use,
                                  ncol = parsed$n_market, dimnames = dn)
      market_intercept_means <- market_intercept_sds <- market_intercepts
    }
    for (row in seq_len(n_use)) {
      draw <- draw_indices[[row]]
      tau_sq <- parsed$tau[[draw]]^2
      if (!is.finite(tau_sq))
        stop("Posterior draw ", draw,
             " has an unusable sigma_b_strat squared.", call. = FALSE)
      if (.posterior_demand_shock_model(fit) && isTRUE(fit@stan_data$likelihood_mode == 3L)) {
        for (block in seq_along(parsed$blocks)) {
          units <- seq.int(parsed$block_starts[[block]], length.out = parsed$blocks[[block]])
          q <- stats::rnorm(length(units), sd = parsed$tau[[draw]])
          quality[row, units] <- q - mean(q)
          sds[row, units] <- parsed$tau[[draw]] * sqrt(1 - 1 / length(units))
        }
      } else if (isTRUE(parsed$relative)) {
        sigma_share <- if (is.null(parsed$sigma_share_market))
          rep(parsed$sigma[[draw]], length(parsed$y)) else
          parsed$sigma_share_market[draw, parsed$market]
        residual <- parsed$y - parsed$share_mean[draw, ]
        variance <- sigma_share^2
        if (!is.null(joint)) {
          idx <- joint$idx
          sigma_margin <- joint$sigma_by_market[row, parsed$market[idx]]
          residual[idx] <- residual[idx] - joint$rho_by_market[row, parsed$market[idx]] * sigma_share[idx] /
            sigma_margin * (joint$y[idx] - log(joint$prediction[row, idx] + 1e-9))
          variance[idx] <- variance[idx] * (1 - joint$rho_by_market[row, parsed$market[idx]]^2)
        }
        arguments <- list(residual = residual, variance = variance,
          quality_index = parsed$quality,
          quality_block = rep.int(seq_along(parsed$blocks), parsed$blocks),
          cell = as.integer(fit@stan_data$market_year_idx %||% rep.int(1L, length(parsed$y))),
          tau = parsed$tau[[draw]])
        moments <- do.call(.relative_share_gaussian,
          c(arguments, list(return_moments = TRUE)))
        means[row, ] <- moments$conditional_mean
        sds[row, ] <- sqrt(pmax(0, moments$conditional_variance))
        quality[row, ] <- as.numeric(do.call(.relative_share_quality_draw, arguments))
      } else if (intercept_collapsed) {
        for (m in seq_len(parsed$n_market)) {
          obs <- which(parsed$market == m)
          units <- which(parsed$quality_lookup$market_index == m)
          if (!length(obs) || !length(units))
            stop("Collapsed market-intercept recovery has an empty market block.",
                 call. = FALSE)
          local_group <- match(parsed$quality[obs], units)
          if (anyNA(local_group))
            stop("Collapsed market-intercept recovery has an invalid quality map.",
                 call. = FALSE)
          sigma_share <- if (is.null(parsed$sigma_share_market))
            parsed$sigma[[draw]] else parsed$sigma_share_market[draw, m]
          if (!is.finite(sigma_share) || sigma_share <= 0)
            stop("Posterior draw ", draw,
                 " has an unusable market share residual scale.", call. = FALSE)
          z <- parsed$y[obs] - parsed$share_mean[draw, obs]
          variance <- rep(sigma_share^2, length(obs))
          if (!is.null(joint)) {
            jpos <- which(joint$idx %in% obs)
            if (length(jpos)) {
              jidx <- joint$idx[jpos]
              local_obs <- match(jidx, obs)
              sigma_margin <- joint$sigma_by_market[row, m]
              slope <- joint$rho[[row]] * sigma_share / sigma_margin
              margin_residual <- joint$y[jidx] -
                log(joint$prediction[row, jidx] + 1e-9)
              if (!is.finite(slope) ||
                  !is.finite(sigma_share^2 * (1 - joint$rho[[row]]^2)) ||
                  any(!is.finite(margin_residual)))
                stop("Unusable market-specific Gaussian rho sufficient statistics.",
                     call. = FALSE)
              z[local_obs] <- z[local_obs] - slope * margin_residual
              variance[local_obs] <- sigma_share^2 *
                (1 - joint$rho[[row]]^2)
            }
          }
          recovered <- .posterior_quality_market_intercept_conditional(
            z = z, variance = variance, group = local_group,
            tau = sqrt(tau_sq), mu_b = parsed$mu_b_market[[draw]],
            sigma_b = parsed$sigma_b_market[[draw]]
          )
          quality[row, units] <- recovered$quality_draw
          means[row, units] <- recovered$quality_mean
          sds[row, units] <- recovered$quality_sd
          market_intercepts[row, m] <- recovered$intercept
          market_intercept_means[row, m] <- recovered$intercept_mean
          market_intercept_sds[row, m] <- recovered$intercept_sd
        }
      } else {
        sigma_sq <- parsed$sigma[[draw]]^2
        if (!is.finite(sigma_sq) || sigma_sq <= 0) {
          stop("Posterior draw ", draw,
               " has an unusable sigma_share squared.", call. = FALSE)
        }
        ## Accumulate observation-level sufficient statistics once. The block
        ## loop below only slices these vectors, keeping market scope linear in
        ## N + Q rather than rescanning all observations for every market.
        a_all <- tabulate(parsed$quality, nbins = parsed$n_quality) / sigma_sq
        c_all <- numeric(parsed$n_quality)
        share_residual <- (parsed$y - parsed$share_mean[draw, ]) / sigma_sq
        for (jj in seq_along(parsed$quality)) {
          qj <- parsed$quality[[jj]]
          c_all[[qj]] <- c_all[[qj]] + share_residual[[jj]]
        }
        if (!is.null(joint)) {
          slope <- joint$rho[[row]] * joint$sigma[[row]] / parsed$sigma[[draw]]
          variance <- joint$sigma[[row]]^2 * (1 - joint$rho[[row]]^2)
          if (!is.finite(variance) || variance <= 0 || !is.finite(slope))
            stop("Unusable Gaussian rho sufficient statistics.", call. = FALSE)
          jidx <- joint$idx
          residual <- joint$y[jidx] - log(joint$prediction[row, jidx] + 1e-9) -
            slope * (parsed$y[jidx] - parsed$share_mean[draw, jidx])
          if (any(!is.finite(residual)))
            stop("Unusable Gaussian joint quality residuals.", call. = FALSE)
          for (jj in seq_along(jidx)) {
            qj <- parsed$quality[[jidx[[jj]]]]
            a_all[[qj]] <- a_all[[qj]] + slope^2 / variance
            c_all[[qj]] <- c_all[[qj]] - slope * residual[[jj]] / variance
          }
        }
        if (any(!is.finite(a_all)) || any(!is.finite(c_all)))
          stop("Posterior draw ", draw, " produced non-finite Gaussian sufficient statistics.", call. = FALSE)
        for (block in seq_along(parsed$blocks)) {
          first <- parsed$block_starts[[block]]
          last <- first + parsed$blocks[[block]] - 1L
          units <- first:last
          a <- a_all[units]
          cstat <- c_all[units]
          d <- if (tau_sq == 0) rep(0, length(units)) else tau_sq / (1 + tau_sq * a)
          if (any(!is.finite(d))) stop("Posterior draw ", draw, " produced non-finite conditional quality variances.", call. = FALSE)
          raw_mean <- d * cstat
          den <- sum(d)
          if (!is.finite(den) || any(!is.finite(raw_mean)))
            stop("Posterior draw ", draw, " produced non-finite conditional quality moments.", call. = FALSE)
          if (den > 0) {
            means[row, units] <- raw_mean - d * (sum(raw_mean) / den)
            sds[row, units] <- sqrt(pmax(0, d - d * d / den))
            u <- stats::rnorm(length(units), mean = raw_mean, sd = sqrt(d))
            quality[row, units] <- u - d * (sum(u) / den)
            quality[row, units[[length(units)]]] <- -sum(quality[row, units[-length(units)]])
          }
        }
      }
    }
    mixture_mean <- colMeans(means)
    mixture_variance <- colMeans(sds^2 + means^2) - mixture_mean^2
    mixture_variance <- pmax(0, mixture_variance)
    summary <- parsed$quality_lookup
    summary$posterior_mean <- as.numeric(mixture_mean)
    summary$posterior_sd <- sqrt(as.numeric(mixture_variance))
    summary <- summary[c("market_label", "group_label", "quality_unit", "posterior_mean", "posterior_sd")]
    names(summary)[names(summary) == "market_label"] <- "market"
    names(summary)[names(summary) == "group_label"] <- "group"
    draws_long <- data.frame(
      draw_index = rep(as.integer(draw_indices), each = parsed$n_quality),
      market = rep(parsed$quality_lookup$market_label, times = n_use),
      group = rep(parsed$quality_lookup$group_label, times = n_use),
      quality_unit = rep(parsed$quality_units, times = n_use),
      quality = as.vector(t(quality)),
      conditional_mean = as.vector(t(means)),
      conditional_sd = as.vector(t(sds)),
      stringsAsFactors = FALSE
    )
    if (!is.null(joint$realized_k)) {
      market_for_unit <- parsed$quality_lookup$market_index
      ## A legacy/global one-market fit has no market label on its shared
      ## quality unit, but its exact state is still the sole market state.
      if (ncol(joint$realized_k) == 1L) market_for_unit[is.na(market_for_unit)] <- 1L
      market_for_unit <- match(market_for_unit, seq_len(ncol(joint$realized_k)))
      draws_long$realized_k <- as.vector(t(joint$realized_k[,
        market_for_unit, drop = FALSE]))
    }
    quality_labels <- parsed$quality_labels
    result <- structure(list(
      draws = quality,
      conditional_means = means,
      conditional_sds = sds,
      group_labels = parsed$labels,
      quality_labels = quality_labels,
      quality_units = parsed$quality_units,
      market_labels = as.character(parsed$quality_lookup$market_label),
      quality_scope = parsed$quality_scope,
      market_intercept_parameterization =
        parsed$modes$market_intercept_parameterization,
      realized_k = if (is.null(joint)) NULL else joint$realized_k,
      market_intercept_draws = market_intercepts,
      market_intercept_conditional_means = market_intercept_means,
      market_intercept_conditional_sds = market_intercept_sds,
      market_intercept_labels = if (intercept_collapsed)
        parsed$market_labels else NULL,
      draw_indices = as.integer(draw_indices),
      draws_long = draws_long,
      summary = summary,
      metadata = list(
        available_draws = as.integer(n_available),
        selected_draws = as.integer(n_use),
        quality_parameterization = parsed$modes$quality_parameterization,
        share_likelihood = parsed$modes$share_likelihood,
        likelihood_mode = parsed$modes$likelihood_mode,
        conditioning = if (isTRUE(parsed$relative)) {
          if (is.null(joint)) "observed relative-share contrasts only" else
            "observed relative-share contrasts and matched Gaussian margins"
        } else if (is.null(joint)) "observed log shares only" else
          "observed log shares and matched Gaussian margins",
        conditional_moments = if (intercept_collapsed) {
          "quality_moments_integrate_market_intercept_conditional_on_observed_state"
        } else if (!is.null(joint$realized_k)) {
          "conditional_on_realized_k"
        } else {
          "conditional_on_observed_quality_state"
        },
        market_intercept_moments = if (intercept_collapsed)
          "conditional_on_observed_shares_margins_and_realized_k" else NULL,
        quality_scope = parsed$quality_scope,
        centered_constraint = if (parsed$quality_scope == "global")
          "sum(group quality) = 0" else "sum(group quality) = 0 within each market"
      )
    ), class = c("posterior_quality", "list"))
    if (isTRUE(parsed$relative)) {
      derived <- .lq_relative_derived_intercepts(fit, draws, draw_indices, quality)
      result$derived_market_year_intercept_draws <- derived$draws
      result$derived_market_year_labels <- derived$labels
      result$metadata$share_basis <- "relative"
      result$metadata$observation_basis <- "ilr_quantity_shares_v1"
    }
    result
  }
  .posterior_quality_with_seed(seed, calculate)
}
