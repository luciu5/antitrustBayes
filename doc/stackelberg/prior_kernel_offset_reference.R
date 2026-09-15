# Constants omitted by the audited version-6 Stan prior statements. Sampling
# targets remain unchanged. Counts use actual parameter coordinates, including
# unidentified but properly normalized raw directions in centered qualities.
.bayes_prior_offset <- function(dat, parameter_dims) {
  n <- function(name) {
    if (!name %in% names(parameter_dims)) return(0)
    d <- parameter_dims[[name]]
    if (!length(d)) return(1)
    if (anyNA(d) || any(d < 0) || any(d != as.integer(d))) {
      stop("Invalid Stan parameter dimensions for ", name, ".")
    }
    prod(d)
  }
  fixed <- function(sd, count) {
    if (!count) return(0)
    if (length(sd) != 1L || !is.finite(sd) || sd <= 0) {
      stop("Proper evidence requires finite positive prior scales.")
    }
    count * (-log(sd) - .5 * log(2 * pi))
  }
  standard <- function(count) fixed(1, count)
  half <- function(sd, count) fixed(sd, count) + count * log(2)
  beta <- function(a, b, count) {
    if (!count) return(0)
    if (length(a) != 1L || length(b) != 1L || !is.finite(a) ||
        !is.finite(b) || a <= 0 || b <= 0) {
      stop("Proper evidence requires positive finite Beta hyperparameters.")
    }
    -count * lbeta(a, b)
  }
  if (isTRUE(dat$use_s0_market_pool == 1L) &&
      isTRUE(dat$use_s0_market_year_pool == 1L)) {
    stop("Both outside-share pooling flags declare coordinates without a prior; proper evidence is undefined.")
  }
  if (isTRUE(dat$use_hmt == 1L) && !isTRUE(dat$likelihood_mode == 2L)) {
    stop("Proper evidence is unsupported with the non-normalized HMT penalty; use bridge_normalization = 'legacy' for raw pseudo-evidence.")
  }
  strat_count <- n("sigma_b_strat")
  strat <- fixed(dat$prior_sigma_meanval_strat, strat_count)
  if (strat_count) strat <- strat - strat_count *
    stats::pnorm(1 / dat$prior_sigma_meanval_strat, log.p = TRUE)
  rho <- 0
  if (n("rho_raw")) {
    eta <- dat$prior_lkj
    if (length(eta) != 1L || !is.finite(eta) || eta <= 0) {
      stop("Proper evidence requires a finite positive correlation concentration.")
    }
    # t=rho^2 gives the exact incomplete-Beta integral of the bounded kernel.
    log_integral <- lbeta(.5, eta) +
      stats::pbeta(.95^2, .5, eta, log.p = TRUE) - log(.95)
    rho <- -n("rho_raw") * log_integral
  }
  df <- 0
  if (n("nu_minus_2")) {
    a <- dat$prior_nu_shape; b <- dat$prior_nu_rate
    if (length(a) != 1L || length(b) != 1L || !is.finite(a) ||
        !is.finite(b) || a <= 0 || b <= 0) {
      stop("Proper evidence requires positive finite Gamma hyperparameters.")
    }
    df <- n("nu_minus_2") * (a * log(b) - lgamma(a))
  }
  components <- c(
    global_alpha = fixed(dat$prior_alpha_sd,
      n("mu_log_a") + n("mu_log_a_free") + n("log_kappa_vec")),
    alpha_scale = half(dat$prior_sigma_alpha, n("sigma_log_a_vec")),
    alpha_raw = standard(n("r_market_a_raw")),
    alpha_centered = standard(n("log_a_market_pool")),
    alpha_cournot_centered = standard(n("log_kappa_market")),
    market_mean = standard(n("mu_b_market_vec")),
    market_scale = half(.5, n("sigma_b_market_vec")),
    market_raw = standard(n("b_market_raw")),
    quality_raw = standard(n("b_group_raw")),
    quality_scale = strat,
    cutoff = beta(dat$prior_cutoff_alpha, dat$prior_cutoff_beta, n("cutoff_share")),
    fringe_quality = fixed(dat$prior_sigma_meanval_fringe, n("mu_b_fringe_vec")),
    outside_mean = fixed(dat$s0_prior_sd, n("mu_s0_pool_vec")),
    outside_scale = half(dat$prior_sigma_s0_pool, n("sigma_s0_pool_vec")),
    outside_market_raw = standard(n("s0_market_raw")),
    outside_period_raw = standard(n("s0_pool_raw")),
    outside_market_centered = standard(n("eta_s0_market")),
    outside_period_centered = standard(n("eta_s0_pool")),
    outside_unpooled = fixed(dat$s0_prior_sd, n("s0_logit")),
    outside_regression = fixed(dat$prior_sigma_beta_s0, n("beta_s0")),
    bargaining = beta(dat$bargaining_prior_alpha, dat$bargaining_prior_beta,
      n("bargaining_unit")),
    share_scale = fixed(dat$prior_log_sigma_share_sd, n("log_sigma_share")),
    margin_scale = fixed(dat$prior_log_sigma_margin_sd, n("log_sigma_margin")),
    share_market_scale_raw = standard(n("log_sigma_share_market_raw")),
    margin_market_scale_raw = standard(n("log_sigma_margin_market_raw")),
    share_market_scale_hyper = half(dat$prior_sigma_resid_market,
      n("sigma_log_sigma_share_market")),
    margin_market_scale_hyper = half(dat$prior_sigma_resid_market,
      n("sigma_log_sigma_margin_market")),
    df = df, correlation = rho
  )
  if (any(!is.finite(components))) stop("The proper evidence offset is not finite.")
  list(offset = sum(components), components = components)
}

.bayes_is_stackelberg <- function(fit) {
  isTRUE(fit@stan_data$structure == 2L) ||
    identical(fit@metadata$structure, "stackelberg")
}

.bayes_evidence_offset <- function(fit) {
  validate_fit_spec(fit)
  sf <- fit@stanfit
  if (!methods::is(sf, "stanfit")) {
    stop("Proper evidence requires a fitted, audited Stan model.")
  }
  canonical <- function(code) trimws(gsub("\r", "", code, fixed = TRUE))
  available <- stanmodels[intersect(names(stanmodels),
    c("logit_supply", "ces_supply", "logit_stackelberg", "ces_stackelberg"))]
  graph <- names(available)[vapply(available, function(model)
    identical(canonical(model@model_code), canonical(sf@stanmodel@model_code)),
    logical(1))]
  if (length(graph) != 1L || !startsWith(graph, paste0(fit@demand_model, "_")) ||
      (.bayes_is_stackelberg(fit) != endsWith(graph, "_stackelberg"))) {
    stop("Proper evidence is unsupported for this unrecognized prior graph.")
  }
  code <- sf@stanmodel@model_code
  block <- regmatches(code, regexec("(?s)\\nparameters\\s*\\{(.*?)\\n\\}",
    code, perl = TRUE))[[1L]]
  if (length(block) != 2L) stop("Cannot audit the Stan parameter block.")
  declarations <- trimws(strsplit(gsub("//[^\n]*", "", block[2L]),
    ";", fixed = TRUE)[[1L]])
  declarations <- declarations[nzchar(declarations)]
  pn <- vapply(declarations, function(z) {
    match <- regmatches(z, regexec("([A-Za-z_][A-Za-z_0-9]*)\\s*$",
      z, perl = TRUE))[[1L]]
    if (length(match) != 2L) stop("Cannot audit a Stan parameter declaration.")
    match[2L]
  }, character(1))
  if (!all(pn %in% names(sf@par_dims))) stop("Stan parameter dimensions are missing.")
  result <- .bayes_prior_offset(fit@stan_data, sf@par_dims[pn])
  result$graph <- graph
  result$interpretation <- "normalized priors and transformed-outcome likelihood"
  result
}

.bayes_comparison_units <- function(fits, loo_unit) {
  for (fit in fits) if (.bayes_is_stackelberg(fit)) {
    dat <- fit@stan_data
    if (loo_unit == "observation") {
      stop("Observation LOO is unsupported for whole-equilibrium Stackelberg mixtures.")
    }
    if (loo_unit == "period" && isTRUE(dat$k_scope_code == 1L)) {
      stop("Period LOO is unsupported when one Stackelberg k is shared across the market's years.")
    }
    if (loo_unit == "cluster") {
      unit <- dat$mky_to_partition_unit[dat$market_year_idx]
      split_clusters <- tapply(dat$cluster, unit, function(x) length(unique(x)))
      if (any(split_clusters != 1L)) {
        stop("LOO clusters must contain whole Stackelberg partition units.")
      }
    }
  }
  invisible(TRUE)
}
