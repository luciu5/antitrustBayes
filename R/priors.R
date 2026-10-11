# Prespecified priors and model-universe weights used for model comparison.

#' Prespecified Parameter-Prior Manifest
#'
#' Returns one of three proper, data-independent prior specifications intended
#' for marginal-likelihood comparisons.  The tight and diffuse manifests keep
#' the baseline centers fixed and multiply the standard deviations of the
#' outside-share, standardized price-coefficient, and margin-noise priors by
#' one half and two, respectively.
#'
#' @details For multiple markets, the market-intercept SD prior can be selected
#'   through \code{priors}: \code{prior_sigma_b_market_family = "half_normal"}
#'   (historical default, canonical code 1) uses
#'   \code{prior_sigma_b_market_scale = 0.5};
#'   \code{prior_sigma_b_market_family = "gamma"} (code 2) uses
#'   \code{prior_sigma_b_market_shape = 2} and
#'   \code{prior_sigma_b_market_rate = 5}. Gamma uses a rate, not a scale.
#'   The prior applies to sampled and collapsed market intercepts in both
#'   demand systems and all conduct families. It is inactive for one market.
#'   These hyperparameters stay fixed across the tight/diffuse manifests.
#' @examples
#' priors <- bayes_prior_manifest()
#' priors$prior_sigma_b_market_family <- "gamma"
#' priors$prior_sigma_b_market_shape <- 2
#' priors$prior_sigma_b_market_rate <- 5
#' # Supply priors to bayes_est(..., priors = priors).
#' @param name Character. One of \code{"baseline"}, \code{"tight"}, or
#'   \code{"diffuse"}.
#' @return A named list suitable for the \code{priors} argument to
#'   \code{\link{bayes_est}}.
#' @export
bayes_prior_manifest <- function(name = c("baseline", "tight", "diffuse")) {
  name <- match.arg(name)
  multiplier <- switch(name, baseline = 1, tight = 0.5, diffuse = 2)

  list(
    prior_sigma_share = 0.5,
    prior_log_sigma_share_sd = 0.5,
    prior_sigma_margin = 0.75,
    prior_log_sigma_margin_sd = 1 * multiplier,
    prior_sigma_meanval_strat = 1,
    prior_sigma_meanval_fringe = 1,
    # alpha is standardized by the empirical price SD before entering Stan.
    # This prior is proper and deliberately does not use observed margins.
    prior_alpha_mean = 0,
    prior_alpha_sd = 1 * multiplier,
    s0_prior_mean = -1.446786403764666,
    s0_prior_sd = 0.456234014550216 * multiplier,
    prior_sigma_alpha = 0.5,
    # Market-intercept SD: 1 = half-normal, 2 = Gamma (shape/rate).
    # Preserve the historical default; callers can select either family.
    prior_sigma_b_market_family = 1L,
    prior_sigma_b_market_scale = 0.5,
    prior_sigma_b_market_shape = 2,
    prior_sigma_b_market_rate = 5,
    prior_sigma_s0_pool = 0.5,
    prior_sigma_beta_s0 = 1,
    prior_sigma_resid_market = 0.35,
    prior_cutoff_alpha = 3,
    prior_cutoff_beta = 100,
    # Fixed default df for the joint bivariate Student-t; 5 is a practical
    # default for sparse information sets (e.g. single-market PNB fits).
    resid_df_fixed = 5,
    # Estimated-df prior (gamma on nu - 2): mean 8, mode 4. Provisional --
    # deferred pending prior-predictive checks, not a settled default.
    prior_nu_shape = 2,
    prior_nu_rate = 0.25,
    # LKJ-style concentration for the share-margin correlation rho.
    prior_lkj = 2.0,
    bargaining_prior_lower = 0,
    bargaining_prior_upper = 0.9999,
    bargaining_prior_alpha = 1,
    bargaining_prior_beta = 1
  )
}

canonical_supply_name <- function(x) {
  key <- tolower(gsub("[^a-z]", "", x))
  out <- c(
    bertrand = "bertrand",
    cournot = "cournot",
    auction = "auction",
    secondscoreauction = "auction",
    monopolistic = "monopolistic",
    monopolisticcompetition = "monopolistic",
    moncom = "monopolistic",
    mc = "monopolistic",
    bargaining = "bargaining",
    nashbargaining = "bargaining"
  )[key]
  if (anyNA(out)) {
    stop("Unknown supply model: ", paste(unique(x[is.na(out)]), collapse = ", "))
  }
  unname(out)
}

#' Prespecified Model Priors for Retail-Deposit Conduct
#'
#' The institutional specification assigns 30 percent prior probability to
#' Bertrand, Cournot, and monopolistic competition and 10 percent to the
#' second-score auction. Equal and core-three specifications are provided as
#' candidate-set sensitivities.
#'
#' @param models Character vector of supply-model names.
#' @param scheme Character. \code{"institutional"}, \code{"equal"}, or
#'   \code{"core_three"}.
#' @return A named numeric vector summing to one in the order of \code{models}.
#' @export
bank_model_priors <- function(
    models = c("bertrand", "cournot", "auction", "monopolistic"),
    scheme = c("institutional", "equal", "core_three")) {
  scheme <- match.arg(scheme)
  canonical <- canonical_supply_name(models)
  if (anyDuplicated(canonical)) {
    stop("'models' contains duplicate conduct models.")
  }

  if (scheme == "equal") {
    weights <- rep(1 / length(models), length(models))
  } else if (scheme == "core_three") {
    if (any(!canonical %in% c("bertrand", "cournot", "monopolistic"))) {
      stop("The core_three scheme only supports Bertrand, Cournot, and MC.")
    }
    weights <- rep(1 / length(models), length(models))
  } else {
    base <- c(
      bertrand = 0.30,
      cournot = 0.30,
      auction = 0.10,
      monopolistic = 0.30
    )
    if (any(!canonical %in% names(base))) {
      stop("The institutional banking prior does not include bargaining.")
    }
    weights <- unname(base[canonical])
    weights <- weights / sum(weights)
  }
  stats::setNames(weights, models)
}

resolve_prior_manifest <- function(prior_manifest, priors) {
  if (!is.list(priors) || (length(priors) && is.null(names(priors)))) {
    stop("'priors' must be a named list.")
  }
  base <- bayes_prior_manifest(prior_manifest)
  if (length(priors)) base[names(priors)] <- priors
  .normalize_market_intercept_sd_prior(base)
}

# Public priors accept readable family names; Stan and saved provenance use
# canonical integer codes, matching the existing numeric hyperparameter list.
.normalize_market_intercept_sd_prior <- function(priors) {
  family <- priors$prior_sigma_b_market_family
  if (is.character(family) && length(family) == 1L && !is.na(family) &&
      family %in% c("half_normal", "gamma")) {
    family <- unname(c(half_normal = 1L, gamma = 2L)[family])
  }
  if (!is.numeric(family) || length(family) != 1L || !is.finite(family) ||
      !family %in% c(1L, 2L)) {
    stop("Prior 'prior_sigma_b_market_family' must be 'half_normal' or 'gamma' ",
         "(or code 1 or 2).", call. = FALSE)
  }
  priors$prior_sigma_b_market_family <- as.integer(family)
  for (name in c("prior_sigma_b_market_scale", "prior_sigma_b_market_shape",
                 "prior_sigma_b_market_rate")) {
    value <- priors[[name]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value <= 0) {
      stop("Prior '", name, "' must be a positive finite numeric scalar.", call. = FALSE)
    }
  }
  priors
}

package_provenance <- function(data, prior_manifest, model_prior = NA_real_) {
  commit <- Sys.getenv("ANTITRUSTBAYES_COMMIT", unset = NA_character_)
  if (is.na(commit)) {
    commit <- tryCatch(
      system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)[1],
      error = function(e) NA_character_,
      warning = function(w) NA_character_
    )
  }
  package_version <- tryCatch(
    as.character(utils::packageVersion("antitrustBayes")),
    error = function(e) {
      description <- tryCatch(
        utils::packageDescription("antitrustBayes"),
        error = function(e) NULL
      )
      if (!is.null(description$Version)) description$Version else "development"
    }
  )
  list(
    data_hash = digest::digest(data, algo = "sha256"),
    prior_manifest = prior_manifest,
    model_prior = model_prior,
    package_version = package_version,
    commit = commit
  )
}
