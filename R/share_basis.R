# Scientific observation basis, distinct from sampler coordinates.
.resolve_share_basis <- function(share_basis = c("relative", "absolute"),
                                 share_likelihood = "gaussian",
                                 margin_likelihood = NULL,
                                 market_intercept_parameterization = NULL,
                                 market_s0_parameterization = "legacy",
                                 priors = list()) {
  basis <- match.arg(share_basis, c("relative", "absolute"))
  if (basis == "relative") {
    if (!identical(share_likelihood[[1L]], "gaussian") ||
        (!is.null(margin_likelihood) &&
         !identical(margin_likelihood[[1L]], "gaussian")))
      stop("share_basis = 'relative' requires Gaussian shares and margins; ",
           "use share_basis = 'absolute' for Student-t or hybrid likelihoods.",
           call. = FALSE)
    if (!is.null(market_intercept_parameterization))
      stop("market_intercept_parameterization applies only to share_basis = 'absolute'; ",
           "relative shares derive the intercept per market-year.", call. = FALSE)
    if (!identical(market_s0_parameterization[[1L]], "legacy"))
      stop("Relative shares require market_s0_parameterization = 'legacy'; ",
           "rotated/joint coordinates contain a sampled intercept.", call. = FALSE)
    intercept_priors <- grep("^prior_sigma_b_market", names(priors), value = TRUE)
    defaults <- bayes_prior_manifest("baseline")
    changed_intercept_prior <- vapply(intercept_priors, function(key) {
      value <- priors[[key]]
      if (key == "prior_sigma_b_market_family" && is.character(value) &&
          length(value) == 1L && !is.na(value))
        value <- switch(value, half_normal = 1L, gamma = 2L, value)
      !isTRUE(all.equal(value, defaults[[key]], check.attributes = FALSE))
    }, logical(1))
    if (any(changed_intercept_prior))
      stop("Market-intercept priors apply only to share_basis = 'absolute'.",
           call. = FALSE)
  }
  list(name = basis, code = if (basis == "relative") 1L else 2L,
       observation_measure = if (basis == "relative")
         "ilr_quantity_shares_v1" else "absolute_log_quantity_shares_v1")
}

.attach_share_basis <- function(stan_data, basis) {
  stan_data$share_basis <- basis$code
  attr(stan_data, "share_basis") <- basis$name
  attr(stan_data, "share_observation_measure") <- basis$observation_measure
  attr(stan_data, "share_dimension") <- if (basis$code == 1L)
    stan_data$N - stan_data$N_market_year else stan_data$N
  stan_data
}

# Missing basis on historical serialized objects denotes their absolute-share
# likelihood. validate_fit_spec() separately rejects an obsolete specification.
.share_basis_from_fit <- function(fit) {
  if (methods::is(fit, "BayesSpecifiedFit")) return("absolute")
  decode <- function(x) {
    if (is.null(x)) return(NULL)
    if (length(x) != 1L || is.na(x))
      stop("Malformed share_basis provenance.", call. = FALSE)
    if (is.numeric(x)) {
      if (!x %in% c(1, 2)) stop("Malformed share_basis code.", call. = FALSE)
      return(if (x == 1) "relative" else "absolute")
    }
    if (!x %in% c("relative", "absolute"))
      stop("Malformed share_basis provenance.", call. = FALSE)
    x
  }
  values <- lapply(list(fit@metadata$share_basis,
                       attr(fit@stan_data, "share_basis", exact = TRUE),
                       fit@stan_data$share_basis), decode)
  values <- values[!vapply(values, is.null, logical(1))]
  if (!length(values)) return("absolute")
  if (length(unique(unlist(values))) != 1L)
    stop("Conflicting share_basis provenance.", call. = FALSE)
  values[[1L]]
}
