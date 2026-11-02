#!/usr/bin/env Rscript

## Bounded production validation for the collapsed Gaussian quality block.
##
## This script deliberately uses only the already installed production Stan
## models and the package's data-preparation helpers.  It creates one-draw
## Fixed_param fits (there is no HMC, compilation, installation, or benchmark
## here), evaluates the Stan model instance directly, and writes a compact
## record of every check before the final strict assertions.

options(stringsAsFactors = FALSE)

cga_script_file <- commandArgs(trailingOnly = FALSE)
cga_script_file <- cga_script_file[startsWith(cga_script_file, "--file=")]
cga_doc_dir <- if (length(cga_script_file)) {
  dirname(normalizePath(sub("^--file=", "", cga_script_file[[1L]]),
                         mustWork = FALSE))
} else {
  file.path(getwd(), "doc")
}
cga_out_file <- file.path(cga_doc_dir,
                          "collapsed_gaussian_production_validation.csv")

if (!requireNamespace("rstan", quietly = TRUE)) {
  stop("The production validation requires the installed rstan package.",
       call. = FALSE)
}
if (!requireNamespace("antitrustBayes", quietly = TRUE)) {
  stop("The production validation requires the installed antitrustBayes package.",
       call. = FALSE)
}
library(rstan)

## The source reference has an intentional global-environment guard around
## its own run_checks() call.  A fresh environment keeps that reference
## independent and prevents its full standalone checks from running here.
cga_reference <- new.env(parent = globalenv())
sys.source(file.path(cga_doc_dir, "collapsed_gaussian_reference.R"),
           envir = cga_reference)

cga_lse <- function(x) {
  x <- as.numeric(x)
  if (!length(x)) return(-Inf)
  m <- max(x)
  if (is.infinite(m)) return(m)
  m + log(sum(exp(x - m)))
}

cga_as_draw <- function(fit, name) {
  z <- rstan::extract(fit, pars = name, permuted = TRUE,
                      inc_warmup = FALSE)[[name]]
  if (is.null(z)) stop("Generated quantity '", name, "' was not returned.",
                       call. = FALSE)
  as.numeric(z)
}

cga_init_one <- function(fit) {
  z <- rstan::get_inits(fit)
  ## get_inits() normally returns one constrained parameter list per chain.
  if (is.list(z) && length(z) == 1L && is.list(z[[1L]]) &&
      !is.null(names(z[[1L]]))) z <- z[[1L]]
  if (!is.list(z) || is.null(names(z))) {
    stop("rstan::get_inits() did not return a named constrained init list.",
         call. = FALSE)
  }
  ## rstan::get_inits() includes transformed parameters and generated
  ## quantities.  Sampling/unconstrain_pars accept the model block only.
  z[intersect(fit@model_pars, names(z))]
}

cga_set_all <- function(init, name, value) {
  if (name %in% names(init)) {
    old <- init[[name]]
    val <- if (length(value) == length(old)) as.numeric(value) else
      rep(as.numeric(value), length(old))
    ## Stan vectors/arrays of length one still require a one-dimensional array
    ## in an init list; rep() alone drops that dimension attribute.
    if (!is.null(dim(old))) dim(val) <- dim(old)
    init[[name]] <- val
  }
  init
}

cga_common_init <- function(fit, demand, supply, beta = NULL,
                            rho = 0, tau = 0.5,
                            log_sigma_share = -1.5,
                            log_sigma_margin = -1.3,
                            s0 = 0.2) {
  init <- cga_init_one(fit)

  ## All validation cases use no outside-share pooling.  Invert the bounded
  ## inverse-logit used by the production graphs so the constrained outside
  ## share is exactly s0 (rather than s0 plus the 1e-6 boundary adjustment).
  s0_eps <- 1e-6
  if (!is.finite(s0) || s0 <= s0_eps || s0 >= 1 - s0_eps) {
    stop("s0 must lie strictly inside the bounded inverse-logit range",
         call. = FALSE)
  }
  init <- cga_set_all(init, "s0_logit",
                      qlogis((s0 - s0_eps) / (1 - 2 * s0_eps)))
  init <- cga_set_all(init, "sigma_b_strat", tau)
  init <- cga_set_all(init, "log_sigma_share", log_sigma_share)
  init <- cga_set_all(init, "log_sigma_margin", log_sigma_margin)
  init <- cga_set_all(init, "rho_raw", atanh(rho / 0.95))
  init <- cga_set_all(init, "mu_b_market_vec", 0)
  init <- cga_set_all(init, "b_market_raw", 0)
  init <- cga_set_all(init, "b_group_centered", if (is.null(beta)) 0 else beta)

  ## Logit's demand scale is exp(mu_log_a_free), except for Cournot where
  ## Stan derives it from log_kappa and the outside share.  CES stores the
  ## common positive scale in mu_log_a directly.
  if (identical(demand, "ces")) {
    init <- cga_set_all(init, "mu_log_a", log(0.8))
  } else if (identical(supply, "cournot")) {
    init <- cga_set_all(init, "log_kappa_vec", log(0.8))
  } else {
    init <- cga_set_all(init, "mu_log_a_free", log(0.8))
  }
  init
}

cga_parameter_labels <- function(fit) {
  ## The order of model_pars and par_dims is the order used by
  ## unconstrain_pars().  All validation models have only scalar/vector
  ## parameter transforms, so this is also the exact unconstrained order.
  out <- character()
  ## stanfit@model_pars also records transformed parameters and generated
  ## quantities.  sigma_share is the first transformed parameter in every
  ## production model, so the preceding roots are exactly the parameters
  ## passed to unconstrain_pars().
  model_end <- match("sigma_share", fit@model_pars)
  if (is.na(model_end)) stop("could not locate transformed-parameter boundary",
                             call. = FALSE)
  model_pars <- fit@model_pars[seq_len(model_end - 1L)]
  for (nm in model_pars) {
    d <- fit@par_dims[[nm]]
    if (length(d) && prod(d) == 0) next
    if (!length(d)) {
      out <- c(out, nm)
    } else {
      ii <- seq_len(prod(d))
      if (length(d) == 1L) {
        out <- c(out, paste0(nm, "[", ii, "]"))
      } else {
        ai <- arrayInd(ii, .dim = d)
        out <- c(out, apply(ai, 1L, function(x) {
          paste0(nm, "[", paste(x, collapse = ","), "]")
        }))
      }
    }
  }
  out
}

cga_model_for <- function(demand, structure = "simultaneous") {
  nm <- if (identical(demand, "ces")) "ces_" else "logit_"
  if (identical(structure, "simultaneous")) nm <- paste0(nm, "supply")
  else nm <- paste0(nm, "stackelberg")
  antitrustBayes:::stanmodels[[nm]]
}

cga_fixture <- function(two_market = FALSE) {
  ## The repeated group labels 1,1 and 4,4 make the centered quality loading
  ## visibly multiproduct.  The two-market split is used only by exact-k.
  list(
    prices = c(1.10, 1.35, 1.60, 1.20, 1.45, 1.80),
    quantities = c(30, 20, 15, 25, 18, 12),
    margins = c(0.30, NA_real_, 0.35, 0.28, NA_real_, 0.32),
    ## Interleave rows so each market has three distinct firms (the repeated
    ## products remain a multiproduct firm in the one-market fixture).
    market = if (two_market) c(1L, 2L, 1L, 2L, 1L, 2L) else rep(1L, 6L),
    year = rep(1L, 6L),
    group = c(1L, 1L, 2L, 3L, 4L, 4L),
    margin_observed = c(TRUE, FALSE, TRUE, TRUE, FALSE, TRUE)
  )
}

## Exact-k supports are configurable for bounded reruns without editing the
## script.  Defaults are the requested fringe {0,1,2,3} and Stackelberg
## {1,2,3}; every configured support must retain at least three states.
cga_support_env <- function(name, default, lower) {
  raw <- Sys.getenv(name, unset = "")
  if (!nzchar(raw)) return(as.integer(default))
  tokens <- trimws(strsplit(raw, ",", fixed = TRUE)[[1L]])
  vals <- unlist(lapply(tokens, function(token) {
    bounds <- trimws(strsplit(token, ":", fixed = TRUE)[[1L]])
    if (length(bounds) == 2L) {
      lo <- suppressWarnings(as.integer(bounds[[1L]]))
      hi <- suppressWarnings(as.integer(bounds[[2L]]))
      if (is.na(lo) || is.na(hi)) return(NA_integer_)
      seq.int(lo, hi)
    } else suppressWarnings(as.integer(token))
  }), use.names = FALSE)
  if (anyNA(vals) || anyDuplicated(vals) || length(vals) < 3L ||
      any(vals < lower)) {
    stop(name, " must be at least three unique integer states at or above ",
         lower, call. = FALSE)
  }
  vals
}
cga_fringe_support <- cga_support_env("CGA_FRINGE_SUPPORT", 0:3, 0L)
cga_stack_support <- cga_support_env("CGA_STACK_SUPPORT", 1:3, 1L)

cga_prepare <- function(demand, supply, fixture, quality = "collapsed",
                        rho = 0, structure = "simultaneous",
                        k_support = NULL, k_prior = NULL,
                        use_rho = NULL) {
  if (is.null(use_rho)) use_rho <- identical(structure, "simultaneous")
  sm <- switch(supply, bertrand = 1L, cournot = 3L,
               stop("unsupported validation supply model", call. = FALSE))
  exact <- !identical(structure, "simultaneous")
  if (exact && is.null(k_support)) {
    stop("exact-k preparation requires k_support", call. = FALSE)
  }
  args <- list(
    prices = fixture$prices,
    quantities = fixture$quantities,
    margins = fixture$margins,
    market = fixture$market,
    year = fixture$year,
    group = fixture$group,
    supply_model = sm,
    use_cutoff = FALSE,
    market_side = "output",
    use_rho = isTRUE(use_rho),
    use_market_resid_scale = FALSE,
    use_hmt = FALSE,
    use_alpha_market_pool = FALSE,
    use_s0_market_year_pool = FALSE,
    use_s0_market_pool = FALSE,
    use_centered_alpha = FALSE,
    use_centered_s0_pool = FALSE,
    priors = list(),
    X_s0 = NULL,
    margin_observed = fixture$margin_observed,
    outside_price = 1,
    likelihood_mode = 1L,
    loo_cluster = fixture$market,
    estimate_df = FALSE,
    resid_df_fixed = 5,
    structure = structure,
    partition_prior = if (exact) "ordered_k" else "cutoff",
    k_support = k_support,
    k_prior = k_prior,
    k_scope = "market",
    share_likelihood = "gaussian",
    quality_parameterization = quality,
    margin_likelihood = "gaussian"
  )
  fn <- if (identical(demand, "ces")) {
    antitrustBayes:::prepare_ces_stan_data
  } else antitrustBayes:::prepare_stan_data
  do.call(fn, args)
}

cga_fixed <- function(model, data, init = 0, seed = 1001L) {
  ## init = 0 is used for the model-instance API fixture.  Target draws use
  ## an explicit constrained list derived from that same init-0 fit.
  rstan::sampling(
    object = model, data = data, chains = 1L, iter = 1L, warmup = 0L,
    algorithm = "Fixed_param", init = init, seed = as.integer(seed),
    refresh = 0
  )
}

cga_bundle <- function(model, data, demand, supply, beta = NULL,
                        rho = 0, seed = 1001L,
                        log_sigma_share = -1.5,
                        log_sigma_margin = -1.3) {
  fit0 <- cga_fixed(model, data, init = 0, seed = seed)
  init <- cga_common_init(
    fit0, demand, supply, beta = beta, rho = rho,
    log_sigma_share = log_sigma_share,
    log_sigma_margin = log_sigma_margin
  )
  fit <- cga_fixed(model, data, init = list(init), seed = seed + 1L)
  list(fit0 = fit0, fit = fit, init = init,
       upars = rstan::unconstrain_pars(fit0, init))
}

cga_log_target <- function(fit, upars) {
  as.numeric(rstan::log_prob(fit, upars, adjust_transform = FALSE))
}

cga_record_rows <- list()
cga_add <- function(check, demand = NA_character_, supply = NA_character_,
                    rho = NA_real_, component = NA_character_,
                    error = NA_real_, tolerance = NA_real_,
                    status = NULL, notes = "") {
  if (is.null(status)) {
    status <- if (is.finite(error) && error <= tolerance) "passed" else "failed"
  }
  cga_record_rows[[length(cga_record_rows) + 1L]] <<- data.frame(
    check = as.character(check), demand = as.character(demand),
    supply = as.character(supply), rho = as.numeric(rho),
    component = as.character(component), max_abs_error = as.numeric(error),
    tolerance = as.numeric(tolerance), status = as.character(status),
    notes = as.character(notes), stringsAsFactors = FALSE
  )
  invisible(NULL)
}

cga_try <- function(check, demand = NA_character_, supply = NA_character_,
                    rho = NA_real_, component = NA_character_, fun) {
  tryCatch(fun(), error = function(e) {
    cga_add(check, demand, supply, rho, component,
            error = Inf, tolerance = 0, status = "failed",
            notes = conditionMessage(e))
    NULL
  })
}

cga_rel_error <- function(ad, fd) {
  abs(fd - ad) / (1 + abs(ad))
}

## A small independent spectral evaluator used only by the tiny-sigma stress
## fixture.  The residual is constructed as H q, so the N-dimensional density
## reduces exactly to the d-dimensional eigenproblem for H'H plus the
## orthogonal residual.  This avoids asking a dense N x N Cholesky to resolve
## 1e-12 eigenvalues against O(1) quality eigenvalues.
cga_share_log_density_spectral <- function(y, mu, H, tau, sigma) {
  r <- as.numeric(y - mu)
  H <- as.matrix(H)
  n <- length(r)
  d <- ncol(H)
  if (!d) return(sum(dnorm(r, 0, sigma, log = TRUE)))
  q <- qr.solve(H, r)
  gram <- crossprod(H)
  eg <- eigen(gram, symmetric = TRUE, only.values = FALSE)
  lam <- pmax(as.numeric(eg$values), 0)
  coeff <- as.numeric(crossprod(eg$vectors, q))
  sigma2 <- sigma * sigma
  tau2 <- tau * tau
  residual_perp <- r - as.numeric(H %*% q)
  logdet <- (n - d) * log(sigma2) +
    sum(log(sigma2 + tau2 * lam))
  quadratic <- sum(residual_perp * residual_perp) / sigma2 +
    sum(lam * coeff * coeff / (sigma2 + tau2 * lam))
  -0.5 * (n * log(2 * pi) + logdet + quadratic)
}

## 1. Production AD gradients versus central finite differences.  Every
## unconstrained coordinate is checked; highlighted coordinates are retained
## in the notes/summary to make alpha/gamma, s0, tau, sigma, and rho visible.
cga_one_gradient <- function(demand, supply, rho) {
  fixture <- cga_fixture(FALSE)
  data <- cga_prepare(demand, supply, fixture, quality = "collapsed", rho = rho)
  bundle <- cga_bundle(cga_model_for(demand), data, demand, supply,
                       rho = rho, seed = 1100L + round(100 * rho))
  fit0 <- bundle$fit0
  u <- bundle$upars
  ad <- as.numeric(rstan::grad_log_prob(
    fit0, u, adjust_transform = FALSE
  ))
  h <- 1e-5
  fd <- vapply(seq_along(u), function(j) {
    up <- u; um <- u
    up[[j]] <- up[[j]] + h
    um[[j]] <- um[[j]] - h
    (cga_log_target(fit0, up) - cga_log_target(fit0, um)) / (2 * h)
  }, numeric(1))
  labels <- cga_parameter_labels(fit0)
  if (length(labels) != length(u)) {
    stop("model_pars/par_dims labels do not cover every unconstrained coordinate",
         call. = FALSE)
  }
  rel <- cga_rel_error(ad, fd)
  highlight <- grepl(
    "mu_log_a|log_kappa|sigma_b_strat|s0_logit|log_sigma_share|log_sigma_margin|rho_raw",
    labels
  )
  worst <- which.max(rel)
  highlighted <- which(highlight)
  highlighted_names <- labels[highlighted]
  notes <- paste0("highlighted=alpha/gamma,s0,tau,sigma,rho; worst=",
                  labels[[worst]], "; highlighted_max=",
                  if (length(highlighted)) {
                    labels[[highlighted[[which.max(rel[highlighted])]]]]
                  } else "none", "; highlighted_coords=",
                  if (length(highlighted_names)) {
                    paste(highlighted_names, collapse = "|")
                  } else "none")
  cga_add("production_ad_gradient", demand, supply, rho, "all_unconstrained",
          max(rel), 1e-6, notes = notes)
  ## Preserve coordinate-level failures in the CSV even when a summary row is
  ## the usual compact output.
  if (any(rel > 1e-6)) for (j in which(rel > 1e-6)) {
    cga_add("production_ad_gradient_coordinate", demand, supply, rho,
            labels[[j]], rel[[j]], 1e-6, notes = paste0("AD=", ad[[j]],
                                                        "; FD=", fd[[j]]))
  }
  invisible(TRUE)
}

for (demand in c("logit", "ces")) for (supply in c("bertrand", "cournot"))
  for (rho in c(0, 0.3)) {
    cga_try("production_ad_gradient", demand, supply, rho,
            "all_unconstrained",
            fun = function() cga_one_gradient(demand, supply, rho))
  }

## Cases below reuse these target fits so that production GQ values and direct
## log_prob values are generated at exactly the same constrained parameters.
cga_one_market <- list()

## 2 and 3. Exact normalized expanded/collapsed target equality and fully
## normalized production GQ densities against the independent dense oracle.
cga_one_density_case <- function(demand, supply, rho,
                                 sigma_share_value = exp(-1.5),
                                 check_name = "production_density_case",
                                 check_target = TRUE,
                                 centered_stress = FALSE) {
  fixture <- cga_fixture(FALSE)
  beta <- c(0.19, -0.11, 0.23)
  posterior <- NULL
  if (!is.finite(sigma_share_value) || sigma_share_value <= 0) {
    stop("sigma_share_value must be positive and finite", call. = FALSE)
  }
  log_sigma_share <- log(sigma_share_value)
  model <- cga_model_for(demand)
  system <- list(
    Z = cga_reference$one_hot(fixture$group),
    C = cga_reference$helmert_basis(4L)
  )
  H <- system$Z %*% system$C
  data_col <- cga_prepare(demand, supply, fixture, quality = "collapsed",
                          rho = rho)
  ## For stress cases, first obtain the actual production mean and then make
  ## the fixed Stan input residual exactly H q.  This is a centered math
  ## fixture, not a simulated data draw and does not invoke HMC.
  if (isTRUE(centered_stress)) {
    preliminary <- cga_bundle(model, data_col, demand, supply, rho = rho,
                               seed = 2050L + round(100 * rho),
                               log_sigma_share = log_sigma_share)
    mu0 <- cga_as_draw(preliminary$fit, "share_mean")
    q0 <- sigma_share_value * c(0.25, -0.10, 0.15)
    data_col$log_shareIn <- array(as.numeric(mu0 + H %*% q0),
                                  dim = data_col$N)
  }
  col <- cga_bundle(model, data_col, demand, supply, rho = rho,
                    seed = 2100L + round(100 * rho),
                    log_sigma_share = log_sigma_share)
  data_exp <- if (isTRUE(check_target)) {
    cga_prepare(demand, supply, fixture, quality = "sampled", rho = rho)
  } else NULL
  if (isTRUE(check_target)) {
    if (isTRUE(centered_stress)) {
      data_exp$log_shareIn <- data_col$log_shareIn
    }
    expd <- cga_bundle(model, data_exp, demand, supply, beta = beta,
                       rho = rho, seed = 2200L + round(100 * rho),
                       log_sigma_share = log_sigma_share)
    lp_difference <- cga_log_target(expd$fit0, rstan::unconstrain_pars(
      expd$fit0, expd$init
    )) - cga_log_target(col$fit0, col$upars)
  }

  share_mean <- cga_as_draw(col$fit, "share_mean")
  pred_margin <- cga_as_draw(col$fit, "pred_marginNormInv")
  margin_idx <- which(fixture$margin_observed)
  y_share <- data_col$log_shareIn
  ## Production's margin likelihood uses log(pred_marginNormInv + 1e-9).
  ## Keep the independent posterior and density references on that exact
  ## convention, including for the exact-k mixture oracle below.
  y_margin <- data_col$log_marginNormInv[margin_idx] -
    log(pred_margin[margin_idx] + 1e-9)
  mu <- share_mean
  nu <- rep(0, length(margin_idx))
  R <- matrix(0, nrow = length(y_share), ncol = length(margin_idx))
  R[cbind(margin_idx, seq_along(margin_idx))] <- 1
  tau <- 0.5
  sigma_share <- sigma_share_value
  sigma_margin <- exp(-1.3)
  if (isTRUE(check_target)) {
    posterior <- cga_reference$joint_quality_posterior(
      y_share, mu, y_margin, nu, H, tau, sigma_share, sigma_margin, rho, R
    )
    rhs <- cga_reference$log_mvnorm(beta, posterior$mean, posterior$covariance)
    cga_add("normalized_expanded_minus_collapsed_target", demand, supply, rho,
            "joint_quality_posterior", abs(lp_difference - rhs), 1e-8,
            notes = paste0("adjust_transform=FALSE; beta arbitrary; sigma_share=",
                           format(sigma_share_value, scientific = TRUE)))
  }

  ref_share <- if (isTRUE(centered_stress)) {
    cga_share_log_density_spectral(y_share, mu, H, tau, sigma_share)
  } else cga_reference$collapsed_share_log_density(
    y_share, mu, system$Z, system$C, tau, sigma_share
  )
  share_gq <- cga_as_draw(col$fit, "log_lik_share_collapsed")
  if (isTRUE(centered_stress)) {
    cga_add(check_name, demand, supply, rho, "centered_share_spectral",
            abs(share_gq - ref_share), 1e-8,
            notes = paste0("independent spectral H'H oracle; centered residual Hq; sigma_share=",
                           format(sigma_share_value, scientific = TRUE)))
  } else {
    ref_joint <- cga_reference$collapsed_joint_log_density(
      y_share, mu, y_margin, nu, system$Z, system$C, tau, sigma_share,
      sigma_margin, rho, R
    )
    joint_gq <- cga_as_draw(col$fit, "log_lik_joint_collapsed")
    cga_add(check_name, demand, supply, rho,
            "target_and_gq", max(abs(c(share_gq - ref_share,
                                       joint_gq - ref_joint))), 1e-8,
            notes = paste0("normalized constants; sigma_share=",
                           format(sigma_share_value, scientific = TRUE)))
    cga_add("production_collapsed_share_gq", demand, supply, rho,
            "normalized_share_density", abs(share_gq - ref_share), 1e-8,
            notes = paste0("normalizing constants included; sigma_share=",
                           format(sigma_share_value, scientific = TRUE)))
    cga_add("production_collapsed_joint_gq", demand, supply, rho,
            "normalized_joint_density", abs(joint_gq - ref_joint), 1e-8,
            notes = paste0("original share-margin error covariance; missing margins masked; sigma_share=",
                           format(sigma_share_value, scientific = TRUE)))
  }
  cga_add("production_prediction_bounds", demand, supply, rho,
          "predictions", if (all(pred_margin > 0) &&
                             all(fixture$prices > 0) &&
                             all(fixture$margins[margin_idx] < 1)) 0 else Inf,
          0, notes = paste0("prices positive and observed output margins below one; sigma_share=",
                            format(sigma_share_value, scientific = TRUE)))
  cga_one_market[[paste(demand, supply, rho, sep = "|")]] <<- list(
    fixture = fixture, data = data_col, bundle = col, share_mean = share_mean,
    pred_margin = pred_margin, beta = beta, posterior = posterior,
    model = model
  )
  invisible(TRUE)
}

for (demand in c("logit", "ces")) for (supply in c("bertrand", "cournot"))
  for (rho in c(0, 0.3)) {
    cga_try("production_density_case", demand, supply, rho,
            "target_and_gq",
            fun = function() cga_one_density_case(demand, supply, rho))
  }

## Small-scale cancellation stress cases for the O(N + G) production
## evaluator.  These use the same fixture and common economic parameters as
## the standard density cases while changing only sigma_share.
for (sigma_share_value in c(1e-4, 1e-6))
  for (demand in c("logit", "ces")) for (supply in c("bertrand", "cournot"))
    for (rho in c(0, 0.3)) {
      cga_try("production_stress_density_case", demand, supply, rho,
              paste0("sigma_share=", format(sigma_share_value,
                                             scientific = TRUE)),
              fun = function() cga_one_density_case(
                demand, supply, rho,
                sigma_share_value = sigma_share_value,
                check_name = "production_stress_density_case",
                check_target = FALSE,
                centered_stress = TRUE
              ))
    }

## Explicitly exercise the documented rejection boundary.  The unsupported
## case is never sent to a Stan sampler.
cga_try("exact_k_multi_market_rho_rejection", component = "prepare_helper",
        fun = function() {
          fixture <- cga_fixture(TRUE)
          got <- tryCatch({
            cga_prepare("logit", "bertrand", fixture, quality = "collapsed",
                        rho = 0.3, structure = "fringe",
                        k_support = cga_fringe_support,
                        k_prior = rep(1 / length(cga_fringe_support),
                                      length(cga_fringe_support)),
                        use_rho = TRUE)
            FALSE
          }, error = function(e) {
            grepl("one market", conditionMessage(e), fixed = TRUE)
          })
          cga_add("exact_k_multi_market_rho_rejection", component = "prepare_helper",
                  error = if (isTRUE(got)) 0 else Inf, tolerance = 0,
                  notes = "collapsed Gaussian exact-k rho is one-market only")
          stopifnot(isTRUE(got))
          invisible(TRUE)
        })

## 4. Exact-k BF/CF/BS/CS market mixtures.  A single common theta and one
## common collapsed share density are used for every market and k.  Singleton
## support fits provide exact production predictions for the independent R
## margin-score oracle; the multi-support fit supplies the production mixture.
cga_exact_case <- function(demand, supply, structure) {
  fixture <- cga_fixture(TRUE)
  ks <- if (identical(structure, "fringe")) cga_fringe_support else
    cga_stack_support
  prior <- rep(1 / length(ks), length(ks))
  model <- cga_model_for(demand, structure)
  multi_data <- cga_prepare(
    demand, supply, fixture, quality = "collapsed", rho = 0,
    structure = structure, k_support = ks, k_prior = prior, use_rho = FALSE
  )
  multi <- cga_bundle(model, multi_data, demand, supply, rho = 0,
                      seed = 3100L + match(structure, c("fringe", "stackelberg")) * 100)
  share_gq <- cga_as_draw(multi$fit, "log_lik_share_collapsed")
  share_mean <- cga_as_draw(multi$fit, "share_mean")
  pred_by_k <- vector("list", length(ks))
  singleton_share <- numeric(length(ks))
  singleton_market_scores <- matrix(NA_real_, nrow = length(ks), ncol = 2L)
  margin_idx <- which(fixture$margin_observed)
  sigma_margin <- exp(-1.3)
  for (kk in seq_along(ks)) {
    singleton_data <- cga_prepare(
      demand, supply, fixture, quality = "collapsed", rho = 0,
      structure = structure, k_support = ks[[kk]], k_prior = 1,
      use_rho = FALSE
    )
    singleton <- cga_bundle(model, singleton_data, demand, supply, rho = 0,
                            seed = 3200L + 10L * kk +
                              if (identical(structure, "stackelberg")) 100L else 0L)
    pred <- cga_as_draw(singleton$fit, "pred_marginNormInv")
    pred_by_k[[kk]] <- pred
    singleton_share[[kk]] <- cga_as_draw(singleton$fit,
                                         "log_lik_share_collapsed")
    for (m in 1:2) {
      ii <- margin_idx[fixture$market[margin_idx] == m]
      singleton_market_scores[kk, m] <- sum(dnorm(
        multi_data$log_marginNormInv[ii], mean = log(pred[ii] + 1e-9),
        sd = sigma_margin, log = TRUE
      ))
    }
  }
  expected_market <- vapply(1:2, function(m) {
    cga_lse(log(prior) + singleton_market_scores[, m])
  }, numeric(1))
  production_market <- cga_as_draw(multi$fit, "log_lik_market")
  cga_add("exact_k_market_complete_mixture", demand, supply, 0,
          structure, max(abs(production_market - expected_market)), 1e-8,
          notes = paste0("support=", paste(ks, collapse = ","),
                         "; one whole-market score per state"))

  ## The share term appears once globally; market mixtures then add once each.
  total_expected <- share_gq + sum(expected_market)
  production_total <- cga_as_draw(multi$fit, "log_lik_joint_collapsed")
  cga_add("exact_k_total_share_plus_market_mixture", demand, supply, 0,
          structure, abs(production_total - total_expected), 1e-8,
          notes = "common theta; share density counted once")
  cga_add("exact_k_share_density_common_across_states", demand, supply, 0,
          structure, max(abs(share_gq - singleton_share)), 1e-8,
          notes = "same theta/share density evaluated once across all k")
  cga_add("exact_k_prediction_bounds", demand, supply, 0, structure,
          if (all(vapply(pred_by_k, function(z) all(z > 0), logical(1)))) 0 else Inf,
          0, notes = "singleton exact economic predictions positive")
  invisible(TRUE)
}

for (demand in c("logit", "ces")) for (supply in c("bertrand", "cournot"))
  for (structure in c("fringe", "stackelberg")) {
    cga_try("exact_k_market_mixture", demand, supply, 0, structure,
            fun = function() cga_exact_case(demand, supply, structure))
  }

## Write failures as well as successes before strict assertions.  If a Stan or
## API detail regresses, the CSV remains a compact diagnostic artifact.
if (!length(cga_record_rows)) stop("No validation checks were recorded.", call. = FALSE)
cga_results <- do.call(rbind, cga_record_rows)
utils::write.csv(cga_results, cga_out_file, row.names = FALSE, quote = TRUE)

## Standalone strict contract: every bounded production check must pass.
stopifnot(all(cga_results$status == "passed"))

cat("Production collapsed-Gaussian validation passed.\n")
cat("Wrote ", normalizePath(cga_out_file, mustWork = FALSE), "\n", sep = "")
invisible(cga_results)
