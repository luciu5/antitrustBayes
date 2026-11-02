#!/usr/bin/env Rscript

## Direct production scaling harness for the collapsed Gaussian quality block.
##
## This is deliberately a small, transparent performance harness.  It loads
## the installed package's precompiled Logit simultaneous-supply model, uses
## prepare_stan_data() to create the exact production data contract, builds one
## one-draw Fixed_param stanfit per data size/layout, and then times only
## log_prob() and grad_log_prob().  Initialization and generated quantities are
## outside the timed regions.  No HMC fit is attempted here.
##
## The quality block is always quality_parameterization = "collapsed".  The
## production Stan function therefore keeps only O(N + G) group statistics;
## this script does not construct an N x N covariance, a G x G contrast matrix,
## or a dense Q.  A dense small-N oracle is intentionally not part of this
## harness; independent equality gates are maintained separately.

options(stringsAsFactors = FALSE)

cga_script_args <- commandArgs(trailingOnly = TRUE)
cga_script_file <- commandArgs(trailingOnly = FALSE)
cga_script_file <- cga_script_file[startsWith(cga_script_file, "--file=")]
cga_doc_dir <- if (length(cga_script_file)) {
  dirname(normalizePath(sub("^--file=", "", cga_script_file[[1L]]),
                         mustWork = FALSE))
} else {
  file.path(getwd(), "doc")
}
cga_output_arg <- cga_script_args[startsWith(cga_script_args, "--output=")]
cga_output_file <- if (length(cga_output_arg)) {
  sub("^--output=", "", cga_output_arg[[1L]])
} else {
  Sys.getenv("CGA_OUTPUT_FILE", unset =
               file.path(cga_doc_dir, "collapsed_gaussian_production_scaling.csv"))
}

if (!requireNamespace("rstan", quietly = TRUE)) {
  stop("The production scaling harness requires the installed rstan package.",
       call. = FALSE)
}
if (!requireNamespace("antitrustBayes", quietly = TRUE)) {
  stop("The production scaling harness requires the installed antitrustBayes package.",
       call. = FALSE)
}
library(rstan)

## RStan's default may write a chain CSV if a caller supplies a default file
## option.  The harness keeps the fit in memory and emits only the compact
## timing table requested by the production run.
options(mc.cores = 1L)

cga_env_integer <- function(name, default, lower = 1L) {
  raw <- Sys.getenv(name, unset = "")
  if (!nzchar(raw)) return(as.integer(default))
  val <- suppressWarnings(as.integer(raw))
  if (is.na(val) || val < lower) {
    stop(name, " must be an integer >= ", lower, call. = FALSE)
  }
  val
}

cga_repetitions_small <- cga_env_integer("CGA_REPETITIONS_SMALL", 7L)
cga_repetitions_large <- cga_env_integer("CGA_REPETITIONS_LARGE", 5L)

cga_lse <- function(x) {
  x <- as.numeric(x)
  if (!length(x)) return(-Inf)
  m <- max(x)
  m + log(sum(exp(x - m)))
}

cga_group_sum <- function(x, group, G) {
  out <- numeric(G)
  for (i in seq_along(x)) out[group[[i]]] <- out[group[[i]]] + x[[i]]
  out
}

## The production transformed-data Helmert Q is never needed to make the
## collapsed fixture.  A centered group residual vector is enough and avoids
## allocating any dense group design in R as well.
cga_centered_group_residual <- function(G) {
  g <- seq_len(G)
  q <- 0.035 * sin(2 * pi * g / G) + 0.015 * cos(4 * pi * g / G)
  q - mean(q)
}

cga_fixture <- function(N, layout, log_a = log(0.8)) {
  N <- as.integer(N)
  if (N < 1L || N %% 4L != 0L) {
    stop("The requested production sizes must be positive multiples of four.",
         call. = FALSE)
  }
  if (identical(layout, "one_product")) {
    G <- N
    group <- seq_len(N)
    market <- rep.int(1L, N)
    year <- rep.int(1L, N)
    product <- seq_len(N)
    layout_note <- "one product per persistent group; G=N; one market-year"
  } else if (identical(layout, "multiproduct_panel")) {
    G <- N %/% 4L
    ## Each persistent group owns two product rows in each of two years.
    ## `product` is retained as fixture documentation; production ownership
    ## is represented by the repeated `group` labels passed to the helper.
    group <- rep(seq_len(G), each = 4L)
    product <- rep(rep(1:2, each = 2L), G)
    year <- rep(rep(1:2, each = 2L), G)
    market <- rep.int(1L, N)
    layout_note <- paste0(
      "two products x two periods per persistent group; G=N/4=", G,
      "; one market and two market-years"
    )
  } else {
    stop("Unknown production scaling layout: ", layout, call. = FALSE)
  }

  ## The price path has a fixed, non-degenerate standard deviation at every
  ## requested N.  Its phase is deterministic, so all timed calls use the
  ## same fixed data and fixed initialization policy.
  row <- seq_len(N)
  price <- 1.0 + 0.25 * sin(2 * pi * row / 19) +
    0.10 * cos(2 * pi * row / 7)
  if (any(price <= 0)) stop("fixture generated a non-positive price")
  price_sd <- stats::sd(price)
  if (!is.finite(price_sd) || price_sd <= 0) stop("fixture price sd is invalid")

  q_group <- cga_centered_group_residual(G)
  mky_key <- paste(market, year, sep = "_")
  mky_levels <- unique(mky_key)
  mky <- match(mky_key, mky_levels)
  mky_n <- tabulate(mky, nbins = length(mky_levels))
  s0_by_mky <- 1 / (mky_n + 1)

  ## Match the package's output-market sign and the fixed common alpha init.
  ## The observed inside shares are the normalized exp(mu + centered group
  ## residual), so share residuals stay on the intended quality scale even at
  ## N=10,000.  This is a compositional fixture, not a post-fit adjustment.
  bounded_s0 <- 1e-6 + (1 - 2e-6) * s0_by_mky
  mu <- log(bounded_s0[mky] / (1 - bounded_s0[mky])) -
    exp(log_a) / price_sd * ((price - mean(price)) / price_sd)
  latent_log_share <- mu + q_group[group]
  log_share <- latent_log_share - vapply(
    seq_along(mky_levels),
    function(mm) cga_lse(latent_log_share[mky == mm]),
    numeric(1)
  )[mky]
  quantities <- exp(log_share)
  if (any(!is.finite(quantities)) || any(quantities <= 0)) {
    stop("fixture generated invalid quantities")
  }

  ## Construct margins from the same fixed Bertrand output-market equations.
  ## This keeps the joint likelihood near its intended target and prevents a
  ## large, irrelevant margin residual from dominating the gradient timing.
  inside_by_group <- numeric(N)
  for (mm in seq_along(mky_n)) {
    ii <- which(mky == mm)
    inside_by_group[ii] <- cga_group_sum(quantities[ii], group[ii], G)[group[ii]]
  }
  alpha <- exp(log_a) / price_sd
  inv_markup <- price * alpha *
    (1 - inside_by_group * bounded_s0[mky])
  if (any(!is.finite(inv_markup)) || any(inv_markup <= 1)) {
    stop("fixture generated invalid output-market inverse markups")
  }
  margins <- 1 / inv_markup
  if (any(margins <= 0) || any(margins >= 1)) {
    stop("fixture generated output margins outside (0, 1)")
  }

  list(
    N = N,
    G = G,
    layout = layout,
    layout_note = layout_note,
    product = product,
    prices = as.numeric(price),
    quantities = as.numeric(quantities),
    margins = as.numeric(margins),
    market = as.integer(market),
    year = as.integer(year),
    group = as.integer(group),
    margin_observed = rep(TRUE, N),
    mky = as.integer(mky),
    mky_n = as.integer(mky_n),
    q_group = q_group,
    target_log_share = as.numeric(log_share),
    fixed_log_a = as.numeric(log_a),
    fixed_s0_logit = qlogis(s0_by_mky)
  )
}

cga_prepare <- function(fixture) {
  ## Every hierarchy switch is explicit in this call.  The no-rho path is
  ## intentional: it exercises the production collapsed Gaussian share block
  ## and the independent Gaussian margin equation without adding a second
  ## dense covariance route.
  data <- antitrustBayes:::prepare_stan_data(
    prices = fixture$prices,
    quantities = fixture$quantities,
    margins = fixture$margins,
    market = fixture$market,
    year = fixture$year,
    group = fixture$group,
    supply_model = 1L,
    use_cutoff = FALSE,
    market_side = "output",
    use_rho = FALSE,
    use_market_resid_scale = FALSE,
    use_hmt = FALSE,
    use_alpha_market_pool = FALSE,
    use_s0_market_year_pool = FALSE,
    use_s0_market_pool = FALSE,
    use_centered_alpha = FALSE,
    use_centered_s0_pool = FALSE,
    priors = list(),
    X_s0 = NULL,
    grainsize = 1L,
    margin_observed = fixture$margin_observed,
    outside_price = 1,
    likelihood_mode = 1L,
    loo_cluster = fixture$market,
    estimate_df = FALSE,
    resid_df_fixed = 5,
    structure = "simultaneous",
    partition_prior = "cutoff",
    k_support = NULL,
    k_prior = NULL,
    k_scope = "market",
    share_likelihood = "gaussian",
    quality_parameterization = "collapsed",
    margin_likelihood = "gaussian"
  )
  alignment_error <- max(abs(as.numeric(data$log_shareIn) -
                             fixture$target_log_share))
  if (!is.finite(alignment_error) || alignment_error > 1e-8) {
    stop("prepared log_shareIn is not aligned with the fixed fixture target")
  }
  stopifnot(
    identical(as.integer(data$quality_parameterization), 2L),
    identical(as.integer(data$use_rho), 0L),
    identical(as.integer(data$use_market_resid_scale), 0L),
    identical(as.integer(data$use_alpha_market_pool), 0L),
    identical(as.integer(data$use_s0_market_year_pool), 0L),
    identical(as.integer(data$use_s0_market_pool), 0L),
    identical(as.integer(data$use_hmt), 0L),
    identical(as.integer(data$likelihood_mode), 1L),
    identical(as.integer(data$share_likelihood), 2L),
    identical(as.integer(data$margin_likelihood), 2L)
  )
  attr(data, "cga_share_alignment_error") <- alignment_error
  data
}

cga_common_init <- function(data) {
  ## These values are reused across layouts and sizes.  Only the number of
  ## explicitly unpooled s0 coordinates follows the number of market-years.
  ## The target values are deterministic and are not sampled from a fit.
  list(
    ## RStan requires zero-length declarations to be present in an explicit
    ## initialization list.  They are data-dependent zero-dimensional Stan
    ## vectors, so they add no parameter coordinates or work to the timing.
    mu_b_market_vec = numeric(0),
    sigma_b_market_vec = numeric(0),
    b_market_raw = numeric(0),
    b_group_centered = numeric(0),
    mu_b_fringe_vec = numeric(0),
    mu_log_a_free = array(log(0.8), dim = 1L),
    sigma_b_strat = 0.5,
    sigma_log_a_vec = numeric(0),
    r_market_a_raw = numeric(0),
    log_a_market_pool = numeric(0),
    log_kappa_vec = numeric(0),
    log_kappa_market = numeric(0),
    s0_logit = array(qlogis(1 / (as.numeric(data$mky_n) + 1)),
                     dim = length(data$mky_n)),
    mu_s0_pool_vec = numeric(0),
    sigma_s0_pool_vec = numeric(0),
    s0_pool_raw = numeric(0),
    eta_s0_pool = numeric(0),
    s0_market_raw = numeric(0),
    eta_s0_market = numeric(0),
    cutoff_share = numeric(0),
    rho_raw = numeric(0),
    log_sigma_share = log(0.2),
    log_sigma_margin = log(0.2),
    log_sigma_share_market_raw = numeric(0),
    log_sigma_margin_market_raw = numeric(0),
    sigma_log_sigma_share_market = numeric(0),
    sigma_log_sigma_margin_market = numeric(0),
    beta_s0 = numeric(0),
    bargaining_unit = numeric(0),
    nu_minus_2 = numeric(0)
  )
}

cga_fixed_param_instance <- function(model, data, init, seed) {
  ## This is the only sampling call.  It creates the stanfit wrapper needed by
  ## rstan::log_prob()/grad_log_prob(); no HMC iterations are run.
  rstan::sampling(
    object = model,
    data = data,
    chains = 1L,
    iter = 1L,
    warmup = 0L,
    algorithm = "Fixed_param",
    init = list(init),
    seed = as.integer(seed),
    refresh = 0L,
    cores = 1L
  )
}

cga_timed_repeated <- function(fun, repetitions, target_seconds = 0.02) {
  ## Warm one direct call outside the recorded timing.  This removes one-time
  ## module-call overhead while retaining the actual production evaluator.
  warm <- tryCatch({
    value <- fun()
    list(ok = TRUE, value = value, message = "")
  }, error = function(e) list(ok = FALSE, value = NULL,
                              message = conditionMessage(e)))
  if (!warm$ok) {
    return(list(median_seconds = NA_real_, repetitions = repetitions,
                inner_iterations = 0L,
                successful_calls = 0L, failure_count = repetitions,
                failure_messages = paste0("warmup: ", warm$message),
                status = "failed"))
  }

  ## A single Stan call at N=100 can be shorter than the elapsed-time clock's
  ## printed precision.  Calibrate a fixed batch once and report seconds per
  ## direct evaluator call; each reported repetition still repeats the same
  ## direct call on the same model instance and upars.
  inner_iterations <- 1L
  repeat {
    probe <- tryCatch({
      elapsed <- system.time({
        for (j in seq_len(inner_iterations)) invisible(fun())
      })[["elapsed"]]
      list(ok = TRUE, elapsed = as.numeric(elapsed), message = "")
    }, error = function(e) list(ok = FALSE, elapsed = NA_real_,
                                message = conditionMessage(e)))
    if (!probe$ok) {
      return(list(median_seconds = NA_real_, repetitions = repetitions,
                  inner_iterations = inner_iterations,
                  successful_calls = 0L, failure_count = repetitions,
                  failure_messages = paste0("calibration: ", probe$message),
                  status = "failed"))
    }
    if (is.finite(probe$elapsed) && probe$elapsed >= target_seconds) break
    if (inner_iterations >= 100000L) break
    inner_iterations <- min(100000L, inner_iterations * 2L)
  }
  elapsed <- numeric(repetitions)
  failures <- character()
  for (rep_id in seq_len(repetitions)) {
    one <- tryCatch({
      value <- NULL
      measured <- system.time({
        for (j in seq_len(inner_iterations)) value <- fun()
      })[["elapsed"]] / inner_iterations
      if (!is.numeric(value) || any(!is.finite(value))) {
        stop("non-finite evaluator output")
      }
      list(ok = TRUE, elapsed = as.numeric(measured), message = "")
    }, error = function(e) list(ok = FALSE, elapsed = NA_real_,
                                message = conditionMessage(e)))
    if (one$ok && is.finite(one$elapsed)) {
      elapsed[[rep_id]] <- one$elapsed
    } else {
      elapsed[[rep_id]] <- NA_real_
      failures <- c(failures, paste0("rep", rep_id, ": ", one$message))
    }
  }
  good <- is.finite(elapsed)
  list(
    median_seconds = if (any(good)) median(elapsed[good]) else NA_real_,
    repetitions = as.integer(repetitions),
    inner_iterations = as.integer(inner_iterations),
    successful_calls = as.integer(sum(good)),
    failure_count = as.integer(sum(!good)),
    failure_messages = if (length(failures)) paste(failures, collapse = " | ") else "",
    status = if (all(good)) "ok" else "failed"
  )
}

cga_parameter_dimension <- function(fit, init) {
  ## rstan's unconstrain_pars() accepts model-block constrained values.  Keep
  ## only explicitly present nonzero parameters; collapsed quality has no
  ## b_group_centered vector, and the all-hierarchy-false fixture has no pool
  ## or market-effect vectors.
  model_end <- match("sigma_share", fit@model_pars)
  if (is.na(model_end)) stop("could not locate parameter boundary")
  model_names <- fit@model_pars[seq_len(model_end - 1L)]
  model_init <- init[intersect(names(init), model_names)]
  upars <- rstan::unconstrain_pars(fit, model_init)
  list(upars = as.numeric(upars), dimension = length(upars))
}

cga_failure_row <- function(fixture, operation, repetitions, message) {
  data.frame(
    N = fixture$N,
    layout = fixture$layout,
    G = fixture$G,
    N_market_year = length(fixture$mky_n),
    operation = operation,
    median_seconds = NA_real_,
    repetitions = as.integer(repetitions),
    inner_iterations = 0L,
    successful_calls = 0L,
    failure_count = as.integer(repetitions),
    parameter_dimension = NA_integer_,
    share_alignment_error = NA_real_,
    status = "failed",
    failure_messages = as.character(message),
    dense_oracle = "prohibited_and_not_run",
    notes = fixture$layout_note,
    stringsAsFactors = FALSE
  )
}

cga_run_case <- function(model, N, layout, seed) {
  fixture <- cga_fixture(N, layout)
  repetitions <- if (N >= 10000L) cga_repetitions_large else cga_repetitions_small
  setup <- tryCatch({
    data <- cga_prepare(fixture)
    init <- cga_common_init(data)
    fit <- cga_fixed_param_instance(model, data, init, seed)
    pars <- cga_parameter_dimension(fit, init)
    if (any(!is.finite(pars$upars))) stop("unconstrained initialization is non-finite")
    list(ok = TRUE, data = data, fit = fit, upars = pars$upars,
         parameter_dimension = pars$dimension,
         share_alignment_error = attr(data, "cga_share_alignment_error"),
         message = "")
  }, error = function(e) list(ok = FALSE, message = conditionMessage(e)))
  if (isTRUE(setup$ok) == FALSE) {
    return(rbind(
      cga_failure_row(fixture, "log_prob", repetitions, setup$message),
      cga_failure_row(fixture, "grad_log_prob", repetitions, setup$message)
    ))
  }

  ## The two closures contain no data construction or generated-quantity work.
  log_fun <- function() rstan::log_prob(
    setup$fit, setup$upars, adjust_transform = FALSE
  )
  grad_fun <- function() rstan::grad_log_prob(
    setup$fit, setup$upars, adjust_transform = FALSE
  )
  timings <- list(
    log_prob = cga_timed_repeated(log_fun, repetitions),
    grad_log_prob = cga_timed_repeated(grad_fun, repetitions)
  )
  do.call(rbind, lapply(names(timings), function(operation) {
    one <- timings[[operation]]
    data.frame(
      N = fixture$N,
      layout = fixture$layout,
      G = fixture$G,
      N_market_year = length(fixture$mky_n),
      operation = operation,
      median_seconds = one$median_seconds,
      repetitions = one$repetitions,
      inner_iterations = one$inner_iterations,
      successful_calls = one$successful_calls,
      failure_count = one$failure_count,
      parameter_dimension = as.integer(setup$parameter_dimension),
      share_alignment_error = as.numeric(setup$share_alignment_error),
      status = one$status,
      failure_messages = one$failure_messages,
      dense_oracle = "prohibited_and_not_run",
      notes = paste0(
        fixture$layout_note,
        "; collapsed quality; all hierarchy pools FALSE; use_rho FALSE; HMT FALSE;",
        " likelihood joint; Fixed_param setup excluded from timings;",
        " adjust_transform FALSE; no HMC runtime claim"
      ),
      stringsAsFactors = FALSE
    )
  }))
}

cga_main <- function() {
  model <- antitrustBayes:::stanmodels[["logit_supply"]]
  sizes <- c(100L, 1000L, 10000L)
  layouts <- c("one_product", "multiproduct_panel")
  rows <- list()
  row_id <- 0L
  for (layout_id in seq_along(layouts)) {
    for (size_id in seq_along(sizes)) {
      row_id <- row_id + 1L
      N <- sizes[[size_id]]
      layout <- layouts[[layout_id]]
      seed <- 271828L + 100L * layout_id + size_id
      cat("Preparing ", layout, " N=", N, " (Fixed_param setup once)...\n",
          sep = "")
      rows[[row_id]] <- cga_run_case(model, N, layout, seed)
      flush.console()
    }
  }
  results <- do.call(rbind, rows)
  output_dir <- dirname(cga_output_file)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE,
                                          showWarnings = FALSE)
  utils::write.csv(results, cga_output_file, row.names = FALSE, quote = TRUE)
  cat("\nWrote production scaling results to ",
      normalizePath(cga_output_file, mustWork = FALSE), "\n", sep = "")
  print(results[, c("N", "layout", "G", "operation", "median_seconds",
                    "repetitions", "parameter_dimension", "status")],
        row.names = FALSE)
  if (any(results$status != "ok")) {
    cat("\nOne or more direct production evaluations failed; see the CSV.\n")
  } else {
    cat("\nAll requested direct production evaluations completed.\n")
  }
  invisible(results)
}

if (sys.nframe() == 0L) cga_main()
