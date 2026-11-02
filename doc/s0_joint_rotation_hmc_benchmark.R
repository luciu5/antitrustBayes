#!/usr/bin/env Rscript

## Equal-budget HMC comparison for market/outside-share coordinates.
##
## By default this runs a small deterministic synthetic diagnostic.  The
## Supreme Court gate must supply an RDS with columns prices, quantities,
## margins, market, year, group (and optional margin_observed) through
## S0_BENCH_DATA.  No fit objects are written; only compact CSV diagnostics.

options(stringsAsFactors = FALSE)
lib <- Sys.getenv("S0_JOINT_LIB", unset = "/tmp/atb-joint-lib")
suppressPackageStartupMessages(library(antitrustBayes, lib.loc = lib))
suppressPackageStartupMessages(library(rstan))
suppressPackageStartupMessages(library(posterior))

env_int <- function(name, default) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, unset = as.character(default))))
  if (length(value) != 1L || is.na(value) || value < 1L) stop(name, " must be positive")
  value
}

env_num <- function(name, default) {
  value <- suppressWarnings(as.numeric(Sys.getenv(name, unset = as.character(default))))
  if (length(value) != 1L || is.na(value) || !is.finite(value)) {
    stop(name, " must be a finite number")
  }
  value
}

truthy <- function(value) tolower(value) %in% c("true", "t", "1", "yes", "y")

make_synthetic <- function() {
  set.seed(20260919L)
  firms <- c(7L, 8L, 6L, 9L, 7L, 8L)
  market <- rep(seq_along(firms), firms)
  group <- unlist(lapply(firms, seq_len), use.names = FALSE)
  N <- length(market)
  prices <- 1 + .18 * sin(seq_len(N) / 3) + .025 * group
  raw_q <- exp(-.85 * prices + .18 * cos(group))
  quantities <- raw_q * rep(seq(75, 115, length.out = length(firms)), firms)
  margins <- .24 + .018 * sin(seq_len(N) / 2)
  data.frame(
    prices = prices, quantities = quantities, margins = margins,
    market = paste0("M", market), year = "Y1",
    group = paste0("G", group), margin_observed = TRUE
  )
}

make_synthetic_panel <- function() {
  base <- make_synthetic()
  periods <- c("Y1", "Y2", "Y3")
  out <- do.call(rbind, lapply(seq_along(periods), function(t) {
    data <- base
    data$year <- periods[t]
    data$prices <- data$prices * exp(c(-.035, 0, .04)[t])
    data$quantities <- data$quantities *
      exp(c(.025, 0, -.03)[t] + .012 * sin(seq_len(nrow(data)) + t))
    data$margins <- pmin(.8, pmax(.03,
      data$margins + c(-.008, 0, .011)[t]))
    data
  }))
  rownames(out) <- NULL
  out
}

data_path <- Sys.getenv("S0_BENCH_DATA", unset = "")
output_path <- Sys.getenv("S0_BENCH_OUTPUT",
                          unset = "doc/s0_joint_rotation_hmc.csv")
posterior_output_path <- Sys.getenv("S0_BENCH_POSTERIOR_OUTPUT", unset = "")
divergence_output_path <- Sys.getenv("S0_BENCH_DIVERGENCE_OUTPUT", unset = "")
divergence_draws_output_path <- Sys.getenv("S0_BENCH_DIVERGENCE_DRAWS_OUTPUT",
                                           unset = "")
if (nzchar(data_path)) {
  dat <- readRDS(data_path)
  fixture <- normalizePath(data_path, mustWork = TRUE)
} else {
  synthetic_fixture <- tolower(Sys.getenv(
    "S0_BENCH_SYNTHETIC", unset = "cross_section"
  ))
  if (!synthetic_fixture %in% c("cross_section", "panel")) {
    stop("S0_BENCH_SYNTHETIC must be cross_section or panel")
  }
  if (identical(synthetic_fixture, "panel")) {
    dat <- make_synthetic_panel()
    fixture <- "synthetic_6_market_3_period_panel"
  } else {
    dat <- make_synthetic()
    fixture <- "synthetic_6_market_cross_section"
  }
}
required <- c("prices", "quantities", "margins", "market", "year", "group")
if (!all(required %in% names(dat))) stop("benchmark data lacks required columns")
if (!"margin_observed" %in% names(dat)) dat$margin_observed <- TRUE

case_table <- data.frame(
  case = c("nopool_rotated", "market_legacy_nc", "market_rotated_nc",
           "market_rotated_centered", "market_joint_nc"),
  pooling = c("none", rep("market", 4L)),
  coordinate = c("rotated", "legacy", "rotated", "rotated", "joint"),
  centered = c(FALSE, FALSE, FALSE, TRUE, FALSE),
  stringsAsFactors = FALSE
)
case_filter <- strsplit(Sys.getenv("S0_BENCH_CASES", unset =
  paste(case_table$case, collapse = ",")), ",", fixed = TRUE)[[1L]]
case_table <- case_table[case_table$case %in% trimws(case_filter), , drop = FALSE]
if (!nrow(case_table)) stop("S0_BENCH_CASES selected no cases")
seeds <- as.integer(strsplit(Sys.getenv("S0_BENCH_SEEDS", unset = "19091,29091"),
                             ",", fixed = TRUE)[[1L]])
if (anyNA(seeds)) stop("S0_BENCH_SEEDS must be comma-separated integers")

chains <- env_int("S0_BENCH_CHAINS", 4L)
iter <- env_int("S0_BENCH_ITER", 800L)
warmup <- env_int("S0_BENCH_WARMUP", 400L)
cores <- env_int("S0_BENCH_CORES", min(chains, 2L))
adapt_delta <- env_num("S0_BENCH_ADAPT_DELTA", .95)
max_treedepth <- env_int("S0_BENCH_MAX_TREEDEPTH", 13L)
supply_model <- tolower(Sys.getenv("S0_BENCH_SUPPLY", unset = "bertrand"))
if (!supply_model %in% c("bertrand", "cournot", "monopolistic", "auction")) {
  stop("S0_BENCH_SUPPLY must be bertrand, cournot, monopolistic, or auction")
}
alpha_market_pool <- truthy(Sys.getenv("S0_BENCH_ALPHA_POOL", unset = "false"))
quality_scale_parameterization <- tolower(Sys.getenv(
  "S0_BENCH_QUALITY_SCALE", unset = "legacy"
))
if (!quality_scale_parameterization %in% c("legacy", "total_fraction")) {
  stop("S0_BENCH_QUALITY_SCALE must be legacy or total_fraction")
}
resume <- truthy(Sys.getenv("S0_BENCH_RESUME", unset = "true"))
if (warmup >= iter) stop("warmup must be less than iter")

diagnostics <- function(fit, elapsed, row, seed) {
  sp <- rstan::get_sampler_params(fit@stanfit, inc_warmup = FALSE)
  divergences <- sum(vapply(sp, function(x) sum(x[, "divergent__"]), numeric(1)))
  td_limit <- max_treedepth
  depth_hits <- sum(vapply(sp, function(x) sum(x[, "treedepth__"] >= td_limit), numeric(1)))
  bfmi <- vapply(sp, function(x) {
    e <- x[, "energy__"]
    if (length(e) < 2L || stats::var(e) <= 0) return(NA_real_)
    mean(diff(e)^2) / stats::var(e)
  }, numeric(1))
  draws <- posterior::as_draws_array(fit@stanfit)
  monitored <- grep("^(mu_log_a|log_kappa|alpha_structural_market|s0_prob|b_market|sigma_s0_pool_vec|sigma_b_market_vec|rho$|sigma_b_strat_economic$|share_variance_fraction$|sigma_share$|sigma_margin$)",
                    posterior::variables(draws), value = TRUE)
  sm <- posterior::summarise_draws(posterior::subset_draws(draws, variable = monitored))
  raw <- rstan::extract(fit@stanfit,
    pars = intersect(c("b_market", "s0_prob", "sigma_s0_pool_vec",
      "sigma_b_market_vec", "b_market_raw", "c_s0_raw", "s0_market_raw"),
      fit@stanfit@model_pars), permuted = TRUE)
  corr_economic <- NA_real_
  if (!is.null(raw$b_market) && !is.null(raw$s0_prob)) {
    eta <- qlogis(pmin(pmax(raw$s0_prob, 1e-10), 1 - 1e-10))
    ## b_market has one column per market, whereas s0_prob has one column per
    ## market-year.  Average eta over periods before flattening so this
    ## diagnostic remains dimensionally correct for panels.
    market_f <- droplevels(as.factor(dat$market))
    year_f <- droplevels(as.factor(dat$year))
    market_int <- as.integer(market_f)
    mky_int <- as.integer(as.factor(paste(market_int, as.integer(year_f), sep = "_")))
    mky_to_market <- vapply(seq_len(max(mky_int)), function(k) {
      unique_market <- unique(market_int[mky_int == k])
      if (length(unique_market) != 1L) stop("market-year maps to multiple markets")
      unique_market
    }, integer(1))
    eta_market <- vapply(seq_len(nlevels(market_f)), function(m) {
      rowMeans(eta[, mky_to_market == m, drop = FALSE])
    }, numeric(nrow(eta)))
    if (is.null(dim(eta_market))) eta_market <- matrix(eta_market, ncol = 1L)
    corr_economic <- suppressWarnings(cor(
      as.numeric(raw$b_market), as.numeric(eta_market)
    ))
  }
  corr_sampling <- NA_real_
  second_coordinate <- if (!is.null(raw$c_s0_raw)) {
    raw$c_s0_raw
  } else if (!is.null(raw$s0_market_raw)) {
    raw$s0_market_raw
  } else {
    NULL
  }
  if (!is.null(raw$b_market_raw) && !is.null(second_coordinate) &&
      length(raw$b_market_raw) == length(second_coordinate)) {
    corr_sampling <- suppressWarnings(cor(
      as.numeric(raw$b_market_raw), as.numeric(second_coordinate)
    ))
  }
  ratio_q95 <- NA_real_
  if (!is.null(raw$b_market) && !is.null(raw$sigma_s0_pool_vec)) {
    sig <- rep(as.numeric(raw$sigma_s0_pool_vec), times = ncol(raw$b_market))
    ratio_q95 <- unname(quantile(abs(as.numeric(raw$b_market) / sig), .95,
                                 na.rm = TRUE))
  }
  result <- data.frame(
    fixture = fixture, N = nrow(dat), N_market = length(unique(dat$market)),
    supply_model = supply_model, alpha_market_pool = alpha_market_pool,
    quality_scale_parameterization = quality_scale_parameterization,
    case = row$case, pooling = row$pooling, coordinate = row$coordinate,
    centered = row$centered, seed = seed, chains = chains, iter = iter,
    warmup = warmup, cores = cores, adapt_delta = adapt_delta,
    max_treedepth = max_treedepth, elapsed_seconds = elapsed,
    divergences = divergences, treedepth_hits = depth_hits,
    min_bfmi = min(bfmi, na.rm = TRUE),
    max_rhat = max(sm$rhat, na.rm = TRUE),
    min_bulk_ess = min(sm$ess_bulk, na.rm = TRUE),
    min_tail_ess = min(sm$ess_tail, na.rm = TRUE),
    mean_stepsize = mean(vapply(sp, function(x) mean(x[, "stepsize__"]), numeric(1))),
    economic_b_eta_correlation = corr_economic,
    sampling_coordinate_correlation = corr_sampling,
    b_over_sigma_s0_abs_q95 = ratio_q95,
    sigma_s0_q025 = if (is.null(raw$sigma_s0_pool_vec)) NA_real_ else
      unname(quantile(raw$sigma_s0_pool_vec, .025)),
    sigma_s0_q50 = if (is.null(raw$sigma_s0_pool_vec)) NA_real_ else
      unname(quantile(raw$sigma_s0_pool_vec, .5)),
    usable = divergences == 0L && depth_hits == 0L &&
      min(bfmi, na.rm = TRUE) >= .30 && max(sm$rhat, na.rm = TRUE) <= 1.01 &&
      min(sm$ess_bulk, na.rm = TRUE) >= 100 &&
      min(sm$ess_tail, na.rm = TRUE) >= 100,
    error = "", stringsAsFactors = FALSE
  )
  if (nzchar(posterior_output_path)) {
    posterior_monitored <- monitored[!grepl("_raw|^c_s0", monitored)]
    posterior_sm <- posterior::summarise_draws(
      posterior::subset_draws(draws, variable = posterior_monitored),
      "mean", "sd", "mcse_mean", "rhat", "ess_bulk", "ess_tail"
    )
    posterior_sm$fixture <- fixture
    posterior_sm$case <- row$case
    posterior_sm$seed <- seed
    attr(result, "posterior_summary") <- posterior_sm[, c(
      "fixture", "case", "seed", "variable", "mean", "sd", "mcse_mean",
      "rhat", "ess_bulk", "ess_tail"
    )]
  }
  if (nzchar(divergence_output_path) || nzchar(divergence_draws_output_path)) {
    pars <- intersect(c("sigma_s0_pool_vec", "sigma_b_market_vec",
                        "mu_s0_pool_vec", "mu_b_market_vec",
                        "sigma_b_strat_economic", "sigma_b_strat",
                        "log_sigma_share", "share_variance_fraction",
                        "sigma_share", "sigma_margin",
                        "rho", "rho_raw", "alpha_structural_market",
                        "s0_prob", "market_s0_joint_scale", "b_market",
                        "b_market_raw", "c_s0_raw", "s0_market_raw"),
                      fit@stanfit@model_pars)
    native <- rstan::extract(fit@stanfit, pars = pars, permuted = FALSE,
                             inc_warmup = FALSE)
    native_names <- dimnames(native)[[3L]]
    take <- function(pattern) {
      idx <- grep(pattern, native_names)
      if (!length(idx)) return(NULL)
      native[, , idx, drop = FALSE]
    }
    sigma_arr <- take("^sigma_s0_pool_vec\\[1\\]$")
    sigma_b_arr <- take("^sigma_b_market_vec\\[1\\]$")
    mu_s_arr <- take("^mu_s0_pool_vec\\[1\\]$")
    mu_b_arr <- take("^mu_b_market_vec\\[1\\]$")
    sigma_quality_arr <- take("^sigma_b_strat_economic$")
    if (is.null(sigma_quality_arr)) sigma_quality_arr <- take("^sigma_b_strat$")
    contrast_total_arr <- if (identical(quality_scale_parameterization,
                                        "total_fraction")) {
      take("^sigma_b_strat$")
    } else NULL
    variance_fraction_arr <- take("^share_variance_fraction$")
    sigma_share_arr <- take("^sigma_share$")
    sigma_margin_arr <- take("^sigma_margin$")
    rho_arr <- take("^rho$")
    rho_raw_arr <- take("^rho_raw\\[1\\]$")
    alpha_arr <- take("^alpha_structural_market\\[")
    s0_arr <- take("^s0_prob\\[")
    joint_scale_arr <- take("^market_s0_joint_scale$")
    b_arr <- take("^b_market\\[")
    b_raw_arr <- take("^b_market_raw\\[")
    c_arr <- take("^c_s0_raw\\[")
    z_arr <- take("^s0_market_raw\\[")
    collapse_scalar <- function(x) if (is.null(x)) NULL else as.numeric(x[, , 1L])
    collapse_max <- function(x) {
      if (is.null(x)) return(NULL)
      apply(abs(x), c(1L, 2L), max)
    }
    sigma_draw <- collapse_scalar(sigma_arr)
    b_over_sigma <- if (is.null(sigma_draw) || is.null(b_arr)) NULL else
      as.numeric(collapse_max(b_arr)) / sigma_draw
    reconstructed_z <- NULL
    if (!is.null(sigma_draw) && !is.null(b_arr) && !is.null(c_arr) &&
        identical(row$coordinate, "rotated") && !isTRUE(row$centered)) {
      z_values <- c_arr
      for (m in seq_len(dim(b_arr)[3L])) {
        z_values[, , m] <- c_arr[, , m] - b_arr[, , m] /
          matrix(sigma_draw, nrow = dim(b_arr)[1L], ncol = dim(b_arr)[2L])
      }
      reconstructed_z <- as.numeric(collapse_max(z_values))
    } else if (!is.null(z_arr)) {
      reconstructed_z <- as.numeric(collapse_max(z_arr))
    }
    divergent <- as.logical(unlist(lapply(sp, function(x) x[, "divergent__"]),
                                   use.names = FALSE))
    sigma_quality_draw <- collapse_scalar(sigma_quality_arr)
    sigma_share_draw <- collapse_scalar(sigma_share_arr)
    rho_draw <- collapse_scalar(rho_arr)
    quality_to_residual <- if (is.null(sigma_quality_draw) ||
                               is.null(sigma_share_draw)) NULL else
      sigma_quality_draw / sigma_share_draw
    quality_to_conditional_residual <- if (is.null(quality_to_residual) ||
                                           is.null(rho_draw)) NULL else
      quality_to_residual / sqrt(1 - rho_draw^2)
    metrics <- list(
      sigma_s0_pool = sigma_draw,
      sigma_b_market = collapse_scalar(sigma_b_arr),
      mu_s0_pool = collapse_scalar(mu_s_arr),
      mu_b_market = collapse_scalar(mu_b_arr),
      sigma_quality = sigma_quality_draw,
      contrast_total_scale = collapse_scalar(contrast_total_arr),
      share_variance_fraction = collapse_scalar(variance_fraction_arr),
      sigma_share = sigma_share_draw,
      sigma_margin = collapse_scalar(sigma_margin_arr),
      rho = rho_draw,
      rho_raw = collapse_scalar(rho_raw_arr),
      quality_to_residual = quality_to_residual,
      quality_to_conditional_residual = quality_to_conditional_residual,
      alpha_structural_max = if (is.null(alpha_arr)) NULL else
        as.numeric(apply(alpha_arr, c(1L, 2L), max)),
      s0_min = if (is.null(s0_arr)) NULL else
        as.numeric(apply(s0_arr, c(1L, 2L), min)),
      s0_max = if (is.null(s0_arr)) NULL else
        as.numeric(apply(s0_arr, c(1L, 2L), max)),
      joint_scale = collapse_scalar(joint_scale_arr),
      max_abs_b_market = if (is.null(b_arr)) NULL else as.numeric(collapse_max(b_arr)),
      max_abs_joint_v = if (is.null(b_raw_arr)) NULL else
        as.numeric(collapse_max(b_raw_arr)),
      max_abs_joint_h = if (is.null(c_arr)) NULL else
        as.numeric(collapse_max(c_arr)),
      max_abs_b_over_sigma_s0 = b_over_sigma,
      max_abs_reconstructed_z_s0 = reconstructed_z
    )
    metrics <- metrics[!vapply(metrics, is.null, logical(1))]
    draw_metrics <- data.frame(
      fixture = fixture, case = row$case, seed = seed,
      draw = seq_along(divergent), divergent = divergent,
      as.data.frame(metrics, check.names = FALSE), stringsAsFactors = FALSE
    )
    summarize_metric <- function(value, name, flag) {
      selected <- value[divergent == flag]
      data.frame(
        fixture = fixture, case = row$case, seed = seed,
        divergent = flag, n_draws = length(selected), metric = name,
        mean = mean(selected), median = median(selected),
        q05 = unname(quantile(selected, .05)),
        q95 = unname(quantile(selected, .95)), stringsAsFactors = FALSE
      )
    }
    divergence_rows <- list()
    for (name in names(metrics)) if (!is.null(metrics[[name]])) {
      for (flag in c(FALSE, TRUE)) if (any(divergent == flag)) {
        divergence_rows[[length(divergence_rows) + 1L]] <-
          summarize_metric(metrics[[name]], name, flag)
      }
    }
    attr(result, "divergence_summary") <- do.call(rbind, divergence_rows)
    attr(result, "divergence_draws") <- draw_metrics
  }
  result
}

`%||%` <- function(x, y) if (is.null(x)) y else x
existing <- NULL
if (resume && file.exists(output_path)) {
  existing <- tryCatch(read.csv(output_path, stringsAsFactors = FALSE),
                       error = function(e) NULL)
  if (!is.null(existing) && nrow(existing) &&
      !"quality_scale_parameterization" %in% names(existing)) {
    existing$quality_scale_parameterization <- "legacy"
  }
  if (!is.null(existing) && nrow(existing) && "error" %in% names(existing)) {
    existing$error[is.na(existing$error)] <- ""
  }
}
rows <- if (is.null(existing) || !nrow(existing)) list() else
  split(existing, seq_len(nrow(existing)))
resume_rows <- function(path) {
  if (!resume || !nzchar(path) || !file.exists(path)) return(list())
  value <- tryCatch(read.csv(path, stringsAsFactors = FALSE),
                    error = function(e) NULL)
  if (is.null(value) || !nrow(value)) list() else
    split(value, seq_len(nrow(value)))
}
posterior_rows <- resume_rows(posterior_output_path)
divergence_rows <- resume_rows(divergence_output_path)
divergence_draw_rows <- resume_rows(divergence_draws_output_path)
for (i in seq_len(nrow(case_table))) for (seed in seeds) {
  spec <- case_table[i, ]
  if (!is.null(existing) && nrow(existing) &&
      all(c("fixture", "case", "seed", "supply_model",
            "alpha_market_pool", "iter", "warmup") %in% names(existing))) {
    complete <- existing$fixture == fixture & existing$case == spec$case &
      existing$seed == seed & existing$supply_model == supply_model &
      as.logical(existing$alpha_market_pool) == alpha_market_pool &
      existing$quality_scale_parameterization ==
        quality_scale_parameterization &
      existing$iter == iter & existing$warmup == warmup &
      is.finite(existing$divergences)
    if (any(complete, na.rm = TRUE)) {
      cat("skipping completed", spec$case, "seed", seed, "\n")
      next
    }
  }
  cat("running", spec$case, "seed", seed, "\n")
  started <- proc.time()[["elapsed"]]
  result <- tryCatch({
    fit <- bayes_est(
      prices = dat$prices, quantities = dat$quantities,
      margins = dat$margins, market = dat$market, year = dat$year,
      group = dat$group, margin_observed = dat$margin_observed,
      demand_model = "logit", supply_model = supply_model,
      structure = "simultaneous", s0_pooling = spec$pooling,
      market_s0_parameterization = spec$coordinate,
      use_centered_s0_pool = spec$centered,
      use_alpha_market_pool = alpha_market_pool, use_rho = TRUE,
      share_likelihood = "gaussian", margin_likelihood = "gaussian",
      quality_parameterization = "collapsed", quality_scope = "market",
      quality_scale_parameterization = quality_scale_parameterization,
      stan_options = list(chains = chains, cores = cores, iter = iter,
        warmup = warmup, refresh = 0,
        seed = seed,
        control = list(adapt_delta = adapt_delta,
                       max_treedepth = max_treedepth)),
      prior_manifest = "baseline"
    )
    diagnostics(fit, proc.time()[["elapsed"]] - started, spec, seed)
  }, error = function(e) data.frame(
    fixture = fixture, N = nrow(dat), N_market = length(unique(dat$market)),
    supply_model = supply_model, alpha_market_pool = alpha_market_pool,
    quality_scale_parameterization = quality_scale_parameterization,
    case = spec$case, pooling = spec$pooling, coordinate = spec$coordinate,
    centered = spec$centered, seed = seed, chains = chains, iter = iter,
    warmup = warmup, cores = cores, adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    elapsed_seconds = proc.time()[["elapsed"]] - started,
    divergences = NA, treedepth_hits = NA, min_bfmi = NA, max_rhat = NA,
    min_bulk_ess = NA, min_tail_ess = NA, mean_stepsize = NA,
    economic_b_eta_correlation = NA, sampling_coordinate_correlation = NA,
    b_over_sigma_s0_abs_q95 = NA,
    sigma_s0_q025 = NA, sigma_s0_q50 = NA, usable = FALSE,
    error = conditionMessage(e), stringsAsFactors = FALSE
  ))
  posterior_result <- attr(result, "posterior_summary")
  if (!is.null(posterior_result)) {
    posterior_rows[[length(posterior_rows) + 1L]] <- posterior_result
    write.csv(do.call(rbind, posterior_rows), posterior_output_path,
              row.names = FALSE)
  }
  divergence_result <- attr(result, "divergence_summary")
  if (!is.null(divergence_result)) {
    divergence_rows[[length(divergence_rows) + 1L]] <- divergence_result
    write.csv(do.call(rbind, divergence_rows), divergence_output_path,
              row.names = FALSE)
  }
  divergence_draw_result <- attr(result, "divergence_draws")
  if (!is.null(divergence_draw_result) &&
      nzchar(divergence_draws_output_path)) {
    divergence_draw_rows[[length(divergence_draw_rows) + 1L]] <-
      divergence_draw_result
    write.csv(do.call(rbind, divergence_draw_rows),
              divergence_draws_output_path, row.names = FALSE)
  }
  rows[[length(rows) + 1L]] <- result
  write.csv(do.call(rbind, rows), output_path, row.names = FALSE)
}

out <- do.call(rbind, rows)
print(out)
if (any(nzchar(out$error))) quit(status = 2L)
