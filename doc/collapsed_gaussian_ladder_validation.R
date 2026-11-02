#!/usr/bin/env Rscript

## Fixed-parameter validation for the information-aware outside-share ladder.
## This harness intentionally exercises the installed package only.  It does
## not compile Stan, run HMC, or source the working-tree R files.

args <- commandArgs(trailingOnly = TRUE)
repo_root <- if (length(args)) normalizePath(args[[1L]], mustWork = TRUE) else getwd()
lib_root <- Sys.getenv("CGALIB", unset = "/tmp/cga-Rlib")
csv_path <- file.path(repo_root, "doc", "collapsed_gaussian_ladder_validation.csv")

if (!dir.exists(lib_root)) stop("Installed package library is missing: ", lib_root)
suppressPackageStartupMessages(library(antitrustBayes, lib.loc = lib_root))
suppressPackageStartupMessages(library(rstan))

results <- list()
record <- function(model, cell, check, status, value = NA_character_, message = "") {
  results[[length(results) + 1L]] <<- data.frame(
    model = as.character(model), cell = as.character(cell),
    check = as.character(check), status = as.character(status),
    value = as.character(value), message = as.character(message),
    stringsAsFactors = FALSE
  )
}

assert_true <- function(ok, model, cell, check, value = NA_character_) {
  if (isTRUE(ok)) {
    record(model, cell, check, "pass", value)
    invisible(TRUE)
  } else {
    stop(sprintf("%s [%s/%s]", check, model, cell))
  }
}

assert_equal <- function(actual, expected, model, cell, check) {
  ok <- isTRUE(all.equal(actual, expected, check.attributes = FALSE))
  assert_true(ok, model, cell, check,
              paste(capture.output(dput(actual)), collapse = ""))
}

run_case <- function(model, cell, expr) {
  tryCatch(
    force(expr),
    error = function(e) {
      record(model, cell, "fit_or_assertions", "fail", message = conditionMessage(e))
      NULL
    }
  )
}

fixture <- function(n_markets, n_years) {
  grid <- expand.grid(
    market_id = seq_len(n_markets), year_id = seq_len(n_years),
    group_id = seq_len(4L), product_id = seq_len(2L)
  )
  prices <- 1 + 0.05 * grid$market_id + 0.03 * grid$year_id +
    0.01 * grid$group_id + 0.002 * grid$product_id
  quantities <- exp(0.05 * grid$group_id + 0.02 * grid$product_id)
  list(
    prices = as.numeric(prices), quantities = as.numeric(quantities),
    margins = as.numeric(rep(0.30, nrow(grid))),
    market = paste0("M", grid$market_id),
    year = paste0("Y", grid$year_id),
    ## Each group has two persistent products in every market-year.
    group = paste0("G", grid$group_id)
  )
}

fit_one <- function(model, d, X_s0 = NULL, s0_pooling = NULL,
                    legacy_none = FALSE, explicit_market_year = FALSE,
                    init = NULL) {
  is_ces <- identical(model, "ces")
  args <- list(
    prices = d$prices, quantities = d$quantities,
    margins = if (is_ces) d$margins else 1 / (80 * 0.7 * d$prices),
    market = d$market, year = d$year, group = d$group,
    demand_model = model,
    supply_model = "bertrand",
    market_side = if (is_ces) "output" else "input",
    use_rho = FALSE,
    share_likelihood = "gaussian",
    quality_parameterization = "collapsed",
    margin_likelihood = "gaussian",
    use_market_resid_scale = FALSE,
    use_hmt = FALSE,
    X_s0 = X_s0,
    stan_options = list(
      algorithm = "Fixed_param", chains = 1L, iter = 1L, warmup = 0L,
      cores = 1L, refresh = 0L, seed = 90417L,
      init = if (is.null(init)) 0 else init
    )
  )
  if (!is.null(s0_pooling)) args$s0_pooling <- s0_pooling
  if (legacy_none) {
    args$use_s0_market_pool <- FALSE
    args$use_s0_market_year_pool <- FALSE
  }
  if (explicit_market_year) args$s0_pooling <- "market_year"
  if (!is.null(init)) args$stan_options$init <- init
  do.call(bayes_est, args)
}

## RStan stores vector declarations of length one with an explicit dim
## attribute.  Preserve those attributes when changing a constrained init;
## passing c(value) drops the declaration and is rejected as a scalar.
init_one <- function(fit) {
  z <- rstan::get_inits(fit@stanfit)
  if (is.list(z) && length(z) == 1L && is.list(z[[1L]]) &&
      !is.null(names(z[[1L]]))) z <- z[[1L]]
  if (!is.list(z) || is.null(names(z))) {
    stop("get_inits() did not return a named init list")
  }
  z[intersect(fit@stanfit@model_pars, names(z))]
}

set_init_values <- function(init, name, value) {
  if (!name %in% names(init)) stop("Missing init parameter: ", name)
  old <- init[[name]]
  value <- as.numeric(value)
  if (length(value) != length(old)) {
    stop("Init length mismatch for ", name, ": expected ", length(old),
         ", found ", length(value))
  }
  if (!is.null(dim(old))) dim(value) <- dim(old)
  init[[name]] <- value
  init
}

par_dim <- function(fit, name) {
  dims <- fit@stanfit@par_dims[[name]]
  if (is.null(dims)) integer() else as.integer(dims)
}

check_dim <- function(fit, name, expected, model, cell) {
  assert_equal(par_dim(fit, name), as.integer(expected), model, cell,
               paste0("parameter_dimension:", name))
}

draw_matrix <- function(fit, name) {
  x <- rstan::extract(fit@stanfit, pars = name, permuted = TRUE)[[name]]
  if (is.null(dim(x))) return(matrix(as.numeric(x), nrow = 1L))
  if (length(dim(x)) == 2L) return(x)
  matrix(x, nrow = dim(x)[1L], ncol = prod(dim(x)[-1L]))
}

check_common_metadata <- function(fit, model, cell, resolved) {
  assert_equal(fit@metadata$s0_pooling_requested, "auto", model, cell,
               "metadata_requested_s0_pooling")
  assert_equal(fit@metadata$s0_pooling_resolved, resolved, model, cell,
               "metadata_resolved_s0_pooling")
  assert_equal(attr(fit@stan_data, "s0_pooling_requested"), "auto", model, cell,
               "stan_data_requested_s0_pooling")
  assert_equal(attr(fit@stan_data, "s0_pooling_resolved"), resolved, model, cell,
               "stan_data_resolved_s0_pooling")
  assert_equal(fit@metadata$prior_normalization$s0_pooling_resolved,
               resolved, model, cell, "prior_normalization_s0_pooling")
  assert_true(is.finite(fit@stan_data$s0_prior_mean) &&
                is.finite(fit@stan_data$s0_prior_sd) &&
                fit@stan_data$s0_prior_sd > 0,
              model, cell, "proper_s0_prior_constants")
  assert_equal(fit@stan_data$quality_parameterization, 2L, model, cell,
               "collapsed_quality_code")
  assert_equal(fit@stan_data$use_rho, 0L, model, cell, "rho_disabled")
  assert_equal(fit@stan_data$use_market_resid_scale, 0L, model, cell,
               "market_residual_scale_disabled")
  assert_equal(fit@stan_data$use_hmt, 0L, model, cell, "hmt_disabled")
  ## Persistent groups are deliberately represented by repeated product rows;
  ## collapsed quality has no separately sampled quality implementation.
  assert_true(anyDuplicated(fit@stan_data$group) > 0L,
              model, cell, "persistent_quality_groups")
  check_dim(fit, "b_group_centered", 0L, model, cell)
}

cells <- list(
  one_market_one_year = c(n_markets = 1L, n_years = 1L, resolved = "none"),
  one_market_three_year = c(n_markets = 1L, n_years = 3L, resolved = "market"),
  three_market_one_year = c(n_markets = 3L, n_years = 1L, resolved = "market"),
  three_market_two_year = c(n_markets = 3L, n_years = 2L, resolved = "market")
)

for (model in c("logit", "ces")) {
  for (cell in names(cells)) {
    spec <- cells[[cell]]
    d <- fixture(as.integer(spec[["n_markets"]]), as.integer(spec[["n_years"]]))
    run_case(model, cell, {
      fit <- fit_one(model, d)
      check_common_metadata(fit, model, cell, unname(spec[["resolved"]]))
      n_market <- as.integer(spec[["n_markets"]])
      n_year <- as.integer(spec[["n_years"]])
      n_mky <- n_market * n_year
      if (n_market == 1L && n_year == 1L) {
        check_dim(fit, "s0_logit", c(n_mky), model, cell)
        check_dim(fit, "mu_s0_pool_vec", 0L, model, cell)
        check_dim(fit, "sigma_s0_pool_vec", 0L, model, cell)
      } else if (n_market == 1L) {
        check_dim(fit, "s0_logit", 0L, model, cell)
        check_dim(fit, "mu_s0_pool_vec", 1L, model, cell)
        check_dim(fit, "sigma_s0_pool_vec", 0L, model, cell)
        check_dim(fit, "s0_market_raw", 0L, model, cell)
        check_dim(fit, "eta_s0_market", 0L, model, cell)
        s0 <- draw_matrix(fit, "s0")
        assert_true(max(abs(s0[, -1L, drop = FALSE] - s0[, 1L])) < 1e-12,
                    model, cell, "single_market_shared_base_s0")
      } else {
        check_dim(fit, "s0_logit", 0L, model, cell)
        check_dim(fit, "mu_s0_pool_vec", 1L, model, cell)
        check_dim(fit, "sigma_s0_pool_vec", 1L, model, cell)
        check_dim(fit, "s0_market_raw", n_market, model, cell)
        check_dim(fit, "eta_s0_market", 0L, model, cell)
        check_dim(fit, "sigma_log_a_vec", 1L, model, cell)
        check_dim(fit, "r_market_a_raw", n_market, model, cell)
      }
    })
  }
}

## X_s0 remains a market-year shift on top of the automatic market base.
for (model in c("logit", "ces")) {
  cell <- "one_market_three_year_X_s0"
  d <- fixture(1L, 3L)
  X_s0 <- matrix(c(-0.5, 0, 0.5), ncol = 1L)
  run_case(model, cell, {
    fit <- fit_one(model, d, X_s0 = X_s0)
    check_common_metadata(fit, model, cell, "market")
    check_dim(fit, "mu_s0_pool_vec", 1L, model, cell)
    check_dim(fit, "sigma_s0_pool_vec", 0L, model, cell)
    check_dim(fit, "beta_s0", 1L, model, cell)
    assert_equal(fit@stan_data$K_s0, 1L, model, cell, "X_s0_column_count")
    assert_equal(as.numeric(fit@stan_data$X_s0), as.numeric(X_s0), model, cell,
                 "X_s0_values")
  })
}

## Explicit market-year pooling keeps the old global market-year hierarchy.
for (model in c("logit", "ces")) {
  cell <- "explicit_market_year"
  d <- fixture(1L, 3L)
  run_case(model, cell, {
    ## Start from an accepted Fixed_param initialization, then retain its
    ## exact Stan array dimensions while separating the three period effects.
    ## This also records the dimensions RStan actually exposes through
    ## get_inits(), which guards the fixture against scalar/vector coercion.
    baseline <- fit_one(model, d, explicit_market_year = TRUE)
    baseline_init <- init_one(baseline)
    assert_equal(dim(baseline_init$s0_pool_raw), 3L, model, cell,
                 "get_inits_s0_pool_raw_dimension")
    assert_equal(dim(baseline_init$sigma_s0_pool_vec), 1L, model, cell,
                 "get_inits_sigma_s0_pool_dimension")
    init <- set_init_values(baseline_init, "s0_pool_raw", c(-1, 0, 1))
    fit <- fit_one(model, d, explicit_market_year = TRUE,
                   init = function() init)
    assert_equal(fit@metadata$s0_pooling_requested, "market_year", model, cell,
                 "explicit_requested_market_year")
    assert_equal(fit@metadata$s0_pooling_resolved, "market_year", model, cell,
                 "explicit_resolved_market_year")
    check_dim(fit, "mu_s0_pool_vec", 1L, model, cell)
    check_dim(fit, "sigma_s0_pool_vec", 1L, model, cell)
    check_dim(fit, "s0_pool_raw", 3L, model, cell)
    check_dim(fit, "s0_logit", 0L, model, cell)
    s0 <- draw_matrix(fit, "s0")
    assert_true(max(s0) - min(s0) > 1e-8, model, cell,
                "explicit_market_year_periods_can_differ")
  })
}

## Explicit legacy FALSE/FALSE maps to independent market-year outside shares.
for (model in c("logit", "ces")) {
  cell <- "legacy_false_false"
  d <- fixture(1L, 3L)
  run_case(model, cell, {
    fit <- fit_one(model, d, legacy_none = TRUE)
    assert_equal(fit@metadata$s0_pooling_requested, "none", model, cell,
                 "legacy_requested_none")
    assert_equal(fit@metadata$s0_pooling_resolved, "none", model, cell,
                 "legacy_resolved_none")
    check_dim(fit, "s0_logit", 3L, model, cell)
    check_dim(fit, "mu_s0_pool_vec", 0L, model, cell)
    check_dim(fit, "sigma_s0_pool_vec", 0L, model, cell)
  })
}

out <- if (length(results)) do.call(rbind, results) else {
  data.frame(model = character(), cell = character(), check = character(),
             status = character(), value = character(), message = character(),
             stringsAsFactors = FALSE)
}
dir.create(dirname(csv_path), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(out, csv_path, row.names = FALSE, na = "")

failed <- out$status == "fail"
if (any(failed)) {
  stop(sum(failed), " validation checks failed; see ", csv_path)
}
cat("Validated ", nrow(out), " checks across ",
    length(unique(out$model)), " demand models; CSV: ", csv_path, "\n", sep = "")
