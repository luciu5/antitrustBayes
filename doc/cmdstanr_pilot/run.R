#!/usr/bin/env Rscript

# Run from the repository root after installing this package and CmdStanR.
# Usage: Rscript doc/cmdstanr_pilot/run.R logit|ces output-directory
# The CSV inputs are frozen subsets of the benchmark-v2 B and CES_B fixtures.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || !args[[1L]] %in% c("logit", "ces")) {
  stop("Usage: run.R logit|ces output-directory", call. = FALSE)
}
demand <- args[[1L]]
out <- normalizePath(args[[2L]], mustWork = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
stopifnot(dir.exists(out))

suppressPackageStartupMessages(library(antitrustBayes))
if ("inference" %in% names(formals(antitrustBayes::bayes_est))) {
  stop("This archived RStan-versus-CmdStanR pilot requires the pre-migration antitrustBayes package. The installed package already uses CmdStanR.",
       call. = FALSE)
}
stopifnot(requireNamespace("cmdstanr", quietly = TRUE),
          requireNamespace("posterior", quietly = TRUE))

input <- file.path("doc", "cmdstanr_pilot", paste0(demand, ".csv"))
dat <- read.csv(input, stringsAsFactors = FALSE)
stopifnot(nrow(dat) == 36L, length(unique(dat$outside_price)) == 1L)
outside_price <- unique(dat$outside_price)
dat$outside_price <- NULL
kernel <- paste0(demand, "_supply")
stan_source <- file.path("inst", "stan", paste0(kernel, ".stan"))
stan_copy <- file.path(out, paste0(kernel, ".stan"))
same_source <- file.exists(stan_copy) &&
  identical(unname(tools::md5sum(stan_source)),
            unname(tools::md5sum(stan_copy)))
if (!same_source && !file.copy(stan_source, stan_copy, overwrite = TRUE)) {
  stop("Could not copy the Stan program to the writable pilot directory.")
}

# Sampled quality provides pointwise log likelihoods for the LOO check.
# Both backends receive the exact Stan data prepared by bayes_est().
seed <- 271828L
chains <- 4L
warmup <- 300L
sampling <- 300L
parallel_chains <- 2L
est_args <- c(as.list(dat), list(
  demand_model = demand, supply_model = "bertrand", structure = "simultaneous",
  market_side = "output", outside_price = outside_price, use_rho = FALSE,
  share_likelihood = "gaussian", margin_likelihood = "gaussian",
  quality_parameterization = "sampled", quality_scope = "market",
  prior_manifest = "baseline", likelihood = "joint",
  stan_options = list(
    chains = chains, iter = warmup + sampling, warmup = warmup,
    cores = parallel_chains, seed = seed, init = 0, refresh = 0,
    control = list(adapt_delta = 0.95, max_treedepth = 13L,
                   metric = "diag_e")
  )
))

clock <- function(expr) {
  elapsed <- system.time(value <- eval.parent(substitute(expr)))[["elapsed"]]
  list(value = value, seconds = unname(elapsed))
}
attempt <- function(expr) {
  tryCatch(expr, error = function(e) structure(
    list(error = conditionMessage(e)), class = "pilot_error"))
}
value_or_na <- function(x, name) {
  if (inherits(x, "pilot_error")) return(NA_real_)
  as.numeric(x[[name]])
}
error_text <- function(x) {
  if (inherits(x, "pilot_error")) x$error else ""
}

cat("RStan public fit:", demand, "\n")
rstan_timed <- clock(do.call(antitrustBayes::bayes_est, est_args))
rfit <- rstan_timed$value
stan_data <- rfit@stan_data
stopifnot(methods::is(rfit@stanfit, "stanfit"),
          identical(as.character(rfit@stanfit@model_name), kernel))
saveRDS(rfit, file.path(out, "rstan-fit.rds"))

cat("CmdStan model compilation:", demand, "\n")
# Avoid building a precompiled header in the shared CmdStan installation.
compile_timed <- clock(cmdstanr::cmdstan_model(
  stan_copy, force_recompile = TRUE,
  cpp_options = list(PRECOMPILED_HEADERS = "false")
))
model <- compile_timed$value
data_names <- names(model$variables()$data)
missing <- setdiff(data_names, names(stan_data))
if (length(missing)) stop("Missing CmdStan data: ", paste(missing, collapse = ", "))
cmd_data <- stan_data[data_names]
stopifnot(!anyDuplicated(names(cmd_data)))
data_path <- file.path(out, "stan-data.json")
cmdstanr::write_stan_json(cmd_data, data_path)
dir.create(file.path(out, "csv"), showWarnings = FALSE)

cat("CmdStan sampling:", demand, "\n")
cmdstan_timed <- clock(model$sample(
  data = data_path, seed = seed, init = 0, refresh = 0,
  chains = chains, parallel_chains = parallel_chains,
  iter_warmup = warmup, iter_sampling = sampling,
  adapt_delta = 0.95, max_treedepth = 13L, metric = "diag_e",
  sig_figs = 18L, output_dir = file.path(out, "csv")
))
cfit <- cmdstan_timed$value
stopifnot(inherits(cfit, "CmdStanMCMC"),
          all(cfit$return_codes() == 0L), cfit$num_chains() == chains)
# All downstream checks run on a reloaded object, including bridge sampling.
cmdstan_fit_path <- file.path(out, "cmdstan-fit.rds")
cfit$save_object(cmdstan_fit_path)
rm(cfit)
invisible(gc())
cfit <- readRDS(cmdstan_fit_path)

key <- c("mu_log_a", "sigma_share", "sigma_margin")
model_variables <- model$variables()
key <- intersect(key, c(names(model_variables$parameters),
                        names(model_variables$transformed_parameters)))
r_draws <- posterior::as_draws_array(rstan::extract(
  rfit@stanfit, pars = key, permuted = FALSE, inc_warmup = FALSE))
c_draws <- cfit$draws(variables = key, format = "draws_array")
r_summary <- posterior::summarise_draws(r_draws)
c_summary <- posterior::summarise_draws(c_draws)
summary <- merge(as.data.frame(r_summary), as.data.frame(c_summary),
                 by = "variable", suffixes = c("_rstan", "_cmdstan"))
write.csv(summary, file.path(out, "posterior-summary.csv"), row.names = FALSE)
r_diag <- rstan::get_sampler_params(rfit@stanfit, inc_warmup = FALSE)
c_diag <- as.array(cfit$sampler_diagnostics())
stopifnot(all(c("divergent__", "treedepth__") %in% dimnames(c_diag)[[3L]]))
r_divergences <- sum(vapply(r_diag, function(x) sum(x[, "divergent__"]), numeric(1)))
c_divergences <- sum(c_diag[, , "divergent__"])
r_depth_hits <- sum(vapply(r_diag, function(x)
  sum(x[, "treedepth__"] >= 13L), numeric(1)))
c_depth_hits <- sum(c_diag[, , "treedepth__"] >= 13L)

cat("LOO and WAIC:", demand, "\n")
loo_r <- attempt({
  ll <- loo::extract_log_lik(rfit@stanfit, parameter_name = "log_lik",
                             merge_chains = FALSE)
  list(looic = loo::loo(ll, r_eff = loo::relative_eff(exp(ll)))$
         estimates["looic", "Estimate"],
       waic = loo::waic(ll)$estimates["waic", "Estimate"])
})
loo_c <- attempt({
  ll <- as.array(cfit$draws(variables = "log_lik", format = "draws_array"))
  list(looic = cfit$loo(variables = "log_lik")$
         estimates["looic", "Estimate"],
       waic = loo::waic(ll)$estimates["waic", "Estimate"])
})

cat("Normalized target check:", demand, "\n")
target <- attempt({
  n <- as.integer(rstan::get_num_upars(rfit@stanfit))
  points <- c(0, 0.2, -0.2)
  params_path <- file.path(out, "target-points.json")
  output_path <- file.path(out, "target-values.csv")
  cmdstanr::write_stan_json(
    list(params_r = matrix(rep(points, each = n), nrow = length(points),
                           byrow = TRUE)), params_path)
  cli <- suppressWarnings(system2(model$exe_file(), c(
    "log_prob", paste0("unconstrained_params=", shQuote(params_path)),
    "jacobian=true", "data", paste0("file=", shQuote(data_path)),
    "output", paste0("file=", shQuote(output_path)), "sig_figs=18"
  ), stdout = TRUE, stderr = TRUE))
  if (!is.null(attr(cli, "status")) && attr(cli, "status") != 0L) {
    stop(paste(cli, collapse = "\n"))
  }
  cmdstan_targets <- read.csv(output_path, comment.char = "#",
                              check.names = FALSE)
  stopifnot(nrow(cmdstan_targets) == length(points))
  targets <- do.call(rbind, lapply(seq_along(points), function(i) {
    point <- points[[i]]
    u <- rep(point, n)
    r <- rstan::log_prob(rfit@stanfit, u, adjust_transform = TRUE)
    c <- cmdstan_targets$lp__[[i]]
    data.frame(point = point, rstan = as.numeric(r), cmdstan = as.numeric(c),
               abs_difference = abs(as.numeric(r) - as.numeric(c)))
  }))
  list(n_upars = n, points = targets,
       abs_difference = max(targets$abs_difference))
})

cat("Bridge evidence:", demand, "\n")
bridge_r <- attempt(clock(bridgesampling::bridge_sampler(
  rfit@stanfit, repetitions = 1L, maxiter = 1000L, silent = TRUE
)))
bridge_c <- attempt(clock(bridgesampling::bridge_sampler(
  cfit, repetitions = 1L, maxiter = 1000L, silent = TRUE
)))

saveRDS(list(rstan = bridge_r, cmdstan = bridge_c, target = target,
             loo_rstan = loo_r, loo_cmdstan = loo_c),
        file.path(out, "comparison-details.rds"))

bridge_r_value <- if (inherits(bridge_r, "pilot_error")) NA_real_ else
  as.numeric(bridge_r$value$logml)
bridge_c_value <- if (inherits(bridge_c, "pilot_error")) NA_real_ else
  as.numeric(bridge_c$value$logml)
row <- data.frame(
  demand = demand, observations = nrow(dat), chains = chains,
  warmup = warmup, sampling = sampling, parallel_chains = parallel_chains,
  rstan_fit_seconds = rstan_timed$seconds,
  cmdstan_compile_seconds = compile_timed$seconds,
  cmdstan_sample_seconds = cmdstan_timed$seconds,
  rstan_min_bulk_ess = min(r_summary$ess_bulk, na.rm = TRUE),
  cmdstan_min_bulk_ess = min(c_summary$ess_bulk, na.rm = TRUE),
  rstan_min_bulk_ess_per_fit_call_second = min(r_summary$ess_bulk, na.rm = TRUE) /
    rstan_timed$seconds,
  cmdstan_min_bulk_ess_per_sample_second = min(c_summary$ess_bulk, na.rm = TRUE) /
    cmdstan_timed$seconds,
  rstan_max_rhat = max(r_summary$rhat, na.rm = TRUE),
  cmdstan_max_rhat = max(c_summary$rhat, na.rm = TRUE),
  rstan_divergences = r_divergences, cmdstan_divergences = c_divergences,
  rstan_treedepth_hits = r_depth_hits, cmdstan_treedepth_hits = c_depth_hits,
  rstan_looic = value_or_na(loo_r, "looic"),
  cmdstan_looic = value_or_na(loo_c, "looic"),
  rstan_waic = value_or_na(loo_r, "waic"),
  cmdstan_waic = value_or_na(loo_c, "waic"),
  rstan_logml = bridge_r_value, cmdstan_logml = bridge_c_value,
  rstan_bridge_seconds = if (inherits(bridge_r, "pilot_error")) NA_real_ else
    bridge_r$seconds,
  cmdstan_bridge_seconds = if (inherits(bridge_c, "pilot_error")) NA_real_ else
    bridge_c$seconds,
  target_abs_difference = if (inherits(target, "pilot_error")) NA_real_ else
    target$abs_difference,
  rstan_bridge_error = error_text(bridge_r),
  cmdstan_bridge_error = error_text(bridge_c),
  target_error = error_text(target),
  loo_error = paste(Filter(nzchar, c(error_text(loo_r), error_text(loo_c))),
                    collapse = " | "),
  stringsAsFactors = FALSE
)
write.csv(row, file.path(out, "pilot-results.csv"), row.names = FALSE)
saveRDS(list(input_sha256 = digest::digest(file = input, algo = "sha256"),
             stan_sha256 = digest::digest(file = stan_source, algo = "sha256"),
             stan_data_sha256 = digest::digest(cmd_data, algo = "sha256"),
             stan_data_json_sha256 = digest::digest(file = data_path, algo = "sha256"),
             r_version = R.version.string,
             rstan_version = as.character(packageVersion("rstan")),
             cmdstanr_version = as.character(packageVersion("cmdstanr")),
             cmdstan_version = as.character(cmdstanr::cmdstan_version()),
             bridgesampling_version = as.character(packageVersion("bridgesampling")),
             seed = seed, settings = est_args$stan_options),
        file.path(out, "manifest.rds"))
print(row)
