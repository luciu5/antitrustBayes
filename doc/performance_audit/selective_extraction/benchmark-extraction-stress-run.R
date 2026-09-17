## Run one extraction-only stress cell.  The same Fixed_param kernel and seed
## are used by full and selected processes; this cell is not posterior
## inference and never runs HMC.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("usage: Rscript ... input.rds full|selected output.rds")
input <- readRDS(args[[1L]]); mode <- args[[2L]]; out <- args[[3L]]
if (!mode %in% c("full", "selected")) stop("mode must be full or selected")
library(antitrustBayes)
dat <- input$data
common <- list(
  prices = dat$prices, quantities = dat$quantities, margins = dat$margins,
  market = dat$market, year = dat$year, group = dat$group,
  demand_model = "logit", market_side = "output", supply_model = "bertrand",
  use_rho = FALSE, df = 5, share_likelihood = "gaussian",
  quality_parameterization = "collapsed", margin_likelihood = "gaussian",
  use_market_resid_scale = FALSE, use_alpha_market_pool = FALSE,
  ## Omit all s0 pooling arguments: passing the new and legacy spellings
  ## together is intentionally rejected by the public API.  With one market
  ## and no legacy flags, the documented auto resolution is the no-pool path.
  market_s0_parameterization = "legacy",
  use_centered_alpha = FALSE, use_centered_s0_pool = FALSE, use_hmt = FALSE,
  prior_manifest = "baseline", likelihood = "joint",
  outside_price = input$outside_price,
  stan_options = list(chains = 1L, iter = 200L, warmup = 0L, cores = 1L,
                      refresh = 0L, algorithm = "Fixed_param", init = 0,
                      seed = 20260916L)
)
fit_time <- system.time({ fit <- do.call(bayes_est, common) })[["elapsed"]]
if (!methods::is(fit@stanfit, "stanfit")) stop("Fixed_param did not return a stanfit")
full_time <- system.time({
  full <- rstan::extract(fit@stanfit, permuted = TRUE, inc_warmup = FALSE)
})[["elapsed"]]
full_bytes <- as.numeric(object.size(full))
full_count <- if (is.null(dim(full$a_market))) length(full$a_market) else dim(full$a_market)[[1L]]
selected_time <- system.time({
  selected <- antitrustBayes:::.posterior_extract(fit, "simulation")
})[["elapsed"]]
selected_bytes <- as.numeric(object.size(selected))
selected_count <- antitrustBayes:::.posterior_draw_count_from_draws(selected, fit)
draws <- if (identical(mode, "full")) full else selected
count <- if (identical(mode, "full")) full_count else selected_count
if (count != 200L) stop("Fixed_param row count was not 200: ", count)
if (!identical(mode, "full") && !identical(names(selected), c("a_market", "s0_prob")))
  stop("unexpected selected contract: ", paste(names(selected), collapse = ","))
if (identical(mode, "full") && !all(c("a_market", "s0_prob") %in% names(full)))
  stop("full extraction lacks required canonical fields")
finite <- all(vapply(draws, function(x) all(is.finite(as.numeric(x))), logical(1)))
canonical_hashes <- vapply(c("a_market", "s0_prob"), function(nm)
  digest::digest(draws[[nm]]), character(1))
names(canonical_hashes) <- c("a_market", "s0_prob")
result <- list(
  N = input$N, G = input$G, T = input$T, mode = mode,
  fixed_param = TRUE, posterior_inference = FALSE, HMC = FALSE,
  kernel = fit@stanfit@stanmodel@model_name,
  fit_seconds = fit_time, full_extract_seconds = full_time,
  selected_extract_seconds = selected_time,
  retained_mode_seconds = if (identical(mode, "full")) full_time else selected_time,
  full_bytes = full_bytes, selected_bytes = selected_bytes,
  full_count = full_count, selected_count = selected_count,
  retained_bytes = if (identical(mode, "full")) full_bytes else selected_bytes,
  retained_count = count, retained_fields = names(draws), finite = finite,
  canonical_hashes = canonical_hashes,
  fixed_rows = 200L,
  note = input$description
)
saveRDS(result, out)
write.csv(as.data.frame(result[setdiff(names(result), c("retained_fields", "note"))]),
          sub("\\.rds$", ".csv", out), row.names = FALSE)
cat("stress N=", input$N, " G=", input$G, " mode=", mode,
    " kernel=", result$kernel, " rows=", count,
    " fit_s=", sprintf("%.3f", fit_time),
    " extract_s=", sprintf("%.3f", result$retained_mode_seconds),
    " bytes=", result$retained_bytes, " finite=", finite, "\n", sep = "")
