## Standalone GQ row amortization.  All draws are constrained raw-parameter
## columns selected by the package schema helper.  This measures rstan::gqs
## calls only; the fitted intercept captures per-call setup/context overhead.

out_dir <- "/tmp/performance-stan-profile"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages({
  library(antitrustBayes)
  library(rstan)
})
ns <- asNamespace("antitrustBayes")
models <- get("stanmodels", envir = ns)

flat_matrix <- function(sf) {
  a <- rstan::extract(sf, permuted = FALSE, inc_warmup = FALSE)
  stopifnot(length(dim(a)) == 3L)
  matrix(a[, 1L, , drop = FALSE], nrow = dim(a)[1L],
         dimnames = list(NULL, sf@sim$fnames_oi))
}
raw_draws <- function(sf, n = 8L) {
  all <- flat_matrix(sf)
  raw <- antitrustBayes:::.bayes_bridge_raw_parameter_names(sf)
  cols <- sf@sim$fnames_oi[
    sub("\\[.*$", "", sf@sim$fnames_oi) %in% raw
  ]
  stopifnot(length(cols) > 0L, nrow(all) >= n)
  all[seq_len(n), cols, drop = FALSE]
}
canonical_data <- function(dat) {
  dat <- dat[!duplicated(names(dat), fromLast = TRUE)]
  if (is.null(dat$market_s0_parameterization))
    dat$market_s0_parameterization <- 1L
  dat
}
mode_text <- function(dat, field, code_map) {
  value <- attr(dat, field, exact = TRUE)
  if (length(value) == 1L && !is.na(value)) return(as.character(value))
  value <- dat[[field]]
  if (length(value) == 1L && !is.na(value)) {
    hit <- names(code_map)[match(as.integer(value), unname(code_map))]
    if (!is.na(hit)) return(hit)
  }
  NA_character_
}

make_data <- function(N, quality = "sampled", G = N) {
  args <- list(
    prices = seq(1, 2, length.out = N), quantities = rep(1, N),
    margins = rep(.2, N), market = rep("M", N), year = rep("Y", N),
    group = paste0("G", ((seq_len(N) - 1L) %% G) + 1L),
    supply_model = 1L, use_cutoff = FALSE, use_rho = FALSE,
    use_market_resid_scale = FALSE, use_hmt = FALSE,
    use_alpha_market_pool = FALSE, use_s0_market_year_pool = FALSE,
    use_s0_market_pool = FALSE, use_centered_alpha = FALSE,
    use_centered_s0_pool = FALSE, market_s0_parameterization = "legacy",
    priors = get("resolve_prior_manifest", ns)("baseline", list()),
    grainsize = 1L, margin_observed = NULL, outside_price = 1,
    likelihood_mode = 1L, loo_cluster = rep("M", N), estimate_df = FALSE,
    resid_df_fixed = 5, structure = "simultaneous", cutoff_scale = NULL,
    partition_prior = "cutoff", k_support = NULL, k_prior = NULL,
    k_scope = "market", share_likelihood = "gaussian",
    margin_likelihood = "gaussian", quality_parameterization = quality
  )
  do.call(get("prepare_stan_data", ns), args)
}

archive_paths <- c(
  logit_BF = "/mnt/d/projects/antitrustbayes/doc/statistical_toolkit_runs/20260916/hybrid/fits/logit_BF_noisy_271828__D.rds",
  ces_simultaneous = "/mnt/d/projects/antitrustbayes/doc/statistical_toolkit_runs/20260916/hybrid/fits/ces_bertrand_noisy_271828__D.rds"
)

fixtures <- list()
for (label in names(archive_paths)) {
  if (!file.exists(archive_paths[[label]])) next
  fit <- readRDS(archive_paths[[label]])
  dat <- canonical_data(fit@stan_data)
  sf <- fit@stanfit
  available <- min(200L, nrow(flat_matrix(sf)))
  fixtures[[label]] <- list(
    label = label, kernel = sf@model_name, data = dat,
    draws = raw_draws(sf, available), source = "archived_posterior"
  )
}

## A large sampled-quality fixture exercises many groups while retaining the
## exact same-family Gaussian/Gaussian likelihood contract.  Fixed_param gives
## one valid constrained raw row; replication supplies a controlled number of
## rows for GQ timing without running HMC.
large_dat <- make_data(1000L, quality = "sampled", G = 1000L)
large_context <- suppressWarnings(rstan::sampling(
  models[["logit_supply"]], data = large_dat, chains = 1L, iter = 1L,
  warmup = 0L, algorithm = "Fixed_param", init = 0, refresh = 0L,
  seed = 48101L
))
large_one <- raw_draws(large_context, 1L)
large_draws <- large_one[rep(1L, 200L), , drop = FALSE]
fixtures[["logit_sampled_N1000G1000_fixed_repeated"]] <- list(
  label = "logit_sampled_N1000G1000_fixed_repeated",
  kernel = "logit_supply", data = large_dat, draws = large_draws,
  source = "fixed_param_row_repeated"
)

counts <- c(1L, 8L, 50L, 200L)
repeats <- 3L
observations <- list()
for (fixture_name in names(fixtures)) {
  f <- fixtures[[fixture_name]]
  for (n in counts) {
    draws <- f$draws[seq_len(n), , drop = FALSE]
    for (replicate_id in seq_len(repeats)) {
      elapsed <- tryCatch({
        tm <- system.time(suppressWarnings(rstan::gqs(
          models[[f$kernel]], data = f$data, draws = draws,
          seed = 48110L + replicate_id + n
        )))
        as.numeric(tm[["elapsed"]])
      }, error = function(e) structure(NA_real_, error = conditionMessage(e)))
      err <- attr(elapsed, "error", exact = TRUE)
      observations[[length(observations) + 1L]] <- data.frame(
        fixture = f$label, kernel = f$kernel, source = f$source,
        N = as.integer(f$data$N), G = as.integer(f$data$N_group),
        share_likelihood = mode_text(f$data, "share_likelihood",
          c(student_t = 1L, gaussian = 2L)),
        margin_likelihood = mode_text(f$data, "margin_likelihood",
          c(student_t = 1L, gaussian = 2L)),
        quality_parameterization = mode_text(f$data, "quality_parameterization",
          c(sampled = 1L, collapsed = 2L)),
        rows = n, replicate = replicate_id, seconds = elapsed,
        error = if (is.null(err)) "" else err, stringsAsFactors = FALSE
      )
    }
  }
}
observations <- do.call(rbind, observations)
write.csv(observations, file.path(out_dir, "gq_amortization_v2_observations.csv"), row.names = FALSE)

fits <- list()
for (fixture_name in unique(observations$fixture)) {
  sub <- observations[observations$fixture == fixture_name &
                        is.finite(observations$seconds), , drop = FALSE]
  if (nrow(sub) < 2L) next
  fit <- lm(seconds ~ rows, data = sub)
  co <- coef(fit)
  f <- fixtures[[fixture_name]]
  fits[[fixture_name]] <- data.frame(
    fixture = fixture_name, kernel = f$kernel, source = f$source,
    N = as.integer(f$data$N), G = as.integer(f$data$N_group),
    successful_calls = nrow(sub), intercept_seconds = unname(co[["(Intercept)"]]),
    seconds_per_row = unname(co[["rows"]]), residual_sd = summary(fit)$sigma,
    r_squared = summary(fit)$r.squared, stringsAsFactors = FALSE
  )
}
write.csv(do.call(rbind, fits), file.path(out_dir, "gq_amortization_v2_fit.csv"), row.names = FALSE)
cat("gq amortization complete\n")
