## Bounded native profile for the baseline compiled antitrustBayes kernels.
## This script only loads an existing package/fit, uses Fixed_param to create
## a model context, and calls log_prob/grad_log_prob/gqs.  It never compiles,
## runs HMC, installs packages, or writes outside its output directory.

out_dir <- "/tmp/performance-stan-profile"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(antitrustBayes)
  library(rstan)
  library(digest)
})

ns <- asNamespace("antitrustBayes")
models <- get("stanmodels", envir = ns)

archive_paths <- c(
  logit_simultaneous = "/mnt/d/projects/antitrustbayes/doc/statistical_toolkit_runs/20260916/hybrid/fits/logit_bertrand_noisy_271828__D.rds",
  logit_BF = "/mnt/d/projects/antitrustbayes/doc/statistical_toolkit_runs/20260916/hybrid/fits/logit_BF_noisy_271828__D.rds",
  logit_BS = "/mnt/d/projects/antitrustbayes/doc/statistical_toolkit_runs/20260916/hybrid/fits/logit_BS_noisy_271828__D.rds",
  ces_simultaneous = "/mnt/d/projects/antitrustbayes/doc/statistical_toolkit_runs/20260916/hybrid/fits/ces_bertrand_noisy_271828__D.rds"
)

scalar <- function(x, default = NA_real_) {
  if (length(x) != 1L || !is.finite(x)) default else as.numeric(x)
}

canonical_data <- function(dat) {
  dat <- dat[!duplicated(names(dat), fromLast = TRUE)]
  if (is.null(dat$market_s0_parameterization))
    dat$market_s0_parameterization <- 1L
  dat
}

flat_matrix <- function(sf) {
  a <- rstan::extract(sf, permuted = FALSE, inc_warmup = FALSE)
  if (length(dim(a)) != 3L)
    stop("stanfit did not return a draw array")
  matrix(a[, 1L, , drop = FALSE], nrow = dim(a)[1L],
         dimnames = list(NULL, sf@sim$fnames_oi))
}

raw_draws <- function(sf, n = 8L) {
  all_draws <- flat_matrix(sf)
  n <- min(as.integer(n), nrow(all_draws))
  raw_names <- antitrustBayes:::.bayes_bridge_raw_parameter_names(sf)
  raw_cols <- sf@sim$fnames_oi[
    sub("\\[.*$", "", sf@sim$fnames_oi) %in% raw_names
  ]
  if (!length(raw_cols)) stop("no raw parameter columns were found")
  all_draws[seq_len(n), raw_cols, drop = FALSE]
}

flat_one <- function(sf) {
  a <- rstan::extract(sf, permuted = FALSE, inc_warmup = FALSE)
  if (length(dim(a)) != 3L)
    stop("stanfit did not return a draw array")
  setNames(as.numeric(a[1L, 1L, ]), sf@sim$fnames_oi)
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

compare_gq <- function(existing, standalone) {
  old <- flat_matrix(existing)
  now <- flat_matrix(standalone)
  common <- intersect(colnames(now), colnames(old))
  if (!length(common))
    stop("standalone GQ has no names in common with existing output")
  n <- min(nrow(now), nrow(old))
  now <- now[seq_len(n), common, drop = FALSE]
  old <- old[seq_len(n), common, drop = FALSE]
  delta <- now - old
  finite <- is.finite(now) & is.finite(old)
  list(
    common = length(common), rows = n,
    finite = all(finite),
    max_abs = if (any(finite)) max(abs(delta[finite])) else NA_real_,
    max_rel = if (any(finite)) max(abs(delta[finite]) /
      (1 + abs(old[finite]))) else NA_real_,
    nonfinite = sum(!finite)
  )
}

base_row <- function(kind, kernel, demand, structure, dat) {
  data.frame(
    kind = kind, kernel = kernel, demand_model = demand,
    structure = structure, N = as.integer(dat$N), G = as.integer(dat$N_group),
    share_likelihood = mode_text(dat, "share_likelihood",
      c(student_t = 1L, gaussian = 2L)),
    margin_likelihood = mode_text(dat, "margin_likelihood",
      c(student_t = 1L, gaussian = 2L)),
    quality_parameterization = mode_text(dat, "quality_parameterization",
      c(sampled = 1L, collapsed = 2L)),
    setup_status = "not_run", setup_seconds = NA_real_,
    target_status = "not_run", target_seconds = NA_real_,
    gradient_status = "not_run", gradient_seconds = NA_real_,
    target_value = NA_real_, gradient_l2 = NA_real_, u_dim = NA_integer_,
    gqs_status = "not_run", gqs_rows = NA_integer_, gqs_seconds = NA_real_,
    gq_rows_compared = NA_integer_, gq_common = NA_integer_, gq_finite = NA, gq_max_abs = NA_real_,
    gq_max_rel = NA_real_, gq_nonfinite = NA_integer_, error = "",
    stringsAsFactors = FALSE
  )
}

profile_one <- function(kind, kernel, demand, structure, dat,
                        archive_sf = NULL, repetitions = 20L,
                        gq_rows = 8L, label = kind) {
  row <- base_row(kind, kernel, demand, structure, dat)
  model <- models[[kernel]]
  context <- NULL
  setup <- tryCatch({
    tm <- system.time({
      context <- suppressWarnings(rstan::sampling(
        model, data = dat, chains = 1L, iter = 1L, warmup = 0L,
        algorithm = "Fixed_param", init = 0, refresh = 0L,
        seed = 48001L
      ))
    })
    row$setup_status <- if (context@mode == 0L) "ok" else "bad_mode"
    row$setup_seconds <- scalar(tm[["elapsed"]])
    row$u_dim <- as.integer(rstan::get_num_upars(context))
    NULL
  }, error = function(e) conditionMessage(e))
  if (!is.null(setup)) {
    row$error <- paste0("setup: ", setup)
    return(row)
  }

  u <- rep(0, row$u_dim)
  lp <- tryCatch(as.numeric(rstan::log_prob(
    context, u, adjust_transform = TRUE, gradient = FALSE
  )), error = function(e) e)
  gr <- tryCatch(as.numeric(rstan::grad_log_prob(
    context, u, adjust_transform = TRUE
  )), error = function(e) e)
  if (inherits(lp, "condition") || inherits(gr, "condition")) {
    row$target_status <- if (inherits(lp, "condition")) "error" else "ok"
    row$gradient_status <- if (inherits(gr, "condition")) "error" else "ok"
    row$error <- paste(
      c(if (inherits(lp, "condition")) paste0("target: ", lp$message),
        if (inherits(gr, "condition")) paste0("gradient: ", gr$message)),
      collapse = " | "
    )
    return(row)
  }
  row$target_value <- lp
  row$gradient_l2 <- sqrt(sum(gr^2))
  row$target_status <- if (is.finite(lp)) "ok" else "nonfinite"
  row$gradient_status <- if (all(is.finite(gr))) "ok" else "nonfinite"

  tlog <- tryCatch({
    tm <- system.time(replicate(repetitions, rstan::log_prob(
      context, u, adjust_transform = TRUE, gradient = FALSE
    )))
    scalar(tm[["elapsed"]]) / repetitions
  }, error = function(e) e)
  tgrad <- tryCatch({
    tm <- system.time(replicate(repetitions, rstan::grad_log_prob(
      context, u, adjust_transform = TRUE
    )))
    scalar(tm[["elapsed"]]) / repetitions
  }, error = function(e) e)
  if (inherits(tlog, "condition") || inherits(tgrad, "condition")) {
    row$error <- paste(
      c(if (inherits(tlog, "condition")) paste0("target timing: ", tlog$message),
        if (inherits(tgrad, "condition")) paste0("gradient timing: ", tgrad$message)),
      collapse = " | "
    )
  } else {
    row$target_seconds <- tlog
    row$gradient_seconds <- tgrad
  }

  fixed_draws <- raw_draws(context, 1L)
  if (!is.null(archive_sf)) {
    fixed_draws <- raw_draws(archive_sf, gq_rows)
  }
  saveRDS(fixed_draws, file.path(out_dir, paste0("v2_", label, "_raw_draws.rds")))
  gq <- tryCatch({
    tm <- system.time({
      gq_fit <- suppressWarnings(rstan::gqs(
        model, data = dat, draws = fixed_draws, seed = 48002L
      ))
    })
    list(fit = gq_fit, seconds = scalar(tm[["elapsed"]]))
  }, error = function(e) e)
  if (inherits(gq, "condition")) {
    row$gqs_status <- "error"
    row$error <- paste(c(row$error, paste0("gqs: ", conditionMessage(gq))),
                       collapse = " | ")
  } else {
    row$gqs_status <- if (gq$fit@mode == 0L) "ok" else "bad_mode"
    row$gqs_rows <- nrow(fixed_draws)
    row$gqs_seconds <- gq$seconds
    cmp <- if (!is.null(archive_sf)) compare_gq(archive_sf, gq$fit)
      else compare_gq(context, gq$fit)
    row$gq_rows_compared <- cmp$rows
    row$gq_common <- cmp$common
    row$gq_finite <- cmp$finite
    row$gq_max_abs <- cmp$max_abs
    row$gq_max_rel <- cmp$max_rel
    row$gq_nonfinite <- cmp$nonfinite
    saveRDS(gq$fit, file.path(out_dir, paste0("v2_", label, "_gqs_fit.rds")))
  }
  row
}

archive_rows <- list()
archive_manifest <- data.frame(
  label = names(archive_paths), path = unname(archive_paths), exists = file.exists(archive_paths),
  stringsAsFactors = FALSE
)
for (label in names(archive_paths)) {
  path <- archive_paths[[label]]
  if (!file.exists(path)) {
    archive_rows[[label]] <- data.frame(kind = "archive", kernel = NA_character_,
      demand_model = NA_character_, structure = NA_character_, N = NA_integer_,
      G = NA_integer_, setup_status = "missing_archive", setup_seconds = NA_real_,
      target_status = "not_run", target_seconds = NA_real_, gradient_status = "not_run",
      gradient_seconds = NA_real_, target_value = NA_real_, gradient_l2 = NA_real_,
      u_dim = NA_integer_, gqs_status = "not_run", gqs_rows = NA_integer_,
      share_likelihood = NA_character_, margin_likelihood = NA_character_,
      quality_parameterization = NA_character_,
      gqs_seconds = NA_real_, gq_rows_compared = NA_integer_, gq_common = NA_integer_, gq_finite = NA,
      gq_max_abs = NA_real_, gq_max_rel = NA_real_, gq_nonfinite = NA_integer_,
      error = "archive missing", stringsAsFactors = FALSE)
    next
  }
  archive <- tryCatch(readRDS(path), error = function(e) e)
  if (inherits(archive, "condition")) {
    archive_rows[[label]] <- data.frame(kind = "archive", kernel = NA_character_,
      demand_model = NA_character_, structure = NA_character_, N = NA_integer_,
      G = NA_integer_, setup_status = "read_error", setup_seconds = NA_real_,
      target_status = "not_run", target_seconds = NA_real_, gradient_status = "not_run",
      gradient_seconds = NA_real_, target_value = NA_real_, gradient_l2 = NA_real_,
      u_dim = NA_integer_, gqs_status = "not_run", gqs_rows = NA_integer_,
      share_likelihood = NA_character_, margin_likelihood = NA_character_,
      quality_parameterization = NA_character_,
      gqs_seconds = NA_real_, gq_rows_compared = NA_integer_, gq_common = NA_integer_, gq_finite = NA,
      gq_max_abs = NA_real_, gq_max_rel = NA_real_, gq_nonfinite = NA_integer_,
      error = conditionMessage(archive), stringsAsFactors = FALSE)
    next
  }
  dat <- canonical_data(archive@stan_data)
  md <- archive@metadata
  archive_rows[[label]] <- profile_one(
    kind = "archive", kernel = archive@stanfit@model_name,
    demand = md$demand_model, structure = md$structure, dat = dat,
    archive_sf = archive@stanfit, repetitions = 20L, gq_rows = 8L,
    label = label
  )
}

## Make compact, valid Gaussian/Gaussian collapsed and sampled simultaneous
## data for scaling. Both modes use G=N through N=1,000 and G=1,000 at
## N=10,000, so the dense sampled path is measured only through G=1,000.
make_data <- function(N, demand, quality, G = if (quality == "sampled") N else 1L) {
  group <- paste0("G", ((seq_len(N) - 1L) %% G) + 1L)
  args <- list(
    prices = seq(1, 2, length.out = N), quantities = rep(1, N),
    margins = rep(.2, N), market = rep("M", N), year = rep("Y", N),
    group = group, supply_model = 1L, use_cutoff = FALSE,
    use_rho = FALSE, use_market_resid_scale = FALSE, use_hmt = FALSE,
    use_alpha_market_pool = FALSE, use_s0_market_year_pool = FALSE,
    use_s0_market_pool = FALSE, use_centered_alpha = FALSE,
    use_centered_s0_pool = FALSE, market_s0_parameterization = "legacy",
    priors = get("resolve_prior_manifest", ns)("baseline", list()),
    grainsize = 1L, margin_observed = NULL, outside_price = 1,
    likelihood_mode = 1L, loo_cluster = rep("M", N), estimate_df = FALSE,
    resid_df_fixed = 5, structure = "simultaneous", cutoff_scale = NULL,
    partition_prior = "cutoff", k_support = NULL, k_prior = NULL,
    k_scope = "market", share_likelihood = "gaussian",
    margin_likelihood = "gaussian",
    quality_parameterization = quality
  )
  fun <- if (demand == "logit") get("prepare_stan_data", ns)
    else get("prepare_ces_stan_data", ns)
  if (demand == "ces") args$market_side <- "output"
  do.call(fun, args)
}

scaling_rows <- list()
for (demand in c("logit", "ces")) {
  kernel <- if (demand == "logit") "logit_supply" else "ces_supply"
  for (N in c(100L, 1000L, 10000L)) {
    for (quality in c("collapsed", "sampled")) {
      G <- if (N <= 1000L) N else 1000L
      label <- paste("scaling", demand, quality, paste0("N", N, "G", G), sep = "_")
      dat <- tryCatch(make_data(N, demand, quality, G), error = function(e) e)
      if (inherits(dat, "condition")) {
        scaling_rows[[label]] <- data.frame(
          kind = "scaling", kernel = kernel, demand_model = demand,
          structure = "simultaneous", N = N, G = G,
          setup_status = "data_error", setup_seconds = NA_real_,
          target_status = "not_run", target_seconds = NA_real_,
          gradient_status = "not_run", gradient_seconds = NA_real_,
          target_value = NA_real_, gradient_l2 = NA_real_, u_dim = NA_integer_,
          gqs_status = "not_run", gqs_rows = NA_integer_, gqs_seconds = NA_real_,
          share_likelihood = NA_character_, margin_likelihood = NA_character_,
          quality_parameterization = NA_character_,
          gq_rows_compared = NA_integer_, gq_common = NA_integer_, gq_finite = NA, gq_max_abs = NA_real_,
          gq_max_rel = NA_real_, gq_nonfinite = NA_integer_,
          error = conditionMessage(dat), stringsAsFactors = FALSE)
        next
      }
      scaling_rows[[label]] <- profile_one(
        kind = "scaling", kernel = kernel, demand = demand,
        structure = "simultaneous", dat = dat, repetitions = if (N >= 10000L) 3L else 10L,
        gq_rows = 1L, label = label
      )
      gc()
    }
  }
}

all_rows <- do.call(rbind, c(archive_rows, scaling_rows))
write.csv(all_rows, file.path(out_dir, "native_profile_v2_results.csv"), row.names = FALSE)
write.csv(archive_manifest, file.path(out_dir, "archive_manifest_v2.csv"), row.names = FALSE)

## The dense sampled-quality Helmert matrix would require G * (G - 1)
## doubles.  Record the N=10,000, G=10,000 storage bound without allocating
## that matrix or asking Stan to validate the corresponding data.
storage_estimates <- do.call(rbind, lapply(c("logit", "ces"), function(demand) {
  data.frame(
    demand_model = demand, quality_parameterization = "sampled",
    N = 10000L, G = 10000L, status = "skipped_storage_cap",
    dense_helmert_bytes = 8 * 10000 * 9999,
    dense_helmert_mib = (8 * 10000 * 9999) / 2^20,
    note = "No native allocation; bound is dense Q only.",
    stringsAsFactors = FALSE
  )
}))
write.csv(storage_estimates, file.path(out_dir, "storage_estimates_v2.csv"), row.names = FALSE)

## Record the RStan threading API and the functional state of the four model
## kernels.  No compile is attempted; this is an API/source readiness probe.
model_code <- vapply(models, slot, character(1), name = "model_code")
thread_old <- rstan::rstan_options("threads_per_chain")
thread_api <- tryCatch({
  rstan::rstan_options(threads_per_chain = 2L)
  got <- rstan::rstan_options("threads_per_chain")
  rstan::rstan_options(threads_per_chain = thread_old)
  identical(as.integer(got), 2L)
}, error = function(e) {
  try(rstan::rstan_options(threads_per_chain = thread_old), silent = TRUE)
  FALSE
})
stan_cxxflags <- paste(capture.output(StanHeaders:::CxxFlags()), collapse = " ")
thread_rows <- data.frame(
  kernel = names(model_code),
  model_code_sha256 = vapply(model_code, digest::digest, character(1), algo = "sha256"),
  has_reduce_sum = grepl("reduce_sum", model_code, fixed = TRUE),
  has_map_rect = grepl("map_rect", model_code, fixed = TRUE),
  has_openmp = grepl("openmp|STAN_THREADS", model_code, ignore.case = TRUE),
  rstan_threads_option_roundtrip = thread_api,
  rstan_threads_option_before = as.integer(thread_old),
  rcppparallel_default_threads = as.integer(RcppParallel::defaultNumThreads()),
  stanheaders_cxxflags = stan_cxxflags,
  stringsAsFactors = FALSE
)
write.csv(thread_rows, file.path(out_dir, "thread_readiness.csv"), row.names = FALSE)

capture.output(sessionInfo(), file = file.path(out_dir, "session_info.txt"))
cat("native profile complete\n")
