#!/usr/bin/env Rscript

## Frozen, resumable HMC comparison for integrating market intercepts.
##
## Safety: the default MIC_BENCH_MODE=smoke only constructs and validates the
## data fixture. It never samples. Set MIC_BENCH_MODE=run explicitly to launch
## the frozen 3-case x 2-seed benchmark. Fit objects are never saved.
##
## Inputs:
##   MIC_BENCH_DATA_RDS=/path/to/ready-cross-section.rds
## or MIC_BENCH_MERGERBAYES_DATA_DIR=/path/to/mergerBayes/data
##
## Outputs:
##   MIC_BENCH_OUTPUT (default doc/market_intercept_collapse_hmc.csv)
##   MIC_BENCH_POSTERIOR_OUTPUT (optional parameter-summary CSV)
##
## Frozen cases (all other model/data/prior/HMC settings match):
##   sampled intercept + joint s0 coordinates + total_fraction quality scale
##   sampled intercept + legacy s0 coordinates + total_fraction control
##   collapsed intercept + legacy s0 coordinates + total_fraction
## All cases start at zero on their unconstrained Stan coordinates. For these
## three representations that maps shared parameters to the same economic
## center; the sampled intercept starts at its integrated conditional center.

options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg)) {
  sub("^--file=", "", script_arg[[1L]])
} else {
  file.path(getwd(), "doc", "market_intercept_collapse_hmc_benchmark.R")
}
script_path <- normalizePath(script_path, mustWork = FALSE)
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)

env_value <- function(name, default = "") {
  value <- Sys.getenv(name, unset = default)
  if (is.na(value) || !nzchar(value)) default else value
}
env_int <- function(name, default, min_value = 1L) {
  value <- suppressWarnings(as.integer(env_value(name, as.character(default))))
  if (length(value) != 1L || is.na(value) || value < min_value) {
    stop(name, " must be an integer >= ", min_value, call. = FALSE)
  }
  value
}
env_num <- function(name, default) {
  value <- suppressWarnings(as.numeric(env_value(name, as.character(default))))
  if (length(value) != 1L || is.na(value) || !is.finite(value)) {
    stop(name, " must be finite", call. = FALSE)
  }
  value
}
truthy <- function(x) tolower(x) %in% c("true", "t", "1", "yes", "y")
`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

hash_files <- function(paths) {
  paths <- normalizePath(paths, mustWork = TRUE)
  paste(paste(basename(paths), unname(tools::md5sum(paths)), sep = "="),
        collapse = ";")
}

git_info <- function(path) {
  if (is.null(path) || !length(path) || !nzchar(path) || !dir.exists(path)) {
    return(list(sha = NA_character_, dirty = NA, diff_md5 = NA_character_))
  }
  sha <- tryCatch(system2("git", c("-C", path, "rev-parse", "HEAD"),
                          stdout = TRUE, stderr = FALSE),
                  error = function(e) character())
  status <- tryCatch(system2("git", c("-C", path, "status", "--porcelain"),
                             stdout = TRUE, stderr = FALSE),
                     error = function(e) character())
  diff_file <- tempfile("mic-bench-git-diff-")
  on.exit(unlink(diff_file), add = TRUE)
  invisible(tryCatch(system2("git", c("-C", path, "diff", "--binary", "HEAD"),
                             stdout = diff_file, stderr = FALSE),
                     error = function(e) 1L))
  diff_md5 <- if (file.exists(diff_file)) {
    unname(tools::md5sum(diff_file))
  } else {
    NA_character_
  }
  list(
    sha = if (length(sha) && !is.na(sha[[1L]])) sha[[1L]] else NA_character_,
    dirty = length(status) > 0L,
    diff_md5 = diff_md5
  )
}

canonical_data_dir <- function() {
  requested <- env_value("MIC_BENCH_MERGERBAYES_DATA_DIR")
  if (nzchar(requested)) return(path.expand(requested))
  candidates <- c(
    file.path(repo_root, "..", "mergerBayes", "data"),
    "~/Projects/mergerBayes/data",
    "~/Projects/mergerbayes/data"
  )
  candidates <- path.expand(candidates)
  present <- candidates[file.exists(file.path(candidates, "supreme_data.RData"))]
  if (!length(present)) return(path.expand(candidates[[1L]]))
  present[[1L]]
}

load_ready_rds <- function(path) {
  path <- normalizePath(path.expand(path), mustWork = TRUE)
  dat <- readRDS(path)
  if (!is.data.frame(dat)) stop("Ready RDS must contain a data.frame.", call. = FALSE)
  required <- c("prices", "quantities", "margins", "market", "year", "group")
  missing <- setdiff(required, names(dat))
  if (length(missing)) {
    stop("Ready RDS lacks required columns: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  if (!"margin_observed" %in% names(dat)) {
    dat$margin_observed <- !is.na(dat$margins)
  }
  if (!is.logical(dat$margin_observed) || anyNA(dat$margin_observed)) {
    stop("margin_observed must be a non-missing logical vector.", call. = FALSE)
  }
  dat$prices <- as.numeric(dat$prices)
  dat$quantities <- as.numeric(dat$quantities)
  dat$margins <- as.numeric(dat$margins)
  dat$market <- as.character(dat$market)
  dat$year <- as.character(dat$year)
  dat$group <- as.character(dat$group)
  if (any(!is.finite(dat$prices) | dat$prices <= 0) ||
      any(!is.finite(dat$quantities) | dat$quantities <= 0)) {
    stop("All ready-RDS prices and quantities must be positive and finite.",
         call. = FALSE)
  }
  if (any(dat$margin_observed &
          (!is.finite(dat$margins) | is.na(dat$margins) | dat$margins <= 0))) {
    stop("Observed margins must be positive and finite.", call. = FALSE)
  }
  if (anyNA(dat$market) || anyNA(dat$year) || anyNA(dat$group) ||
      any(!nzchar(dat$market)) || any(!nzchar(dat$year)) || any(!nzchar(dat$group))) {
    stop("market, year, and group labels must be nonempty.", call. = FALSE)
  }
  list(data = dat, source = path, source_hash = hash_files(path), canonical = FALSE)
}

load_canonical_sc_fixture <- function(data_dir) {
  data_file <- file.path(path.expand(data_dir), "supreme_data.RData")
  events_file <- file.path(path.expand(data_dir), "event_overlaps.RData")
  if (!file.exists(data_file)) stop("Missing canonical input: ", data_file, call. = FALSE)
  if (!file.exists(events_file)) stop("Missing canonical event metadata: ", events_file,
                                      call. = FALSE)
  data_env <- new.env(parent = baseenv())
  event_env <- new.env(parent = baseenv())
  load(data_file, envir = data_env)
  load(events_file, envir = event_env)
  if (!exists("simdata", envir = data_env, inherits = FALSE) ||
      !exists("events", envir = event_env, inherits = FALSE)) {
    stop("Canonical files must define simdata and events respectively.", call. = FALSE)
  }
  simdata <- data_env$simdata
  events <- event_env$events
  required <- c("event_id", "year", "tophold", "rate_deposits",
                "total_deposits", "selected_input_markdown_ratio")
  missing <- setdiff(required, names(simdata))
  if (length(missing)) {
    stop("supreme_data.RData lacks required columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  if (!all(c("event_id", "year") %in% names(events))) {
    stop("event_overlaps.RData events object must include event_id and base year.",
         call. = FALSE)
  }
  event_id <- as.character(events$event_id)
  base_year <- as.character(events$year)
  if (anyNA(event_id) || anyNA(base_year) || anyDuplicated(event_id)) {
    stop("Event ids and base years must be complete and unique.", call. = FALSE)
  }
  event_index <- match(as.character(simdata$event_id), event_id)
  data_year <- as.character(simdata$year)
  price <- as.numeric(simdata$rate_deposits)
  quantity <- as.numeric(simdata$total_deposits)
  margin <- as.numeric(simdata$selected_input_markdown_ratio)
  use <- !is.na(event_index) & data_year == base_year[event_index] &
    is.finite(price) & price > 0 & is.finite(quantity) & quantity > 0 &
    is.finite(margin) & margin > 0
  dat <- data.frame(
    prices = price[use], quantities = quantity[use], margins = margin[use],
    market = as.character(simdata$event_id[use]),
    year = data_year[use], group = as.character(simdata$tophold[use]),
    margin_observed = TRUE, stringsAsFactors = FALSE
  )
  dat <- dat[order(as.integer(dat$market), as.integer(dat$year),
                   -dat$quantities), , drop = FALSE]
  rownames(dat) <- NULL
  n_markets <- length(unique(dat$market))
  if (nrow(dat) != 176L || n_markets != 6L ||
      !setequal(unique(dat$market), event_id)) {
    stop("Canonical year-cross-section validation failed: expected 176 rows " ,
         "across all 6 event markets; found N=", nrow(dat), ", M=", n_markets,
         ". Check the input vintage and margin filters.", call. = FALSE)
  }
  list(
    data = dat,
    source = paste(data_file, events_file, sep = ";"),
    source_hash = hash_files(c(data_file, events_file)),
    canonical = TRUE
  )
}

load_fixture <- function() {
  ready <- env_value("MIC_BENCH_DATA_RDS")
  if (nzchar(ready)) return(load_ready_rds(ready))
  load_canonical_sc_fixture(canonical_data_dir())
}

case_table <- data.frame(
  case = c("sampled_joint_s0", "sampled_legacy_s0_control",
           "collapsed_legacy_s0"),
  intercept = c("sampled", "sampled", "collapsed"),
  s0_coordinates = c("joint", "legacy", "legacy"),
  stringsAsFactors = FALSE
)

fixture <- load_fixture()
dat <- fixture$data
if (length(unique(dat$market)) < 2L) {
  stop("This experiment requires at least two markets.", call. = FALSE)
}

mode <- tolower(env_value("MIC_BENCH_MODE", "smoke"))
if (!mode %in% c("smoke", "run")) {
  stop("MIC_BENCH_MODE must be 'smoke' or 'run'.", call. = FALSE)
}
chains <- env_int("MIC_BENCH_CHAINS", 4L)
iter <- env_int("MIC_BENCH_ITER", 2000L)
warmup <- env_int("MIC_BENCH_WARMUP", 1000L)
cores <- env_int("MIC_BENCH_CORES", min(chains, 2L))
adapt_delta <- env_num("MIC_BENCH_ADAPT_DELTA", 0.99)
max_treedepth <- env_int("MIC_BENCH_MAX_TREEDEPTH", 13L)
seeds <- suppressWarnings(as.integer(strsplit(
  env_value("MIC_BENCH_SEEDS", "19091,29091"), ",", fixed = TRUE
)[[1L]]))
if (!length(seeds) || anyNA(seeds)) {
  stop("MIC_BENCH_SEEDS must be comma-separated integers.", call. = FALSE)
}
if (warmup >= iter) stop("warmup must be less than iter.", call. = FALSE)
if (cores > chains) stop("cores must not exceed chains.", call. = FALSE)

cat("Fixture:", fixture$source, "\n")
cat("Rows:", nrow(dat), " Markets:", length(unique(dat$market)),
    " Years:", length(unique(dat$year)), " Observed margins:",
    sum(dat$margin_observed), "\n")
cat("Fixture hash:", fixture$source_hash, "\n")
cat("Mode:", mode, "(HMC sampling only runs when MIC_BENCH_MODE=run)\n")
cat("Frozen sampler: chains=", chains, " iter=", iter, " warmup=", warmup,
    " cores=", cores, " adapt_delta=", adapt_delta,
    " max_treedepth=", max_treedepth, " seeds=", paste(seeds, collapse = ","),
    "\n", sep = "")
print(case_table, row.names = FALSE)

if (mode == "smoke") {
  if (fixture$canonical &&
      (nrow(dat) != 176L || length(unique(dat$market)) != 6L)) {
    stop("Canonical smoke fixture did not meet its 176 x 6 contract.", call. = FALSE)
  }
  cat("Smoke fixture validation passed. No Stan fit was launched.\n")
  quit(save = "no", status = 0L)
}

lib_path <- env_value("MIC_BENCH_LIB")
if (nzchar(lib_path)) .libPaths(c(path.expand(lib_path), .libPaths()))
if (!requireNamespace("antitrustBayes", quietly = TRUE) ||
    !requireNamespace("rstan", quietly = TRUE) ||
    !requireNamespace("posterior", quietly = TRUE) ||
    !requireNamespace("digest", quietly = TRUE)) {
  stop("Run mode requires installed antitrustBayes, rstan, posterior, and digest packages.",
       call. = FALSE)
}

expected_kernel_hashes <- c(
  ces_stackelberg =
    "89700d43789a4c4109633858f551dc537fd20f7c7a74e50bb107b0c4b7cd9716",
  ces_supply =
    "4eaf21ee78113b8d934a7c1bd123be22dee01e1a5594aa7da9e863a196d22821",
  logit_stackelberg =
    "371dd99d72684f7e8f5ead2b51a8bfc505a9e3b5df6c0623d288e386bf7f5074",
  logit_supply =
    "cb8ba797461801cfba6bbe4e6b2e9d76c554d75de48030f0801160f8a16507fe"
)
installed_models <- get("stanmodels", envir = asNamespace("antitrustBayes"))
if (!all(names(expected_kernel_hashes) %in% names(installed_models))) {
  stop("Installed antitrustBayes lacks one or more frozen benchmark kernels.",
       call. = FALSE)
}
installed_kernel_hashes <- vapply(names(expected_kernel_hashes), function(kernel) {
  digest::digest(installed_models[[kernel]]@model_code,
                 algo = "sha256", serialize = FALSE)
}, character(1))
if (!identical(unname(installed_kernel_hashes),
               unname(expected_kernel_hashes))) {
  mismatch <- names(expected_kernel_hashes)[
    installed_kernel_hashes != expected_kernel_hashes
  ]
  stop("Installed antitrustBayes Stan source does not match the frozen ",
       "market-intercept-collapse stage for: ", paste(mismatch, collapse = ", "),
       ". Reinstall the checked-out package before sampling.", call. = FALSE)
}
installed_kernel_hash_manifest <- paste(
  paste(names(installed_kernel_hashes), installed_kernel_hashes, sep = "="),
  collapse = ";"
)

repo_candidates <- list(
  mergerBayes = env_value("MIC_BENCH_MERGERBAYES_ROOT",
                           file.path(repo_root, "..", "mergerBayes")),
  antitrustBayes = repo_root,
  antitrust = env_value("ANTITRUST_PATH", file.path(repo_root, "..", "antitrust")),
  coordination = env_value("COORDINATION_PATH",
                            file.path(repo_root, "..", "coordination"))
)
repo_prov <- lapply(repo_candidates, git_info)
package_version <- as.character(utils::packageVersion("antitrustBayes"))
package_path <- normalizePath(system.file(package = "antitrustBayes"),
                              mustWork = FALSE)

prior_spec <- list(
  prior_manifest = "baseline",
  prior_sigma_resid_market = 0.35,
  prior_sigma_alpha = 0.5,
  prior_sigma_s0_pool = 0.5,
  s0_prior_sd = 0.456234014550216
)
model_spec <- list(
  demand_model = "logit", supply_model = "bertrand", structure = "simultaneous",
  market_side = "input", likelihood = "joint", use_rho = TRUE,
  share_likelihood = "gaussian", margin_likelihood = "gaussian",
  quality_parameterization = "collapsed", quality_scope = "market",
  quality_scale_parameterization = "total_fraction", s0_pooling = "market",
  use_alpha_market_pool = FALSE, use_centered_s0_pool = FALSE,
  use_market_resid_scale = FALSE
)

run_key <- function(case_row, seed) {
  fields <- c(
    fixture$source_hash, case_row$case, seed, chains, iter, warmup, cores,
    adapt_delta, max_treedepth, package_version,
    repo_prov$antitrustBayes$sha, repo_prov$antitrustBayes$diff_md5,
    installed_kernel_hash_manifest, "init=zero_unconstrained_common",
    repo_prov$antitrust$sha, repo_prov$coordination$sha,
    paste(names(model_spec), unlist(model_spec), sep = "=", collapse = ";"),
    paste(names(prior_spec), unlist(prior_spec), sep = "=", collapse = ";")
  )
  hash_file <- tempfile("mic-bench-run-key-")
  on.exit(unlink(hash_file), add = TRUE)
  writeLines(paste(fields, collapse = "\n"), hash_file, useBytes = TRUE)
  unname(tools::md5sum(hash_file))
}

output_path <- path.expand(env_value(
  "MIC_BENCH_OUTPUT", file.path(repo_root, "doc", "market_intercept_collapse_hmc.csv")
))
summary_path <- path.expand(env_value("MIC_BENCH_POSTERIOR_OUTPUT"))
resume <- truthy(env_value("MIC_BENCH_RESUME", "true"))
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
if (nzchar(summary_path)) {
  dir.create(dirname(summary_path), recursive = TRUE, showWarnings = FALSE)
}

atomic_write_csv <- function(x, path) {
  temp <- tempfile(pattern = ".mic-bench-", tmpdir = dirname(path), fileext = ".csv")
  on.exit(unlink(temp), add = TRUE)
  utils::write.csv(x, temp, row.names = FALSE, na = "")
  if (!file.rename(temp, path)) stop("Could not atomically replace ", path, call. = FALSE)
  invisible(path)
}
read_existing <- function(path) {
  if (!resume || !file.exists(path)) return(NULL)
  x <- tryCatch(utils::read.csv(path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(x) || !nrow(x)) return(NULL)
  if (!"benchmark_version" %in% names(x) || any(x$benchmark_version != 1L)) {
    stop("Existing output has an incompatible schema; choose a new MIC_BENCH_OUTPUT.",
         call. = FALSE)
  }
  x
}
results <- read_existing(output_path)
posterior_results <- if (nzchar(summary_path)) read_existing(summary_path) else NULL

economic_parameters <- c(
  "alpha_structural_market", "mu_log_a", "log_kappa", "s0_prob",
  "mu_b_market_vec", "sigma_b_market_vec", "mu_s0_pool_vec",
  "sigma_s0_pool_vec", "sigma_b_strat_economic", "sigma_b_strat",
  "sigma_share", "sigma_margin", "rho", "share_variance_fraction"
)

summarise_fit <- function(fit, elapsed_seconds, case_row, seed, key) {
  stanfit <- fit@stanfit
  sampling_start <- proc.time()[["elapsed"]]
  sampler <- rstan::get_sampler_params(stanfit, inc_warmup = FALSE)
  divergences <- sum(vapply(sampler, function(x) sum(x[, "divergent__"]), numeric(1)))
  depth_hits <- sum(vapply(sampler, function(x) {
    sum(x[, "treedepth__"] >= max_treedepth)
  }, numeric(1)))
  bfmi <- vapply(sampler, function(x) {
    energy <- x[, "energy__"]
    if (length(energy) < 2L || !is.finite(stats::var(energy)) ||
        stats::var(energy) <= 0) return(NA_real_)
    mean(diff(energy)^2) / stats::var(energy)
  }, numeric(1))
  fnames <- stanfit@sim$fnames_oi
  available <- unique(sub("\\[.*$", "", fnames))
  pars <- intersect(economic_parameters, available)
  if (!length(pars)) stop("Could not find monitored economic parameters in stanfit.")
  raw <- rstan::extract(stanfit, pars = pars, permuted = FALSE,
                        inc_warmup = FALSE)
  draws <- posterior::as_draws_array(raw)
  parameter_summary <- posterior::summarise_draws(
    draws, "mean", "sd", "mcse_mean", "rhat", "ess_bulk", "ess_tail"
  )
  parameter_summary$ess_bulk_per_second <- parameter_summary$ess_bulk / elapsed_seconds
  parameter_summary$ess_tail_per_second <- parameter_summary$ess_tail / elapsed_seconds
  finite_rhat <- parameter_summary$rhat[is.finite(parameter_summary$rhat)]
  finite_bulk <- parameter_summary$ess_bulk[is.finite(parameter_summary$ess_bulk)]
  finite_tail <- parameter_summary$ess_tail[is.finite(parameter_summary$ess_tail)]
  finite_bulk_sec <- parameter_summary$ess_bulk_per_second[
    is.finite(parameter_summary$ess_bulk_per_second)]
  finite_tail_sec <- parameter_summary$ess_tail_per_second[
    is.finite(parameter_summary$ess_tail_per_second)]
  max_rhat <- if (length(finite_rhat)) max(finite_rhat) else NA_real_
  min_bulk <- if (length(finite_bulk)) min(finite_bulk) else NA_real_
  min_tail <- if (length(finite_tail)) min(finite_tail) else NA_real_
  min_bulk_sec <- if (length(finite_bulk_sec)) min(finite_bulk_sec) else NA_real_
  min_tail_sec <- if (length(finite_tail_sec)) min(finite_tail_sec) else NA_real_
  min_bulk_var <- if (length(finite_bulk_sec)) {
    parameter_summary$variable[which.min(parameter_summary$ess_bulk_per_second)]
  } else NA_character_
  min_tail_var <- if (length(finite_tail_sec)) {
    parameter_summary$variable[which.min(parameter_summary$ess_tail_per_second)]
  } else NA_character_
  native_time <- tryCatch(rstan::get_elapsed_time(stanfit), error = function(e) NULL)
  mean_warmup <- if (!is.null(native_time) && "warmup" %in% colnames(native_time)) {
    mean(native_time[, "warmup"])
  } else NA_real_
  mean_sampling <- if (!is.null(native_time) && "sampling" %in% colnames(native_time)) {
    mean(native_time[, "sampling"])
  } else NA_real_
  postfit_processing_seconds <- proc.time()[["elapsed"]] - sampling_start
  usable <- isTRUE(divergences == 0L && depth_hits == 0L &&
    length(bfmi) == chains && all(is.finite(bfmi)) && min(bfmi) >= 0.30 &&
    is.finite(max_rhat) && max_rhat <= 1.01 &&
    is.finite(min_bulk) && min_bulk >= 100 &&
    is.finite(min_tail) && min_tail >= 100)
  row <- data.frame(
    benchmark_version = 1L, run_key = key, status = "completed",
    fixture_source = fixture$source, fixture_hash = fixture$source_hash,
    canonical_fixture = fixture$canonical, N = nrow(dat),
    N_market = length(unique(dat$market)), N_year = length(unique(dat$year)),
    case = case_row$case, market_intercept_parameterization = case_row$intercept,
    market_s0_parameterization = case_row$s0_coordinates,
    supply_model = "bertrand", likelihood = "joint_gaussian_rho",
    quality_parameterization = "collapsed", quality_scope = "market",
    quality_scale_parameterization = "total_fraction", s0_pooling = "market",
    alpha_market_pool = FALSE, prior_manifest = "baseline",
    prior_sigma_resid_market = prior_spec$prior_sigma_resid_market,
    prior_sigma_alpha = prior_spec$prior_sigma_alpha,
    prior_sigma_s0_pool = prior_spec$prior_sigma_s0_pool,
    s0_prior_sd = prior_spec$s0_prior_sd,
    seed = seed, chains = chains, iter = iter, warmup = warmup, cores = cores,
    adapt_delta = adapt_delta, max_treedepth = max_treedepth,
    wall_seconds = elapsed_seconds,
    postfit_processing_seconds = postfit_processing_seconds,
    native_warmup_seconds_mean = mean_warmup,
    native_sampling_seconds_mean = mean_sampling,
    divergences = divergences, treedepth_hits = depth_hits,
    min_bfmi = if (all(is.na(bfmi))) NA_real_ else min(bfmi, na.rm = TRUE),
    max_rank_rhat = max_rhat, min_bulk_ess = min_bulk, min_tail_ess = min_tail,
    min_bulk_ess_per_second = min_bulk_sec,
    min_tail_ess_per_second = min_tail_sec,
    min_bulk_ess_variable = min_bulk_var, min_tail_ess_variable = min_tail_var,
    usable = usable, error = "",
    antitrustBayes_version = package_version,
    antitrustBayes_git = repo_prov$antitrustBayes$sha,
    antitrustBayes_dirty = repo_prov$antitrustBayes$dirty,
    antitrustBayes_diff_md5 = repo_prov$antitrustBayes$diff_md5,
    installed_kernel_hashes = installed_kernel_hash_manifest,
    init_mode = "zero_unconstrained_common",
    antitrust_git = repo_prov$antitrust$sha,
    coordination_git = repo_prov$coordination$sha,
    mergerBayes_git = repo_prov$mergerBayes$sha,
    R_version = R.version.string, package_path = package_path,
    stringsAsFactors = FALSE
  )
  for (i in seq_along(bfmi)) row[[paste0("bfmi_chain", i)]] <- bfmi[[i]]
  if (nzchar(summary_path)) {
    parameter_summary$benchmark_version <- 1L
    parameter_summary$run_key <- key
    parameter_summary$case <- case_row$case
    parameter_summary$seed <- seed
    parameter_summary$wall_seconds <- elapsed_seconds
    parameter_summary <- parameter_summary[, c(
      "benchmark_version", "run_key", "case", "seed", "variable", "mean", "sd",
      "mcse_mean", "rhat", "ess_bulk", "ess_tail", "ess_bulk_per_second",
      "ess_tail_per_second", "wall_seconds"
    )]
    attr(row, "posterior_summary") <- parameter_summary
  }
  row
}

append_and_write <- function(existing, addition, path) {
  if (is.null(existing) || !nrow(existing)) {
    combined <- addition
  } else {
    all_names <- union(names(existing), names(addition))
    for (name in setdiff(all_names, names(existing))) existing[[name]] <- NA
    for (name in setdiff(all_names, names(addition))) addition[[name]] <- NA
    combined <- rbind(existing[, all_names, drop = FALSE],
                      addition[, all_names, drop = FALSE])
  }
  atomic_write_csv(combined, path)
  combined
}

for (i in seq_len(nrow(case_table))) {
  spec <- case_table[i, , drop = FALSE]
  for (seed in seeds) {
    key <- run_key(spec, seed)
    completed <- !is.null(results) &&
      any(results$run_key == key & results$status == "completed")
    summary_present <- !nzchar(summary_path) ||
      (!is.null(posterior_results) && nrow(posterior_results) &&
         any(posterior_results$run_key == key))
    if (completed && summary_present) {
      cat("Skipping completed", spec$case, "seed", seed, "(resume).\n")
      next
    }
    cat("Running", spec$case, "seed", seed, "\n")
    started <- proc.time()[["elapsed"]]
    row <- tryCatch({
      fit <- antitrustBayes::bayes_est(
        prices = dat$prices, quantities = dat$quantities, margins = dat$margins,
        market = dat$market, year = dat$year, group = dat$group,
        margin_observed = dat$margin_observed,
        demand_model = "logit", supply_model = "bertrand", structure = "simultaneous",
        market_side = "input", likelihood = "joint", use_rho = TRUE,
        share_likelihood = "gaussian", margin_likelihood = "gaussian",
        quality_parameterization = "collapsed", quality_scope = "market",
        quality_scale_parameterization = "total_fraction",
        market_intercept_parameterization = spec$intercept,
        s0_pooling = "market", market_s0_parameterization = spec$s0_coordinates,
        use_centered_s0_pool = FALSE, use_alpha_market_pool = FALSE,
        use_market_resid_scale = FALSE,
        prior_manifest = "baseline",
        priors = list(
          prior_sigma_resid_market = prior_spec$prior_sigma_resid_market,
          prior_sigma_alpha = prior_spec$prior_sigma_alpha,
          prior_sigma_s0_pool = prior_spec$prior_sigma_s0_pool,
          s0_prior_sd = prior_spec$s0_prior_sd
        ),
        stan_options = list(
          chains = chains, iter = iter, warmup = warmup, cores = cores,
          refresh = 0, seed = seed, init = 0,
          control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth)
        )
      )
      elapsed <- proc.time()[["elapsed"]] - started
      summarise_fit(fit, elapsed, spec, seed, key)
    }, error = function(e) {
      data.frame(
        benchmark_version = 1L, run_key = key, status = "error",
        fixture_source = fixture$source, fixture_hash = fixture$source_hash,
        canonical_fixture = fixture$canonical, N = nrow(dat),
        N_market = length(unique(dat$market)), N_year = length(unique(dat$year)),
        case = spec$case, market_intercept_parameterization = spec$intercept,
        market_s0_parameterization = spec$s0_coordinates,
        seed = seed, chains = chains, iter = iter, warmup = warmup, cores = cores,
        adapt_delta = adapt_delta, max_treedepth = max_treedepth,
        wall_seconds = proc.time()[["elapsed"]] - started,
        divergences = NA_integer_, treedepth_hits = NA_integer_,
        min_bfmi = NA_real_, max_rank_rhat = NA_real_, min_bulk_ess = NA_real_,
        min_tail_ess = NA_real_, min_bulk_ess_per_second = NA_real_,
        min_tail_ess_per_second = NA_real_, usable = FALSE,
        error = conditionMessage(e), antitrustBayes_version = package_version,
        antitrustBayes_git = repo_prov$antitrustBayes$sha,
        antitrustBayes_dirty = repo_prov$antitrustBayes$dirty,
        antitrustBayes_diff_md5 = repo_prov$antitrustBayes$diff_md5,
        installed_kernel_hashes = installed_kernel_hash_manifest,
        init_mode = "zero_unconstrained_common",
        antitrust_git = repo_prov$antitrust$sha,
        coordination_git = repo_prov$coordination$sha,
        mergerBayes_git = repo_prov$mergerBayes$sha,
        R_version = R.version.string, package_path = package_path,
        stringsAsFactors = FALSE
      )
    })
    posterior_row <- attr(row, "posterior_summary")
    row <- as.data.frame(row, stringsAsFactors = FALSE)
    results <- append_and_write(results, row, output_path)
    if (nzchar(summary_path) && !is.null(posterior_row)) {
      if (is.null(posterior_results) || !nrow(posterior_results)) {
        posterior_results <- posterior_row
      } else {
        posterior_results <- rbind(posterior_results, posterior_row)
      }
      posterior_results <- posterior_results[
        !duplicated(paste(posterior_results$run_key,
                          posterior_results$variable, sep = "::")), , drop = FALSE
      ]
      atomic_write_csv(posterior_results, summary_path)
    }
    if (exists("fit")) rm(fit)
    invisible(gc())
    cat("  status=", row$status, " usable=", row$usable %||% FALSE,
        " divergences=", row$divergences %||% NA,
        " max_rhat=", row$max_rank_rhat %||% NA, "\n", sep = "")
  }
}

results <- read_existing(output_path)
print(results[, intersect(c("case", "seed", "status", "wall_seconds",
                            "divergences", "treedepth_hits", "min_bfmi",
                            "max_rank_rhat", "min_bulk_ess", "min_tail_ess",
                            "min_bulk_ess_per_second", "usable", "error"),
                          names(results)), drop = FALSE], row.names = FALSE)
