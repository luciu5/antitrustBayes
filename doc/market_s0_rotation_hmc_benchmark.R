#!/usr/bin/env Rscript

## Frozen paired HMC runner for the market/outside-share coordinate rotation.
##
## This script prepares its plan and constrained starts without sampling by
## default.  HMC requires an explicit --run and the independent native gate
## marker ANTITRUSTBAYES_ROTATION_NATIVE_GATE=approved.  Legacy and rotated
## workers run in separate R processes so they cannot accidentally share a
## loaded package namespace or native DLL.

options(stringsAsFactors = FALSE)

rotation_script_file <- commandArgs(trailingOnly = FALSE)
rotation_script_file <- rotation_script_file[startsWith(rotation_script_file, "--file=")]
rotation_script_dir <- if (length(rotation_script_file)) {
  dirname(normalizePath(sub("^--file=", "", rotation_script_file[[1L]]),
                         mustWork = FALSE))
} else file.path(getwd(), "doc")
rotation_repo_dir <- dirname(rotation_script_dir)
rotation_archive_root <- "/mnt/d/projects/antitrustbayes/doc/market_s0_rotation_runs/20260916"
rotation_default_archive <- file.path(
  rotation_archive_root, "..", "..", "collapsed_gaussian_runs",
  "20260916", "cga-ladder-benchmark-20260916", "stanfit",
  "logit_bertrand_noisy_271828__B__ladder_multi_market.rds"
)
rotation_local_archive <- file.path(
  rotation_repo_dir, "doc", "collapsed_gaussian_runs", "20260916",
  "cga-ladder-benchmark-20260916", "stanfit",
  "logit_bertrand_noisy_271828__B__ladder_multi_market.rds"
)
if (file.exists(rotation_local_archive)) rotation_default_archive <- rotation_local_archive

rotation_parse_option <- function(args, name, default = NULL) {
  key <- paste0("--", name, "=")
  hit <- args[startsWith(args, key)]
  if (!length(hit)) return(default)
  sub(key, "", hit[[1L]], fixed = TRUE)
}

rotation_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (any(args %in% c("--help", "-h"))) return(list(action = "help"))
  known <- c("--prepare", "--run", "--resume", "--worker", "--scale",
             "--help", "-h")
  unknown <- args[startsWith(args, "--") &
    !args %in% known &
    !startsWith(args, "--out=") &
    !startsWith(args, "--archive=") &
    !startsWith(args, "--mode=") &
    !startsWith(args, "--seed=") &
    !startsWith(args, "--scale-m=")]
  if (length(unknown)) stop("Unknown option: ", unknown[[1L]], call. = FALSE)
  if (sum(args %in% c("--prepare", "--run", "--worker")) > 1L)
    stop("Use exactly one of --prepare, --run, or --worker.", call. = FALSE)
  action <- if ("--worker" %in% args) "worker" else
    if ("--run" %in% args) "run" else "prepare"
  list(
    action = action,
    resume = "--resume" %in% args,
    scale = "--scale" %in% args,
    out = rotation_parse_option(args, "out"),
    archive = rotation_parse_option(args, "archive", rotation_default_archive),
    mode = rotation_parse_option(args, "mode"),
    seed = rotation_parse_option(args, "seed"),
    scale_m = rotation_parse_option(args, "scale-m")
  )
}

rotation_require_archive <- function(path) {
  path <- normalizePath(path, mustWork = FALSE)
  if (!file.exists(path)) stop("Archived legacy fit is missing: ", path,
                              call. = FALSE)
  path
}

rotation_require_output <- function(path, resume = FALSE) {
  if (is.null(path) || !nzchar(path)) {
    path <- file.path(rotation_archive_root, "paired-hmc-20260916")
  }
  path <- normalizePath(path, mustWork = FALSE)
  root <- normalizePath(rotation_archive_root, mustWork = FALSE)
  if (!identical(path, root) && !startsWith(path, paste0(root, "/")))
    stop("--out must be under the ignored Windows benchmark archive: ", root,
         call. = FALSE)
  if (dir.exists(path) && length(list.files(path, all.files = TRUE,
                                           no.. = TRUE)) && !resume)
    stop("Existing benchmark output must be preserved; use a new --out or --resume.",
         call. = FALSE)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(path)) stop("Unable to create benchmark output: ", path,
                              call. = FALSE)
  path
}

rotation_require_digest <- function() {
  if (!requireNamespace("digest", quietly = TRUE))
    stop("The frozen runner requires the digest package.", call. = FALSE)
  invisible(TRUE)
}

rotation_source_benchmark <- function() {
  path <- file.path(rotation_repo_dir, "doc", "collapsed_gaussian_benchmark.R")
  if (!file.exists(path)) stop("Cannot locate collapsed Gaussian fixture source: ",
                              path, call. = FALSE)
  e <- new.env(parent = globalenv())
  sys.source(path, e)
  e
}

rotation_fixture_hash <- function(panel) {
  digest::digest(panel$data, algo = "sha256")
}

rotation_archived_fit <- function(path) {
  if (!requireNamespace("antitrustBayes", quietly = TRUE))
    stop("Load the preserved baseline antitrustBayes library before preparing starts.",
         call. = FALSE)
  value <- readRDS(path)
  if (!methods::is(value, "BayesFit"))
    stop("Archived object is not a BayesFit.", call. = FALSE)
  value
}

rotation_primary_model_parameters <- c(
  "mu_b_market_vec", "sigma_b_market_vec", "b_market_raw",
  "b_group_centered", "sigma_b_strat", "mu_b_fringe_vec",
  "mu_log_a_free", "sigma_log_a_vec", "r_market_a_raw",
  "log_a_market_pool", "log_kappa_vec", "log_kappa_market", "s0_logit",
  "mu_s0_pool_vec", "sigma_s0_pool_vec", "s0_pool_raw", "eta_s0_pool",
  "s0_market_raw", "eta_s0_market", "cutoff_share", "rho_raw",
  "log_sigma_share", "log_sigma_margin", "log_sigma_share_market_raw",
  "log_sigma_margin_market_raw", "sigma_log_sigma_share_market",
  "sigma_log_sigma_margin_market", "beta_s0", "bargaining_unit", "nu_minus_2"
)

rotation_map_primary_starts <- function(archived_fit) {
  st <- archived_fit@stanfit
  if (length(st@inits) != 4L)
    stop("Archived legacy fit must contain exactly four chain starts.", call. = FALSE)
  ## The BayesFit wrapper owns stan_data; the embedded stanfit does not.
  old_mky_market <- as.integer(archived_fit@stan_data$mky_to_market)
  legacy <- vector("list", 4L)
  rotated <- vector("list", 4L)
  rows <- vector("list", 4L)
  for (chain in seq_len(4L)) {
    old <- st@inits[[chain]]
    if (!all(rotation_primary_model_parameters %in% names(old)))
      stop("Archived chain ", chain, " is missing a required constrained parameter.",
           call. = FALSE)
    old_model <- old[rotation_primary_model_parameters]
    q <- as.numeric(old$s0_logit)
    b <- as.numeric(old$mu_b_market_vec[[1L]] +
      old$sigma_b_market_vec[[1L]] * old$b_market_raw)
    if (length(q) != length(old_mky_market) || length(b) < max(old_mky_market))
      stop("Archived chain ", chain, " has incompatible market/outside dimensions.",
           call. = FALSE)
    c_s0 <- q + b[old_mky_market]
    new_model <- old_model
    new_model$s0_logit <- NULL
    new_model$c_s0 <- array(c_s0, dim = length(c_s0))
    shared <- setdiff(names(old_model), "s0_logit")
    shared_start_exact <- isTRUE(all(vapply(shared, function(nm) {
      isTRUE(all.equal(old_model[[nm]], new_model[[nm]],
                       check.attributes = FALSE))
    }, logical(1))))
    legacy[[chain]] <- old_model
    rotated[[chain]] <- new_model
    rows[[chain]] <- data.frame(
      chain = chain,
      legacy_start_sha256 = digest::digest(old_model, algo = "sha256"),
      rotated_start_sha256 = digest::digest(new_model, algo = "sha256"),
      n_market_year = length(q), n_market = length(b),
      max_reconstruction_error = max(abs((c_s0 - b[old_mky_market]) - q)),
      shared_start_exact = shared_start_exact,
      rotation_map = "c_s0 = s0_logit + b_market[mky_to_market]",
      stringsAsFactors = FALSE
    )
  }
  list(legacy = legacy, rotated = rotated, manifest = do.call(rbind, rows),
       market_map = old_mky_market)
}

rotation_primary_plan <- function() {
  seeds <- c(271828L, 314159L)
  rows <- do.call(rbind, lapply(seeds, function(seed) {
    data.frame(
      profile = "primary",
      case_id = paste0("primary_", c("legacy", "rotated"), "_seed", seed),
      mode = c("legacy", "rotated"), seed = seed,
      fixture = "cga_ladder_fixtures()$multi_market",
      N = 96L, G = 48L, N_market = 2L, N_year = 2L,
      share_likelihood = "gaussian", quality_parameterization = "collapsed",
      margin_likelihood = "gaussian", use_rho = FALSE,
      use_market_resid_scale = FALSE, use_alpha_market_pool = FALSE,
      use_s0_market_pool = FALSE, use_s0_market_year_pool = FALSE,
      use_centered_alpha = FALSE, use_centered_s0_pool = FALSE,
      likelihood = "joint", prior_manifest = "baseline",
      chains = 4L, iter = 1000L, warmup = 500L, adapt_delta = .95,
      max_treedepth = 13L, metric = "diag_e", cores = 1L,
      stringsAsFactors = FALSE
    )
  }))
  rownames(rows) <- NULL
  rows
}

rotation_scaling_plan <- function() {
  M <- c(10L, 100L, 500L, 2000L)
  seed <- c(271828L, 314159L, 161803L, 141421L)
  plan <- rbind(
    data.frame(profile = "scaling", case_id = "scaling_legacy_M10",
               mode = "legacy", M = 10L, N = 120L, G = 60L,
               seed = seed[[1L]], condition = "both primary rotated fits usable",
               stringsAsFactors = FALSE),
    data.frame(profile = "scaling", case_id = paste0("scaling_rotated_M", M),
               mode = "rotated", M = M, N = 12L * M, G = 6L * M,
               seed = seed, condition = "both primary rotated fits usable",
               stringsAsFactors = FALSE)
  )
  plan$fixture <- "rotation_scaling_fixture(M, seed)"
  plan$fixture_sha256 <- vapply(seq_len(nrow(plan)), function(i) {
    rotation_fixture_hash(rotation_scaling_fixture(plan$M[[i]], plan$seed[[i]]))
  }, character(1))
  plan
}

rotation_validate_primary_archive <- function(panel, archived_fit) {
  dat <- panel$data
  if (nrow(dat) != 96L || length(unique(dat$market)) != 2L ||
      length(unique(dat$year)) != 2L || length(unique(dat$group)) != 48L)
    stop("Frozen primary fixture dimensions changed.", call. = FALSE)
  expected_hash <- "047aee8bf159b4104d5f08f7e795219aa21398ea9aa1d224185678fda69d0084"
  actual_hash <- rotation_fixture_hash(panel)
  if (!identical(actual_hash, expected_hash))
    stop("Frozen primary fixture hash changed: ", actual_hash, call. = FALSE)
  old_data <- archived_fit@data
  for (nm in c("prices", "quantities", "margins")) {
    if (max(abs(as.numeric(dat[[nm]]) - as.numeric(old_data[[nm]]))) > 0)
      stop("Archived legacy fit and fixture differ in ", nm, call. = FALSE)
  }
  for (nm in c("market", "year", "group")) {
    if (!identical(as.character(dat[[nm]]), as.character(old_data[[nm]])))
      stop("Archived legacy fit and fixture differ in ", nm, call. = FALSE)
  }
  sd <- archived_fit@stan_data
  if (!identical(as.integer(sd$share_likelihood), 2L) ||
      !identical(as.integer(sd$quality_parameterization), 2L) ||
      !identical(as.integer(sd$margin_likelihood), 2L) ||
      !identical(as.integer(sd$use_rho), 0L) ||
      !identical(as.integer(sd$use_s0_market_pool), 0L) ||
      !identical(as.integer(sd$use_s0_market_year_pool), 0L))
    stop("Archived fit does not have the frozen collapsed/unpooled specification.",
         call. = FALSE)
  list(status = "passed", fixture_sha256 = actual_hash,
       n_observations = nrow(dat), n_market = length(unique(dat$market)),
       n_year = length(unique(dat$year)), n_group = length(unique(dat$group)))
}

rotation_make_manifest <- function(panel, archived_path, starts, checks,
                                   scaling_start_manifest = NULL,
                                   out_dir = NULL) {
  data.frame(
    schema = "market_s0_rotation_hmc_v1",
    frozen_date = "2026-09-16",
    protocol = "doc/market_s0_rotation_benchmark_protocol.md",
    fixture = "cga_ladder_fixtures()$multi_market",
    fixture_sha256 = rotation_fixture_hash(panel),
    archived_legacy_fit = normalizePath(archived_path),
    archived_legacy_fit_sha256 = digest::digest(file = archived_path, algo = "sha256"),
    primary_seeds = "271828,314159",
    chains = 4L, iter = 1000L, warmup = 500L, adapt_delta = .95,
    max_treedepth = 13L, metric = "diag_e", cores = 1L,
    share_likelihood = "gaussian", quality_parameterization = "collapsed",
    margin_likelihood = "gaussian", df = 5, use_rho = FALSE,
    use_market_resid_scale = FALSE, use_alpha_market_pool = FALSE,
    use_s0_market_pool = FALSE, use_s0_market_year_pool = FALSE,
    use_centered_alpha = FALSE, use_centered_s0_pool = FALSE,
    likelihood = "joint", prior_manifest = "baseline",
    outside_market = "output", primary_rotation = "legacy versus rotated",
    legacy_argument = "market_s0_parameterization omitted",
    rotated_argument = "market_s0_parameterization = 'rotated'",
    gates = "divergences=0;depth_hits=0;BFMI>=0.30;max_Rhat<=1.01;bulk/tail_ESS>=100",
    scaling_condition = "execute only when both primary rotated fits are usable",
    scaling_design = "M=10,100,500,2000; two periods; six firms per market; N=12M",
    native_gate_required = "ANTITRUSTBAYES_ROTATION_NATIVE_GATE=approved",
    no_hmc_in_prepare = TRUE,
    start_map = "archived legacy constrained s0_logit; rotated c_s0=s0_logit+b_market",
    runner_sha256 = digest::digest(file = file.path(
      rotation_script_dir, "market_s0_rotation_hmc_benchmark.R"), algo = "sha256"),
    protocol_sha256 = digest::digest(file = file.path(
      rotation_repo_dir, "doc", "market_s0_rotation_benchmark_protocol.md"),
      algo = "sha256"),
    frozen_attempt = if (is.null(out_dir)) NA_character_ else basename(out_dir),
    diagnostics_revision = "strict nonfinite gates; within-market correlations",
    primary_plan_sha256 = digest::digest(rotation_primary_plan(), algo = "sha256"),
    scaling_plan_sha256 = digest::digest(rotation_scaling_plan(), algo = "sha256"),
    primary_start_manifest_sha256 = digest::digest(starts$manifest, algo = "sha256"),
    scaling_start_manifest_sha256 = if (is.null(scaling_start_manifest)) NA_character_
      else digest::digest(scaling_start_manifest, algo = "sha256"),
    start_max_error = max(starts$manifest$max_reconstruction_error),
    preparation_checks = paste(checks$status, collapse = ";"),
    stringsAsFactors = FALSE
  )
}

rotation_write_prepare <- function(out_dir, archive_path) {
  rotation_require_digest()
  env <- rotation_source_benchmark()
  panel <- env$cga_ladder_fixtures()$multi_market
  archived <- rotation_archived_fit(archive_path)
  checks <- list()
  add_check <- function(name, status, detail) {
    checks[[length(checks) + 1L]] <<- data.frame(
      check = name, status = status, detail = detail, stringsAsFactors = FALSE)
  }
  validation <- tryCatch(rotation_validate_primary_archive(panel, archived),
    error = function(e) e)
  if (inherits(validation, "condition")) stop(validation)
  add_check("fixture_archive_identity", "passed",
            paste0("sha256=", validation$fixture_sha256))
  starts <- rotation_map_primary_starts(archived)
  if (max(starts$manifest$max_reconstruction_error) > 1e-12)
    stop("Frozen start map is not numerically exact.", call. = FALSE)
  add_check("four_chain_start_map", "passed", "max error <= 1e-12")
  if (!all(starts$manifest$shared_start_exact))
    stop("Mapped rotated starts do not preserve every shared economic coordinate.",
         call. = FALSE)
  add_check("paired_economic_starts", "passed",
            "all non-outside coordinates identical by chain")
  if (!grepl("market_s0_parameterization", paste(readLines(
      file.path(rotation_repo_dir, "R", "estimate.R"), warn = FALSE), collapse = "\n"),
      fixed = TRUE))
    stop("Current source does not expose market_s0_parameterization.", call. = FALSE)
  add_check("production_api_source", "passed", "mode argument present")
  legacy_args <- rotation_worker_args(panel, "legacy", starts$legacy,
                                      seed = 271828L)
  rotated_args <- rotation_worker_args(panel, "rotated", starts$rotated,
                                       seed = 271828L)
  if ("market_s0_parameterization" %in% names(legacy_args) ||
      !identical(rotated_args$market_s0_parameterization, "rotated"))
    stop("Legacy/rotated worker argument contract is not frozen.", call. = FALSE)
  add_check("worker_argument_contract", "passed",
            "legacy omits mode flag; rotated passes mode flag")
  if (!identical(nrow(rotation_primary_plan()), 4L))
    stop("Primary plan is not exactly four paired fits.", call. = FALSE)
  add_check("primary_plan", "passed", "two seeds x legacy/rotated")

  scaling <- rotation_scaling_plan()
  scaling_starts <- lapply(seq_len(nrow(scaling)), function(i) {
    panel_i <- rotation_scaling_fixture(scaling$M[[i]], scaling$seed[[i]])
    contract <- list(mky_n = panel_i$mky_n, N_market = panel_i$M,
                     N_year = panel_i$T)
    rotation_scaling_initials(contract, scaling$mode[[i]])
  })
  names(scaling_starts) <- paste0(scaling$mode, "_M", scaling$M)
  scaling_start_manifest <- data.frame(
    case_id = scaling$case_id,
    mode = scaling$mode,
    M = scaling$M,
    fixture_sha256 = scaling$fixture_sha256,
    start_sha256 = vapply(scaling_starts, digest::digest, character(1),
                          algo = "sha256"),
    stringsAsFactors = FALSE
  )
  add_check("scaling_fixture_hashes", "passed",
            paste0("frozen custom O(N) fixtures for M=",
                   paste(scaling$M, collapse = ",")))
  add_check("scaling_start_freeze", "passed",
            "deterministic starts serialized for every planned scaling case")

  utils::write.csv(rotation_primary_plan(), file.path(out_dir, "primary_plan.csv"),
                   row.names = FALSE)
  utils::write.csv(rotation_scaling_plan(), file.path(out_dir, "scaling_plan.csv"),
                   row.names = FALSE)
  utils::write.csv(scaling_start_manifest,
                   file.path(out_dir, "scaling_initial_states_manifest.csv"),
                   row.names = FALSE)
  utils::write.csv(starts$manifest, file.path(out_dir, "initial_states_manifest.csv"),
                   row.names = FALSE)
  check_df <- do.call(rbind, checks)
  utils::write.csv(check_df, file.path(out_dir, "preparation_checks.csv"),
                   row.names = FALSE)
  utils::write.csv(rotation_make_manifest(panel, archive_path, starts, check_df,
                                         scaling_start_manifest, out_dir),
                   file.path(out_dir, "benchmark_manifest.csv"), row.names = FALSE)
  writeLines(c(
    "Frozen market/outside-share paired HMC benchmark",
    "Preparation completed without HMC execution.",
    "Run requires ANTITRUSTBAYES_ROTATION_NATIVE_GATE=approved.",
    "Primary: cga_ladder_fixtures()$multi_market, legacy versus rotated,",
    "Gaussian shares + collapsed quality + Gaussian margins + rho=0,",
    "two seeds (271828, 314159), four chains, 1000/500, delta=.95, depth=13.",
    "Primary starts are mapped from archived four-chain states; scaling fixtures",
    "and starts are deterministic and serialized during preparation.",
    "Scaling is conditional on both rotated primary fits passing every gate.",
    "All fit RDS, logs, and sampler artifacts belong below this directory."
  ), file.path(out_dir, "README.txt"))
  starts$scaling <- scaling_starts
  starts$scaling_manifest <- scaling_start_manifest
  saveRDS(starts, file.path(out_dir, "initial_states.rds"), version = 3L)
  invisible(list(output = out_dir, checks = check_df, starts = starts$manifest))
}

rotation_with_seed <- function(seed, expr) {
  had <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = .GlobalEnv)
  set.seed(as.integer(seed))
  on.exit({
    if (had) assign(".Random.seed", old, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  force(expr)
}

rotation_lse <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

rotation_scaling_fixture <- function(M, seed) {
  M <- as.integer(M)
  if (length(M) != 1L || M < 1L) stop("M must be positive", call. = FALSE)
  F <- 6L; T <- 2L; N <- M * F * T
  market <- rep(seq_len(M), each = F * T)
  year <- rep(rep(seq_len(T), each = F), M)
  firm <- rep(seq_len(F), times = M * T)
  group <- (market - 1L) * F + firm
  mky_key <- paste(market, year, sep = "_")
  mky_levels <- unique(mky_key)
  mky <- match(mky_key, mky_levels)
  mky_n <- tabulate(mky, nbins = length(mky_levels))
  price <- 9.8 + 0.015 * market + 0.09 * year + 0.025 * firm +
    0.08 * sin(seq_len(N) / 11)
  price_sd <- stats::sd(price)
  price_std <- (price - mean(price)) / price_sd
  s0_mky <- stats::plogis(-2.2 + 0.08 * sin(seq_along(mky_n) / 5))
  quality <- 0.04 * sin(group / 7) + 0.02 * cos(group / 13)
  quality <- quality - mean(quality)
  latent <- stats::qlogis(s0_mky[mky]) - 6 * price_std + quality
  noisy <- rotation_with_seed(seed, latent + stats::rnorm(N, 0, 0.04))
  log_share <- numeric(N)
  for (j in seq_along(mky_n)) {
    ii <- which(mky == j)
    log_share[ii] <- noisy[ii] - rotation_lse(noisy[ii])
  }
  quantities <- exp(log_share)
  inside <- quantities
  for (j in seq_along(mky_n)) {
    ii <- which(mky == j)
    inside[ii] <- quantities[ii] / sum(quantities[ii])
  }
  inv_markup <- price * 6 * (1 - inside * (1 - s0_mky[mky]))
  if (any(!is.finite(inv_markup)) || any(inv_markup <= 1))
    stop("scaling fixture produced invalid inverse markups", call. = FALSE)
  data <- data.frame(
    prices = price, quantities = quantities, margins = 1 / inv_markup,
    market = paste0("market_", market), year = as.character(year),
    group = paste0("market_", market, "_firm_", firm),
    stringsAsFactors = FALSE
  )
  list(data = data, parameters = list(outside_price = rep(1.3, length(mky_n))),
       N = N, G = M * F, M = M, T = T, mky = mky, mky_n = mky_n,
       seed = as.integer(seed), fixture = "controlled_two_period_six_firm")
}

rotation_scaling_initials <- function(data, mode) {
  mky_n <- as.integer(data$mky_n)
  s0 <- stats::qlogis(1 / (mky_n + 1))
  init <- list(
    mu_b_market_vec = 0,
    sigma_b_market_vec = 0.5,
    b_market_raw = rep(0, data$N_market),
    b_group_centered = numeric(0), sigma_b_strat = 0.5,
    mu_b_fringe_vec = numeric(0), mu_log_a_free = log(0.8),
    sigma_log_a_vec = numeric(0), r_market_a_raw = numeric(0),
    log_a_market_pool = numeric(0), log_kappa_vec = numeric(0),
    log_kappa_market = numeric(0), s0_logit = s0,
    mu_s0_pool_vec = numeric(0), sigma_s0_pool_vec = numeric(0),
    s0_pool_raw = numeric(0), eta_s0_pool = numeric(0),
    s0_market_raw = numeric(0), eta_s0_market = numeric(0),
    cutoff_share = numeric(0), rho_raw = numeric(0),
    log_sigma_share = log(0.2), log_sigma_margin = log(0.2),
    log_sigma_share_market_raw = numeric(0),
    log_sigma_margin_market_raw = numeric(0),
    sigma_log_sigma_share_market = numeric(0),
    sigma_log_sigma_margin_market = numeric(0), beta_s0 = numeric(0),
    bargaining_unit = numeric(0), nu_minus_2 = numeric(0)
  )
  if (identical(mode, "rotated")) {
    mky_market <- match(seq_along(data$mky_n),
                        rep(seq_len(data$N_market), each = data$N_year))
    init$c_s0 <- array(s0 + init$mu_b_market_vec +
      init$sigma_b_market_vec * init$b_market_raw[mky_market], dim = length(s0))
    init$s0_logit <- NULL
  }
  init
}

rotation_unwrap <- function(value) {
  if (is.null(value)) return(NULL)
  if (methods::is(value, "BayesFit")) return(value@stanfit)
  if (methods::is(value, "stanfit")) return(value)
  NULL
}

rotation_matrix <- function(x) {
  if (is.null(x)) return(NULL)
  d <- dim(x)
  if (is.null(d)) return(matrix(as.numeric(x), ncol = 1L))
  matrix(as.numeric(x), nrow = d[[1L]], ncol = prod(d[-1L]))
}

rotation_cor <- function(x, y) {
  if (is.null(x) || is.null(y)) return(NA_real_)
  x <- as.numeric(x); y <- as.numeric(y)
  n <- min(length(x), length(y)); if (n < 3L) return(NA_real_)
  x <- x[seq_len(n)]; y <- y[seq_len(n)]
  if (stats::sd(x) <= 0 || stats::sd(y) <= 0) return(NA_real_)
  as.numeric(stats::cor(x, y))
}

rotation_sampler_diagnostics <- function(st, max_treedepth = 13L) {
  sp <- rstan::get_sampler_params(st, inc_warmup = FALSE)
  div <- sum(vapply(sp, function(x) sum(x[, "divergent__"] > 0), numeric(1)))
  depth <- sum(vapply(sp, function(x) sum(x[, "treedepth__"] >= max_treedepth), numeric(1)))
  bfmi <- vapply(sp, function(x) {
    e <- as.numeric(x[, "energy__"])
    if (length(e) < 2L || stats::var(e) <= 0) return(NA_real_)
    mean(diff(e)^2) / stats::var(e)
  }, numeric(1))
  step <- vapply(sp, function(x) mean(x[, "stepsize__"]), numeric(1))
  c(divergences = div, depth_hits = depth,
    bfmi_nonfinite = !all(is.finite(bfmi)),
    min_bfmi = if (all(is.finite(bfmi))) min(bfmi) else NA_real_,
    step_size_mean = mean(step), step_size_min = min(step), step_size_max = max(step))
}

rotation_fit_diagnostics <- function(fit, mode, seed, elapsed, start_hash,
                                     case_id, failure = "", market_map = NULL) {
  base <- data.frame(case_id = case_id, mode = mode, seed = seed,
                     status = if (nzchar(failure)) "fit_error" else "diagnostic_failure",
                     error = failure, elapsed_seconds = elapsed,
                     native_stan_seconds = NA_real_, divergences = NA_real_,
                     depth_hits = NA_real_, min_bfmi = NA_real_, max_rhat = NA_real_,
                     min_bulk_ess = NA_real_, min_tail_ess = NA_real_,
                     monitor_nonfinite = NA, constant_columns_ignored = NA_integer_,
                     monitored_columns = NA_integer_,
                     structural_ess_per_stan_second = NA_real_,
                     step_size_mean = NA_real_, step_size_min = NA_real_,
                     step_size_max = NA_real_, sigma_b_strat_sigma_share_corr = NA_real_,
                     sigma_b_strat_sigma_margin_corr = NA_real_,
                     sigma_share_sigma_margin_corr = NA_real_,
                     economic_intercept_s0_corr = NA_real_,
                     economic_b_market_s0_corr = NA_real_,
                     economic_intercept_s0_corr_summary = "max_abs_per_market",
                     economic_b_market_s0_corr_summary = "max_abs_per_market",
                     economic_within_market_max_abs_corr = NA_real_,
                     economic_within_mky_max_abs_corr = NA_real_,
                     coordinate_braw_s0_corr = NA_real_,
                     coordinate_braw_c_s0_corr = NA_real_,
                     coordinate_braw_s0_corr_summary = "max_abs_per_cell",
                     coordinate_braw_c_s0_corr_summary = "max_abs_per_cell",
                     coordinate_braw_s0_max_abs_corr = NA_real_,
                     coordinate_braw_c_s0_max_abs_corr = NA_real_,
                     usable = FALSE, gate_failures = "", stringsAsFactors = FALSE)
  st <- rotation_unwrap(fit)
  if (is.null(st) || nzchar(failure)) return(base)
  smp <- rotation_sampler_diagnostics(st)
  base[names(smp)] <- smp
  base$native_stan_seconds <- sum(rstan::get_elapsed_time(st))
  ## Keep scaling diagnostics bounded in memory.  In particular, do not
  ## extract observation-level generated quantities for M=2000; the gate is
  ## based on sampled parameters plus the small structural summaries below.
  diagnostic_pars <- unique(c(rotation_primary_model_parameters, "c_s0",
                              "sigma_share", "sigma_margin", "b_market",
                              "gamma_market", "s0_prob", "s0",
                              "alpha_structural_market"))
  diagnostic_pars <- intersect(diagnostic_pars, st@sim$pars_oi)
  raw <- rstan::extract(st, pars = diagnostic_pars, permuted = TRUE,
                        inc_warmup = FALSE)
  arr <- rstan::extract(st, pars = diagnostic_pars, permuted = FALSE,
                        inc_warmup = FALSE)
  draws <- posterior::as_draws_array(arr)
  tab <- posterior::summarise_draws(draws)
  draw_matrix <- as.matrix(posterior::as_draws_matrix(draws))
  draw_variables <- posterior::variables(draws)
  constant_draw <- vapply(seq_len(ncol(draw_matrix)), function(j) {
    x <- draw_matrix[, j]
    ## A nonfinite column is never considered constant for gate purposes.
    all(is.finite(x)) && length(x) > 0L && all(x == x[[1L]])
  }, logical(1))
  tab_constant <- constant_draw[match(tab$variable, draw_variables)]
  tab_constant[is.na(tab_constant)] <- FALSE
  monitored <- !tab_constant
  tab_finite <- is.finite(tab$rhat) & is.finite(tab$ess_bulk) &
    is.finite(tab$ess_tail)
  base$monitor_nonfinite <- any(!tab_finite[monitored])
  if (!any(monitored)) base$monitor_nonfinite <- TRUE
  base$constant_columns_ignored <- sum(tab_constant)
  base$monitored_columns <- sum(monitored)
  if (any(monitored)) {
    base$max_rhat <- max(tab$rhat[monitored])
    base$min_bulk_ess <- min(tab$ess_bulk[monitored])
    base$min_tail_ess <- min(tab$ess_tail[monitored])
  } else {
    base$max_rhat <- Inf
    base$min_bulk_ess <- Inf
    base$min_tail_ess <- Inf
  }
  structural <- grepl(
    "^(mu_b_market_vec|sigma_b_market_vec|b_market_raw|c_s0|s0_logit|s0_prob|gamma_market|b_market|alpha_structural_market|sigma_b_strat|sigma_share|sigma_margin|log_sigma_share|log_sigma_margin)",
    tab$variable
  )
  if (!any(structural)) structural <- rep(TRUE, nrow(tab))
  structural <- structural & monitored
  if (any(structural)) {
    base$min_bulk_ess <- min(tab$ess_bulk[structural])
    base$min_tail_ess <- min(tab$ess_tail[structural])
  }
  if (is.finite(base$native_stan_seconds) && base$native_stan_seconds > 0)
    base$structural_ess_per_stan_second <- base$min_bulk_ess / base$native_stan_seconds
  get <- function(nm) if (!is.null(raw[[nm]])) rotation_matrix(raw[[nm]]) else NULL
  sb <- get("sigma_b_strat"); ss <- get("sigma_share"); smg <- get("sigma_margin")
  base$sigma_b_strat_sigma_share_corr <- rotation_cor(sb, ss)
  base$sigma_b_strat_sigma_margin_corr <- rotation_cor(sb, smg)
  base$sigma_share_sigma_margin_corr <- rotation_cor(ss, smg)
  s0 <- get("s0_prob"); gm <- get("gamma_market"); bm <- get("b_market")
  if (!is.null(s0) && !is.null(gm) && !is.null(bm) && length(market_map)) {
    map <- as.integer(market_map)
    n_market <- max(map)
    ## The economic ridge is measured within each market, against the mean
    ## logit outside level.  Flattening market columns would mix market means.
    s0_logit_mky <- qlogis(pmin(pmax(s0, .Machine$double.eps),
                                1 - .Machine$double.eps))
    s0m <- do.call(cbind, lapply(seq_len(n_market), function(m)
      rowMeans(s0_logit_mky[, map == m, drop = FALSE])))
    gm_corr <- vapply(seq_len(n_market), function(m)
      rotation_cor(gm[, m], s0m[, m]), numeric(1))
    bm_corr <- vapply(seq_len(n_market), function(m)
      rotation_cor(bm[, m], s0m[, m]), numeric(1))
    base$economic_intercept_s0_corr <- if (any(is.finite(gm_corr)))
      max(abs(gm_corr), na.rm = TRUE) else NA_real_
    base$economic_b_market_s0_corr <- if (any(is.finite(bm_corr)))
      max(abs(bm_corr), na.rm = TRUE) else NA_real_
    base$economic_within_market_max_abs_corr <- base$economic_b_market_s0_corr
    for (m in seq_len(n_market)) {
      base[[paste0("economic_intercept_market_", m, "_corr")]] <- gm_corr[[m]]
      base[[paste0("economic_b_market_", m, "_corr")]] <- bm_corr[[m]]
    }
    mky_corr <- vapply(seq_along(map), function(j)
      rotation_cor(bm[, map[[j]]], s0_logit_mky[, j]), numeric(1))
    base$economic_within_mky_max_abs_corr <- if (any(is.finite(mky_corr)))
      max(abs(mky_corr), na.rm = TRUE) else NA_real_
    for (j in seq_along(map))
      base[[paste0("economic_mky_", j, "_b_market_logit_s0_corr")]] <- mky_corr[[j]]

    br <- get("b_market_raw"); cs <- get("c_s0"); q <- get("s0_logit")
    if (!is.null(br)) {
      bmap <- br[, map, drop = FALSE]
      c_corr <- if (!is.null(cs)) vapply(seq_along(map), function(j)
        rotation_cor(bmap[, j], cs[, j]), numeric(1)) else rep(NA_real_, length(map))
      q_corr <- if (!is.null(q)) vapply(seq_along(map), function(j)
        rotation_cor(bmap[, j], q[, j]), numeric(1)) else rep(NA_real_, length(map))
      base$coordinate_braw_c_s0_corr <- if (any(is.finite(c_corr)))
        max(abs(c_corr), na.rm = TRUE) else NA_real_
      base$coordinate_braw_s0_corr <- if (any(is.finite(q_corr)))
        max(abs(q_corr), na.rm = TRUE) else NA_real_
      base$coordinate_braw_c_s0_max_abs_corr <- base$coordinate_braw_c_s0_corr
      base$coordinate_braw_s0_max_abs_corr <- base$coordinate_braw_s0_corr
      for (j in seq_along(map)) {
        base[[paste0("coordinate_braw_c_s0_cell_", j, "_corr")]] <- c_corr[[j]]
        base[[paste0("coordinate_braw_s0_cell_", j, "_corr")]] <- q_corr[[j]]
      }
    }
  }
  gates <- c(
    divergences = base$divergences == 0,
    depth_hits = base$depth_hits == 0,
    bfmi_finite = isTRUE(base$bfmi_nonfinite == 0),
    min_bfmi = base$min_bfmi >= .30,
    max_rhat = base$max_rhat <= 1.01,
    min_bulk_ess = base$min_bulk_ess >= 100,
    min_tail_ess = base$min_tail_ess >= 100,
    monitored_finite = !isTRUE(base$monitor_nonfinite)
  )
  base$usable <- all(is.finite(gates)) && all(gates)
  base$gate_failures <- if (base$usable) "" else
    paste(names(gates)[!gates | is.na(gates)], collapse = ";")
  base$status <- if (base$usable) "usable" else "diagnostic_failure"
  base
}

rotation_rbind_fill <- function(frames) {
  if (!length(frames)) return(data.frame())
  columns <- unique(unlist(lapply(frames, names), use.names = FALSE))
  frames <- lapply(frames, function(x) {
    missing <- setdiff(columns, names(x))
    for (nm in missing) x[[nm]] <- NA
    x[columns]
  })
  do.call(rbind, frames)
}

rotation_worker_args <- function(panel, mode, init, seed) {
  args <- list(
    prices = panel$data$prices, quantities = panel$data$quantities,
    margins = panel$data$margins, market = panel$data$market,
    year = panel$data$year, group = panel$data$group,
    demand_model = "logit", market_side = "output", supply_model = "bertrand",
    structure = "simultaneous", use_cutoff = FALSE, use_rho = FALSE, df = 5,
    share_likelihood = "gaussian", quality_parameterization = "collapsed",
    margin_likelihood = "gaussian", use_market_resid_scale = FALSE,
    use_alpha_market_pool = FALSE, use_s0_market_year_pool = FALSE,
    use_s0_market_pool = FALSE, use_centered_alpha = FALSE,
    use_centered_s0_pool = FALSE, use_hmt = FALSE, prior_manifest = "baseline",
    likelihood = "joint", margin_observed = rep(TRUE, nrow(panel$data)),
    outside_price = as.numeric(panel$parameters$outside_price),
    stan_options = list(chains = 4L, iter = 1000L, warmup = 500L, cores = 1L,
                        refresh = 0L, init = init, seed = as.integer(seed),
                        control = list(adapt_delta = .95, max_treedepth = 13L,
                                       metric = "diag_e"))
  )
  if (identical(mode, "rotated")) args$market_s0_parameterization <- "rotated"
  args
}

rotation_run_worker <- function(mode, seed, out_dir, archive_path,
                                scale_m = NULL) {
  mode <- match.arg(mode, c("legacy", "rotated"))
  if (!identical(Sys.getenv("ANTITRUSTBAYES_ROTATION_NATIVE_GATE"), "approved"))
    stop("HMC is locked until ANTITRUSTBAYES_ROTATION_NATIVE_GATE=approved.",
         call. = FALSE)
  if (!requireNamespace("antitrustBayes", quietly = TRUE) ||
      !requireNamespace("rstan", quietly = TRUE) ||
      !requireNamespace("posterior", quietly = TRUE) ||
      !requireNamespace("digest", quietly = TRUE))
    stop("Worker requires antitrustBayes, rstan, posterior, and digest.", call. = FALSE)
  env <- rotation_source_benchmark()
  panel <- if (is.null(scale_m)) env$cga_ladder_fixtures()$multi_market else
    rotation_scaling_fixture(as.integer(scale_m), as.integer(seed))
  if (is.null(scale_m)) {
    archived <- rotation_archived_fit(archive_path)
    mapped <- rotation_map_primary_starts(archived)
    init <- mapped[[mode]][[1L]]
    ## The primary starts are chain-specific; pass all four below.
    starts <- mapped[[mode]]
    case_id <- paste0("primary_", mode, "_seed", seed)
  } else {
    ## Scaling fixtures and starts were frozen by --prepare.  Load exactly
    ## those starts rather than rebuilding a package-dependent contract here.
    frozen_path <- file.path(out_dir, "initial_states.rds")
    if (!file.exists(frozen_path))
      stop("Frozen initial_states.rds is missing; run --prepare first.", call. = FALSE)
    frozen <- readRDS(frozen_path)
    key <- paste0(mode, "_M", as.integer(scale_m))
    if (is.null(frozen$scaling[[key]]))
      stop("Frozen scaling starts are missing for ", key, call. = FALSE)
    starts <- replicate(4L, frozen$scaling[[key]], simplify = FALSE)
    scaling_plan <- rotation_scaling_plan()
    expected <- scaling_plan$fixture_sha256[
      scaling_plan$mode == mode & scaling_plan$M == as.integer(scale_m)]
    if (length(expected) != 1L ||
        !identical(rotation_fixture_hash(panel), unname(expected[[1L]])))
      stop("Scaling fixture does not match the frozen preparation manifest.",
           call. = FALSE)
    case_id <- paste0("scaling_", mode, "_M", scale_m)
  }
  start_hash <- digest::digest(starts, algo = "sha256")
  log_dir <- file.path(out_dir, "logs"); fit_dir <- file.path(out_dir, "stanfit")
  result_dir <- file.path(out_dir, "results")
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  log_path <- file.path(log_dir, paste0(case_id, ".log"))
  fit_path <- file.path(fit_dir, paste0(case_id, ".rds"))
  result_path <- file.path(result_dir, paste0(case_id, ".csv"))
  started <- proc.time()[["elapsed"]]
  fit <- NULL; fit_error <- ""; warns <- character(); output <- character()
  output <- capture.output({
    value <- tryCatch(withCallingHandlers({
      fn <- getExportedValue("antitrustBayes", "bayes_est")
      do.call(fn, rotation_worker_args(panel, mode, starts, seed))
    }, warning = function(w) {
      warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
    }), error = function(e) {
      fit_error <<- conditionMessage(e); NULL
    })
    fit <- value
  }, type = "output")
  elapsed <- proc.time()[["elapsed"]] - started
  if (!is.null(fit) && !nzchar(fit_error)) saveRDS(fit, fit_path)
  writeLines(c(
    paste0("case_id=", case_id), paste0("mode=", mode), paste0("seed=", seed),
    paste0("start_sha256=", start_hash), output,
    if (length(warns)) paste0("WARNING: ", warns),
    if (nzchar(fit_error)) paste0("ERROR: ", fit_error)
  ), log_path)
  market_map <- if (is.null(scale_m)) mapped$market_map else {
    match(seq_along(panel$mky_n), rep(seq_len(panel$M), each = panel$T))
  }
  d <- rotation_fit_diagnostics(fit, mode, seed, elapsed, start_hash, case_id,
                                fit_error, market_map = market_map)
  d$profile <- if (is.null(scale_m)) "primary" else "scaling"
  d$M <- if (is.null(scale_m)) 2L else as.integer(scale_m)
  d$N <- nrow(panel$data)
  d$G <- length(unique(panel$data$group))
  d$start_sha256 <- start_hash
  d$log_path <- normalizePath(log_path)
  utils::write.csv(d, result_path, row.names = FALSE, na = "")
  invisible(d)
}

rotation_library_env <- function(mode) {
  baseline <- Sys.getenv("ROTATION_BASELINE_RLIB", "/tmp/rotation-baseline-Rlib")
  fresh <- Sys.getenv("ROTATION_FRESH_RLIB", "/tmp/rotation-fresh-Rlib")
  deps <- Sys.getenv("ROTATION_DEP_LIB",
                     "/tmp/cga-Rlib:/home/ctaragin/R/x86_64-pc-linux-gnu-library/4.6:/usr/local/lib/R/site-library:/usr/lib/R/site-library")
  paste0("R_LIBS=", if (identical(mode, "legacy")) baseline else fresh,
         ":", deps)
}

rotation_invoke_worker <- function(mode, seed, out_dir, archive_path,
                                   scale_m = NULL) {
  case_id <- if (is.null(scale_m)) paste0("primary_", mode, "_seed", seed) else
    paste0("scaling_", mode, "_M", as.integer(scale_m))
  args <- c("--vanilla", normalizePath(file.path(
    rotation_script_dir, "market_s0_rotation_hmc_benchmark.R"), mustWork = FALSE),
    "--worker", paste0("--mode=", mode), paste0("--seed=", seed),
    paste0("--out=", out_dir), paste0("--archive=", archive_path))
  if (!is.null(scale_m)) args <- c(args, paste0("--scale-m=", scale_m))
  output <- system2(file.path(R.home("bin"), "Rscript"), args,
                    env = rotation_library_env(mode), stdout = TRUE, stderr = TRUE)
  code <- attr(output, "status", exact = TRUE)
  if (is.null(code)) code <- 0L
  dir.create(file.path(out_dir, "logs"), recursive = TRUE,
             showWarnings = FALSE)
  writeLines(c(paste0("case_id=", case_id), paste0("exit_status=", code),
               if (length(output)) output),
             file.path(out_dir, "logs", paste0(case_id, ".launcher.log")))
  as.integer(code)
}

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

rotation_run <- function(out_dir, archive_path, resume = FALSE) {
  if (!identical(Sys.getenv("ANTITRUSTBAYES_ROTATION_NATIVE_GATE"), "approved"))
    stop("HMC is locked until ANTITRUSTBAYES_ROTATION_NATIVE_GATE=approved.",
         call. = FALSE)
  out_dir <- rotation_require_output(out_dir, resume)
  plan <- rotation_primary_plan()
  if (!file.exists(file.path(out_dir, "benchmark_manifest.csv")))
    rotation_write_prepare(out_dir, archive_path)
  status <- list()
  for (i in seq_len(nrow(plan))) {
    row <- plan[i, , drop = FALSE]
    result_path <- file.path(out_dir, "results", paste0(row$case_id, ".csv"))
    if (resume && file.exists(result_path)) {
      status[[length(status) + 1L]] <- data.frame(case_id = row$case_id,
        worker = "reused", exit_status = 0L, stringsAsFactors = FALSE)
      next
    }
    code <- rotation_invoke_worker(row$mode, row$seed, out_dir, archive_path)
    status[[length(status) + 1L]] <- data.frame(case_id = row$case_id,
      worker = row$mode, exit_status = as.integer(code), stringsAsFactors = FALSE)
  }
  status_df <- do.call(rbind, status)
  utils::write.csv(status_df, file.path(out_dir, "worker_status.csv"), row.names = FALSE)
  result_files <- list.files(file.path(out_dir, "results"), pattern = "^primary_.*\\.csv$",
                             full.names = TRUE)
  results <- if (length(result_files)) do.call(rbind, lapply(result_files,
    utils::read.csv, stringsAsFactors = FALSE)) else data.frame()
  if (nrow(results)) utils::write.csv(results, file.path(out_dir, "primary_results.csv"),
                                      row.names = FALSE, na = "")
  rotated <- results[results$mode == "rotated", , drop = FALSE]
  primary_gate <- nrow(rotated) == 2L && all(rotated$usable)
  if (!primary_gate) {
    reason <- if (!nrow(rotated)) "both rotated primary results missing" else
      paste(rotated$case_id[!rotated$usable], collapse = ";")
    utils::write.csv(data.frame(status = "skipped", reason = reason,
                                condition = "both primary rotated fits usable"),
                     file.path(out_dir, "scaling_status.csv"), row.names = FALSE)
    return(invisible(list(output = out_dir, primary = results,
                          scaling = "skipped", worker_status = status_df)))
  }
  scaling <- rotation_scaling_plan()
  scale_status <- list()
  for (i in seq_len(nrow(scaling))) {
    row <- scaling[i, , drop = FALSE]
    result_path <- file.path(out_dir, "results", paste0(row$case_id, ".csv"))
    if (resume && file.exists(result_path)) next
    code <- rotation_invoke_worker(row$mode, row$seed, out_dir, archive_path,
                                   scale_m = row$M)
    scale_status[[length(scale_status) + 1L]] <- data.frame(
      case_id = row$case_id, worker = row$mode, exit_status = as.integer(code),
      stringsAsFactors = FALSE)
  }
  scale_df <- if (length(scale_status)) do.call(rbind, scale_status) else
    data.frame(case_id = character(), worker = character(), exit_status = integer())
  utils::write.csv(data.frame(status = "executed",
                              condition = "both primary rotated fits usable"),
                   file.path(out_dir, "scaling_status.csv"), row.names = FALSE)
  utils::write.csv(scale_df, file.path(out_dir, "scaling_worker_status.csv"), row.names = FALSE)
  scale_files <- list.files(file.path(out_dir, "results"), pattern = "^scaling_.*\\.csv$",
                            full.names = TRUE)
  if (length(scale_files)) utils::write.csv(rotation_rbind_fill(lapply(scale_files,
    utils::read.csv, stringsAsFactors = FALSE)),
    file.path(out_dir, "scaling_results.csv"), row.names = FALSE, na = "")
  invisible(list(output = out_dir, primary = results, scaling = scale_df,
                 worker_status = status_df))
}

rotation_main <- function() {
  cli <- rotation_cli()
  if (cli$action == "help") {
    cat("Usage: Rscript market_s0_rotation_hmc_benchmark.R [--prepare|--run --out=...]\n",
        "Default is --prepare; --run requires the native gate marker.\n",
        "Primary is four paired fits; scaling is conditional on both rotated gates.\n",
        sep = "")
    return(invisible(NULL))
  }
  archive <- rotation_require_archive(cli$archive)
  if (cli$action == "worker") {
    if (is.null(cli$mode) || is.null(cli$seed) || is.null(cli$out))
      stop("--worker requires --mode, --seed, and --out.", call. = FALSE)
    out <- rotation_require_output(cli$out, resume = TRUE)
    rotation_run_worker(cli$mode, as.integer(cli$seed), out, archive,
                        scale_m = if (is.null(cli$scale_m)) NULL else as.integer(cli$scale_m))
    return(invisible(NULL))
  }
  out <- rotation_require_output(cli$out, cli$resume)
  if (cli$action == "prepare") {
    rotation_write_prepare(out, archive)
    cat("PREPARED: ", normalizePath(out), " (no HMC executed)\n", sep = "")
  } else {
    rotation_run(out, archive, resume = cli$resume)
    cat("RUN COMPLETE: ", normalizePath(out), "\n", sep = "")
  }
  invisible(NULL)
}

if (sys.nframe() == 0L) rotation_main()
