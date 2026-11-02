## Shared, noninteractive helpers for the enterprise market-intercept benchmark.
##
## This file never launches Stan.  It owns the frozen task matrix, deterministic
## fixtures, provenance, task isolation, and result schemas used by the runner,
## aggregator, and Slurm array entry point.

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

.mic_enterprise_common_path <- tryCatch(
  normalizePath(sys.frame(1L)$ofile, mustWork = FALSE),
  error = function(e) NA_character_
)
.mic_enterprise_repo_hint <- if (is.character(.mic_enterprise_common_path) &&
                                   length(.mic_enterprise_common_path) == 1L &&
                                   !is.na(.mic_enterprise_common_path)) {
  normalizePath(file.path(dirname(.mic_enterprise_common_path), "..", "..", ".."),
                mustWork = FALSE)
} else NA_character_

mic_enterprise_root <- function() {
  requested <- Sys.getenv("MIC_ENTERPRISE_REPO", unset = "")
  if (nzchar(requested) && file.exists(file.path(requested, "DESCRIPTION")))
    return(normalizePath(requested, mustWork = TRUE))
  if (!is.na(.mic_enterprise_repo_hint) &&
      file.exists(file.path(.mic_enterprise_repo_hint, "DESCRIPTION")))
    return(.mic_enterprise_repo_hint)
  script <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  here <- if (length(script)) dirname(normalizePath(sub("^--file=", "", script[[1L]]),
                                                     mustWork = FALSE)) else getwd()
  candidate <- normalizePath(file.path(here, "..", ".."), mustWork = FALSE)
  if (file.exists(file.path(candidate, "DESCRIPTION"))) return(candidate)
  normalizePath(getwd(), mustWork = TRUE)
}

mic_env <- function(name, default = "") {
  value <- Sys.getenv(name, unset = default)
  if (!length(value) || is.na(value) || !nzchar(value)) default else value
}

mic_truthy <- function(x) tolower(x) %in% c("true", "t", "1", "yes", "y")

mic_hash <- function(x) {
  if (!requireNamespace("digest", quietly = TRUE))
    stop("The enterprise benchmark requires digest for SHA-256 provenance.",
         call. = FALSE)
  digest::digest(x, algo = "sha256")
}

mic_hash_files <- function(paths) {
  paths <- normalizePath(paths, mustWork = TRUE)
  if (!requireNamespace("digest", quietly = TRUE))
    stop("The enterprise benchmark requires digest for SHA-256 provenance.",
         call. = FALSE)
  hashes <- vapply(paths, digest::digest, character(1), algo = "sha256", file = TRUE)
  paste(paste(basename(paths), hashes, sep = "="), collapse = ";")
}

mic_expected_kernel_hashes <- function() c(
  ces_stackelberg = "89700d43789a4c4109633858f551dc537fd20f7c7a74e50bb107b0c4b7cd9716",
  ces_supply = "4eaf21ee78113b8d934a7c1bd123be22dee01e1a5594aa7da9e863a196d22821",
  logit_stackelberg = "371dd99d72684f7e8f5ead2b51a8bfc505a9e3b5df6c0623d288e386bf7f5074",
  logit_supply = "cb8ba797461801cfba6bbe4e6b2e9d76c554d75de48030f0801160f8a16507fe"
)

mic_atomic_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temp <- tempfile(".mic-enterprise-", tmpdir = dirname(path), fileext = ".csv")
  on.exit(unlink(temp), add = TRUE)
  utils::write.csv(x, temp, row.names = FALSE, na = "")
  if (!file.rename(temp, path)) stop("Could not atomically write ", path, call. = FALSE)
  invisible(path)
}

mic_git_info <- function(path) {
  if (!dir.exists(path)) return(list(sha = NA_character_, dirty = NA, diff = NA_character_))
  sha <- tryCatch(system2("git", c("-C", path, "rev-parse", "HEAD"), stdout = TRUE,
                          stderr = FALSE), error = function(e) character())
  status <- tryCatch(system2("git", c("-C", path, "status", "--porcelain"), stdout = TRUE,
                             stderr = FALSE), error = function(e) character())
  list(sha = if (length(sha)) sha[[1L]] else NA_character_,
       dirty = length(status) > 0L, diff = mic_hash(status))
}

mic_model_contract <- function(model) {
  contracts <- list(
    B  = list(demand_role = "atomic", supply = "bertrand",     structure = "simultaneous"),
    C  = list(demand_role = "atomic", supply = "cournot",      structure = "simultaneous"),
    MC = list(demand_role = "atomic", supply = "monopolistic", structure = "simultaneous"),
    A  = list(demand_role = "atomic", supply = "auction",       structure = "simultaneous"),
    BF = list(demand_role = "rich",   supply = "bertrand",     structure = "fringe"),
    CF = list(demand_role = "rich",   supply = "cournot",      structure = "fringe"),
    BS = list(demand_role = "rich",   supply = "bertrand",     structure = "stackelberg"),
    CS = list(demand_role = "rich",   supply = "cournot",      structure = "stackelberg")
  )
  out <- contracts[[as.character(model)]]
  if (is.null(out)) stop("Unknown benchmark model: ", model, call. = FALSE)
  out
}

mic_conduct_schedule <- function(q, model, demand, s0, gamma = 3.5) {
  q <- pmin(pmax(as.numeric(q), 1e-12), 1 - 1e-12)
  if (demand == "logit") {
    return(switch(model,
      B = 1 - q,
      C = s0 / (s0 + q),
      MC = rep(1, length(q)),
      A = q / -log1p(-q),
      stop("Atomic Logit schedule is undefined for ", model, call. = FALSE)))
  }
  switch(model,
    B = gamma - (gamma - 1) * q,
    C = 1 / (1 / gamma + (gamma - 1) * q /
      (gamma * (1 + (gamma - 1) * s0))),
    MC = rep(gamma, length(q)),
    A = 1 / (1 - (1 - q)^(1 / (gamma - 1))),
    stop("Atomic CES schedule is undefined for ", model, call. = FALSE))
}

mic_fixture_engine <- local({
  cache <- NULL
  function(repo = mic_enterprise_root()) {
    if (!is.null(cache)) return(cache)
    env <- new.env(parent = globalenv())
    old <- Sys.getenv(c("MQ_REPO", "MQ_BENCH_ROOT"), unset = NA_character_)
    on.exit({
      for (name in names(old)) {
        if (is.na(old[[name]])) Sys.unsetenv(name) else do.call(Sys.setenv,
          stats::setNames(list(old[[name]]), name))
      }
    }, add = TRUE)
    Sys.setenv(MQ_REPO = repo, MQ_BENCH_ROOT = file.path(tempdir(), "mic-mq-unused"))
    sys.source(file.path(repo, "doc", "market_quality_benchmark.R"), envir = env)
    cache <<- env
    cache
  }
})

mic_rewrite_atomic_margins <- function(fixture, model) {
  if (!model %in% c("MC", "A")) return(fixture)
  dat <- fixture$data
  M <- fixture$M; T <- fixture$T; F <- fixture$F
  old_noise <- numeric(nrow(dat))
  new_inv <- numeric(nrow(dat))
  for (m in seq_len(M)) for (tt in seq_len(T)) {
    jj <- (m - 1L) * T + tt
    ii <- (jj - 1L) * F + seq_len(F)
    p <- dat$prices[ii]
    r <- dat$quantities[ii] / sum(dat$quantities[ii])
    if (fixture$demand == "logit") {
      q <- (1 - fixture$truth$s0[[jj]]) * r
      old_inv <- p * fixture$truth$alpha * (1 - q)
      new_inv[ii] <- p * fixture$truth$alpha *
        mic_conduct_schedule(q, model, "logit", fixture$truth$s0[[jj]])
    } else {
      rev <- p * r / sum(p * r)
      R0 <- stats::plogis(stats::qlogis(fixture$truth$s0[[jj]]) +
                            log(fixture$outside_price) - log(sum(p * r)))
      R <- (1 - R0) * rev
      old_inv <- mic_conduct_schedule(R, "B", "ces", R0, fixture$truth$gamma)
      new_inv[ii] <- mic_conduct_schedule(R, model, "ces", R0,
                                           fixture$truth$gamma)
    }
    old_noise[ii] <- log(1 / dat$margins[ii]) - log(old_inv)
  }
  dat$margins <- exp(-(log(new_inv) + old_noise))
  fixture$data <- dat
  fixture$conduct <- if (model == "MC") "monopolistic" else "auction"
  fixture$structure <- "simultaneous"
  fixture$truth$inverse_margin <- new_inv
  fixture$truth$margin_noise <- old_noise
  fixture
}

mic_rebuild_margins <- function(fixture, model, rho = .30, sigma_margin = .06,
                                seed = 8675309L) {
  dat <- fixture$data
  M <- fixture$M; T <- fixture$T; F <- fixture$F
  if (nrow(dat) != M * T * F)
    stop("Margin rebuilding requires a balanced synthetic fixture.", call. = FALSE)
  spec <- mic_model_contract(model)
  engine <- mic_fixture_engine()
  inverse <- numeric(nrow(dat)); share_residual <- numeric(nrow(dat))
  rank <- fixture$truth$fixed_rank
  for (m in seq_len(M)) for (tt in seq_len(T)) {
    jj <- (m - 1L) * T + tt
    ii <- (jj - 1L) * F + seq_len(F)
    p <- dat$prices[ii]; r <- dat$quantities[ii] / sum(dat$quantities[ii])
    s0 <- fixture$truth$s0[[jj]]
    share_residual[ii] <- log(r) - fixture$truth$latent_log_share[ii]
    if (fixture$demand == "logit") {
      if (spec$structure == "simultaneous" && model %in% c("B", "C", "MC", "A")) {
        q <- (1 - s0) * r
        conduct <- mic_conduct_schedule(q, model, "logit", s0)
      } else {
        conduct <- engine$mq_foc$cga_role_conduct_logit(
          r, s0, spec$supply, spec$structure, rank[m, ], 2L)
      }
      inverse[ii] <- p * fixture$truth$alpha * conduct
    } else {
      rev <- p * r / sum(p * r)
      R0 <- stats::plogis(stats::qlogis(s0) + log(fixture$outside_price) -
                            log(sum(p * r)))
      if (spec$structure == "simultaneous" && model %in% c("B", "C", "MC", "A")) {
        R <- (1 - R0) * rev
        inverse[ii] <- mic_conduct_schedule(R, model, "ces", R0,
                                             fixture$truth$gamma)
      } else {
        inverse[ii] <- engine$mq_foc$cga_role_inverse_lerner_ces(
          rev, R0, fixture$truth$gamma, spec$supply, spec$structure,
          rank[m, ], 2L)
      }
    }
  }
  sigma_share <- stats::sd(share_residual)
  if (!is.finite(sigma_share) || sigma_share <= 1e-10)
    stop("Share residual scale is degenerate in stress fixture.", call. = FALSE)
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(as.integer(seed))
  margin_noise <- rho * sigma_margin / sigma_share * share_residual +
    stats::rnorm(nrow(dat), 0, sigma_margin * sqrt(1 - rho^2))
  dat$margins <- exp(-(log(inverse) + margin_noise))
  fixture$data <- dat
  fixture$truth$inverse_margin <- inverse
  fixture$truth$margin_noise <- margin_noise
  fixture$truth$share_residual_sd <- sigma_share
  fixture$truth$rho <- rho
  fixture$truth$sigma_margin <- sigma_margin
  fixture
}

mic_apply_stress <- function(fixture, scenario) {
  scenario <- as.character(scenario)
  if (scenario %in% c("standard", "compositional")) return(fixture)
  dat <- fixture$data
  if (scenario == "partial_margins") {
    first_market <- unique(dat$market)[[1L]]
    dat$margin_observed <- dat$market != first_market & seq_len(nrow(dat)) %% 3L != 0L
    dat$margins[!dat$margin_observed] <- NA_real_
  } else if (scenario == "unbalanced") {
    keep <- !(dat$market == unique(dat$market)[[1L]] &
      as.integer(sub("^firm", "", dat$group)) > max(2L, fixture$F - 2L))
    dat <- dat[keep, , drop = FALSE]
  } else if (scenario == "low_share_noise") {
    ## Shrink the normalized compositional innovation toward its latent share.
    ## Persistent quality then dominates transitory share noise without
    ## altering the quality prior or conduct equation.
    cell <- interaction(dat$market, dat$year, drop = TRUE, lex.order = TRUE)
    for (lev in levels(cell)) {
      ii <- which(cell == lev)
      r <- dat$quantities[ii] / sum(dat$quantities[ii])
      latent <- fixture$truth$latent_log_share[ii]
      score <- latent + .20 * (log(r) - latent)
      score <- score - max(score)
      r_new <- exp(score) / sum(exp(score))
      dat$quantities[ii] <- 1000 * r_new
    }
  } else if (scenario == "high_rho") {
    ## Shares are unchanged; margins are rebuilt below with rho=.9.
  } else {
    stop("Unknown synthetic scenario: ", scenario, call. = FALSE)
  }
  fixture$data <- dat
  fixture$scenario <- scenario
  fixture
}

mic_apply_market_intercept_truth <- function(fixture, sigma_b_market = 0,
                                             seed = 272729L) {
  sigma_b_market <- as.numeric(sigma_b_market)
  if (length(sigma_b_market) != 1L || !is.finite(sigma_b_market) ||
      sigma_b_market < 0)
    stop("sigma_b_market must be one finite nonnegative number.", call. = FALSE)
  if (fixture$M < 2L && sigma_b_market > 0)
    stop("A nonzero market-intercept dispersion requires multiple markets.",
         call. = FALSE)

  b_market <- numeric(fixture$M)
  if (sigma_b_market > 0) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(as.integer(seed))
    raw <- stats::rnorm(fixture$M)
    raw <- raw - mean(raw)
    if (!is.finite(stats::sd(raw)) || stats::sd(raw) <= 1e-12)
      stop("Could not generate a nondegenerate market-intercept fixture.",
           call. = FALSE)
    b_market <- sigma_b_market * raw / stats::sd(raw)
  }

  ## The observed outcomes are conditional inside shares.  A common market
  ## intercept therefore cancels from relative inside shares and instead
  ## changes the outside share.  If eta_old generated the original fixture,
  ## eta_new = eta_old - b_m makes eta_new + b_m exactly equal eta_old.
  old_eta <- stats::qlogis(fixture$truth$s0)
  b_mky <- rep(b_market, each = fixture$T)
  fixture$truth$s0 <- stats::plogis(old_eta - b_mky)
  fixture$truth$market_intercept <- b_market
  fixture$truth$mu_b_market <- mean(b_market)
  fixture$truth$sigma_b_market <- if (length(b_market) > 1L)
    stats::sd(b_market) else 0
  fixture$truth$market_intercept_identity_error <- max(abs(
    stats::qlogis(fixture$truth$s0) + b_mky - old_eta
  ))
  fixture
}

mic_make_synthetic <- function(model, demand = "logit", M = 10L, T = 1L,
                               F = 6L, seed = 271828L,
                               scenario = "standard",
                               sigma_b_market = 0) {
  spec <- mic_model_contract(model)
  engine <- mic_fixture_engine()
  base_supply <- if (spec$supply %in% c("bertrand", "cournot")) spec$supply else "bertrand"
  fixture <- engine$mq_fixture(M = as.integer(M), T = as.integer(T), F = as.integer(F),
    demand = demand, conduct = base_supply, structure = spec$structure,
    seed = as.integer(seed))
  fixture <- mic_rewrite_atomic_margins(fixture, model)
  fixture <- mic_apply_stress(fixture, scenario)
  fixture <- mic_apply_market_intercept_truth(
    fixture, sigma_b_market = sigma_b_market,
    seed = as.integer(seed) + 901L
  )
  if (sigma_b_market > 0 || scenario %in% c("low_share_noise", "high_rho")) {
    fixture <- mic_rebuild_margins(fixture, model,
      rho = if (scenario == "high_rho") .90 else .30,
      seed = as.integer(seed) + 7301L)
  }
  fixture$model <- model
  fixture$scenario <- scenario
  fixture$source <- paste0("synthetic:", model, ":", demand, ":M", M, ":T", T,
                           ":F", F, ":", scenario, ":sigma_b", sigma_b_market,
                           ":seed", seed)
  fixture$source_hash <- mic_hash(fixture$data)
  fixture
}

mic_sc_data_dir <- function() {
  requested <- mic_env("MIC_ENTERPRISE_SC_DATA")
  candidates <- unique(path.expand(c(requested,
    "~/Projects/mergerBayes/data", "~/Projects/mergerbayes/data",
    file.path(mic_enterprise_root(), "..", "mergerBayes", "data"))))
  candidates <- candidates[nzchar(candidates)]
  hit <- candidates[file.exists(file.path(candidates, "supreme_data.RData"))]
  if (!length(hit)) stop("Could not find mergerBayes Supreme Court data.", call. = FALSE)
  hit[[1L]]
}

mic_load_sc <- function(fixture = c("sc_cross_section", "sc_panel"), data_dir = mic_sc_data_dir()) {
  fixture <- match.arg(fixture)
  data_file <- file.path(data_dir, "supreme_data.RData")
  event_file <- file.path(data_dir, "event_overlaps.RData")
  de <- new.env(parent = baseenv()); ee <- new.env(parent = baseenv())
  load(data_file, envir = de); load(event_file, envir = ee)
  if (!exists("simdata", de, inherits = FALSE) || !exists("events", ee, inherits = FALSE))
    stop("SC inputs must define simdata and events.", call. = FALSE)
  x <- de$simdata; events <- ee$events
  required <- c("event_id", "year", "tophold", "rate_deposits", "total_deposits",
                "selected_input_markdown_ratio")
  missing <- setdiff(required, names(x))
  if (length(missing)) stop("SC data lacks: ", paste(missing, collapse = ", "), call. = FALSE)
  event_id <- as.character(events$event_id); base_year <- as.integer(events$year)
  pos <- match(as.character(x$event_id), event_id)
  year <- as.integer(as.character(x$year))
  use_window <- if (fixture == "sc_cross_section") year == base_year[pos] else
    year >= base_year[pos] - 1L & year <= base_year[pos] + 2L
  p <- as.numeric(x$rate_deposits); q <- as.numeric(x$total_deposits)
  margin <- as.numeric(x$selected_input_markdown_ratio)
  keep <- !is.na(pos) & !is.na(use_window) & use_window & is.finite(p) & p > 0 &
    is.finite(q) & q > 0 & is.finite(margin) & margin > 0 & !is.na(x$tophold)
  dat <- data.frame(prices = p[keep], quantities = q[keep], margins = margin[keep],
    market = as.character(x$event_id[keep]), year = as.character(year[keep]),
    group = as.character(x$tophold[keep]), margin_observed = TRUE,
    stringsAsFactors = FALSE)
  dat <- dat[order(as.integer(dat$market), as.integer(dat$year), -dat$quantities),,
             drop = FALSE]
  rownames(dat) <- NULL
  if (fixture == "sc_cross_section" &&
      (nrow(dat) != 176L || length(unique(dat$market)) != 6L)) {
    stop("Canonical SC cross-section must contain 176 rows in 6 markets; found ",
         nrow(dat), " rows in ", length(unique(dat$market)), " markets.", call. = FALSE)
  }
  if (fixture == "sc_panel" &&
      (nrow(dat) != 633L || length(unique(dat$market)) != 6L ||
       length(unique(interaction(dat$market, dat$year, drop = TRUE))) != 22L)) {
    stop("Canonical SC panel must contain 633 rows in 6 markets and 22 ",
         "market-years; found ", nrow(dat), " rows in ",
         length(unique(dat$market)), " markets and ",
         length(unique(interaction(dat$market, dat$year, drop = TRUE))),
         " market-years.", call. = FALSE)
  }
  list(data = dat, source = paste(data_file, event_file, sep = ";"),
       source_hash = mic_hash_files(c(data_file, event_file)), model = NA_character_,
       demand = "logit", outside_price = 1, scenario = fixture)
}

mic_make_manifest <- function() {
  rows <- list(); add <- function(...) rows[[length(rows) + 1L]] <<- data.frame(...,
    stringsAsFactors = FALSE)
  models <- c("B", "C", "MC", "A", "BF", "CF", "BS", "CS")
  for (demand in c("logit", "ces")) for (T in c(1L, 4L)) for (model in models) {
    add(tier = "gate", source = "synthetic", fixture = "compositional", demand = demand,
        model = model, intercept = "collapsed", M = 10L, T = T, F = 6L,
        dgp_seed = 271828L, sampler_seed = 19091L, scenario = "compositional")
  }
  gate_copy <- do.call(rbind, rows)
  for (i in seq_len(nrow(gate_copy))) {
    x <- gate_copy[i, , drop = FALSE]
    add(tier = "replicate", source = x$source, fixture = x$fixture,
        demand = x$demand, model = x$model, intercept = x$intercept,
        M = x$M, T = x$T, F = x$F, dgp_seed = x$dgp_seed,
        sampler_seed = 29091L, scenario = x$scenario)
  }
  for (demand in c("logit", "ces")) for (T in c(1L, 4L))
    for (model in c("B", "C", "BS", "CS")) {
      add(tier = "comparison", source = "synthetic", fixture = "compositional",
          demand = demand, model = model, intercept = "sampled", M = 10L, T = T,
          F = 6L, dgp_seed = 271828L, sampler_seed = 19091L,
          scenario = "compositional")
    }
  for (fixture in c("sc_cross_section", "sc_panel")) for (model in models) {
    add(tier = "sc", source = "supreme_court", fixture = fixture, demand = "logit",
        model = model, intercept = "collapsed", M = NA_integer_, T = NA_integer_,
        F = NA_integer_, dgp_seed = NA_integer_, sampler_seed = 19091L,
        scenario = "observed")
  }
  for (fixture in c("sc_cross_section", "sc_panel")) for (model in c("B", "C")) {
    add(tier = "sc_comparison", source = "supreme_court", fixture = fixture,
        demand = "logit", model = model, intercept = "sampled", M = NA_integer_,
        T = NA_integer_, F = NA_integer_, dgp_seed = NA_integer_,
        sampler_seed = 19091L, scenario = "observed")
  }
  stress <- data.frame(
    demand = rep(c("logit", "ces"), 4L),
    model = c("B", "C", "BS", "CS", "BF", "CF", "B", "BS"),
    scenario = c("partial_margins", "partial_margins", "low_share_noise",
                 "low_share_noise", "high_rho", "high_rho",
                 "unbalanced", "unbalanced"),
    stringsAsFactors = FALSE)
  for (i in seq_len(nrow(stress))) add(tier = "stress", source = "synthetic",
    fixture = stress$scenario[[i]], demand = stress$demand[[i]], model = stress$model[[i]],
    intercept = "collapsed", M = 10L, T = 4L, F = 6L, dgp_seed = 314159L,
    sampler_seed = 29091L, scenario = stress$scenario[[i]])
  ## Optional timing probes.  Their short, single-chain budget measures
  ## scaling only and is deliberately excluded from posterior usability gates.
  for (demand in c("logit", "ces")) for (model in c("B", "BS"))
    for (M in c(10L, 100L, 500L, 2000L)) {
      add(tier = "scale", source = "synthetic", fixture = "compositional",
          demand = demand, model = model, intercept = "collapsed", M = M,
          T = 2L, F = 4L, dgp_seed = 161803L, sampler_seed = 39091L,
          scenario = "compositional")
    }
  out <- do.call(rbind, rows)
  out$task_id <- seq_len(nrow(out))
  out$chains <- 4L; out$iter <- 2000L; out$warmup <- 1000L; out$cores <- 4L
  out$adapt_delta <- .99; out$max_treedepth <- 13L
  out$quality <- "collapsed_market_total_fraction"
  out$true_k <- ifelse(out$source == "synthetic" &
                         out$model %in% c("BF", "CF", "BS", "CS"), 2L, NA_integer_)
  out$true_role_scope <- ifelse(is.na(out$true_k), "none", "market")
  scale <- out$tier == "scale"
  out$chains[scale] <- 1L; out$iter[scale] <- 400L
  out$warmup[scale] <- 200L; out$cores[scale] <- 1L
  out$benchmark_kind <- ifelse(scale, "timing", "hmc_gate")
  out$rho <- TRUE; out$s0_pooling <- "market"
  out$s0_coordinates <- ifelse(out$intercept == "sampled", "joint", "legacy")
  out$run_key <- vapply(seq_len(nrow(out)), function(i) mic_hash(as.list(out[i, ])),
                        character(1))
  out <- out[, c("task_id", "run_key", setdiff(names(out), c("task_id", "run_key")))]
  attr(out, "manifest_hash") <- mic_hash(out)
  out
}

mic_validate_manifest <- function(x) {
  required <- c("task_id", "run_key", "tier", "source", "fixture", "demand", "model",
                "intercept", "sampler_seed", "chains", "iter", "warmup", "cores",
                "benchmark_kind")
  if (!all(required %in% names(x))) stop("Manifest lacks required columns.", call. = FALSE)
  if (!identical(as.integer(x$task_id), seq_len(nrow(x))) || anyDuplicated(x$run_key))
    stop("Manifest task IDs or run keys are not unique and ordered.", call. = FALSE)
  if (any(x$warmup >= x$iter) || any(x$cores > x$chains))
    stop("Manifest has invalid sampler resources.", call. = FALSE)
  invisible(TRUE)
}

mic_manifest_hash <- function(x) {
  x$manifest_hash <- NULL
  text <- paste(capture.output(utils::write.csv(x, row.names = FALSE, na = "")),
                collapse = "\n")
  mic_hash(text)
}

mic_task_result_path <- function(run_dir, task_id, attempt = 1L) {
  file.path(run_dir, "tasks", sprintf("task-%04d", as.integer(task_id)),
            sprintf("attempt-%03d.csv", as.integer(attempt)))
}

mic_load_fixture <- function(task) {
  if (task$source == "synthetic") {
    mic_make_synthetic(task$model, task$demand, task$M, task$T, task$F,
      task$dgp_seed, task$scenario)
  } else if (task$source == "supreme_court") {
    mic_load_sc(task$fixture)
  } else stop("Unknown fixture source: ", task$source, call. = FALSE)
}

mic_fit_args <- function(task, fixture) {
  spec <- mic_model_contract(task$model)
  dat <- fixture$data
  args <- list(prices = dat$prices, quantities = dat$quantities, margins = dat$margins,
    market = dat$market, year = dat$year, group = dat$group,
    margin_observed = dat$margin_observed,
    demand_model = task$demand, supply_model = spec$supply, structure = spec$structure,
    market_side = if (task$source == "supreme_court") "input" else "output",
    likelihood = "joint", use_rho = TRUE,
    share_likelihood = "gaussian", margin_likelihood = "gaussian",
    quality_parameterization = "collapsed", quality_scope = "market",
    quality_scale_parameterization = "total_fraction",
    market_intercept_parameterization = task$intercept,
    s0_pooling = "market", market_s0_parameterization = task$s0_coordinates,
    use_centered_s0_pool = FALSE, use_alpha_market_pool = FALSE,
    use_market_resid_scale = FALSE, prior_manifest = "baseline",
    outside_price = if (task$demand == "ces") fixture$outside_price else 1,
    stan_options = list(chains = task$chains, iter = task$iter, warmup = task$warmup,
      cores = task$cores, seed = task$sampler_seed, refresh = 0, init = 0,
      control = list(adapt_delta = task$adapt_delta,
                     max_treedepth = task$max_treedepth)))
  if (spec$structure != "simultaneous") {
    ## k counts firms, not observation rows.  Repeated periods and
    ## multiproduct rows must never inflate the admissible support.
    Fmin <- min(vapply(split(as.character(dat$group), as.character(dat$market)),
                       function(x) length(unique(x)), integer(1)))
    if (Fmin < 2L) stop("Ordered-k fixtures require at least two firms per market.")
    args$partition_prior <- "ordered_k"
    args$k_support <- if (spec$structure == "fringe") 0:(Fmin - 1L) else
      seq_len(Fmin - 1L)
    args$k_prior <- rep(1 / length(args$k_support), length(args$k_support))
    args$k_scope <- "market"
  }
  args
}

mic_rank_hash <- function(fixture) {
  rank <- fixture$truth$fixed_rank %||% NULL
  if (is.null(rank)) return(NA_character_)
  mic_hash(rank)
}

mic_smoke_task <- function(task) {
  fixture <- mic_load_fixture(task)
  dat <- fixture$data
  required <- c("prices", "quantities", "margins", "market", "year", "group",
                "margin_observed")
  if (!all(required %in% names(dat)) || !nrow(dat)) stop("Fixture contract failed.")
  if (any(!is.finite(dat$prices) | dat$prices <= 0) ||
      any(!is.finite(dat$quantities) | dat$quantities <= 0)) stop("Invalid price/quantity.")
  if (any(dat$margin_observed & (!is.finite(dat$margins) | dat$margins <= 0)))
    stop("Invalid observed margin.")
  args <- mic_fit_args(task, fixture)
  if (task$intercept == "collapsed" && args$market_s0_parameterization != "legacy")
    stop("Collapsed intercept must use the nonredundant legacy s0 coordinates.")
  if (task$intercept == "sampled" && args$market_s0_parameterization != "joint")
    stop("Sampled comparison must use the validated joint s0 coordinates.")
  data.frame(task_id = task$task_id, run_key = task$run_key, source = task$source,
    fixture = task$fixture, demand = task$demand, model = task$model,
    intercept = task$intercept, N = nrow(dat), N_market = length(unique(dat$market)),
    N_period = length(unique(interaction(dat$market, dat$year, drop = TRUE))),
    observed_margins = sum(dat$margin_observed), fixture_hash = fixture$source_hash,
    true_k = task$true_k, true_role_scope = task$true_role_scope,
    fixed_rank_hash = mic_rank_hash(fixture),
    status = "smoke_passed", stringsAsFactors = FALSE)
}
