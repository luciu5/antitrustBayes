#!/usr/bin/env Rscript
# Sourcing defines functions only. --prepare seals fresh fixtures; --run samples.
# Old Nutpie/CmdStan fits target specification11 and are never reused here.

smi_hash <- function(path) digest::digest(file = path, algo = "sha256")

smi_geometry <- function(r, s0, prices, alpha) {
  stopifnot(length(r) == length(prices), all(r > 0), abs(sum(r) - 1) < 1e-12,
    s0 > 0, s0 < 1, alpha > 0)
  location <- mean(log(r)) - stats::qlogis(s0)
  centered_quality <- log(r) - mean(log(r)) + alpha * (prices - mean(prices))
  delta <- location + centered_quality - alpha * (prices - mean(prices))
  unconditional <- exp(delta) / (1 + sum(exp(delta)))
  list(intercept = location, quality = centered_quality, delta = delta,
    shares = unconditional, outside_share = 1 - sum(unconditional),
    conditional_shares = unconditional / sum(unconditional),
    old_implied_outside = stats::plogis(mean(log(r))))
}

smi_margin_profile <- function(conduct, r, inverse_markup, s0) {
  multiplier <- if (conduct == "MC") rep(1, length(r)) else if (conduct == "C")
    s0 / (s0 + (1 - s0) * r) else stop("Only the MC/C identification identity is used.")
  alpha <- sum(multiplier * inverse_markup) / sum(multiplier^2)
  c(alpha = alpha, squared_error = sum((inverse_markup - alpha * multiplier)^2))
}

smi_mc_limit <- function(r, prior_eta_mean, prior_eta_sd) {
  # Zero-share-noise conditional level likelihood c = a + g(eta), a ~ N(0,1).
  # MC margins contain alpha, but no s0. This conditional limiting calculation
  # diagnoses prior sensitivity; it is not a fitted posterior/recovery claim.
  effective <- function(eta) stats::qlogis(1e-6 + (1 - 2e-6) * stats::plogis(eta))
  density <- function(eta) exp(stats::dnorm(mean(log(r)) - effective(eta), log = TRUE) +
    stats::dnorm(eta, prior_eta_mean, prior_eta_sd, log = TRUE))
  moments <- vapply(0:2, function(k) stats::integrate(function(eta) eta^k * density(eta),
    -Inf, Inf, rel.tol = 1e-10)$value, numeric(1))
  c(eta_mean = moments[2L] / moments[1L],
    eta_sd = sqrt(moments[3L] / moments[1L] - (moments[2L] / moments[1L])^2))
}

smi_design <- function() {
  fixtures <- data.frame(id = c("C_M1T1", "C_M1T4", "C_M4T1", "C_M4T4", "CF_M1T1", "CS_M1T1"),
    conduct = c(rep("C", 4L), "CF", "CS"), markets = c(1L, 1L, 4L, 4L, 1L, 1L),
    periods = c(1L, 4L, 1L, 4L, 1L, 1L), stringsAsFactors = FALSE)
  tasks <- list()
  add <- function(...) tasks[[length(tasks) + 1L]] <<- list(...)
  for (id in fixtures$id) for (backend in c("cmdstanr", "nutpie"))
    add(id = paste("fit", id, backend, sep = "_"), stage = "fit", fixture = id, backend = backend)
  for (id in fixtures$id) add(id = paste0("pair_", id), stage = "pair", fixture = id)
  for (id in c("C_M1T1", "CF_M1T1", "CS_M1T1")) for (backend in c("cmdstanr", "nutpie"))
    add(id = paste("bridge", id, backend, sep = "_"), stage = "bridge", fixture = id, backend = backend)
  add(id = "MC_intercept_prior", stage = "prior", fixture = "C_M1T1")
  list(fixtures = fixtures, tasks = tasks, conditional_shares = c(.25, .20, .16, .13, .11, .09, .06),
    outside_share = .2, k = 3L, alpha = 2.4, costs = c(.70, .82, .96, 1.10, 1.27, 1.45, 1.64),
    controls = list(chains = 4L, iter = 1000L, warmup = 500L, cores = 2L, init = 0,
      refresh = 250L, control = list(adapt_delta = .95, max_treedepth = 13L)),
    fit_seed = 452001L, bridge_seed = 452701L, bridge_repetitions = 3L,
    scope = "Specification12 fresh-target checks; no old-fit reuse, retuning, failed-state replacement or production grid.")
}

smi_fixture <- function(row, design) {
  rows <- truths <- gates <- list(); common <- NULL
  for (m in seq_len(row$markets)) {
    arguments <- list(conduct = row$conduct, J = 7L, mode = "parameter",
      alpha = design$alpha * exp(.06 * (m - 1)), costs = design$costs,
      conditional_shares = design$conditional_shares, outside_share = design$outside_share,
      seed = 451901L + m, inside_size = 1000)
    if (row$conduct %in% c("CF", "CS")) arguments$k <- design$k
    base <- do.call(antitrustBayes::bayes_calibrated_market, arguments)
    if (is.null(common)) common <- base$estimation$inputs
    for (tt in seq_len(row$periods)) {
      truth <- if (tt == 1L) base else antitrustBayes::bayes_resolve_market(base,
        costs = base$truth$costs * exp(.035 * sin(seq_len(7L) + tt)))
      d <- truth$truth$data; d$market <- paste0("m", m); d$year <- paste0("t", tt)
      # Small fixed observation errors, identical in both sampler routes.
      q <- d$quantities * exp(.01 * sin(seq_len(7L) + 2L * m + tt))
      d$quantities <- 1000 * q / sum(q)
      d$margins <- d$margins * exp(.01 * cos(seq_len(7L) + m + 2L * tt))
      stopifnot(all(d$margins > 0 & d$margins < 1), isTRUE(truth$fresh_oracle$passed))
      rows[[length(rows) + 1L]] <- d[c("prices", "quantities", "margins", "market", "year", "group")]
      truths[[length(truths) + 1L]] <- data.frame(market = d$market, year = d$year, product = d$group,
        alpha = truth$truth$alpha, s0 = truth$truth$s0, price = truth$truth$data$prices,
        share = truth$truth$shares, margin = truth$truth$margins, cost = truth$truth$costs,
        meanval = truth$truth$meanval)
      gates[[length(gates) + 1L]] <- data.frame(market = m, period = tt,
        price_error = truth$fresh_oracle$max_relative_error, share_error = truth$fresh_oracle$max_share_error,
        cost_error = truth$fresh_oracle$max_cost_error, foc_error = truth$fresh_oracle$max_foc_residual,
        passed = truth$fresh_oracle$passed)
    }
  }
  data <- do.call(rbind, rows); truth <- do.call(rbind, truths)
  common[names(data)] <- as.list(data)
  common$margin_observed <- rep(TRUE, nrow(data)); common$use_rho <- TRUE
  common$share_likelihood <- common$margin_likelihood <- "gaussian"
  common$quality_parameterization <- "collapsed"; common$quality_scope <- "market"
  common$market_intercept_parameterization <- "sampled"
  common$use_alpha_market_pool <- TRUE; common$s0_pooling <- "none"
  common$priors <- antitrustBayes::bayes_prior_manifest(); common$inference <- "hmc"
  if (row$conduct %in% c("CF", "CS")) { common$k_support <- design$k; common$k_prior <- 1 }
  truth$intercept <- ave(truth$meanval, truth$market, FUN = mean) - truth$alpha * mean(data$prices)
  list(id = row$id, arguments = common, truth = truth, equilibrium = do.call(rbind, gates))
}

smi_runtime <- function(repository) {
  files <- c(sort(list.files(file.path(repository, "R"), "[.]R$", full.names = TRUE)),
    file.path(repository, "DESCRIPTION"), file.path(repository, "doc/validation",
      c("single_market_intercept_validation.R", "nutpie_backend_validation.R", "nutpie_evidence_validation.R", "nutpie_closeout_validation.R")),
    vapply(c("logit_supply", "logit_stackelberg"), antitrustBayes:::.cmdstan_kernel_file, character(1)))
  list(files = stats::setNames(vapply(files, smi_hash, character(1)), files),
    cmdstan = list(path = cmdstanr::cmdstan_path(), version = as.character(cmdstanr::cmdstan_version())),
    libraries = .libPaths(), packages = nutpie_validation_runtime(repository)$packages,
    r_version = R.version.string, bridgestan_source = Sys.getenv("BRIDGESTAN"),
    caches = list(cmdstan = getOption("antitrustBayes.cmdstan_cache"), nutpie = getOption("antitrustBayes.nutpie_cache")))
}

smi_assert <- function(manifest) {
  if (!identical(manifest$runtime, smi_runtime(manifest$repository)) ||
      !identical(manifest$design, smi_design()) ||
      !identical(manifest$fixtures_hash, digest::digest(manifest$fixtures, algo = "sha256")))
    stop("Fresh-target runtime, source, design or fixtures changed after sealing; retain this run.")
  if (!all(vapply(names(manifest$legacy$hashes), function(path)
    identical(smi_hash(path), manifest$legacy$hashes[[path]]), logical(1))))
    stop("The read-only specification11 prior-normalization reference changed.")
  invisible(TRUE)
}

smi_legacy_reference <- function(repository) {
  manifest_path <- file.path(repository, ".scratch/nutpie-closeout/manifest.rds")
  old <- readRDS(manifest_path); path <- old$original_paths[["MC_nutpie"]]
  if (!identical(smi_hash(path), old$original_hashes[[path]])) stop("Original MC fit seal changed.")
  native <- readRDS(path)@stanfit$metadata
  diagnostic <- file.path(repository, ".scratch/single-market-intercept/prechange-diagnostic.json")
  expected <- jsonlite::read_json(diagnostic)$old_source_sha256$logit_supply
  source <- file.path(dirname(native$lib_path), basename(native$source_file))
  stopifnot(identical(native$source_hash, expected), identical(smi_hash(source), expected),
    identical(smi_hash(native$lib_path), native$lib_hash))
  files <- c(manifest_path, diagnostic, path, source, native$lib_path)
  list(metadata = native, hashes = stats::setNames(vapply(files, smi_hash, character(1)), files),
    scope = "Verified specification11 compiled target only; no legacy posterior draws or scientific fit are reused.")
}

smi_prepare <- function(repository, root) {
  path <- file.path(root, "manifest.rds")
  if (file.exists(path)) { m <- readRDS(path); smi_assert(m); return(m) }
  design <- smi_design()
  if (!"outside_share" %in% names(formals(antitrustBayes::bayes_calibrated_market)))
    stop("The explicit-outside-share fixture API is required before preparation.")
  fixtures <- lapply(seq_len(nrow(design$fixtures)), function(i) {
    warnings <- character()
    fixture <- withCallingHandlers(smi_fixture(design$fixtures[i, ], design), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    })
    fixture$warnings <- warnings
    fixture
  })
  names(fixtures) <- design$fixtures$id
  m <- list(repository = repository, root = root, prepared_at = Sys.time(), design = design, fixtures = fixtures,
    fixtures_hash = digest::digest(fixtures, algo = "sha256"), runtime = smi_runtime(repository),
    legacy = smi_legacy_reference(repository),
    starting_git_commit = trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)))
  saveRDS(m, path)
  m
}

smi_fit_path <- function(manifest, fixture, backend) file.path(manifest$root,
  paste("fit", fixture, backend, sep = "_"), "fit.rds")

smi_fit <- function(manifest, fixture, backend) {
  path <- smi_fit_path(manifest, fixture, backend)
  seal <- readRDS(file.path(dirname(path), "fit-seal.rds"))
  if (!identical(seal$hash, smi_hash(path))) stop("Saved fresh-target fit changed.")
  fit <- readRDS(path)
  if (!identical(seal$identity, nutpie_validation_identity(fit))) stop("Saved fresh-target fit identity changed.")
  if (!identical(fit@metadata$stan_source_hash, smi_hash(antitrustBayes:::.cmdstan_kernel_file(fit@metadata$stan_kernel))))
    stop("Saved fit targets an older Stan source.")
  fit
}

smi_prior_check <- function(fit, directory, legacy) {
  # Auxiliary MC data only: no fitted object or scientific result is replaced.
  data <- antitrustBayes:::.bayes_bridge_data(fit); data$supply_model <- 2L; data$use_rho <- 0L
  data$likelihood_mode <- 3L
  json <- file.path(directory, "MC-margin-only-data.json"); cmdstanr::write_stan_json(data, json)
  native <- fit@stanfit$metadata
  loaded <- antitrustBayes:::.bayes_nutpie_bridge_load(native$lib_path, json, 453191L,
    as.character(utils::packageVersion("bridgestan")), native$lib_hash)
  model <- loaded$model
  names <- as.character(unlist(model$param_unc_names(), use.names = FALSE))
  index <- which(names %in% c("mu_b_market_vec.1", "mu_b_market_vec[1]"))
  stopifnot(length(index) == 1L, !any(grepl("sigma_b_market", names, fixed = TRUE)))
  q <- rep(0, length(names)); origin <- model$log_density(q, propto = FALSE, jacobian = TRUE)
  points <- c(-1.1, -.4, 0, .7, 1.3)
  observed <- vapply(points, function(a) { changed <- q; changed[index] <- a
    model$log_density(changed, propto = FALSE, jacobian = TRUE) - origin }, numeric(1))
  expected <- stats::dnorm(points, log = TRUE) - stats::dnorm(0, log = TRUE)
  stopifnot(all(is.finite(observed)), max(abs(observed - expected)) < 1e-8)
  old <- antitrustBayes:::.bayes_nutpie_bridge_load(legacy$metadata$lib_path, json, 453191L,
    as.character(utils::packageVersion("bridgestan")), legacy$metadata$lib_hash)$model
  old_names <- as.character(unlist(old$param_unc_names(), use.names = FALSE))
  stopifnot(identical(names[-index], old_names), length(old_names) + 1L == length(names))
  old_density <- old$log_density(rep(0, length(old_names)), propto = FALSE, jacobian = TRUE)
  absolute_difference <- observed + origin - old_density
  absolute_expected <- stats::dnorm(points, log = TRUE)
  stopifnot(all(is.finite(absolute_difference)), max(abs(absolute_difference - absolute_expected)) < 1e-8)
  list(points = points, observed = observed, expected = expected,
    absolute_new_minus_old = absolute_difference, absolute_expected = absolute_expected,
    absolute_maximum_error = max(abs(absolute_difference - absolute_expected)),
    maximum_error = max(abs(observed - expected)), auxiliary_data_hash = digest::digest(data, algo = "sha256"),
    normalization_identity = stats::integrate(function(a) exp(stats::dnorm(a, log = TRUE)), -Inf, Inf),
    scope = "MC margin-only likelihood excludes a: absolute specification12-minus11 normalized target is dnorm(a,0,1,log=TRUE), including its constant. Verified old binary is a reference only; no old posterior or scientific result is reused.")
}

smi_task <- function(task, manifest, directory) {
  if (task$stage == "fit") {
    path <- smi_fit_path(manifest, task$fixture, task$backend)
    if (file.exists(path)) {
      fit <- smi_fit(manifest, task$fixture, task$backend)
      return(list(reused = TRUE, identity = nutpie_validation_identity(fit),
        summary = nutpie_validation_parameters(fit)$summary,
        diagnostics = nutpie_validation_diagnostics(fit)$summary))
    }
    arguments <- manifest$fixtures[[task$fixture]]$arguments
    arguments$backend <- task$backend; arguments$stan_options <- manifest$design$controls
    arguments$stan_options$seed <- manifest$design$fit_seed + match(task$fixture, names(manifest$fixtures))
    if (task$backend == "cmdstanr") {
      arguments$stan_options$output_dir <- file.path(directory, "cmdstan-output")
      dir.create(arguments$stan_options$output_dir)
    }
    else arguments$nutpie_options <- list(adaptation = "diag")
    fit <- do.call(antitrustBayes::bayes_est, arguments)
    stopifnot(fit@stan_data$is_single_market == as.integer(grepl("M1T", task$fixture)))
    saveRDS(fit, path)
    identity <- nutpie_validation_identity(fit)
    saveRDS(list(hash = smi_hash(path), identity = identity), file.path(dirname(path), "fit-seal.rds"))
    return(list(identity = identity, summary = nutpie_validation_parameters(fit)$summary,
      diagnostics = nutpie_validation_diagnostics(fit)$summary))
  }
  if (task$stage == "pair") {
    left <- smi_fit(manifest, task$fixture, "cmdstanr"); right <- smi_fit(manifest, task$fixture, "nutpie")
    stopifnot(identical(nutpie_validation_identity(left), nutpie_validation_identity(right)))
    paths <- stats::setNames(c(smi_fit_path(manifest, task$fixture, "cmdstanr"),
      smi_fit_path(manifest, task$fixture, "nutpie")), c("cmdstanr", "nutpie"))
    # Reuse the fixed-point mathematical checker, never its original fit paths.
    target_function <- nutpie_evidence_target
    environment(target_function) <- new.env(parent = environment(nutpie_evidence_target))
    environment(target_function)$nutpie_evidence_fit_path <- function(root, model, backend, shared = FALSE) paths[[backend]]
    target <- target_function(list(model = task$fixture, kernel = left@metadata$stan_kernel), manifest$root, directory)
    stopifnot(identical(target$status, "passed"))
    a <- nutpie_validation_parameters(left)$summary; b <- nutpie_validation_parameters(right)$summary
    comparison <- nutpie_validation_compare(a, b, identity_equal = TRUE)
    parsed <- antitrustBayes:::.bayes_backend_extract(right@stanfit, c("mu_b_market_vec", "s0_prob", "alpha_structural_market"))
    stopifnot(length(parsed$mu_b_market_vec) >= 2000L, all(is.finite(parsed$mu_b_market_vec)))
    quality <- antitrustBayes::posterior_quality(right, n_draws = 4L,
      seed = 453901L + match(task$fixture, names(manifest$fixtures)))
    stopifnot(all(is.finite(quality$draws)))
    return(list(target = target, parameters = comparison, extracted_dimensions = lapply(parsed, dim),
      posterior_quality = quality))
  }
  if (task$stage == "bridge") {
    fit <- smi_fit(manifest, task$fixture, task$backend)
    return(nutpie_closeout_public(list(component = fit), directory,
      seed = manifest$design$bridge_seed + 1000L * match(task$id, vapply(manifest$design$tasks, `[[`, character(1), "id"))))
  }
  smi_prior_check(smi_fit(manifest, task$fixture, "nutpie"), directory, manifest$legacy)
}

smi_collect <- function(manifest) {
  report <- file.path(manifest$repository, "doc/validation")
  rows <- list(); parameters <- targets <- diagnostics <- bridges <- prior <- list()
  for (task in manifest$design$tasks) {
    path <- file.path(manifest$root, task$id, "result.rds")
    result <- if (file.exists(path)) readRDS(path) else list(status = "not_run", reason = "Not executed.")
    rows[[length(rows) + 1L]] <- data.frame(task = task$id, stage = task$stage,
      status = result$status, reason = result$reason, warnings = paste(result$warnings, collapse = " | "))
    if (result$status != "completed") next
    v <- result$value
    if (task$stage == "fit" && !is.null(v$diagnostics)) { d <- v$diagnostics; d$task <- task$id; diagnostics[[task$id]] <- d }
    if (task$stage == "pair") {
      p <- v$parameters; p$fixture <- task$fixture; parameters[[task$id]] <- p
      t <- v$target$checks; t$fixture <- task$fixture; targets[[task$id]] <- t
    }
    if (task$stage == "bridge") {
      s <- v$summaries$component
      bridges[[task$id]] <- data.frame(task = task$id, fixture = task$fixture, backend = task$backend,
        logml = s$logml, uncertainty = s$uncertainty, requested = s$requested, completed = s$completed)
    }
    if (task$stage == "prior") prior[[task$id]] <- data.frame(intercept = v$points,
      observed_log_density_ratio = v$observed, expected_log_density_ratio = v$expected,
      difference = v$observed - v$expected, normalized_reference_mass = v$normalization_identity$value,
      absolute_new_minus_old = v$absolute_new_minus_old, absolute_expected = v$absolute_expected,
      absolute_error = v$absolute_new_minus_old - v$absolute_expected)
  }
  evidence_checks <- list()
  if (length(bridges)) {
    evidence <- do.call(rbind, bridges)
    for (fixture in unique(evidence$fixture)) {
      a <- evidence[evidence$fixture == fixture & evidence$backend == "cmdstanr", ]
      b <- evidence[evidence$fixture == fixture & evidence$backend == "nutpie", ]
      if (nrow(a) == 1L && nrow(b) == 1L) {
        bound <- 4 * sqrt(a$uncertainty^2 + b$uncertainty^2)
        evidence_checks[[fixture]] <- data.frame(fixture, cmdstanr = a$logml, nutpie = b$logml,
          difference = b$logml - a$logml, four_mcse = bound,
          agreement = abs(b$logml - a$logml) <= bound,
          precision = a$completed == 3L && b$completed == 3L && is.finite(bound) && bound <= .25)
      }
    }
  }
  outputs <- list(status = rows, parameters = parameters, targets = targets, diagnostics = diagnostics,
    evidence = bridges, evidence_comparisons = evidence_checks, prior = prior)
  for (name in names(outputs)) if (length(outputs[[name]])) utils::write.csv(do.call(rbind, outputs[[name]]),
    file.path(report, paste0("single_market_intercept_", name, ".csv")), row.names = FALSE)
  fixture_rows <- lapply(manifest$fixtures, function(x) { d <- x$equilibrium; d$fixture <- x$id
    d$warnings <- paste(x$warnings, collapse = " | "); d })
  utils::write.csv(do.call(rbind, fixture_rows), file.path(report, "single_market_intercept_equilibria.csv"), row.names = FALSE)
  status <- do.call(rbind, rows)
  writeLines(c("# Single-market Logit intercept validation", "",
    paste("Collected", format(Sys.time(), tz = "UTC"), "UTC."), "",
    paste(sum(status$status == "completed"), "of", nrow(status), "tasks completed;", sum(status$status == "failed"), "failed."),
    "Fresh specification12 source/data/runtime and fits are sealed under .scratch/single-market-intercept; specification11 fits are not reused.",
    "Four chains,500 warmup plus500 retained draws per chain are a bounded posterior/extraction/evidence smoke, not broad recovery certification.",
    "Posterior means, medians, SDs and central intervals are compared with four combined MCSEs. Diagnostics, unavailable MCSEs and disagreements remain visible and do not gate scientific calculations.",
    "Six fixed-point paired target checks cover both corrected kernels and all four market/time layouts. Six public bridge calls request three repetitions each, including fixed-k CF/CS.",
    "The MC auxiliary margin-only target independently isolates the single N(0,1) intercept factor. This is not a new fitted scientific model.",
    "See the single_market_intercept_status/parameters/targets/diagnostics/evidence/prior/equilibria CSVs beside this report. Raw fits, complete warnings and immutable attempts remain in the ignored output root."),
    file.path(report, "single_market_intercept_validation_results.md"))
  invisible(status)
}

smi_main <- function(action = "--prepare") {
  if (length(action) != 1L || is.na(action) || !action %in% c("--prepare", "--run", "--collect"))
    stop("Use --prepare, --run or --collect.")
  repository <- normalizePath(".", mustWork = TRUE)
  old <- readRDS(file.path(repository, ".scratch/nutpie-closeout/manifest.rds")); .libPaths(old$libraries)
  pkgload::load_all(repository, quiet = TRUE)
  cmdstanr::set_cmdstan_path(old$runtime$cmdstan$path)
  root <- file.path(repository, ".scratch/single-market-intercept"); dir.create(root, recursive = TRUE, showWarnings = FALSE)
  options(antitrustBayes.cmdstan_cache = root, antitrustBayes.nutpie_cache = root)
  for (file in c("nutpie_backend_validation.R", "nutpie_evidence_validation.R", "nutpie_closeout_validation.R"))
    sys.source(file.path(repository, "doc/validation", file), envir = environment(smi_task))
  manifest <- smi_prepare(repository, root)
  if (action == "--run") for (task in manifest$design$tasks) {
    smi_assert(manifest)
    parent <- file.path(root, task$id); dir.create(parent, recursive = TRUE, showWarnings = FALSE)
    result_path <- file.path(parent, "result.rds")
    if (file.exists(result_path)) {
      previous <- readRDS(result_path)
      if (!identical(previous$manifest_hash, smi_hash(file.path(root, "manifest.rds")))) stop("Task belongs to a different sealed run.")
      if (previous$status == "completed") next
      stop("A failed attempt is retained; inspect it before an explicitly reviewed retry.")
    }
    directory <- tempfile("attempt-", tmpdir = parent); dir.create(directory)
    warnings <- character(); cat("Starting", task$id, "\n"); flush.console()
    value <- tryCatch(withCallingHandlers(smi_task(task, manifest, directory), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }), error = identity)
    result <- list(task = task, status = if (inherits(value, "error")) "failed" else "completed",
      reason = if (inherits(value, "error")) conditionMessage(value) else "", warnings = warnings,
      value = if (inherits(value, "error")) NULL else value,
      manifest_hash = smi_hash(file.path(root, "manifest.rds")))
    saveRDS(result, file.path(directory, "result.rds")); saveRDS(result, result_path)
    smi_collect(manifest); smi_assert(manifest)
    if (result$status == "failed") stop("Task failed; retained attempt: ", directory)
  }
  smi_collect(manifest)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  smi_main(if (length(args)) args[[1L]] else "--prepare")
}
