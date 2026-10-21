#!/usr/bin/env Rscript
# Saved-fit audit only. Sourcing defines functions and never compiles or samples.
# AB_NUTPIE_SERIALIZATION_EXECUTE=true Rscript --vanilla doc/validation/nutpie_serialization_validation.R

nutpie_serialization_hash <- function(value) digest::digest(value, algo = "sha256")

nutpie_serialization_libraries <- function(paths) {
  paths <- unique(normalizePath(paths[dir.exists(paths)], mustWork = TRUE))
  paths[!vapply(paths, function(path)
    any(file.exists(file.path(path, c("nutpieR", "bridgestan"), "DESCRIPTION"))), logical(1))]
}

nutpie_serialization_no_native <- function() {
  for (package in c("nutpieR", "bridgestan")) {
    if (package %in% loadedNamespaces() || requireNamespace(package, quietly = TRUE))
      stop("Optional native package must be unavailable: ", package)
  }
  invisible(TRUE)
}

nutpie_serialization_results <- function(fit) {
  stopifnot(inherits(fit@stanfit, "NutpieFit"))
  draws <- antitrustBayes:::.bayes_backend_draws(fit@stanfit)
  parameters <- intersect(c("a_market", "s0_prob", "sigma_share", "sigma_b_strat"),
    antitrustBayes:::.bayes_backend_variable_names(fit@stanfit))
  extracted <- antitrustBayes:::.bayes_backend_extract(fit@stanfit, parameters)
  summary <- antitrustBayes:::.bayes_backend_summary(fit@stanfit, parameters)
  diagnostics <- antitrustBayes:::.bayes_backend_diagnostic_summary(fit@stanfit)
  quality <- antitrustBayes::posterior_quality(fit, n_draws = 2L, seed = 721091L)
  plan <- antitrustBayes::posterior_policy_draws(fit, n_draws = 2L, seed = 721092L)
  prepared <- suppressWarnings(antitrustBayes::posterior_policy_models(fit,
    draw_plan = plan, market = "m1", year = "t1", insideSize = 1000,
    demand_anchoring = "latent_quality"))
  if (nrow(prepared$failures)) stop("Retained baseline realization failure.")
  data <- fit@data[fit@data$market == "m1" & fit@data$year == "t1", ]
  firm_size <- aggregate(data$quantities, list(firm = data$group), sum)
  pair <- as.character(firm_size$firm[order(-firm_size$x)[1:2]])
  owner <- as.character(data$group)
  owner[owner == pair[[2L]]] <- pair[[1L]]
  tariff <- ifelse(data$group %in% pair, 0, .25 / 1.25)
  policy <- list()
  baseline <- list()
  for (position in seq_along(plan$draw_indices)) {
    source <- prepared$models[[1L]][[position]]
    if (is.null(source)) stop("Missing retained baseline model.")
    identity <- attr(source, "bayes_policy_draw", exact = TRUE)
    stopifnot(identical(identity$draw_index, plan$draw_indices[[position]]))
    promoted <- trade::as_trade_fit(source, policy = "tariff", tariffPre = 0,
      cost_basis = "effective", margin_basis = "net_revenue")
    stopifnot(identical(identity, attr(promoted, "bayes_policy_draw", exact = TRUE)))
    baseline[[position]] <- list(class = class(source), identity = identity,
      accounts = trade::tariff_accounts(promoted, preMerger = TRUE))
    for (regime in c("merger", "tariff25")) {
      remap <- function(roles) {
        if (regime == "merger") roles[roles == pair[[2L]]] <- pair[[1L]]
        unique(roles)
      }
      roles <- if (methods::is(source, "AntitrustFit")) list() else
        if (methods::is(source, "StackelbergLogit") || methods::is(source, "StackelbergCES"))
          list(leadersPost = remap(promoted@tariff_state$leadersPre)) else
          list(corePost = remap(promoted@tariff_state$corePre))
      solved <- suppressWarnings(do.call(trade::simulate, c(list(promoted,
        ownerPost = if (regime == "merger") owner else as.character(data$group),
        tariffPost = if (regime == "tariff25") tariff else 0), roles)))
      before <- trade::tariff_accounts(solved, preMerger = TRUE)
      after <- trade::tariff_accounts(solved, preMerger = FALSE)
      cost_error <- max(abs(before$physical_cost - after$physical_cost))
      stopifnot(all(is.finite(after$price)), is.finite(cost_error), cost_error < 1e-8)
      policy[[length(policy) + 1L]] <- list(regime = regime,
        draw_index = plan$draw_indices[[position]], price = after$price,
        physical_cost_error = cost_error)
    }
  }
  list(draws_hash = nutpie_serialization_hash(draws),
    metadata_hash = nutpie_serialization_hash(fit@metadata),
    native_metadata_hash = nutpie_serialization_hash(fit@stanfit$metadata),
    extraction = extracted, summary = summary, diagnostics = diagnostics,
    quality = quality, plan = plan, baseline = baseline, policy = policy)
}

nutpie_serialization_expect_error <- function(expression, pattern) {
  message <- tryCatch({ force(expression); stop("Expected artifact rejection did not occur.") },
    error = conditionMessage)
  if (!grepl(pattern, message)) stop("Unexpected artifact rejection: ", message)
  message
}

nutpie_serialization_worker <- function(request_path) {
  request <- readRDS(request_path)
  .libPaths(request$libraries)
  nutpie_serialization_no_native()
  pkgload::load_all(request$repository, quiet = TRUE)
  nutpie_serialization_no_native()
  original <- readRDS(request$fit_path)
  reloaded_path <- file.path(request$directory, "reloaded-fit.rds")
  saveRDS(original, reloaded_path)
  fit <- readRDS(reloaded_path)
  # do.call() can retain a closure in @call; its source-reference environments
  # acquire new pointer identities on reload. Compare every other slot exactly
  # and compare the serialized call bytes rather than environment pointers.
  stable_slots <- setdiff(methods::slotNames(original), "call")
  stopifnot(all(vapply(stable_slots, function(slot)
      identical(methods::slot(fit, slot), methods::slot(original, slot)), logical(1))),
    identical(nutpie_serialization_hash(fit@call), nutpie_serialization_hash(original@call)),
    identical(digest::digest(file = request$fit_path, algo = "sha256"), request$fit_hash))
  actual <- nutpie_serialization_results(fit)
  expected <- readRDS(request$reference_path)
  checks <- vapply(names(expected), function(field) identical(actual[[field]], expected[[field]]), logical(1))
  if (!all(checks)) stop("Fresh-process results changed: ", paste(names(checks)[!checks], collapse = ", "))
  # The serialized path is an evidence-restoration dependency, not a posterior dependency.
  missing <- fit
  missing@stanfit$metadata$lib_path <- file.path(request$directory, "absent-native-library.so")
  stopifnot(!file.exists(missing@stanfit$metadata$lib_path))
  missing_result <- nutpie_serialization_results(missing)
  posterior_fields <- setdiff(names(actual), "native_metadata_hash")
  stopifnot(identical(actual[posterior_fields], missing_result[posterior_fields]))
  native_error <- nutpie_serialization_expect_error(
    antitrustBayes:::.bayes_nutpie_bridge_context(fit), "requires the optional 'bridgestan'")
  nutpie_serialization_no_native()
  stopifnot(identical(digest::digest(file = request$fit_path, algo = "sha256"), request$fit_hash))
  output <- list(status = "passed", checks = checks, optional_dependencies_available = FALSE,
    missing_library_posterior_equal = TRUE, native_error = native_error,
    libraries = .libPaths(), namespaces = loadedNamespaces(), results = actual)
  saveRDS(output, request$output_path)
  invisible(output)
}

nutpie_serialization_main <- function(repository = ".") {
  if (!identical(Sys.getenv("AB_NUTPIE_SERIALIZATION_EXECUTE", "false"), "true"))
    stop("Saved-fit execution requires AB_NUTPIE_SERIALIZATION_EXECUTE=true; no sampling is performed.")
  repository <- normalizePath(repository, mustWork = TRUE)
  root <- file.path(repository, ".scratch", "nutpie-certification", "serialization")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  optional <- file.path(repository, ".scratch", "nutpie-library")
  compatible <- Sys.getenv("AB_NUTPIE_DEPENDENCY_LIBRARY",
    "/mnt/d/Projects/mergerBayes/.scratch/paper-library")
  native_libraries <- unique(c(optional, compatible, .libPaths()))
  .libPaths(native_libraries)
  pkgload::load_all(repository, quiet = TRUE)
  clean_libraries <- nutpie_serialization_libraries(c(compatible, .libPaths()))
  stopifnot(length(clean_libraries) > 0L, !any(file.exists(file.path(clean_libraries, "nutpieR", "DESCRIPTION"))),
    !any(file.exists(file.path(clean_libraries, "bridgestan", "DESCRIPTION"))))
  fixtures <- c("B", "BS", "CES_B", "CES_CS", "CF_pooled")
  script <- file.path(repository, "doc", "validation", "nutpie_serialization_validation.R")
  r_source <- sort(list.files(file.path(repository, "R"), pattern = "[.]R$", full.names = TRUE))
  seal_files <- c(script, r_source, unlist(lapply(fixtures, function(fixture)
    file.path(repository, ".scratch", "nutpie-validation", paste0(fixture, "__nutpie_diag"), "fit.rds"))))
  seals <- vapply(seal_files, digest::digest, character(1), file = TRUE, algo = "sha256")
  saveRDS(list(started_at = Sys.time(), files = seals, R = R.version.string,
    libraries = native_libraries, child_libraries = clean_libraries,
    source_revision = system2("git", c("-C", shQuote(repository), "rev-parse", "HEAD"), stdout = TRUE)),
    file.path(root, "manifest.rds"))
  rows <- list()
  for (fixture in fixtures) {
    directory <- file.path(root, fixture)
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    fit_path <- file.path(repository, ".scratch", "nutpie-validation", paste0(fixture, "__nutpie_diag"), "fit.rds")
    fit <- readRDS(fit_path)
    reference_path <- file.path(directory, "reference.rds")
    reference <- nutpie_serialization_results(fit)
    saveRDS(reference, reference_path)
    request <- list(repository = repository, libraries = clean_libraries, directory = directory,
      fit_path = fit_path, fit_hash = digest::digest(file = fit_path, algo = "sha256"),
      reference_path = reference_path, output_path = file.path(directory, "result.rds"))
    request_path <- file.path(directory, "request.rds")
    saveRDS(request, request_path)
    log_path <- file.path(directory, "worker.log")
    # R reads these before startup; .libPaths alone cannot undo already loaded namespaces.
    environment <- c(paste0("R_LIBS=", shQuote(paste(clean_libraries, collapse = .Platform$path.sep))),
      paste0("R_LIBS_USER=", shQuote(paste(clean_libraries, collapse = .Platform$path.sep))),
      paste0("R_LIBS_SITE=", shQuote(paste(clean_libraries, collapse = .Platform$path.sep))))
    status <- system2(file.path(R.home("bin"), "Rscript"),
      c("--vanilla", shQuote(script), "--worker", shQuote(request_path)),
      stdout = log_path, stderr = log_path, env = environment)
    if (status != 0L) stop("Fresh-process audit failed for ", fixture, "; see ", log_path)
    output <- readRDS(request$output_path)
    rows[[fixture]] <- data.frame(fixture, kernel = fit@metadata$stan_kernel,
      backend = "nutpie_diag", status = output$status,
      retained_draws = prod(dim(fit@stanfit$draws_array)[1:2]), checked_policy_draws = 2L,
      exact_result_fields = sum(output$checks), optional_packages_available = FALSE,
      missing_native_path_posterior_equal = output$missing_library_posterior_equal,
      merger_solves = sum(vapply(output$results$policy, function(x) x$regime == "merger", logical(1))),
      tariff_solves = sum(vapply(output$results$policy, function(x) x$regime == "tariff25", logical(1))),
      maximum_physical_cost_error = max(vapply(output$results$policy, `[[`, numeric(1), "physical_cost_error")),
      fit_sha256 = request$fit_hash, stringsAsFactors = FALSE)
    message("Fresh-process optional-free audit passed: ", fixture)
  }
  # Native-present restoration negatives mutate copies only and never load a changed binary.
  fit <- readRDS(file.path(repository, ".scratch", "nutpie-validation", "B__nutpie_diag", "fit.rds"))
  original_library <- fit@stanfit$metadata$lib_path
  library_hash <- digest::digest(file = original_library, algo = "sha256")
  missing <- fit
  missing@stanfit$metadata$lib_path <- file.path(root, "absent-native-library.so")
  absent_error <- nutpie_serialization_expect_error(
    antitrustBayes:::.bayes_nutpie_bridge_context(missing), "shared library is unavailable")
  changed <- fit
  changed_library <- file.path(root, paste0("changed-", basename(original_library)))
  stopifnot(file.copy(original_library, changed_library, overwrite = TRUE))
  connection <- file(changed_library, "ab")
  writeBin(as.raw(0), connection)
  close(connection)
  changed@stanfit$metadata$lib_path <- changed_library
  changed_error <- nutpie_serialization_expect_error(
    antitrustBayes:::.bayes_nutpie_bridge_context(changed), "shared library hash differs")
  stopifnot(identical(library_hash, digest::digest(file = original_library, algo = "sha256")),
    identical(seals, vapply(seal_files, digest::digest, character(1), file = TRUE, algo = "sha256")))
  saveRDS(list(status = "passed", absent_error = absent_error, changed_error = changed_error,
    original_library_unchanged = TRUE, source_and_fit_seals_unchanged = TRUE), file.path(root, "artifact-negatives.rds"))
  table <- do.call(rbind, rows)
  rownames(table) <- NULL
  utils::write.csv(table, file.path(repository, "doc", "validation", "nutpie_serialization_results.csv"), row.names = FALSE)
  saveRDS(table, file.path(root, "results.rds"))
  print(table)
  invisible(table)
}

if (sys.nframe() == 0L) {
  arguments <- commandArgs(trailingOnly = TRUE)
  if (length(arguments) == 2L && identical(arguments[[1L]], "--worker")) {
    nutpie_serialization_worker(arguments[[2L]])
  } else if (!length(arguments)) {
    nutpie_serialization_main()
  } else stop("Use no arguments or --worker <request.rds>.")
}
