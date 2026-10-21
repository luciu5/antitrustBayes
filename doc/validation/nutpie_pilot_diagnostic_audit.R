#!/usr/bin/env Rscript
# No-HMC audit of constant-field availability after the sealed nine-fit pilot.
# AB_NUTPIE_DIAGNOSTIC_AUDIT_EXECUTE=true Rscript --vanilla doc/validation/nutpie_pilot_diagnostic_audit.R

nutpie_pilot_audit_reference <- function(fixture, variant, manifest) {
  if (identical(variant, "cmdstanr")) manifest$reference_fits[[fixture]]$path else
    file.path(manifest$original_root, paste0(fixture, "__nutpie_diag"), "fit.rds")
}

nutpie_pilot_audit_reason <- function(draws, availability) {
  vapply(seq_len(nrow(availability)), function(i) {
    if (availability$diagnostics_available[i]) return("")
    x <- as.numeric(draws[, , availability$variable[i]])
    if (any(!is.finite(x))) "nonfinite retained values" else
      if (availability$constant[i]) "constant absent from reference or not exactly equal" else
        "nonconstant field has unavailable R-hat/bulk ESS/tail ESS"
  }, character(1))
}

nutpie_pilot_diagnostic_audit <- function(repository = ".") {
  if (!identical(Sys.getenv("AB_NUTPIE_DIAGNOSTIC_AUDIT_EXECUTE", "false"), "true"))
    stop("No-HMC audit requires AB_NUTPIE_DIAGNOSTIC_AUDIT_EXECUTE=true.")
  repository <- normalizePath(repository, mustWork = TRUE)
  .libPaths(c(file.path(repository, ".scratch/nutpie-library"),
    Sys.getenv("AB_NUTPIE_DEPENDENCY_LIBRARY", "/mnt/d/Projects/mergerBayes/.scratch/paper-library"), .libPaths()))
  pkgload::load_all(repository, quiet = TRUE)
  root <- normalizePath(Sys.getenv("AB_NUTPIE_PILOT_OUTPUT",
    file.path(repository, ".scratch/nutpie-certification/pilot-accept99")), mustWork = TRUE)
  manifest_path <- file.path(root, "manifest.rds")
  manifest <- readRDS(manifest_path)
  frozen <- file.path(root, "frozen-runner.R")
  helper <- new.env(parent = globalenv())
  sys.source(file.path(repository, "doc/validation/nutpie_backend_validation.R"), envir = helper)
  if (!identical(helper$nutpie_validation_runtime(repository), manifest$runtime$runtime))
    stop("Production source/runtime differs from the original pilot manifest.")
  if (!identical(digest::digest(file = frozen, algo = "sha256"), manifest$runtime$runner))
    stop("Frozen pilot helper differs from the pilot manifest.")
  sys.source(frozen, envir = helper)
  if (!identical(manifest$design, helper$nutpie_pilot_design())) stop("Frozen pilot design differs.")
  tasks <- manifest$design$tasks
  if (nrow(tasks) != 9L || anyDuplicated(tasks$id)) stop("Audit requires all nine planned fits.")
  result_paths <- file.path(root, tasks$id, "result.rds")
  if (!all(file.exists(result_paths))) stop("All nine pilot results must exist before auditing.")
  results <- lapply(result_paths, readRDS)
  if (!all(vapply(results, function(x) identical(x$status, "completed"), logical(1))))
    stop("All nine pilot tasks must be completed before auditing.")
  fit_paths <- file.path(root, tasks$id, "fit.rds")
  native_references <- vapply(unique(tasks$fixture), nutpie_pilot_audit_reference, character(1),
    variant = "nutpie_diag", manifest = manifest)
  original_references <- vapply(manifest$reference_fits, `[[`, character(1), "path")
  production <- file.path(repository, names(manifest$runtime$runtime$source_hashes))
  stan_files <- vapply(names(manifest$stan_source_hash), antitrustBayes:::.cmdstan_kernel_file, character(1))
  if (!identical(vapply(stan_files, digest::digest, character(1), file = TRUE, algo = "sha256"),
      manifest$stan_source_hash)) stop("Installed Stan sources differ from the pilot manifest.")
  raw_reports <- file.path(repository, "doc/validation",
    c("nutpie_certification_pilot_status.csv", "nutpie_certification_pilot_diagnostics.csv"))
  files <- unique(c(manifest_path, frozen, result_paths, fit_paths, native_references,
    original_references, manifest$original_manifest, production, stan_files, raw_reports,
    file.path(repository, "doc/validation/nutpie_pilot_diagnostic_audit.R")))
  seals <- vapply(files, digest::digest, character(1), file = TRUE, algo = "sha256")
  if (!identical(unname(seals[manifest$original_manifest]), manifest$original_manifest_hash) ||
      any(vapply(names(manifest$reference_fits), function(name)
        !identical(unname(seals[manifest$reference_fits[[name]]$path]),
          manifest$reference_fits[[name]]$sha256), logical(1)))) stop("Original manifest/reference seal changed.")
  output <- file.path(root, paste0("diagnostic-audit-", format(Sys.time(), "%Y%m%dT%H%M%S")))
  if (!dir.create(output, recursive = FALSE)) stop("Choose a new audit output directory.")
  saveRDS(list(created_at = Sys.time(), files = seals,
    rule = "Same frozen helper and exact equality; same-backend constants only.",
    public_gate_changed = FALSE), file.path(output, "manifest.rds"))
  chains <- variables <- list()
  for (i in seq_len(nrow(tasks))) {
    task <- tasks[i, ]; raw <- results[[i]]; fit <- readRDS(fit_paths[i])
    rownames(task) <- NULL
    old_reference <- readRDS(manifest$reference_fits[[task$fixture]]$path)
    reference_path <- nutpie_pilot_audit_reference(task$fixture, task$variant, manifest)
    reference <- readRDS(reference_path)
    identity <- helper$nutpie_validation_identity(fit)
    stopifnot(identical(identity, helper$nutpie_validation_identity(old_reference)),
      identical(identity, helper$nutpie_validation_identity(reference)), isTRUE(raw$identity_equal),
      identical(raw$fit_sha256, unname(seals[fit_paths[i]])))
    draws <- antitrustBayes:::.bayes_backend_draws(fit@stanfit)
    stopifnot(identical(raw$posterior_sha256, digest::digest(draws, algo = "sha256")))
    summary <- raw$diagnostics$common$posterior
    constants <- function(x) helper$nutpie_pilot_constant_reference(
      antitrustBayes:::.bayes_backend_draws(x@stanfit))
    original <- helper$nutpie_pilot_posterior_availability(draws, summary, constants(old_reference))
    audited <- helper$nutpie_pilot_posterior_availability(draws, summary, constants(reference))
    stopifnot(identical(original, raw$posterior_diagnostic_availability))
    diagnostics <- raw$diagnostics$summary
    stopifnot(all(diagnostics$unavailable_posterior_diagnostics == sum(!original$diagnostics_available)))
    diagnostics$raw_unavailable_posterior_diagnostics <- diagnostics$unavailable_posterior_diagnostics
    diagnostics$unavailable_posterior_diagnostics <- sum(!audited$diagnostics_available)
    raw_pass <- helper$nutpie_pilot_diagnostic_pass(raw$diagnostics$summary, manifest$design$diagnostic_thresholds)
    audited_pass <- helper$nutpie_pilot_diagnostic_pass(diagnostics, manifest$design$diagnostic_thresholds)
    stopifnot(identical(raw_pass, raw$diagnostics_passed))
    chains[[i]] <- cbind(task, diagnostics, raw_pass, audited_pass,
      reference_variant = if (task$variant == "cmdstanr") "cmdstanr" else "nutpie_diag")
    variables[[i]] <- data.frame(task, variable = audited$variable,
      raw_available = original$diagnostics_available, audited_available = audited$diagnostics_available,
      constant = audited$constant, raw_constant_match = original$constant_matches_reference,
      audited_constant_match = audited$constant_matches_reference,
      raw_reason = nutpie_pilot_audit_reason(draws, original),
      audited_reason = nutpie_pilot_audit_reason(draws, audited), stringsAsFactors = FALSE)
  }
  stopifnot(identical(seals, vapply(files, digest::digest, character(1), file = TRUE, algo = "sha256")))
  chains <- do.call(rbind, chains); variables <- do.call(rbind, variables)
  result <- list(chains = chains, variables = variables, original_seals_unchanged = TRUE,
    output_directory = output, public_gate_changed = FALSE)
  saveRDS(result, file.path(output, "result.rds"))
  for (name in c("chains", "variables")) utils::write.csv(result[[name]],
    file.path(repository, "doc/validation", paste0("nutpie_pilot_diagnostic_audit_", name, ".csv")), row.names = FALSE)
  print(unique(chains[c("fixture", "variant", "raw_pass", "audited_pass")]))
  message("No-HMC audit preserved all ", length(seals), " original seals: ", output)
  invisible(result)
}

if (sys.nframe() == 0L) nutpie_pilot_diagnostic_audit()
