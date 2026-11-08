#!/usr/bin/env Rscript
# Extended, sampling-free replay of sealed specification-14 fits. The historical
# runner and manifests remain unchanged; outputs belong to a fresh directory.
saved_gate_hash <- function(path) digest::digest(file = path, algo = "sha256")

saved_gate_check_sources <- function(hashes, allow_description_change = FALSE) {
  required <- names(hashes)
  if (is.null(required) || anyDuplicated(required) || !all(file.exists(required)))
    stop("Missing or ambiguous sealed source files.")
  production <- c(sort(list.files("R", "[.]R$", full.names = TRUE)),
                  sort(list.files("inst/stan", "[.]stan$", full.names = TRUE)))
  if (!setequal(production, required[grepl("^(R/|inst/stan/)", required)]))
    stop("Production source file set differs from the sealed inference source.")
  current <- vapply(required, saved_gate_hash, character(1))
  changed <- required[current != hashes]
  if (length(setdiff(changed, "DESCRIPTION")))
    stop("Sealed inference source changed: ", paste(setdiff(changed, "DESCRIPTION"), collapse = ", "))
  if ("DESCRIPTION" %in% changed && !isTRUE(allow_description_change))
    stop("DESCRIPTION changed; review packaging/dependency changes and explicitly acknowledge them.")
  data.frame(path = required, inference_sha256 = unname(hashes),
             current_sha256 = unname(current), identical = current == hashes)
}

saved_gate_check_calls <- function(manifest) {
  if (!identical(manifest$calls_hash, digest::digest(manifest$calls, algo = "sha256")))
    stop("Saved fixture/prior/sampler call identity changed.")
  if (!identical(manifest$status, "sampling_complete") ||
      !identical(manifest$model_specification, 14L) ||
      !identical(names(manifest$calls), c("logit_C_M4T4", "ces_C_M4T4_finite_noise")))
    stop("Expected the complete two-fit specification-14 gate; incomplete fits cannot be omitted.")
  invisible(TRUE)
}

saved_gate_check_artifact <- function(path, expected) {
  if (!file.exists(path) || !identical(saved_gate_hash(path), unname(expected)))
    stop("Saved artifact seal differs: ", path)
  invisible(TRUE)
}

saved_gate_output <- function(saved_root, output) {
  saved_root <- normalizePath(saved_root, mustWork = TRUE)
  output <- normalizePath(output, mustWork = FALSE)
  if (identical(output, saved_root) || startsWith(output, paste0(saved_root, "/")) ||
      startsWith(saved_root, paste0(output, "/")))
    stop("Replay output must be separate from the read-only saved gate.")
  if (dir.exists(output) && length(list.files(output, all.files = TRUE, no.. = TRUE)))
    stop("Preserve existing replay output; select a fresh directory.")
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  normalizePath(output, mustWork = TRUE)
}

saved_gate_fit <- function(root, manifest, id) {
  path <- file.path(root, id, "fit.rds")
  saved_gate_check_artifact(path, manifest$attempts[[id]]$artifact_sha256)
  if (!identical(manifest$attempts[[id]]$status, "fit_complete"))
    stop("Saved fit is incomplete: ", id)
  fit <- readRDS(path)
  getFromNamespace("validate_fit_spec", "antitrustBayes")(fit)
  md <- fit@metadata
  args <- manifest$calls[[id]]
  kernel <- paste0(args$demand_model, "_supply")
  stopifnot(identical(as.integer(md$model_spec_version), 14L),
            identical(md$share_basis, "relative"), identical(md$backend, "cmdstanr"),
            identical(md$stan_kernel, kernel),
            identical(md$stan_source_hash, unname(manifest$source_hashes[paste0("inst/stan/", kernel, ".stan")])),
            identical(md$stan_data_hash, digest::digest(fit@stan_data[md$stan_data_fields], algo = "sha256")),
            identical(dim(fit@stanfit$draws())[1:2], c(500L, 4L)))
  controls <- md$sampler_controls
  stopifnot(controls$seed == args$stan_options$seed, controls$chains == 4L,
            controls$warmup == 500L, controls$retained_draws == 500L,
            controls$target_accept == .99, controls$max_depth == 13L,
            identical(controls$adaptation, "diag_e"))
  for (field in c("prices", "quantities", "margins", "market", "year", "group"))
    stopifnot(identical(as.character(fit@data[[field]]), as.character(args[[field]])))
  fit
}

saved_gate_policy_replay <- function(saved_root, output, allow_description_change = FALSE) {
  saved_root <- normalizePath(saved_root, mustWork = TRUE)
  seals <- read.csv("doc/validation/coordinated_release_gate/artifact_hashes.csv", stringsAsFactors = FALSE)
  before <- vapply(file.path(saved_root, seals$artifact), saved_gate_hash, character(1))
  for (i in seq_len(nrow(seals)))
    saved_gate_check_artifact(file.path(saved_root, seals$artifact[i]), seals$sha256[i])
  m <- readRDS(file.path(saved_root, "manifest.rds"))
  saved_gate_check_calls(m)
  sources <- saved_gate_check_sources(m$source_hashes, allow_description_change)
  sealed_sources <- read.csv("doc/validation/coordinated_release_gate/source_hashes.csv", stringsAsFactors = FALSE)
  stopifnot(identical(sealed_sources$path, names(m$source_hashes)),
            identical(sealed_sources$sha256, unname(m$source_hashes)))
  historical <- readRDS(file.path(saved_root, "analysis.rds"))
  stopifnot(nrow(historical$gates) == 20L, all(historical$gates$passed))
  # These helpers remain byte-identical to the reviewed sealed oracle. We do
  # not replace its strict source guard or call its output-writing entrypoints.
  oracle <- new.env(parent = environment())
  sys.source("doc/validation/coordinated_release_gate.R", envir = oracle)
  fits <- lapply(names(m$calls), function(id) saved_gate_fit(saved_root, m, id))
  names(fits) <- names(m$calls)
  output <- saved_gate_output(saved_root, output)
  rows <- list(); diagnostics <- list(); warnings <- character()
  quiet <- function(expr) withCallingHandlers(expr, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  })
  for (id in names(fits)) {
    fit <- fits[[id]]; args <- m$calls[[id]]; truth <- m$truth[[id]]
    diag <- as.data.frame(posterior::summarise_draws(fit@stanfit$draws(oracle$release_gate_active(fit))))
    diagnostics[[id]] <- data.frame(case = id, max_rhat = max(diag$rhat),
      min_bulk_ess = min(diag$ess_bulk), divergences = sum(fit@stanfit$sampler_diagnostics()[,,"divergent__"]))
    plan <- posterior_policy_draws(fit, n_draws = 1L, seed = m$thresholds$policy_plan_seed)
    stopifnot(identical(plan$draw_indices, m$thresholds$fixed_policy_draw))
    prepared <- quiet(posterior_policy_models(fit, draw_plan = plan, year = "t1",
      insideSize = 1000, demand_anchoring = "latent_quality"))
    for (market in paste0("m", 1:4)) for (state in c("truth", "posterior_fixed_draw")) {
      rows[[length(rows) + 1L]] <- tryCatch({
        if (state == "truth") {
          ix <- which(as.character(truth$market) == market & as.character(truth$year) == "t1")
          if (args$demand_model == "logit") {
            p <- truth$price[ix]; q <- truth$share[ix]; margin <- truth$margin[ix]
            owner <- as.character(truth$product[ix])
            params <- bayes_params(alpha = unique(truth$alpha[ix]), s0 = unique(truth$s0[ix]), supply_model = 3L)
          } else {
            p <- truth$prices[ix]; q <- truth$quantities[ix]; margin <- truth$margins[ix]
            owner <- as.character(truth$group[ix]); params <- bayes_params(gamma = 6, s0 = .2, supply_model = 3L)
          }
          market_data <- bayes_market(p, q, margin, group = owner, insideSize = 1000,
            outside_price = if (is.null(args$outside_price)) 1 else args$outside_price, market_side = "output")
          source_object <- quiet(posterior_policy_models(bayes_specify(params, market_data), n_draws = 1L))$models[[1L]][[1L]]
        } else {
          uid <- as.character(prepared$units$unit_id[prepared$units$market == market])
          stopifnot(length(uid) == 1L)
          source_object <- prepared$models[[uid]][[1L]]
          owner <- as.character(args$group[as.character(args$market) == market & as.character(args$year) == "t1"])
        }
        object <- if (methods::is(source_object, "AntitrustFit")) source_object@model else source_object
        stopifnot(!is.null(object))
        original <- serialize(source_object, NULL)
        baseline <- max(abs(oracle$release_gate_raw_foc(object, TRUE)))
        owner_post <- owner; firms <- unique(owner); owner_post[owner_post == firms[2L]] <- firms[1L]
        solved <- quiet(antitrust::simulate(source_object, antitrust::counterfactual(ownership = owner_post)))
        stopifnot(identical(original, serialize(source_object, NULL)))
        if (methods::is(solved, "AntitrustFit")) solved <- solved@model
        stopifnot(any(rowSums(solved@ownerPost) > 1L))
        merger <- max(abs(oracle$release_gate_raw_foc(solved, FALSE)))
        data.frame(case = id, market, state, regime = c("baseline", "merger"), residual = c(baseline, merger), reason = "")
      }, error = function(e) data.frame(case = id, market, state, regime = c("baseline", "merger"), residual = Inf, reason = conditionMessage(e)))
    }
  }
  table <- do.call(rbind, rows); table$threshold <- m$thresholds$maximum_foc_residual
  table$passed <- is.finite(table$residual) & table$residual <= table$threshold
  diagnostic <- do.call(rbind, diagnostics)
  diagnostic$passed <- with(diagnostic, is.finite(max_rhat) & max_rhat <= m$thresholds$maximum_rhat &
    min_bulk_ess >= m$thresholds$minimum_bulk_ess & divergences == m$thresholds$divergences)
  after <- vapply(file.path(saved_root, seals$artifact), saved_gate_hash, character(1))
  stopifnot(identical(before, after))
  for (id in names(m$calls)) saved_gate_check_artifact(file.path(saved_root, id, "fit.rds"), m$attempts[[id]]$artifact_sha256)
  write.csv(table, file.path(output, "policy_focs.csv"), row.names = FALSE)
  write.csv(diagnostic, file.path(output, "saved_fit_diagnostics.csv"), row.names = FALSE)
  write.csv(sources, file.path(output, "inference_source_identity.csv"), row.names = FALSE)
  writeLines(warnings, file.path(output, "warnings.txt"))
  runtime <- oracle$release_gate_packages()
  saveRDS(list(created = Sys.time(), scope = "Saved-fit identity and current-policy replay only; no new inference",
    saved_root = saved_root, original_manifest_sha256 = saved_gate_hash(file.path(saved_root, "manifest.rds")),
    historical_inference_gates = historical$gates, current_policy = table, diagnostics = diagnostic,
    sources = sources, runtime = runtime, warnings = warnings,
    saved_artifacts_unchanged = identical(before, after), sampled = FALSE), file.path(output, "replay.rds"))
  print(table, row.names = FALSE); print(diagnostic, row.names = FALSE)
  if (!all(table$passed) || !all(diagnostic$passed)) stop("Saved-fit/current-policy replay failed; inspect retained output.")
  invisible(table)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE)
  if (length(args) < 2L || length(args) > 3L ||
      (length(args) == 3L && args[3L] != "--allow-description-change"))
    stop("Usage: Rscript doc/validation/saved_gate_policy_replay.R SAVED_GATE FRESH_OUTPUT [--allow-description-change]")
  pkgload::load_all(".", quiet = TRUE)
  saved_gate_policy_replay(args[1L], args[2L], length(args) == 3L)
}
