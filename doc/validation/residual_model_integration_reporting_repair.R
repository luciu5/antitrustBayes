# Reporting-only completion after the six original workers finished. It never
# runs HMC, changes a fit, or rewrites the original execution manifest/statuses.
rmi_reporting_repair <- function(output) {
  manifest_path <- file.path(output, "manifest.rds")
  manifest <- readRDS(manifest_path); rmi_assert(manifest, output)
  stopifnot(length(manifest$attempts) == 6L,
    all(vapply(manifest$attempts, function(x) x$status != "started", logical(1))))
  receipt_path <- file.path(output, "reporting_repair_receipt.rds")
  if (file.exists(receipt_path)) stop("Completed reporting repair is preserved.")
  original_manifest_hash <- rmi_hash(manifest_path, TRUE)
  fixtures <- readRDS(file.path(output, "fixtures.rds")); receipts <- list()
  for (id in names(manifest$calls)) {
    directory <- file.path(output, id); fit_path <- file.path(directory, "fit.rds")
    fit_seal <- readRDS(file.path(directory, "seal.rds"))
    stopifnot(identical(rmi_hash(fit_path, TRUE), fit_seal$fit_sha256))
    fit <- readRDS(fit_path); fixture <- fixtures[[sub("__.*", "", id)]]
    stopifnot(identical(as.integer(fit@metadata$model_spec_version), 15L),
      identical(dim(fit@stanfit$draws())[1:2], c(500L, 4L)),
      identical(fit@metadata$stan_source_hash,
        unname(manifest$source_hashes[paste0("inst/stan/", fit@metadata$stan_kernel, ".stan")])))
    frozen <- manifest$prepared_stan_data[[id]]; actual <- fit@stan_data
    differences <- names(frozen)[vapply(names(frozen), function(field)
      !isTRUE(all.equal(frozen[[field]], actual[[field]])), logical(1))]
    # bayes_est deterministically resolves reduce_sum grain size from N/cores.
    # It is computational scheduling; all scientific data must match exactly.
    stopifnot(all(differences %in% "grainsize"))
    inverse <- as.matrix(fit@stanfit$draws("pred_marginNormInv", format = "draws_matrix"))
    stopifnot(all(is.finite(inverse)), all(inverse > 0), ncol(inverse) == nrow(fixture$data))
    margin <- 1 / inverse
    margin_path <- file.path(directory, "margin_fit.csv")
    stopifnot(!file.exists(margin_path))
    write.csv(data.frame(row = seq_len(nrow(fixture$data)),
      observed_margin = fixture$data$margins, posterior_mean_margin = colMeans(margin),
      error = colMeans(margin) - fixture$data$margins,
      interpretation = "Posterior mean structural FOC margin, inverse of pred_marginNormInv; excludes replicated observation error."),
      margin_path, row.names = FALSE)
    policy <- read.csv(file.path(directory, "policy_all_rows.csv"))
    prepared <- readRDS(file.path(directory, "baseline_models.rds"))
    stopifnot(nrow(policy) == 100L * nrow(prepared$units))
    counts <- data.frame(unit = prepared$units$unit_id,
      attempted = as.integer(table(factor(policy$unit, levels = prepared$units$unit_id))),
      invalid = as.integer(table(factor(policy$unit[policy$status != "valid"], levels = prepared$units$unit_id))))
    invalid_path <- file.path(directory, "invalid_counts.csv")
    stopifnot(!file.exists(invalid_path)); write.csv(counts, invalid_path, row.names = FALSE)
    receipts[[id]] <- list(original_worker_status = manifest$attempts[[id]]$status,
      inference_completed = TRUE, original_error_stage = "final margin report extraction after policy checkpoints",
      original_error = "Can't find the following variable(s) in the output: pred_margin",
      fixed_report_field = "1 / pred_marginNormInv", fit_sha256 = rmi_hash(fit_path, TRUE),
      fit_seal_sha256 = rmi_hash(file.path(directory, "seal.rds"), TRUE),
      actual_stan_data_sha256 = rmi_hash(actual), prepared_stan_data_sha256 = rmi_hash(frozen),
      computational_difference = differences, prepared_grainsize = frozen$grainsize,
      actual_grainsize = actual$grainsize, scientific_data_identical = TRUE,
      margin_csv_sha256 = rmi_hash(margin_path, TRUE),
      policy_rows_sha256 = rmi_hash(file.path(directory, "policy_all_rows.csv"), TRUE),
      invalid_counts_sha256 = rmi_hash(invalid_path, TRUE))
    stopifnot(identical(receipts[[id]]$fit_sha256, fit_seal$fit_sha256))
  }
  stopifnot(identical(original_manifest_hash, rmi_hash(manifest_path, TRUE)))
  receipt <- list(at = Sys.time(), original_manifest_sha256 = original_manifest_hash,
    script_sha256 = rmi_hash("doc/validation/residual_model_integration_reporting_repair.R", TRUE),
    fit_receipts = receipts, new_fits = 0L, new_policy_draws = 0L,
    explanation = "All six HMC fits and all 2200 policy unit/draw states had completed before a nonexistent reporting field caused worker exit 1. This repair only writes missing margin and invalid-count tables. Original execution statuses, logs, fits, plans and seals remain unchanged.")
  saveRDS(receipt, receipt_path); invisible(receipt)
}
