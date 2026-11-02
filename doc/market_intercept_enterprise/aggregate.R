#!/usr/bin/env Rscript

args <- commandArgs(FALSE)
script <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[[1L]])
source(file.path(dirname(normalizePath(script)), "..", "..", "inst", "benchmarks",
                 "market_intercept_enterprise", "common.R"))

run_dir <- path.expand(mic_env("MIC_ENTERPRISE_RUN_DIR",
  file.path(tempdir(), "market-intercept-enterprise")))
manifest_path <- path.expand(mic_env("MIC_ENTERPRISE_MANIFEST",
                                     file.path(run_dir, "manifest.csv")))
if (!file.exists(manifest_path)) stop("Missing manifest: ", manifest_path, call. = FALSE)
manifest_all <- utils::read.csv(manifest_path, stringsAsFactors = FALSE,
                                na.strings = c("", "NA"))
mic_validate_manifest(manifest_all)
manifest_hash <- unique(manifest_all$manifest_hash)
if (length(manifest_hash) != 1L ||
    !identical(as.character(mic_manifest_hash(manifest_all)), as.character(manifest_hash))) {
  stop("Manifest content does not match its frozen hash.", call. = FALSE)
}
tiers_text <- mic_env("MIC_ENTERPRISE_TIERS")
tiers <- if (nzchar(tiers_text)) trimws(strsplit(tiers_text, ",", fixed = TRUE)[[1L]]) else
  character()
if (length(tiers) && any(!tiers %in% unique(manifest_all$tier)))
  stop("Unknown MIC_ENTERPRISE_TIERS value.", call. = FALSE)
manifest <- if (length(tiers)) manifest_all[manifest_all$tier %in% tiers, , drop = FALSE] else
  manifest_all
suffix <- if (length(tiers)) paste0("-", paste(sort(tiers), collapse = "-")) else ""

selected <- vector("list", nrow(manifest)); audit <- vector("list", nrow(manifest))
for (i in seq_len(nrow(manifest))) {
  task <- manifest[i, , drop = FALSE]
  task_dir <- file.path(run_dir, "tasks", sprintf("task-%04d", task$task_id))
  paths <- sort(list.files(task_dir, pattern = "^attempt-[0-9]{3}\\.csv$", full.names = TRUE))
  rows <- lapply(paths, function(path) tryCatch(utils::read.csv(path,
    stringsAsFactors = FALSE), error = function(e) NULL))
  rows <- rows[!vapply(rows, is.null, logical(1))]
  valid <- rows[vapply(rows, function(x) nrow(x) == 1L &&
    identical(as.character(x$run_key), as.character(task$run_key)) &&
    identical(as.character(x$manifest_hash), as.character(manifest_hash)), logical(1))]
  chosen <- if (length(valid)) valid[[length(valid)]] else NULL
  selected[[i]] <- chosen
  audit[[i]] <- data.frame(task_id = task$task_id, run_key = task$run_key,
    tier = task$tier, attempts = length(paths), valid_attempts = length(valid),
    status = if (is.null(chosen)) "missing" else as.character(chosen$status),
    selected_attempt = if (is.null(chosen)) NA_integer_ else chosen$attempt,
    stringsAsFactors = FALSE)
}

bind_rows <- function(xs) {
  xs <- xs[!vapply(xs, is.null, logical(1))]
  if (!length(xs)) return(data.frame())
  columns <- unique(unlist(lapply(xs, names), use.names = FALSE))
  xs <- lapply(xs, function(x) {
    for (name in setdiff(columns, names(x))) x[[name]] <- NA
    x[, columns, drop = FALSE]
  })
  do.call(rbind, xs)
}
results <- bind_rows(selected)
audit <- do.call(rbind, audit)
if (nrow(results)) results <- results[order(results$task_id), , drop = FALSE]
mic_atomic_csv(results, file.path(run_dir, paste0("results", suffix, ".csv")))
mic_atomic_csv(audit, file.path(run_dir, paste0("aggregation_audit", suffix, ".csv")))
completed <- results[results$status == "completed", , drop = FALSE]
provenance <- if (!nrow(completed)) {
  data.frame(field = character(), distinct_values = integer(), consistent = logical())
} else {
  fields <- intersect(c("antitrustBayes_git", "antitrustBayes_version",
                        "model_specification_version", "installed_kernel_hashes",
                        "R_version"), names(completed))
  data.frame(field = fields,
    distinct_values = vapply(fields, function(name)
      length(unique(completed[[name]][!is.na(completed[[name]])])), integer(1)),
    consistent = vapply(fields, function(name)
      length(unique(completed[[name]][!is.na(completed[[name]])])) <= 1L, logical(1)),
    stringsAsFactors = FALSE)
}
mic_atomic_csv(provenance,
  file.path(run_dir, paste0("aggregation_provenance", suffix, ".csv")))

cat("Aggregated", nrow(results), "of", nrow(manifest), "selected tasks.\n")
print(table(audit$status, useNA = "ifany"))
strict <- mic_truthy(mic_env("MIC_ENTERPRISE_STRICT", "false"))
unusable_gate <- if (!nrow(results) || !all(c("benchmark_kind", "usable") %in%
                                            names(results))) FALSE else {
  gate <- results$status == "completed" & grepl("^hmc_", results$benchmark_kind)
  any(gate & (is.na(results$usable) | !as.logical(results$usable)))
}
if (strict && (any(audit$status != "completed") || unusable_gate ||
               (nrow(provenance) && any(!provenance$consistent)))) {
  stop("Strict aggregation failed: incomplete, unusable, or provenance-inconsistent tasks remain.",
       call. = FALSE)
}
