#!/usr/bin/env Rscript

args <- commandArgs(FALSE)
script <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[[1L]])
source(file.path(dirname(normalizePath(script)), "..", "..", "inst", "benchmarks",
                 "market_intercept_enterprise", "common.R"))

run_dir <- path.expand(mic_env("MIC_ENTERPRISE_RUN_DIR",
  file.path(tempdir(), "market-intercept-enterprise")))
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
manifest_path <- file.path(run_dir, "manifest.csv")
metadata_path <- file.path(run_dir, "manifest_metadata.csv")
tier_path <- file.path(run_dir, "tier_ranges.csv")

manifest <- mic_make_manifest()
mic_validate_manifest(manifest)
manifest_hash <- mic_manifest_hash(manifest)
manifest$manifest_hash <- manifest_hash

if (file.exists(manifest_path)) {
  old <- utils::read.csv(manifest_path, stringsAsFactors = FALSE,
                         colClasses = "character")
  new <- utils::read.csv(text = paste(capture.output(utils::write.csv(
    manifest, row.names = FALSE, na = "")), collapse = "\n"),
    stringsAsFactors = FALSE, colClasses = "character")
  if (!identical(old, new)) {
    stop("Existing manifest differs. Preserve it and choose a new run directory.",
         call. = FALSE)
  }
  cat("Verified existing immutable manifest:", manifest_path, "\n")
  quit(save = "no", status = 0L)
}

repo <- mic_enterprise_root()
prov <- mic_git_info(repo)
metadata <- data.frame(
  schema_version = 1L,
  manifest_hash = manifest_hash,
  tasks = nrow(manifest),
  generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  antitrustBayes_git = prov$sha,
  antitrustBayes_dirty = prov$dirty,
  expected_kernel_hashes = paste(mic_expected_kernel_hashes(), collapse = ";"),
  generator = normalizePath(script),
  R_version = R.version.string,
  stringsAsFactors = FALSE
)
mic_atomic_csv(manifest, manifest_path)
mic_atomic_csv(metadata, metadata_path)
tier_ranges <- do.call(rbind, lapply(split(manifest$task_id, manifest$tier), function(x)
  data.frame(first_task = min(x), last_task = max(x), tasks = length(x))))
tier_ranges$tier <- rownames(tier_ranges); rownames(tier_ranges) <- NULL
tier_ranges <- tier_ranges[, c("tier", "first_task", "last_task", "tasks")]
mic_atomic_csv(tier_ranges, tier_path)
dir.create(file.path(run_dir, "tasks"), recursive = TRUE, showWarnings = FALSE)
cat("Wrote", nrow(manifest), "tasks to", manifest_path, "\n")
print(table(manifest$tier))
print(tier_ranges, row.names = FALSE)
