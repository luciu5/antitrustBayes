#!/usr/bin/env Rscript

args <- commandArgs(FALSE)
script <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[[1L]])
source(file.path(dirname(normalizePath(script)), "..", "..", "inst", "benchmarks",
                 "market_intercept_enterprise", "common.R"))

run_dir <- path.expand(mic_env("MIC_ENTERPRISE_RUN_DIR",
  file.path(mic_enterprise_root(), "doc", "statistical_toolkit_runs",
            "20260922", "market-intercept-sparse-panel-pilot")))
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
manifest_path <- file.path(run_dir, "manifest.csv")
if (file.exists(manifest_path)) {
  stop("Sparse-panel pilot manifest already exists; preserve it or choose a new run directory.",
       call. = FALSE)
}

full <- mic_make_manifest()
design <- data.frame(
  demand = c("logit", "logit", "ces", "ces"),
  model = c("B", "BS", "C", "CS"),
  stringsAsFactors = FALSE
)
rows <- list()
for (i in seq_len(nrow(design))) {
  for (intercept in c("collapsed", "sampled")) {
    take <- full$demand == design$demand[[i]] &
      full$model == design$model[[i]] &
      full$intercept == intercept &
      full$source == "synthetic" & full$tier %in% c("gate", "comparison") &
      full$M == 10L & full$T == 4L
    template <- full[which(take)[[1L]], , drop = FALSE]
    template$M <- 2L
    template$T <- 2L
    template$tier <- "sparse_panel"
    template$pilot_case <- paste(design$demand[[i]], design$model[[i]],
                                 "2x2", intercept, sep = "_")
    rows[[length(rows) + 1L]] <- template
  }
}

## Ordinary Logit panel controls distinguish sparse-hierarchy behavior from
## behavior specific to the two-market/two-period information rung.
for (intercept in c("collapsed", "sampled")) {
  take <- full$demand == "logit" & full$model == "B" &
    full$intercept == intercept & full$source == "synthetic" &
    full$tier %in% c("gate", "comparison") & full$M == 10L & full$T == 4L
  template <- full[which(take)[[1L]], , drop = FALSE]
  template$tier <- "panel_control"
  template$pilot_case <- paste("logit_B_10x4", intercept, sep = "_")
  rows[[length(rows) + 1L]] <- template
}

pilot <- do.call(rbind, rows)
pilot$benchmark_kind <- "hmc_pilot"
pilot$chains <- 2L
pilot$iter <- 1200L
pilot$warmup <- 600L
pilot$cores <- 2L
pilot$task_id <- seq_len(nrow(pilot))
pilot$run_key <- vapply(seq_len(nrow(pilot)), function(i) {
  row <- pilot[i, , drop = FALSE]
  row$run_key <- NULL
  mic_hash(as.list(row))
}, character(1))
pilot <- pilot[, c("task_id", "run_key", setdiff(names(pilot),
                                                   c("task_id", "run_key")))]
mic_validate_manifest(pilot)
manifest_hash <- mic_manifest_hash(pilot)
pilot$manifest_hash <- manifest_hash

repo <- mic_enterprise_root()
prov <- mic_git_info(repo)
metadata <- data.frame(
  schema_version = 1L,
  manifest_hash = manifest_hash,
  tasks = nrow(pilot),
  generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  antitrustBayes_git = prov$sha,
  antitrustBayes_dirty = prov$dirty,
  expected_kernel_hashes = paste(mic_expected_kernel_hashes(), collapse = ";"),
  generator = normalizePath(script),
  R_version = R.version.string,
  stringsAsFactors = FALSE
)
mic_atomic_csv(pilot, manifest_path)
mic_atomic_csv(metadata, file.path(run_dir, "manifest_metadata.csv"))
mic_atomic_csv(data.frame(tier = "sparse_panel", first_task = 1L,
                          last_task = nrow(pilot), tasks = nrow(pilot)),
               file.path(run_dir, "tier_ranges.csv"))
dir.create(file.path(run_dir, "tasks"), recursive = TRUE, showWarnings = FALSE)
cat("Wrote sparse-panel pilot with", nrow(pilot), "fits to", manifest_path, "\n")
print(pilot[, c("task_id", "pilot_case", "model", "demand", "intercept",
                "M", "T", "chains", "iter", "warmup")], row.names = FALSE)
