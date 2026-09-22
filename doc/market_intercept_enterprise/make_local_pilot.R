#!/usr/bin/env Rscript

args <- commandArgs(FALSE)
script <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[[1L]])
source(file.path(dirname(normalizePath(script)), "..", "..", "inst", "benchmarks",
                 "market_intercept_enterprise", "common.R"))

run_dir <- path.expand(mic_env("MIC_ENTERPRISE_RUN_DIR",
  file.path(mic_enterprise_root(), "doc", "statistical_toolkit_runs",
            "20260922", "market-intercept-local-pilot")))
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
manifest_path <- file.path(run_dir, "manifest.csv")
if (file.exists(manifest_path)) {
  stop("Local pilot manifest already exists; preserve it or choose a new run directory.",
       call. = FALSE)
}

full <- mic_make_manifest()
selected_ids <- c(1L, 15L, 26L, 32L, 65L, 71L, 103L, 105L, 81L, 90L, 47L)
pilot <- full[match(selected_ids, full$task_id), , drop = FALSE]
pilot$pilot_case <- c(
  "logit_B_cross_collapsed_seed1",
  "logit_BS_panel_collapsed_seed1",
  "ces_C_panel_collapsed_seed1",
  "ces_CS_panel_collapsed_seed1",
  "logit_B_cross_sampled_seed1",
  "logit_BS_panel_sampled_seed1",
  "logit_BS_low_noise_collapsed",
  "logit_BF_high_rho_collapsed",
  "sc_B_cross_collapsed",
  "sc_C_panel_collapsed_seed1",
  "logit_BS_panel_collapsed_seed2"
)
sc_repeat <- pilot[pilot$pilot_case == "sc_C_panel_collapsed_seed1", , drop = FALSE]
sc_repeat$sampler_seed <- 29091L
sc_repeat$pilot_case <- "sc_C_panel_collapsed_seed2"
pilot <- rbind(pilot, sc_repeat)
pilot$tier <- "local_pilot"
pilot$benchmark_kind <- "hmc_pilot"
pilot$chains <- 2L; pilot$iter <- 1200L; pilot$warmup <- 600L; pilot$cores <- 2L
pilot$task_id <- seq_len(nrow(pilot))
pilot$run_key <- vapply(seq_len(nrow(pilot)), function(i) {
  row <- pilot[i, , drop = FALSE]
  row$run_key <- NULL
  mic_hash(as.list(row))
}, character(1))
pilot <- pilot[, c("task_id", "run_key", setdiff(names(pilot), c("task_id", "run_key")))]
mic_validate_manifest(pilot)
manifest_hash <- mic_manifest_hash(pilot)
pilot$manifest_hash <- manifest_hash

repo <- mic_enterprise_root(); prov <- mic_git_info(repo)
metadata <- data.frame(schema_version = 1L, manifest_hash = manifest_hash,
  tasks = nrow(pilot), generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  antitrustBayes_git = prov$sha, antitrustBayes_dirty = prov$dirty,
  expected_kernel_hashes = paste(mic_expected_kernel_hashes(), collapse = ";"),
  generator = normalizePath(script), R_version = R.version.string,
  stringsAsFactors = FALSE)
mic_atomic_csv(pilot, manifest_path)
mic_atomic_csv(metadata, file.path(run_dir, "manifest_metadata.csv"))
mic_atomic_csv(data.frame(tier = "local_pilot", first_task = 1L,
                          last_task = nrow(pilot), tasks = nrow(pilot)),
               file.path(run_dir, "tier_ranges.csv"))
dir.create(file.path(run_dir, "tasks"), recursive = TRUE, showWarnings = FALSE)
cat("Wrote local pilot with", nrow(pilot), "fits to", manifest_path, "\n")
print(pilot[, c("task_id", "pilot_case", "model", "demand", "intercept",
                "scenario", "chains", "iter", "warmup")], row.names = FALSE)

