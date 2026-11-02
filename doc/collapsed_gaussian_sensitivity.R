#!/usr/bin/env Rscript
## Only the three sensitivities frozen in benchmark_protocol.md. No rescue.
cga_sensitivity <- function(out_dir) {
  e <- new.env()
  sys.source("doc/collapsed_gaussian_benchmark.R", e)
  out_dir <- e$cga_require_tmp(out_dir)
  if (file.exists(file.path(out_dir, "sensitivity_manifest.csv"))) stop("Preserve the existing sensitivity run.")
  base <- subset(e$cga_make_plan(), design == "logit_bertrand_noisy_271828" & mode == "B")
  panel <- e$cga_make_dgp(seed = base$seed, demand_model = "logit", supply_model = "bertrand",
    structure = "simultaneous", share_scenario = "noisy")
  sampler_seed <- e$cga_fit_args(panel, base)$stan_options$seed
  settings <- data.frame(experiment = c("quality_prior_width_half", "quality_prior_width_double", "rho_enabled"),
    width = c(.5, 2, 1), use_rho = c(FALSE, FALSE, TRUE), stringsAsFactors = FALSE)
  write.csv(settings, file.path(out_dir, "sensitivity_manifest.csv"), row.names = FALSE)
  rows <- list()
  for (i in seq_len(nrow(settings))) {
    setting <- settings[i, ]
    case <- base
    case$case_id <- paste0(base$case_id, "__", setting$experiment)
    case$bridge_requested <- FALSE
    case$use_rho <- setting$use_rho
    estimator <- function(...) {
      args <- list(...)
      args$priors <- list(prior_sigma_meanval_strat = setting$width)
      args$use_rho <- setting$use_rho
      args$stan_options$seed <- sampler_seed
      do.call(antitrustBayes::bayes_est, args)
    }
    cat("Starting prespecified sensitivity: ", setting$experiment, "\n", sep = "")
    value <- e$cga_run_one(panel, case, out_dir, estimator)
    row <- value$diagnostics
    row$experiment <- setting$experiment
    row$quality_prior_width <- setting$width
    row$use_rho <- setting$use_rho
    row$predictive_quality_conditioning <- if (setting$use_rho)
      "conditional Gaussian quality given shares and matched margins" else
      "conditional Gaussian quality given shares"
    rows[[i]] <- row
    write.csv(do.call(rbind, rows), file.path(out_dir, "sensitivity.csv"), row.names = FALSE, na = "")
    cat("Finished ", setting$experiment, ": ", row$status, "\n", sep = "")
  }
  invisible(do.call(rbind, rows))
}
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2L || args[1] != "--run")
    stop("Usage: Rscript doc/collapsed_gaussian_sensitivity.R --run /tmp/new-cga-sensitivity")
  cga_sensitivity(args[2])
}
