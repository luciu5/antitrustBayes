#!/usr/bin/env Rscript
## Two mode-B information-ladder checks prespecified in the primary protocol.
## These use its explicit unpooled flags, not the separate public auto preset.
cga_ladder_benchmark <- function(out_dir) {
  e <- new.env(); sys.source("doc/collapsed_gaussian_benchmark.R", e)
  out_dir <- e$cga_require_tmp(out_dir)
  if (file.exists(file.path(out_dir, "ladder_manifest.csv"))) stop("Preserve the existing ladder run.")
  base <- subset(e$cga_make_plan(), design == "logit_bertrand_noisy_271828" & mode == "B")
  fixtures <- e$cga_ladder_fixtures()
  write.csv(data.frame(rung = names(fixtures), profile = "prespecified_B_only_unpooled"),
    file.path(out_dir, "ladder_manifest.csv"), row.names = FALSE)
  e$cga_api_ready()
  rows <- list()
  for (rung in names(fixtures)) {
    case <- base; case$case_id <- paste0(base$case_id, "__ladder_", rung)
    case$bridge_requested <- FALSE
    cat("Starting prespecified ladder check: ", rung, "\n", sep = "")
    value <- e$cga_run_one(fixtures[[rung]], case, out_dir)
    row <- value$diagnostics
    row$rung <- rung
    row$N <- nrow(fixtures[[rung]]$data)
    row$G <- length(unique(fixtures[[rung]]$data$group))
    rows[[rung]] <- row
    write.csv(do.call(rbind, rows), file.path(out_dir, "ladder.csv"), row.names = FALSE, na = "")
    cat("Finished ", rung, ": ", row$status, "\n", sep = "")
  }
  invisible(do.call(rbind, rows))
}
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2L || args[1] != "--run")
    stop("Usage: Rscript doc/collapsed_gaussian_ladder_benchmark.R --run /tmp/new-cga-ladder")
  cga_ladder_benchmark(args[2])
}
