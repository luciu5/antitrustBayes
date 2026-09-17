#!/usr/bin/env Rscript

# Exercise the complete runner's retention/bookkeeping without sampling.
args <- commandArgs(trailingOnly = TRUE)
repo <- if (length(args)) normalizePath(args[[1L]]) else normalizePath(".")
e <- new.env(parent = globalenv())
sys.source(file.path(repo, "doc/market_s0_rotation_hmc_benchmark.R"), e)
sentinel <- list(retained_fit = "mock posterior, not an HMC draw")
panel <- list(data = data.frame(group = c("a", "b")), mky_n = c(1, 1))
fixtures <- new.env(parent = emptyenv())
fixtures$cga_ladder_fixtures <- function() list(multi_market = panel)
e$rotation_source_benchmark <- function() fixtures
e$rotation_archived_fit <- function(...) NULL
e$rotation_map_primary_starts <- function(...) {
  starts <- replicate(4, list(marker = 1), simplify = FALSE)
  list(legacy = starts, rotated = starts, market_map = c(1L, 2L))
}
e$rotation_worker_args <- function(...) list()
e$getExportedValue <- function(...) function() sentinel
e$rotation_fit_diagnostics <- function(fit, ...) {
  stopifnot(identical(fit, sentinel))
  data.frame(case_id = "mock", usable = TRUE)
}
old_gate <- Sys.getenv("ANTITRUSTBAYES_ROTATION_NATIVE_GATE", unset = NA_character_)
Sys.setenv(ANTITRUSTBAYES_ROTATION_NATIVE_GATE = "approved")
out <- tempfile("rotation-retention-", tmpdir = "/tmp")
dir.create(out)
checks <- lapply(c("legacy", "rotated"), function(mode) {
  result <- e$rotation_run_worker(mode, 271828L, out, "mock archive")
  id <- paste0("primary_", mode, "_seed271828")
  path <- file.path(out, "stanfit", paste0(id, ".rds"))
  stopifnot(file.exists(path), identical(readRDS(path), sentinel),
            isTRUE(result$usable),
            file.exists(file.path(out, "results", paste0(id, ".csv"))))
  data.frame(mode = mode, retained_fit = TRUE, diagnostic_received_fit = TRUE,
             result_written = TRUE, hmc_executed = FALSE)
})
if (is.na(old_gate)) Sys.unsetenv("ANTITRUSTBAYES_ROTATION_NATIVE_GATE") else
  Sys.setenv(ANTITRUSTBAYES_ROTATION_NATIVE_GATE = old_gate)
print(do.call(rbind, checks), row.names = FALSE)
