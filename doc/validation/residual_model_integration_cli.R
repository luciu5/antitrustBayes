args <- commandArgs(TRUE)
stopifnot(length(args) >= 2L)
.libPaths(c("/mnt/d/Projects/mergerBayes/.scratch/paper-library", .libPaths()))
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1", STAN_NUM_THREADS = "1", MAKEFLAGS = "-j1")
options(antitrustBayes.cmdstan_cache = normalizePath(".scratch", mustWork = TRUE),
  antitrustBayes.nutpie_cache = normalizePath(".scratch", mustWork = TRUE))
pkgload::load_all(".", quiet = TRUE)
source("doc/validation/residual_model_integration.R")
output <- normalizePath(args[[2L]], mustWork = TRUE)
switch(args[[1L]], prepare = rmi_prepare(output), seal = rmi_seal(output),
  worker = rmi_worker(output, args[[3L]]), sample = rmi_run(output),
  stop("Unknown mode"))
