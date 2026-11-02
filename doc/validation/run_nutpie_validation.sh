#!/usr/bin/env bash
# Run from the repository root. No arguments prepares/collects without sampling.
set -euo pipefail

mode="${1:---prepare}"
case "$mode" in
  --prepare|--compile|--run|--evidence) ;;
  *) printf 'Usage: bash doc/validation/run_nutpie_validation.sh [--prepare|--compile|--run|--evidence]\n' >&2; exit 2 ;;
esac
test -f R/estimate.R && test -f doc/validation/nutpie_backend_validation.R || {
  printf 'Run this command from the antitrustBayes repository root.\n' >&2
  exit 2
}

repository="$(pwd -P)"
export R_LIBS="${AB_NUTPIE_R_LIBS:-$repository/.scratch/nutpie-library:$repository/../mergerBayes/.scratch/paper-library}${R_LIBS:+:$R_LIBS}"
export AB_NUTPIE_OUTPUT="${AB_NUTPIE_OUTPUT:-$repository/.scratch/nutpie-validation}"
export AB_NUTPIE_MAX_SECONDS="${AB_NUTPIE_MAX_SECONDS:-14400}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export CARGO_BUILD_JOBS=2 MAKEFLAGS=-j2 CMAKE_BUILD_PARALLEL_LEVEL=2
mkdir -p "$AB_NUTPIE_OUTPUT"
command -v flock >/dev/null || { printf 'GNU flock is required to prevent duplicate validation jobs.\n' >&2; exit 2; }
exec 9>"$AB_NUTPIE_OUTPUT/run.lock"
flock -n 9 || { printf 'Another validation invocation holds %s/run.lock\n' "$AB_NUTPIE_OUTPUT" >&2; exit 2; }

run_validation() {
case "$mode" in
  --prepare)
    export AB_NUTPIE_EXECUTE=false
    Rscript --vanilla doc/validation/nutpie_backend_validation.R
    ;;
  --compile)
    Rscript --vanilla -e '
      pkgload::load_all(".", quiet = TRUE)
      root <- normalizePath(Sys.getenv("AB_NUTPIE_OUTPUT"), mustWork = TRUE)
      options(antitrustBayes.cmdstan_cache = file.path(root, "cmdstan-cache"),
              antitrustBayes.nutpie_cache = file.path(root, "nutpie-cache"))
      values <- lapply(c("logit_supply", "logit_stackelberg", "ces_supply", "ces_stackelberg"), function(kernel) {
        model <- antitrustBayes:::.nutpie_model_for_kernel(kernel)
        cat(kernel, model$metadata$source_hash, model$metadata$compile_time, "seconds\n")
        model$metadata
      })
      saveRDS(values, file.path(root, "precompiled.rds"))
      cat("Four unchanged kernels compiled; no inference ran.\n")'
    ;;
  --run)
    export AB_NUTPIE_EXECUTE=true
    Rscript --vanilla doc/validation/nutpie_backend_validation.R
    ;;
  --evidence)
    export AB_NUTPIE_EVIDENCE_EXECUTE=true AB_NUTPIE_EVIDENCE_SAMPLE=true
    Rscript --vanilla doc/validation/nutpie_evidence_validation.R
    ;;
esac
}

log_path="$AB_NUTPIE_OUTPUT/invocation-$(date +%Y%m%dT%H%M%S)-${mode#--}.log"
printf 'Validation log: %s\n' "$log_path"
run_validation 2>&1 | tee "$log_path"
