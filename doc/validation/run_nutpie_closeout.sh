#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
case "${1:---prepare}" in
  --prepare|--collect) export AB_NUTPIE_CLOSEOUT_EXECUTE=false ;;
  --run) export AB_NUTPIE_CLOSEOUT_EXECUTE=true AB_NUTPIE_CLOSEOUT_SAMPLE=true ;;
  *) printf '%s\n' 'Use --prepare, --run, or --collect.' >&2; exit 2 ;;
esac
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export CARGO_BUILD_JOBS=2 MAKEFLAGS=-j2 CMAKE_BUILD_PARALLEL_LEVEL=2
export AB_NUTPIE_CLOSEOUT_OUTPUT="${AB_NUTPIE_CLOSEOUT_OUTPUT:-$PWD/.scratch/nutpie-closeout}"
export AB_NUTPIE_CLOSEOUT_ACTION="${1:---prepare}"
mkdir -p "$AB_NUTPIE_CLOSEOUT_OUTPUT"
exec 9>"$AB_NUTPIE_CLOSEOUT_OUTPUT/run.lock"
flock -n 9 || { printf '%s\n' 'Another closeout invocation holds the output lock.' >&2; exit 1; }
Rscript --vanilla -e '
  optional <- Sys.getenv("AB_NUTPIE_LIBRARY", ".scratch/nutpie-library")
  structural <- Sys.getenv("AB_NUTPIE_STRUCTURAL_LIBRARY", "/mnt/d/Projects/mergerBayes/.scratch/paper-library")
  .libPaths(c(optional, structural, .libPaths()))
  source("doc/validation/nutpie_closeout_validation.R")
  nutpie_closeout_main(Sys.getenv("AB_NUTPIE_CLOSEOUT_ACTION"))
'
