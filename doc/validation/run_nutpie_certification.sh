#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
case "${1:---prepare}" in
  --prepare|--collect) export AB_NUTPIE_PILOT_EXECUTE=false ;;
  --pilot) export AB_NUTPIE_PILOT_EXECUTE=true ;;
  *) printf '%s\n' 'Use --prepare, --pilot, or --collect.' >&2; exit 2 ;;
esac
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export CARGO_BUILD_JOBS=2 MAKEFLAGS=-j2 CMAKE_BUILD_PARALLEL_LEVEL=2
export AB_NUTPIE_PILOT_OUTPUT="${AB_NUTPIE_PILOT_OUTPUT:-$PWD/.scratch/nutpie-certification/pilot-accept99}"
mkdir -p "$AB_NUTPIE_PILOT_OUTPUT"
exec 9>"$AB_NUTPIE_PILOT_OUTPUT/run.lock"
flock -n 9 || { printf '%s\n' 'Another pilot invocation holds this output lock.' >&2; exit 1; }
export AB_NUTPIE_PILOT_ACTION="${1:---prepare}"
Rscript --vanilla -e '
  optional <- Sys.getenv("AB_NUTPIE_LIBRARY", ".scratch/nutpie-library")
  structural <- Sys.getenv("AB_NUTPIE_STRUCTURAL_LIBRARY", "/mnt/d/Projects/mergerBayes/.scratch/paper-library")
  .libPaths(c(optional, structural, .libPaths()))
  source("doc/validation/nutpie_certification_pilot.R")
  nutpie_pilot_main(Sys.getenv("AB_NUTPIE_PILOT_ACTION"))
'
