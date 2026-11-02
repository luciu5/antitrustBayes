#!/usr/bin/env bash
#SBATCH --job-name=mic-enterprise
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=mic-enterprise-%A_%a.out
#SBATCH --error=mic-enterprise-%A_%a.err

set -euo pipefail

: "${MIC_ENTERPRISE_REPO:?Set MIC_ENTERPRISE_REPO to the antitrustBayes checkout}"
: "${MIC_ENTERPRISE_RUN_DIR:?Set MIC_ENTERPRISE_RUN_DIR to the frozen run directory}"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MIC_ENTERPRISE_MODE=run

cd "${MIC_ENTERPRISE_REPO}"
Rscript --vanilla doc/market_intercept_enterprise/run_task.R

