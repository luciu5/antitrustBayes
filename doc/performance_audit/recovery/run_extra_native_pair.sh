#!/usr/bin/env bash
set -euo pipefail

audit_repo=${PERF_RECOVERY_REPO_ROOT:-/mnt/d/projects/antitrustbayes}
audit_recovery="$audit_repo/doc/statistical_toolkit_runs/20260916/performance-audit/recovery"
audit_code="$audit_repo/doc/performance_audit/recovery/native_review"
audit_userlib=/home/ctaragin/R/x86_64-pc-linux-gnu-library/4.6
audit_artifacts="$audit_recovery/extra-native"
mkdir -p "$audit_artifacts" "$audit_recovery/logs"
export HELMERT_REVIEW_ROOT="$audit_code"
export HELMERT_ARTIFACT_DIR="$audit_artifacts"

# Separate processes prevent namespace/DSO reuse between kernel generations.
# An old/old comparison cannot qualify as an old/candidate comparison.
for audit_role in candidate baseline; do
  export R_LIBS_USER="$audit_recovery/Rlib/$audit_role:$audit_recovery/Rlib:$audit_userlib"
  Rscript --vanilla "$audit_repo/doc/performance_audit/recovery/verify_native_generation.R" \
    "$audit_role" > "$audit_artifacts/verify-$audit_role.txt" 2>&1
done

export R_LIBS_USER="$audit_recovery/Rlib/candidate:$audit_recovery/Rlib:$audit_userlib"
HELMERT_ROLE=candidate Rscript --vanilla "$audit_code/run_extra_helmert_gate.R" \
  > "$audit_recovery/logs/extra-native-candidate.log" 2>&1
export R_LIBS_USER="$audit_recovery/Rlib/baseline:$audit_recovery/Rlib:$audit_userlib"
HELMERT_ROLE=old Rscript --vanilla "$audit_code/run_extra_helmert_gate.R" \
  > "$audit_recovery/logs/extra-native-old.log" 2>&1
Rscript --vanilla "$audit_code/compare_extra_helmert_gate.R" \
  > "$audit_recovery/logs/extra-native-comparison.log" 2>&1
cat "$audit_recovery/logs/extra-native-comparison.log"
