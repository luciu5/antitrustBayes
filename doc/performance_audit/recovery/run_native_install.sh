#!/usr/bin/env bash
set -euo pipefail

repo_root=${PERF_RECOVERY_REPO_ROOT:-/mnt/d/projects/antitrustbayes}
recovery="$repo_root/doc/statistical_toolkit_runs/20260916/performance-audit/recovery"
worktrees="$recovery/worktrees"
rlib="$recovery/Rlib"
user_rlib=/home/ctaragin/R/x86_64-pc-linux-gnu-library/4.6
mkdir -p "$rlib/baseline" "$rlib/candidate"

export MAKEFLAGS=-j1
export R_LIBS_USER="$rlib:$user_rlib"

compiler=$(readlink -f "$(command -v g++)")
compiler_sha256=$(sha256sum "$compiler" | awk '{print $1}')
flags_hash=$(printf '%s\n' "$MAKEFLAGS" "$(R CMD config CXX17)" \
  "$(R CMD config CXXFLAGS)" | sha256sum | awk '{print $1}')

dependency_manifest="$recovery/dependency-manifest.txt"
if [[ ! -s "$dependency_manifest" ]]; then
  dependency_tmp="$dependency_manifest.tmp.$$"
  {
    sha256sum "$recovery/dependency-source-hashes.sha256" \
      "$recovery/antitrust-promoted-source-hashes.sha256"
    R_LIBS_USER="$rlib:$user_rlib" Rscript --vanilla -e \
      'for (p in c("antitrust", "coordination", "rstan", "StanHeaders", "Rcpp", "RcppParallel")) cat(p, "=", as.character(packageVersion(p)), "\n", sep = "")'
  } > "$dependency_tmp"
  mv -f "$dependency_tmp" "$dependency_manifest"
fi
dependency_manifest_sha256=$(sha256sum "$dependency_manifest" | awk '{print $1}')

write_record() {
  target=$1
  shift
  tmp="$target.tmp.$$"
  {
    printf '%s\n' "$@"
  } > "$tmp"
  mv -f "$tmp" "$target"
}

record_value() {
  sed -n "s/^$2=//p" "$1" | head -1
}

source_hash() {
  python3 "$repo_root/doc/performance_audit/recovery/package_source_hash.py" "$1"
}

manifest="$repo_root/doc/performance_audit/recovery/native-build-manifest.txt"
{
  printf 'started_utc=%s\n' "$(date -u +%FT%TZ)"
  printf 'repo_root=%s\n' "$repo_root"
  printf 'recovery=%s\n' "$recovery"
  printf 'R=%s\n' "$(command -v R)"
  R --version | head -2
  printf 'R_LIBS_USER=%s\n' "$R_LIBS_USER"
  printf 'MAKEFLAGS=%s\n' "$MAKEFLAGS"
  printf 'CC=%s\n' "$(R CMD config CC)"
  printf 'CXX=%s\n' "$(R CMD config CXX)"
  printf 'CXX17=%s\n' "$(R CMD config CXX17)"
  printf 'CXXFLAGS=%s\n' "$(R CMD config CXXFLAGS)"
  for tool in "$(command -v gcc)" "$(command -v g++)"; do
    printf 'compiler=%s\n' "$tool"
    readlink -f "$tool"
    sha256sum "$(readlink -f "$tool")"
    "$tool" --version | head -1
  done
  printf 'selected_cxx_path=%s\n' "$compiler"
  printf 'selected_cxx_sha256=%s\n' "$compiler_sha256"
  printf 'flags_hash=%s\n' "$flags_hash"
  printf 'dependency_manifest=%s\n' "$dependency_manifest"
  printf 'dependency_manifest_sha256=%s\n' "$dependency_manifest_sha256"
  printf '\ncommands (identical flags; role-specific source/library only):\n'
  printf 'generate headers before fingerprinting; R CMD INSTALL --no-configure --library="$ROLELIB" --no-multiarch "$SOURCE"\n'
} > "$manifest"

for role in ${PERF_RECOVERY_ROLES:-candidate}; do
  if [[ "$role" != baseline && "$role" != candidate ]]; then
    printf 'invalid recovery role: %s\n' "$role" >&2
    exit 2
  fi
  if [[ "$role" = candidate ]]; then
    source="$worktrees/bayes-helmert-candidate"
  else
    source="$worktrees/bayes-baseline"
  fi
  rolelib="$rlib/$role"
  # Freeze generated source before recording the build inputs. INSTALL must
  # not regenerate headers behind the fingerprint's back.
  mkdir -p "$recovery/logs"
  (cd "$source" && R_LIBS_USER="$rolelib:$rlib:$user_rlib" \
    Rscript --vanilla -e 'rstantools::rstan_config()') \
    > "$recovery/logs/generate-source-$role.log" 2>&1
  source_sha256=$(source_hash "$source")
  started_record="$recovery/stage-$role.started"
  complete_record="$recovery/stage-$role.complete"
  dll="$rolelib/antitrustBayes/libs/antitrustBayes.so"

  if [[ -s "$complete_record" &&
        "$(record_value "$complete_record" status)" = COMPLETE &&
        "$(record_value "$complete_record" source_manifest_version)" = 2 &&
        "$(record_value "$complete_record" source_sha256)" = "$source_sha256" &&
        "$(record_value "$complete_record" compiler_sha256)" = "$compiler_sha256" &&
        "$(record_value "$complete_record" flags_hash)" = "$flags_hash" &&
        "$(record_value "$complete_record" dependency_manifest_sha256)" = "$dependency_manifest_sha256" &&
        -f "$dll" &&
        "$(sha256sum "$dll" | awk '{print $1}')" = "$(record_value "$complete_record" dll_sha256)" ]]; then
    if R_LIBS_USER="$rolelib:$rlib:$user_rlib" Rscript --vanilla \
      "$repo_root/doc/performance_audit/recovery/verify_native_generation.R" \
      "$role" >/dev/null 2>&1; then
      printf 'role=%s action=SKIP reason=verified-complete\n' "$role"
      continue
    fi
  fi

  log="$recovery/logs/install-bayes-$role.log"
  # A known interrupted single-writer build can retain its compiled objects
  # only when all recorded build inputs are identical. Otherwise preserve the
  # objects elsewhere so make cannot silently link a stale header generation.
  reuse_objects=false
  if [[ -s "$started_record" &&
        "$(record_value "$started_record" source_manifest_version)" = 2 &&
        "$(record_value "$started_record" source_sha256)" = "$source_sha256" &&
        "$(record_value "$started_record" compiler_sha256)" = "$compiler_sha256" &&
        "$(record_value "$started_record" flags_hash)" = "$flags_hash" &&
        "$(record_value "$started_record" dependency_manifest_sha256)" = "$dependency_manifest_sha256" ]]; then
    if [[ "${PERF_RECOVERY_RESUME_INCOMPLETE:-false}" != true ]]; then
      printf 'matching unfinished build requires explicit PERF_RECOVERY_RESUME_INCOMPLETE=true after confirming no live writer\n' >&2
      exit 2
    fi
    reuse_objects=true
  fi
  retained="$recovery/retained-builds/$(date -u +%Y%m%dT%H%M%SZ)-$role"
  mkdir -p "$retained"
  if [[ "$reuse_objects" = false ]]; then
    while IFS= read -r -d '' object; do
      mv "$object" "$retained/$(basename "$object")"
    done < <(find "$source/src" -maxdepth 1 -type f \
      \( -name '*.o' -o -name '*.so' \) -print0)
  else
    # A killed compiler can leave a zero-byte object that make considers
    # current. Preserve it and rebuild it instead of linking a partial output.
    while IFS= read -r -d '' object; do
      if [[ ! -s "$object" ]] || ! readelf -h -S "$object" >/dev/null 2>&1; then
        mv "$object" "$retained/$(basename "$object")"
      fi
    done < <(find "$source/src" -maxdepth 1 -type f -name '*.o' -print0)
  fi
  for record in "$started_record" "$complete_record"; do
    if [[ -f "$record" ]]; then cp "$record" "$retained/$(basename "$record")"; fi
  done
  lock="$rolelib/00LOCK-$(basename "$source")"
  if [[ -d "$lock" ]]; then mv "$lock" "$retained/$(basename "$lock")"; fi
  printf 'role=%s reuse_matching_objects=%s retained=%s\n' "$role" "$reuse_objects" "$retained"
  write_record "$started_record" \
    "status=STARTED" "role=$role" "source_manifest_version=2" "source_sha256=$source_sha256" \
    "compiler_sha256=$compiler_sha256" "flags_hash=$flags_hash" \
    "dependency_manifest_sha256=$dependency_manifest_sha256" \
    "started_utc=$(date -u +%FT%TZ)"
  set +e
  {
    printf 'role=%s\n' "$role"
    printf 'source=%s\n' "$source"
    printf 'rolelib=%s\n' "$rolelib"
    printf 'command=R_LIBS_USER=%s MAKEFLAGS=-j1 R CMD INSTALL --no-configure --library=%s --no-multiarch %s\n' "$rolelib:$rlib:$user_rlib" "$rolelib" "$source"
    sha256sum "$source/inst/stan/"*.stan "$source/src/"*.h
    R_LIBS_USER="$rolelib:$rlib:$user_rlib" MAKEFLAGS=-j1 R CMD INSTALL \
      --no-configure --library="$rolelib" --no-multiarch "$source" || exit $?
    R_LIBS_USER="$rolelib:$rlib:$user_rlib" Rscript -e \
      'cat("antitrustBayes ", as.character(packageVersion("antitrustBayes")), "\n", sep=""); print(.libPaths())' || exit $?
    R_LIBS_USER="$rolelib:$rlib:$user_rlib" Rscript --vanilla \
      "$repo_root/doc/performance_audit/recovery/verify_native_generation.R" "$role" || exit $?
  } 2>&1 | tee -a "$log"
  install_status=${PIPESTATUS[0]}
  set -e
  if [[ "$install_status" -ne 0 ]]; then
    write_record "$recovery/stage-$role.failed" \
      "status=FAILED" "role=$role" "exit_code=$install_status" \
      "source_sha256=$source_sha256" "compiler_sha256=$compiler_sha256" \
      "flags_hash=$flags_hash" "dependency_manifest_sha256=$dependency_manifest_sha256" \
      "failed_utc=$(date -u +%FT%TZ)"
    exit "$install_status"
  fi
  dll_sha256=$(sha256sum "$dll" | awk '{print $1}')
  package_version=$(R_LIBS_USER="$rolelib:$rlib:$user_rlib" Rscript --vanilla -e \
    'cat(as.character(packageVersion("antitrustBayes")))')
  write_record "$complete_record" \
    "status=COMPLETE" "role=$role" "exit_code=0" "source_manifest_version=2" \
    "source_sha256=$source_sha256" "compiler_sha256=$compiler_sha256" \
    "flags_hash=$flags_hash" "dependency_manifest_sha256=$dependency_manifest_sha256" \
    "dll_sha256=$dll_sha256" "package_version=$package_version" \
    "completed_utc=$(date -u +%FT%TZ)"
done

printf 'completed_utc=%s\n' "$(date -u +%FT%TZ)" >> "$manifest"
