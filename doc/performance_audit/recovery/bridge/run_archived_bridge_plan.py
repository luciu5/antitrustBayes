"""Sequential benchmark ownership; R persists each immutable attempt."""
from pathlib import Path
import os
import subprocess
import sys

compact = Path("doc/performance_audit/recovery/bridge")
recovery = Path("doc/statistical_toolkit_runs/20260916/performance-audit/recovery")
env = os.environ.copy()
env["R_LIBS_USER"] = ":".join(str(p.resolve()) for p in (
    recovery / "Rlib/candidate", recovery / "Rlib",
    Path("/home/ctaragin/R/x86_64-pc-linux-gnu-library/4.6")))
# The first seed's two defaults were already completed by the same runner.
jobs = [("candidate_seeded_serial", 707101), ("candidate_seeded_parallel", 707101),
        ("candidate_seeded_parallel", 707102), ("candidate_seeded_serial", 707102),
        ("candidate_default", 707102), ("legacy_default", 707102)]
failed = []
for mode, seed in jobs:
    print("START", mode, seed, flush=True)
    log = recovery / "logs" / f"bridge-{mode}-{seed}-v2.log"
    if log.exists():
        raise RuntimeError(f"Log already exists; retain it and inspect checkpoint before resuming: {log}")
    command = [sys.executable, "doc/performance_audit/process_memory.py",
               str(compact / f"{mode}-{seed}-memory.csv"), "Rscript", "--vanilla",
               str(compact / "run_archived_bridge_benchmark.R"), mode, str(seed)]
    with log.open("x") as out:
        result = subprocess.run(command, env=env, stdout=out, stderr=subprocess.STDOUT)
    print("FINISHED", mode, seed, "exit", result.returncode, flush=True)
    if result.returncode:
        failed.append((mode, seed, result.returncode))
print("FAILED JOBS", failed, flush=True)
sys.exit(bool(failed))
