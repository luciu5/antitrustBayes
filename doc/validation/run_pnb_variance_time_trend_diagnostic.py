"""Enforce the absolute diagnostic cap around preparation, compilation and fits."""
from pathlib import Path
import argparse
import atexit
import csv
import datetime as dt
import json
import os
import shutil
import subprocess
import sys

root = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser()
mode = parser.add_mutually_exclusive_group()
mode.add_argument("--plan", action="store_true")
mode.add_argument("--execute", action="store_true")
parser.add_argument("--output", default=str(root / ".scratch/variance_time_trend_diagnostic_20261109"))
parser.add_argument("--paper-root", default=str(root.parent / "mergerBayes"))
args = parser.parse_args()
command = ["Rscript", "--vanilla", str(Path(__file__).with_suffix(".R")),
           f"--output={args.output}", f"--paper-root={args.paper_root}"]
if not args.execute:
    sys.exit(subprocess.call(command + ["--plan"]))
if not shutil.which("timeout"):
    parser.error("Execution needs GNU timeout to enforce the whole-process cap.")
output = Path(args.output).resolve()
output.mkdir(parents=True, exist_ok=True)
lock = Path(str(output) + ".lock")
try:
    lock.mkdir()
except FileExistsError:
    parser.error(f"An existing controller lock is preserved: {lock}")
atexit.register(lock.rmdir)
clock_file = output / "clock.json"
now = dt.datetime.now(dt.timezone.utc)
if not clock_file.exists():
    clock_file.write_text(json.dumps({"started_utc": now.isoformat(),
        "deadline_utc": (now + dt.timedelta(seconds=5400)).isoformat(),
        "budget_seconds": 5400}, indent=2))
clock = json.loads(clock_file.read_text())
deadline = dt.datetime.fromisoformat(clock["deadline_utc"].replace("Z", "+00:00"))
remaining = int((deadline - now).total_seconds())
if remaining <= 15:
    parser.error("The original diagnostic deadline has expired; its clock is preserved.")
# Reserve cleanup/reporting inside the same absolute deadline. Inner workers
# also use that deadline; this outer timeout covers stalled parent preparation.
environment = os.environ.copy()
environment["IOVRSE_VARIANCE_TIME_PILOT_WRAPPED"] = "1"
result = subprocess.run(["timeout", "--signal=TERM", "--kill-after=5s",
    f"{remaining - 10}s"] + command + ["--execute"], cwd=root, env=environment)
if result.returncode in (124, 137, 143):
    status = output / "execution_status.csv"
    if status.exists():
        with status.open(newline="") as stream:
            reader = csv.DictReader(stream)
            fields, rows = reader.fieldnames, list(reader)
        for row in rows:
            if row["Status"] == "running":
                row["Status"], row["Exit"] = "incomplete_time_cap", str(result.returncode)
            elif row["Status"] == "not_run":
                row["Status"] = "not_run_budget"
        with status.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
    (output / "budget_stop.json").write_text(json.dumps({
        "status": "incomplete_time_cap", "deadline_utc": clock["deadline_utc"],
        "exit": result.returncode}, indent=2))
sys.exit(result.returncode)
