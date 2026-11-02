"""Require all prespecified calls and comparisons; never drop failed jobs."""
from pathlib import Path
import csv
import math

root = Path("doc/performance_audit/recovery/bridge")

def records(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))

rows = []
for seed in [707101, 707102]:
    for title, old, new in [
        ("default serial regression", "legacy_default", "candidate_default"),
        ("two-worker bridge candidate", "candidate_seeded_serial", "candidate_seeded_parallel"),
    ]:
        timing = [records(root / f"{mode}-{seed}-timing.csv")[0] for mode in [old, new]]
        memory = [records(root / f"{mode}-{seed}-memory.csv")[0] for mode in [old, new]]
        evidence = [records(root / f"{mode}-{seed}-evidence.csv") for mode in [old, new]]
        expected = {(model, str(rep)) for model in ["BF", "BS"] for rep in [1, 2, 3]}
        for result in evidence:
            assert len(result) == len(expected)
            assert {(row["model"], row["repetition"]) for row in result} == expected
        mapping = [{(r["model"], r["repetition"]): r for r in result} for result in evidence]
        differences = [abs(float(mapping[0][key]["logml"]) - float(mapping[1][key]["logml"]))
                       for key in sorted(expected)]
        assert all(math.isfinite(d) for d in differences)
        iterations_same = all(mapping[0][key]["niter"] == mapping[1][key]["niter"] for key in expected)
        immutable = all(row["fit_immutable"] == "TRUE" for result in evidence for row in result)
        for t, m in zip(timing, memory):
            assert t["status"] == "COMPLETE" and m["exit_status"] == "0"
            assert t["cached"] == "FALSE", "Cached calls cannot supply fresh memory timing"
        old_seconds, new_seconds = [float(t["elapsed_seconds"]) for t in timing]
        rows.append(dict(change=title, benchmark="archived N24 exact BF and BS, separate calls",
            seed=seed, old_wall_seconds=old_seconds, new_wall_seconds=new_seconds,
            speed_ratio=old_seconds/new_seconds,
            old_peak_pss_MiB=float(memory[0]["process_tree_peak_pss_bytes"])/1024**2,
            new_peak_pss_MiB=float(memory[1]["process_tree_peak_pss_bytes"])/1024**2,
            memory_ratio=float(memory[1]["process_tree_peak_pss_bytes"])/float(memory[0]["process_tree_peak_pss_bytes"]),
            max_logml_difference=max(differences), iterations_identical=iterations_same,
            input_fits_immutable=immutable, old_warnings=int(timing[0]["warning_count"]),
            new_warnings=int(timing[1]["warning_count"]),
            old_failures=int(timing[0]["failure_count"]), new_failures=int(timing[1]["failure_count"]),
            recommendation="hold candidate; small speed gain and substantial memory increase" if "worker" in title else "default numerical regression passes",
            scope="no new HMC or cross-DGP Bayes factor"))
with (root / "archived_benchmark_summary.csv").open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
for row in rows:
    print(row["change"], row["seed"], "speed", round(row["speed_ratio"], 3),
          "PSS ratio", round(row["memory_ratio"], 3), "difference", row["max_logml_difference"])
