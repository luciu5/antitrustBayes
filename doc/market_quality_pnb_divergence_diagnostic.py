#!/usr/bin/env python3
"""Read-only, prespecified-case diagnostic of retained native PNB CSV files.

This does not fit, tune or change the model. CSV precision limits numerical
interpretation; ranks of two divergent states do not establish their cause.
"""
import bisect
import csv
import math
import os
from pathlib import Path
import statistics

repo = Path(os.environ.get("MQ_REPO", "/mnt/d/projects/antitrustbayes"))
base = repo / "doc/statistical_toolkit_runs/20260916/performance-audit/recovery/market-quality"
source = base / "pnb-archived"
out = base / "pnb-divergence-diagnostic"
out.mkdir(parents=True, exist_ok=True)
parameters = ["sigma_b_strat", "sigma_share", "sigma_margin", "rho",
              "s0_prob.1", "alpha_structural_market.1"]
values = {name: [] for name in parameters}
divergent = []
chains = []
for chain in range(1, 5):
    path = source / f"B-stan_{chain}.csv"
    with path.open() as handle:
        reader = csv.reader(line for line in handle if not line.startswith("#"))
        header = next(reader)
        positions = {name: header.index(name) for name in parameters + ["divergent__"]}
        count = 0
        chain_divergences = 0
        for row in reader:
            if not row:
                continue
            count += 1
            if count <= 2000:
                continue
            record = {name: float(row[positions[name]]) for name in parameters}
            if not all(math.isfinite(value) for value in record.values()):
                raise ValueError(f"Nonfinite retained parameter in {path}, row {count}")
            for name in parameters:
                values[name].append(record[name])
            if float(row[positions["divergent__"]]) != 0:
                chain_divergences += 1
                divergent.append({"chain": chain, "postwarmup_iteration": count - 2000,
                                  **record})
        if count != 8000:
            raise ValueError(f"Native CSV differs from frozen PNB budget: {path}, {count}")
        chains.append({"chain": chain, "iterations": count, "retained": count - 2000,
                       "divergences": chain_divergences})

def write(name, rows):
    if not rows:
        raise ValueError(f"Empty diagnostic output: {name}")
    with (out / name).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

summary = []
for name in parameters:
    ordered = sorted(values[name])
    for record in divergent:
        record[f"{name}_posterior_percentile"] = 100 * (
            bisect.bisect_left(ordered, record[name]) +
            bisect.bisect_right(ordered, record[name])) / (2 * len(ordered))
    quantiles = statistics.quantiles(ordered, n=40, method="inclusive")
    summary.append({"parameter": name, "mean": statistics.mean(ordered),
                    "sd": statistics.stdev(ordered), "q025": quantiles[0],
                    "median": statistics.median(ordered), "q975": quantiles[-1],
                    "draws": len(ordered)})
write("parameter-summary.csv", summary)
write("chain-divergences.csv", chains)
if divergent:
    write("divergent-state-percentiles.csv", divergent)
correlations = []
for first, second in [("sigma_b_strat", "sigma_share"),
                      ("rho", "sigma_share"), ("rho", "sigma_margin"),
                      ("alpha_structural_market.1", "s0_prob.1")]:
    x, y = values[first], values[second]
    mx, my = statistics.mean(x), statistics.mean(y)
    cross = sum((a - mx) * (b - my) for a, b in zip(x, y))
    denom = math.sqrt(sum((a - mx) ** 2 for a in x) * sum((b - my) ** 2 for b in y))
    correlations.append({"first": first, "second": second, "correlation": cross / denom})
write("parameter-correlations.csv", correlations)
print("Retained draws:", len(values[parameters[0]]), "divergences:", len(divergent))
for record in divergent:
    print(record)
print("Correlations:", correlations)
