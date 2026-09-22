# Enterprise market-intercept benchmark

This harness tests whether analytically integrating the Gaussian market
intercept is a robust sampler improvement across demand systems, conduct
models, panels, exact ordered-k mixtures, missing margins, and the Supreme
Court banking data. It does not change priors, conduct equations, likelihoods,
or sampler budgets.

The frozen manifest contains 124 independent tasks:

| Tier | Tasks | Purpose |
|---|---:|---|
| `gate` | 32 | Collapsed-intercept Logit/CES × B/C/MC/A/BF/CF/BS/CS, each as a 10-market cross-section and four-period panel |
| `replicate` | 32 | A prespecified second sampler seed for every synthetic gate cell |
| `comparison` | 16 | Matched sampled-intercept controls for B/C/BS/CS |
| `sc` | 16 | Collapsed-intercept eight-model Supreme Court cross-section and panel regressions |
| `sc_comparison` | 4 | Sampled-intercept B/C controls on both Supreme Court fixtures |
| `stress` | 8 | Missing-margin, low-share-noise, high-rho, and unbalanced synthetic panels |
| `scale` | 16 | Optional one-chain timing probes for Logit/CES × B/BS at 10, 100, 500, and 2,000 markets |

The synthetic cells use persistent centered firm-market quality, active
share-margin correlation, common residual scales, market-level outside-share
pooling, and fixed market-level ordered-k ranks. They apply a normalized
compositional share perturbation. Because the production likelihood is on log
conditional shares rather than a simplex density, these are deliberately
misspecified robustness fixtures, not exact iid-Gaussian observations. The
ordinary MC and auction margins use the same public conduct schedules as the
package. Rich margins satisfy the relevant FOC schedule at fixed observed
prices; they do not constitute a full marginal-cost/equilibrium round trip and
therefore do not certify counterfactual economics. Fixed-target oracle tests
provide the exact density-equivalence evidence. The Supreme Court cells are
real-data regression and weak-identification stress tests.

Every array element reads one immutable manifest row and writes to its own
task directory. Attempts are never overwritten. Fit objects and posterior
draws are never saved. The aggregator distinguishes missing, failed, and
completed tasks and validates run keys and the manifest hash.

The primary HMC gate is unchanged:

* zero divergences;
* zero maximum-treedepth hits;
* every-chain BFMI at least 0.30;
* maximum rank-normalized R-hat at most 1.01;
* minimum bulk and tail ESS at least 100;
* exact-k probabilities sum to one within `1e-10`.

## Local smoke check

The default runner mode validates fixtures and argument construction without
loading Stan or sampling. Generate a run manifest, then run any task:

```text
MIC_ENTERPRISE_RUN_DIR=/tmp/mic-enterprise \
  Rscript --vanilla doc/market_intercept_enterprise/make_manifest.R

MIC_ENTERPRISE_RUN_DIR=/tmp/mic-enterprise MIC_ENTERPRISE_TASK_ID=1 \
  Rscript --vanilla doc/market_intercept_enterprise/run_task.R

MIC_ENTERPRISE_RUN_DIR=/tmp/mic-enterprise \
  Rscript --vanilla doc/market_intercept_enterprise/aggregate.R
```

## Cluster execution

Create the manifest once on shared storage. Submit tiers separately so a
failed early gate does not spend the later budget:

```text
sbatch --array=1-32%8 \
  --export=ALL,MIC_ENTERPRISE_REPO=/path/to/antitrustBayes,MIC_ENTERPRISE_RUN_DIR=/path/to/run,MIC_ENTERPRISE_LIB=/path/to/R/library \
  doc/market_intercept_enterprise/run_array.sbatch
```

The current task ranges are `1-32` (gate), `33-64` (second-seed replicate),
`65-80` (comparison), `81-96` (Supreme Court), `97-100` (Supreme Court
comparisons), `101-108` (stress), and `109-124` (optional timing scale probes).
Recreate the manifest and inspect `tier_ranges.csv` rather than relying on
these ranges after any versioned design change.

After the array finishes, aggregate with `MIC_ENTERPRISE_STRICT=true` to make
missing or failed tasks an error. A retry requires `MIC_ENTERPRISE_RETRY=true`
and creates a new numbered attempt, preserving the failure.

For a staged gate, restrict strict aggregation to the submitted tier so later
tiers are not treated as missing:

```text
MIC_ENTERPRISE_RUN_DIR=/path/to/run MIC_ENTERPRISE_TIERS=gate \
  MIC_ENTERPRISE_STRICT=true \
  Rscript --vanilla doc/market_intercept_enterprise/aggregate.R
```

This writes `results-gate.csv`, `aggregation_audit-gate.csv`, and the matching
provenance audit. Omit `MIC_ENTERPRISE_TIERS` only after the complete manifest
has run.
