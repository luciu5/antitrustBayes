# Bounded institution-assets variance test (2026-11-09 UTC)

The requested PNB snapshot/history Bertrand comparison completed. This particular
heteroskedastic specification **does not solve baseline shrinkage**: it substantially
worsens snapshot shares and provides little improvement for the three large banks
with history. It is a separate research prototype; production kernels, priors,
simulation defaults, model specification 14 and manuscript results are unchanged.

## Specification and institution assets

Only residual scale changes:

```
log(sigma_share_i) = log(sigma_share) + gamma * institution_log_assets_z_i
gamma ~ Normal(0, .35)
```

`total_assets` is the **contemporaneous institution-wide bank-year total** already
present in the corrected banking input, sourced from `p.assets_avg / 1000` in
`mergerBayes/code/banking/supreme.R`. It is not local-market assets, deposits,
market share or lagged assets. Repeated bank-years across cases must have exactly
the same assets and are deduplicated for scaling. The resulting 1,090 unique
bank-years across all cases determine one mean log assets (3.292770) and SD
(1.535360). The same transformation is used in snapshot and history. Monetary
unit changes leave the standardized covariate unchanged. **Scaling uses all
cases; the coefficient is estimated using PNB only**, as requested.

The two-sided gamma prior allows either size relationship; gamma zero nests the
current model. Its scale is a declared diagnostic regularization, not a calibrated
bank-size precision law. The existing quality prior, base residual-scale prior
(now the scale at global standardized size zero), outside-share prior, total/fraction
Jacobian, margin likelihood and rho prior remain unchanged. Quality stays collapsed.
Heteroskedastic scales enter both the normalized relative-contrast covariance and
the share-on-margin conditional mean `rho * sigma_share_i / sigma_margin`.
There is no added period demand shock or measurement-error decomposition. Policies
retain the recovered persistent quality and discard the share residual as before.

## Execution, provenance and validity

The whole-process clock ran from **05:04:28 to 06:34:28 UTC**. The controller
completed at **05:08:12 UTC (3m44s)**; compilation was included. Three new
fits ran sequentially (fresh snapshot reference, snapshot candidate, history
candidate), each with four chains, 500 warmup plus 500 retained draws per chain,
one numerical thread and at most four cores. The existing four-chain/500-retained
history reference was reused only after exact data attributes/values, unchanged
production kernel and original RDS/CSV seals passed. Nothing was resealed to make
an incompatible fit acceptable. The original shorter snapshot reference was not
used in this comparison. No bridge sampling, averaging, cluster job, package
installation or remote write occurred.

Source revision is `8afbd535d891a4385e69a30d6482bd15528378db` plus the preserved
dirty development checkout, **not a clean Git build**. The manifest records the
actual experiment/input/critical R/dependency source hashes and source paths.
The existing shared library provides antitrust 0.99.43.9000, coordination 0.1.0,
CmdStanR 0.9.0.9002 and CmdStan 2.39.0; antitrustBayes 0.3.0.9001 is loaded from
development source. No production fit is relabeled as this prototype.

| Cell | Maximum R-hat | Minimum bulk ESS | Post-warmup divergences | Depth-limit hits |
| --- | ---: | ---: | ---: | ---: |
| Snapshot reference | 1.00635 | 802 | 0 | 0 |
| Snapshot assets | 1.01043 | 606 | 0 | 0 |
| History reference | 1.00760 | 900 | 0 | 0 |
| History assets | 1.00573 | 943 | 0 | 0 |

These summarize all finite reported-quantity diagnostics, including `lp__`;
constant generated quantities have undefined R-hat and are counted separately.
The snapshot candidate's `lp__` R-hat exceeds the 1.01 pilot criterion slightly:
its convergence classification remains **inconclusive**, despite zero divergences.
History logs retain sporadic rejected warmup proposals from zero-scale Gaussian
evaluations. Stanc syntax passed; pedantic analysis emitted control-flow and
data-dependent-prior warnings from this general-purpose kernel. No stronger
universal sampler-stability claim is made.

Each cell uses one reproducible 100-draw parameter/quality plan. All **400/400**
baseline/merger states, market indices and party indices are valid; no draws were
dropped. All no-op effects are exactly zero. The largest production independent
input-Bertrand FOC check is `1.487219e-8`, below `1e-6`. Paasche denominators use
each draw's solved pre-merger rates. Signed rates and existing validity rules are
preserved. There are no estimated conduct probabilities in a Bertrand-only pilot.

## Results

All shares below are posterior means in percent for the fixed PNB 1960 target.
The complete 16-row comparison CSV includes 95% intervals, shares at observed
rates and after solving equilibrium, bank margin errors and policy intervals.

| Bank | Observed | Snapshot reference | Snapshot assets | History reference | History assets |
| --- | ---: | ---: | ---: | ---: | ---: |
| First Pennsylvania | 22.96 | 18.67 | 2.08 | 22.02 | 22.02 |
| PNB | 21.05 | 17.37 | 2.20 | 19.10 | 19.18 |
| Girard | 14.23 | 9.12 | 1.94 | 14.11 | 14.13 |
| Three banks combined | 58.24 | 45.17 | 6.22 | 55.23 | 55.33 |

The inferred gamma is **positive** in both settings: snapshot mean `.7771`,
95% interval `[.5335, 1.0347]`; history `.1938`, interval `[.0899, .2965]`.
No retained gamma draws are negative. Thus this fit allocates **more** residual
variance to large institutions. Snapshot posterior quality SD drops from 1.533
to .616. Mean residual SDs for First Pennsylvania/PNB/Girard become approximately
6.91/6.65/5.42 in log units. Their persistent quality consequently shrinks heavily
when the counterfactual excludes the residual. This is a conditional-model
allocation, not independent evidence that large-bank accounting measurement is
less accurate. A sign constraint chosen to force closer shares was not tested.

| Metric | Snapshot reference | Snapshot assets | History reference | History assets |
| --- | ---: | ---: | ---: | ---: |
| All-bank log-inverse-margin point RMSE | .42543 | .43501 | .42031 | .41995 |
| All-bank margin predictive RMSE | .59838 | .61069 | .61921 | .61647 |
| All-bank solved-share draw RMSE (pp) | 2.9869 | 5.8311 | .5212 | .6140 |
| Slope per original rate unit | .73804 | .82512 | .73007 | .73057 |
| Whole-market Paasche mean (%) | -8.18 | -.09 | -14.35 | -14.84 |
| Merging-party Paasche mean (%) | -37.60 | -2.86 | -50.33 | -51.94 |

The draw RMSE integrates posterior dispersion; it is not RMSE of posterior mean
shares. Large-bank assessments above use actual posterior mean errors. History's
near-unchanged large-bank means coexist with somewhat worse all-bank share errors.
The saved history reference uses a new fixed policy draw plan in this experiment,
so its 100-draw policy summaries need not equal the previous pilot's summaries.

The snapshot specification fails the requested large-bank 2–3 pp tolerance and
materially changes policy effects. History gives no compelling improvement.
These results do not support adopting institution-assets-dependent **residual**
variance in this form as the shrinkage fix. They do not rule out other observation
models, and do not resolve how much of the existing residual belongs in structural
demand. No primary results or defaults were changed.

## Verification and replay

Seventeen deterministic expectations passed with no failures, warnings or skips.
They cover global institution identity/scaling, unit and market-duplication
invariance, a separately constructed dense normalized Gaussian contrast density
and conditional quality with partial margins/nonzero rho, exact same-seed recovery
at gamma zero, and correct covariance/extraction routes. Independent review also
verified a dense conditional mean to `5.55e-17` and exact constant-scale parity.
The subsequent independent saved-state audit covered every candidate state (200
total): asset-scale identity errors are at most `6.74e-15`, exact seeded dense
conditional-quality replay errors at most `2.40e-13`, and direct demand-share
errors at most `5.56e-17`. Sixteen independently solved pre/post equilibrium
paths reproduce Paasche ratios within `1.97e-8`. The adverse snapshot result is
therefore reproduced independently; no extraction or state-dimension defect was
found. The audit script and receipt are preserved in the output's
`independent_audit/` directory.
This research-only change did not trigger a repeat of the ordinary package suite.

```sh
# Read-only plan:
python3 doc/validation/run_pnb_asset_variance_diagnostic.py --plan
# Deterministic checks, with the existing compatible library selected:
Rscript --vanilla doc/validation/asset_variance_prototype_checks.R
# Saved-output collection; no inference:
Rscript --vanilla doc/validation/report_pnb_asset_variance_diagnostic.R .scratch/asset_variance_diagnostic_20261109
```

Configuration, kernel builder, conditional-quality adapter, runner and collector
live alongside this report. Durable artifacts are under
`.scratch/asset_variance_diagnostic_20261109/`: original clock and manifest,
institution asset index/scaling, input data, generated kernel and executed policy
evaluator, three native fit/CSV sets, four fixed policy draw plans, per-cell
diagnostics/scales/parameter summaries, full comparison, fit metrics, validity
summary and artifact seals. Historical demand-shock outputs remain intact.
