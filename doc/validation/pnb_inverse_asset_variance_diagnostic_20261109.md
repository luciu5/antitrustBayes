# Fixed inverse institution-assets variance pilot (2026-11-09 UTC)

The requested fixed inverse-variance comparison completed with two new Bertrand
fits and two compatible saved references. It **substantially improves snapshot
baseline shares without new sampler parameters or divergences**, but worsens the
fit to historical share variation. This is a research diagnostic, not an adopted
package default or replacement manuscript result.

## Declared specification

```
sigma_share_i^2 = sigma_share^2 * A_reference / total_assets_i
A_reference = exp(mean(log(total_assets))) over unique bank-years across all cases
```

Assets are contemporaneous institution-wide totals from the existing corrected
banking input, not local-market assets, deposits or shares. The same 1,090 unique
bank-years and global transformation used in the earlier unrestricted asset pilot
are reused: reference assets `26.9173267` in the input's units; log-assets SD
`1.5353604`. The corresponding coefficient on standardized log assets is the
**fixed** `-.7676802`, not `-.5`. Doubling assets halves variance and reduces SD
by `sqrt(2)`. The coefficient is not estimated and has no prior.

The original Stan sampler parameter block is unchanged. Existing priors for
quality, the residual scale at reference assets, margins, outside share and rho,
the total/fraction Jacobian and collapsed-quality treatment are retained. The
known asset rule changes the implied residual-scale prior for each institution.
Vector scales enter the normalized relative-share covariance, the joint
share-on-margin conditioning and conditional quality recovery consistently.
Counterfactuals retain recovered quality and omit the residual, as in the current
production model. There is no new period demand shock, extra measurement component
or observed-share anchoring. The fixed target is PNB 1960; history remains the
existing retrospective 1960–1965 information set, not an out-of-sample forecast.

## Execution and provenance

Two candidate fits ran sequentially, each with four chains, 500 warmup and 500
retained draws per chain, one thread per chain/numerical library and at most four
cores. Both original four-chain/500-retained reference fits passed exact base-data
and attribute comparisons, source/kernel checks and original RDS/CSV seals.
Reference policy draw files are **byte-identical** to the preceding asset pilot;
their values were not re-simulated. Each candidate uses the same reproducible
100-row draw-index seed and retains its realized parameters/quality through both
equilibria. Every planned cell remains in the result table.

The fixed absolute 90-minute budget began **05:28:46 UTC**; compilation and all
workers finished **05:32:50 UTC (4m04s)**. There were no retries, new reference fits,
bridge sampling, model averaging, cluster jobs, dependency refreshes or remote
writes. No updated model probabilities were calculated. Historical probabilities
are preserved records, not weights estimated for this changed likelihood.

The development package is antitrustBayes 0.3.0.9001 at base revision
`8afbd535d891a4385e69a30d6482bd15528378db` plus preserved local changes. This is
not a clean Git build. The manifest records 15 actual source/dependency/input
hashes, package paths/versions, CmdStan path/version and original reference SHA
values. The existing shared library provides antitrust 0.99.43.9000,
coordination 0.1.0, CmdStanR 0.9.0.9002 and CmdStan 2.39.0. Production R/Stan,
model specification 14, public defaults and manuscript results are unchanged by
this experiment.

| Cell | Maximum R-hat | Minimum bulk ESS | Post-warmup divergences | Depth-limit hits |
| --- | ---: | ---: | ---: | ---: |
| Snapshot reference | 1.00635 | 802 | 0 | 0 |
| Snapshot inverse assets | 1.00461 | 873 | 0 | 0 |
| History reference | 1.00760 | 900 | 0 | 0 |
| History inverse assets | 1.00330 | 855 | 0 | 0 |

These cover all finite reported-quantity diagnostics; constant generated quantities
and the fixed coefficient have undefined R-hat and are counted separately. The
candidate logs preserve two snapshot and four history rejected warmup proposals
from zero-scale Gaussian evaluations. There were no failed fits or invalid
retained policy states. The bounded results support stability for these two fits,
not a guarantee about every market or model.

All **400/400** baseline/merger states and both Paasche indices are valid. All
no-op effects are exactly zero; maximum independent production FOC residual is
`1.487219e-8` against the existing `1e-6` tolerance. Costs/value-of-funds inversion,
signed-rate policy, outside normalization and index validity conventions are
unchanged. Paasche denominators use each draw's solved pre-merger rates.

## Share and policy results

Solved baseline posterior mean shares (%):

| Bank | Observed | Snapshot current | Snapshot inverse assets | History current | History inverse assets |
| --- | ---: | ---: | ---: | ---: | ---: |
| First Pennsylvania | 22.96 | 18.67 | 22.35 | 22.02 | 21.66 |
| PNB | 21.05 | 17.37 | 20.46 | 19.10 | 18.76 |
| Girard | 14.23 | 9.12 | 13.27 | 14.11 | 13.90 |
| Three banks combined | 58.24 | 45.17 | 56.08 | 55.23 | 54.32 |

Snapshot large-bank mean errors become `-.61/-.60/-.96` percentage points, within
the requested 2–3 pp tolerance. The combined error improves from `-13.07` to
`-2.16` pp. At observed rates the corresponding means are `22.27/20.40/13.19`,
combined `55.87`; the solved baseline improvement therefore starts in demand
reconstruction rather than being created by the equilibrium solver.

History mean errors are `-1.29/-2.29/-.33` pp, versus `-.94/-1.95/-.12` in the
reference. The combined error worsens from `-3.01` to `-3.92` pp. History meets a
3 pp individual-bank tolerance but does not improve its existing baseline.
Snapshot inverse combined-share 95% interval is `[49.65,59.20]`; history is
`[53.08,55.29]`, below the observed combined share. These are conditional posterior
intervals under the imposed precision rule, not evidence validating that rule.

| Metric | Snapshot current | Snapshot inverse assets | History current | History inverse assets |
| --- | ---: | ---: | ---: | ---: |
| All-bank log-inverse-margin point RMSE | .42543 | .42664 | .42031 | .42153 |
| Margin predictive RMSE | .59838 | .60946 | .61921 | .62174 |
| All-bank solved-share draw RMSE (pp) | 2.9869 | .8738 | .5212 | .5541 |
| Slope per original rate unit | .73804 | .71372 | .73007 | .71505 |
| Whole-market Paasche mean (%) | -8.18 | -14.33 | -14.35 | -13.83 |
| Merging-party Paasche mean (%) | -37.60 | -52.28 | -50.33 | -49.97 |

The draw RMSE includes posterior dispersion and is not RMSE of posterior means.
No small-bank accuracy threshold was introduced. Margin point RMSE increases
about .3% in each information set; predictive RMSE rises 1.9% snapshot and .4%
history. Slopes move about -3.3% and -2.1%, respectively. Snapshot policy effects
change materially; history policy means move little. Inverse candidate 95%
Paasche intervals are market `[-20.03,-8.05]`, parties `[-70.52,-36.15]` in the
snapshot, and market `[-16.16,-11.03]`, parties `[-57.81,-40.38]` with history.

## Statistical-fit limitation

Posterior mean normalized collapsed log likelihoods, using all retained draws:

| Quantity | Snapshot current | Snapshot inverse assets | History current | History inverse assets |
| --- | ---: | ---: | ---: | ---: |
| Share equation | -78.69 | -79.25 | -56.78 | -178.54 |
| Joint shares and margins | -97.76 | -99.73 | -183.86 | -306.63 |

These are **in-sample posterior-average likelihood diagnostics**, not marginal
likelihoods, Bayes factors or model probabilities. They nevertheless expose a
substantial historical fit penalty under the fixed inverse pattern, consistent
with the earlier unrestricted pilot's positive size coefficient. Closer snapshot
baseline means do not establish better statistical prediction or validate this
variance law. The snapshot benefit and historical penalty must both be retained
when judging this candidate. A weaker or estimated inverse relationship is
untested. This result does not support making the fixed rule a common default
across snapshot and history on the basis of baseline shares alone.

## Verification and reproducibility

Fourteen deterministic expectations passed without failures, warnings or skips:
inverse variance rather than inverse SD, currency-unit invariance, unchanged Stan
parameter block/no extra prior, and independent dense conditional quality with
partial margins and nonzero rho. Stanc syntax passed. The ordinary package suite
was not repeated for this isolated research-only change.

Independent saved-state review covered all 200 candidate states: inverse-scale
identity error at most `8.17e-15`, exact seeded dense conditional-quality replay
error at most `3.78e-13`, and direct share error at most `8.33e-17`. Sixteen
independent pre/post roots reproduce Paasche ratios within `3.72e-9`. Both reused
draw files, 15 source guards and reference hashes pass. Audit scripts/receipts
are preserved in the output's `independent_audit/` directory.
An additional independent dense likelihood audit checks three posterior draws
from each history fit, using all 204 observations, 198 share contrasts and 202
usable margins. It matches Stan share/joint log likelihoods within `3.87e-12`;
the historical fit penalty is not a covariance or normalization defect.

```sh
# No inference or writes:
python3 doc/validation/run_pnb_inverse_asset_variance_diagnostic.py --plan
# Deterministic checks using the existing compatible library:
R_LIBS_USER=../mergerBayes/.scratch/paper-library Rscript --vanilla doc/validation/inverse_asset_variance_prototype_checks.R
# Collect saved outputs; no inference:
R_LIBS_USER=../mergerBayes/.scratch/paper-library Rscript --vanilla doc/validation/report_pnb_inverse_asset_variance_diagnostic.R .scratch/inverse_asset_variance_diagnostic_20261109
```

Configuration, kernel adapter, runner, checks and collector live alongside this
report. Raw artifacts are in `.scratch/inverse_asset_variance_diagnostic_20261109/`:
clock/completed receipt, source/data manifest, asset index/scaling, generated
kernel and executed evaluator, two new native fits/CSV sets, immutable reference
pointers/copied draw plans, all four policy outputs/diagnostics, full 16-row
comparison with uncertainty intervals, fit metrics, likelihood diagnostic,
validity summary and artifact seals. Every earlier diagnostic remains preserved.
