# Residual-model integration: executed bounded validation

Executed 2026-11-09 against development specification 15. Exactly six new HMC
fits completed: `combined` and `demand_shock` for each of three frozen synthetic
cases. All six passed the declared sampling criteria; all 400 final native
target checks and all 2,200 policy unit/draw checks passed. No bridge sampling,
BMA, conduct selection, cluster jobs, rescue fits, or paper estimates were run.

The new model improves prediction of the observed shares in these cases.
That is not uniform recovery of true shares: snapshot truth-share error became
larger while history and CES truth-share error became smaller. The variance
split remains prior dependent and is not separately identified by this check.

The absolute budget began at **13:51:24 UTC**, the scientific target was sealed
at **14:09:37**, and the final worker ended at **14:13:27**, before the fixed
**15:21:24** deadline. Preparation and compilation count toward that budget.
Each fit used four chains, 500 warmup and 500 retained iterations per chain,
four parallel chains, one thread per chain and numerical library,
`adapt_delta=.999`, depth 13, and the same paired seed. Each fitted model used
100 fixed posterior policy draws; panel checks apply those draws to every unit.

The [configuration](residual_model_integration_config.R),
[runner](residual_model_integration.R),
[native oracle](residual_model_integration_native.R), and
[policy report](residual_model_integration_report.R) retain the execution
contract. Raw artifacts remain under
`/mnt/d/Projects/antitrustBayes/.scratch/residual_model_integration_20261109`.

The snapshot has one market/year, history one market/four distinct years, and
CES two markets/three distinct years. Every unit has six inside products with
ownership `a,a,b,b,c,d` and an outside reference product. The fixtures reuse
`antitrust::synthetic_market(mode="primitives")`; native equilibrium prices
and shares, the generator's Logit FOC diagnostics, and independent CES margin
identities agree before observation error.
Persistent firm quality and period/product demand shocks are nonzero;
their generating SDs are .45 and .25. Measurement SD is .10, margin-error SD
.12, and demand/margin correlation .35. Gaussian log inverse-margin noise is
matched to the implemented likelihood at observed firm shares. The CES case
retains separate market quality blocks. All outcomes and truth were frozen
before sampling. These are bounded working-likelihood fixtures, not a
repeated-sample calibration study.

| Case / residual model | Maximum R-hat | Minimum bulk ESS | Divergences / depth hits |
| --- | ---: | ---: | ---: |
| Input Logit snapshot / combined | 1.00599 | 840 | 0 / 0 |
| Input Logit snapshot / demand shock | 1.00435 | 935 | 0 / 0 |
| Input Logit history / combined | 1.00327 | 840 | 0 / 0 |
| Input Logit history / demand shock | 1.00602 | 1127 | 0 / 0 |
| Output CES panel / combined | 1.00803 | 573 | 0 / 0 |
| Output CES panel / demand shock | 1.00933 | 709 | 0 / 0 |

The gate is zero divergences/depth hits, maximum R-hat at most 1.01, and
minimum bulk ESS at least 400 across all active sampled parameters.
[Diagnostics](residual_model_integration_diagnostics.csv) retain the exact values.
Structural truth (Logit alpha 3; CES gamma 6) lies in all eight reported
market-level 95% intervals. Outside-share truth .20 lies in all 22 reported
cell intervals. Full [parameter and covariance intervals](residual_model_integration_parameters.csv)
include the generating SD, total variance, and share/margin covariance.

The table below uses posterior **mean solved-pre conditional firm shares**,
aggregating products by owner before comparison. Maxima include every firm
and every market/year; they do not select favorable banks or units. Errors are
percentage points. The CSVs also retain predictions at the observed prices.

| Case | Max observed-share error, combined → new | Mean absolute observed-share error | Max true-share error | Mean absolute true-share error |
| --- | ---: | ---: | ---: | ---: |
| Snapshot | 2.297 → .213 | 1.262 → .106 | 1.712 → 2.557 | .972 → 1.279 |
| History | 6.310 → 1.334 | 2.130 → .592 | 5.423 → 3.205 | 2.570 → 1.440 |
| CES panel | 8.327 → 1.332 | 1.817 → .249 | 5.590 → 2.021 | 1.441 → .642 |

Structural margin RMSE changes little: snapshot .03142 → .03127; history
.02040 → .02040; CES .02956 → .02961. These predictions average the structural
FOC margin over posterior draws and exclude replicated observation noise.
See [product shares](residual_model_integration_shares.csv),
[firm shares](residual_model_integration_firm_shares.csv), and
[margin predictions](residual_model_integration_margins.csv).

New-model measurement-SD medians and 95% intervals are .0988 [.0390,.2542],
.0897 [.0365,.1877], and .0982 [.0379,.2244]. The prior median is .10 with log
SD .5. These broad intervals do not establish measurement-versus-demand
identification. Total residual variance (.0725) and share/margin covariance
(.0105) lie in all six reported 95% intervals. For `combined`, the report maps
the sole residual SD and rho to their equivalent total-variance/effective-rho
truth; its fixed measurement SD is exactly zero.

The independent dense Gaussian oracle projects the full observation covariance
and includes collapsed persistent-quality covariance. It checks all four native
kernels, sampled/collapsed quality, positive/negative/zero rho, two measurement
scales, homogeneous/heterogeneous market scales, all five simultaneous conduct
models, and fringe/Stackelberg conditional role weights. Combined mode has no
measurement coordinate and exactly zero measurement variance. Maximum density
error is 1.60e-12; maximum role-probability error is 1.06e-14.
[All 400 checks](residual_model_integration_native_checks.csv) are retained.

The policies merge owners a and b with no efficiency change. All 2,200
unit/draw rows are retained; **zero are invalid**. Independent primitive-demand
FOCs have maximum absolute pre/post residual 1.12e-8, and every no-op price
error is zero. Paasche, Laspeyres and Fisher calculations use each draw's solved
pre-policy prices as their denominator. Posterior indices are descriptive
conditional policies, not estimated model probabilities.
[Policy counts and effects](residual_model_integration_policies.csv) summarize
the complete per-row files and saved 100-draw plans in the raw output.
An independent reviewer also finite-differenced each firm's joint multiproduct
profit for saved new-model draws 1, 50 and 100 in the history and CES cases;
the maximum gradient-versus-oracle discrepancy was 4.69e-11. That reviewer
independently counted 2,200 raw policy rows with zero invalid states and found
no substantive issue in the DGP, FOCs, CES normalization, index denominators
or reporting repair.

Two implementation mistakes in the validation infrastructure remain visible:

- The initial native validator had a helper-scope error in 32 cases. Its
  [failed rows](residual_model_integration_initial_native_failures.csv) and full
  raw initial table were preserved. Correcting the validator and replaying the
  unchanged compiled libraries produced the final 400 passes. Forty-four
  repeated-library-loading warnings (four unique messages) are retained in
  `native/warnings.txt`; source and library hashes were checked throughout.
- All six workers exited with status 1 after saving their completed HMC fits,
  diagnostics, plans, share predictions and all policy rows. The final margin
  report requested nonexistent `pred_margin`. The separate
  [reporting-only repair](residual_model_integration_reporting_repair.R) writes
  the missing tables using `1 / pred_marginNormInv`; it runs no fits or new
  policy draws. Original [execution statuses](residual_model_integration_status.csv),
  logs, source seal, fit bytes and plans remain unchanged. The repair receipt is
  `reporting_repair_receipt.rds` in the raw output. The exact originally executed
  reporter is archived as `reporting-sources/executed_residual_model_integration_report.R`
  there, with its original sealed SHA-256
  `91c1ad6b31c31549040ff9d73217de899dd133f945b599a436cf8596a47d1f23`.
  The repository reporter now uses the corrected inverse-margin expression for
  future separately authorized runs. `reporting_source_correction_receipt.rds`
  records the old/archive/current hashes and exact one-expression difference;
  the current reporter hash is
  `e8045a724e3670872a4a21081c5f7d0fc59da3ee88be72d3da02fe18544a9e22`.
  This correction changes no inference source, source guard, saved fit or draw.

The original execution-manifest SHA-256 is
`b08e7b302b31883bc6107deff978723ce13637d5a49323fdf8e09244c020d798`;
the repair-receipt SHA-256 is
`40661ddcdd77aefc3c99e080a5ec8c80fbc5b03228d3cf1d891d1f12283a1e45`.
The [source seal](residual_model_integration_source_seal.csv),
[fixture seal](residual_model_integration_fixture_seal.csv), and
[individual fit/data seals](residual_model_integration_fit_seals.csv) record the
dirty development checkout precisely. Prepared versus actual scientific Stan
data agree; the sole difference is deterministic computational grain size
(prepared 1; public API resolves 2/6/9 from N and cores). Both hashes are retained.
The compatible shared library used antitrust 0.99.43.9000, coordination .1.0,
CmdStanR .9.0.9002 and CmdStan 2.39.0; native evaluation used BridgeStan 2.8.0.

To inspect the saved evidence without HMC, from the repository root with that
library selected:

```r
pkgload::load_all(".")
source("doc/validation/residual_model_integration_collect.R")
rmi_collect(".scratch/residual_model_integration_20261109")
```

The execution CLI separates `prepare`, `seal`, `sample`, and `worker`; existing
fixtures, seals and attempted fits cannot be overwritten or silently retried.
The recorded run is finished. Re-running sampling, extending the deadline or
creating additional fits requires a separately authorized experiment.

## Package integration and compatibility checks

The production implementation is development version **0.3.0.9002**, model
specification **15**. New relative Gaussian Logit/CES estimation selects the
demand-shock model by default. `share_residual_model = "combined"` selects the
original model explicitly; absolute-share fits retain that original model.
The [statistical and policy contract](../share_residual_demand_state.md)
describes the new priors, complete demand state, identification limits and
relationship to Berry (1994).

The complete development-source non-HMC suite passed with 15 opt-in sampling
skips and six retained warnings. A complete installed-package `R CMD check
--no-manual`, including examples and built/rebuilt vignettes, finished with
**Status: OK** (zero check errors, warnings or NOTEs). Its testthat results
were **8,214 passing expectations, zero failures, six warnings and 33 skips**.
The additional installed skips cover repository-only research/documentation
checks and CRAN-disabled sampling. Warnings concern existing numerical-solver
messages and the documented positive-cost-shock convention; they were not
suppressed. Logs and the built source archive are in
`.scratch/residual_model_package_check/complete/`.

The first ordinary run exposed stale specification-14 expectations and
unversioned test fixtures; tests now distinguish new specification-15
construction from explicit legacy-14 compatibility. The first installed check
also exposed a repository-README extraction assertion in a source package that
intentionally excludes README. That assertion is now repository-only, while
independent public dispatch checks still run in installed packages. Initial
vignette warnings from the deliberately unbuilt-vignette archive were resolved
by building both vignettes. The direct workspace build was interrupted during
local-artifact traversal; the final build used a clean staging copy of package
directories. Staged runtime, installed-source, test and vignette files were verified
byte-identical to the current package files. A final generated-help EOF cleanup
changes no documentation content or executable behavior.

Fresh independent review found and resolved two compatibility defects: restored
Nutpie binaries must match the verified installed source hash while retaining
the original LF/CRLF seal, and historical version-2/3 total-fraction plans must
retain their pre-alias fingerprints. The latter accepts only an exact redundant
quality-SD alias; version-4 state validation remains strict. Four frozen legacy
kernels were verified byte-identical to their original sources. An actual
sealed specification-14 banking fit passed native source/data/evidence checks
and replayed its historical version-3 plan: all ten unchanged-ownership policies
were valid with zero effects. No old posterior, provenance or model probability
was rewritten.

An independent 20,000-draw joint Gaussian conditioning check compared production
quality-plus-demand draws with a separate full-covariance oracle: maximum mean
and covariance discrepancies were 2.125 and 1.636 Monte Carlo standard errors,
respectively, below the declared five-standard-error bound. No HMC was involved.
Focused policy checks also cover input/output Logit, CES, auction, bargaining,
native fringe/Stackelberg roles, missing margins, sampled/collapsed quality,
market covariance and sequential primitive shocks.

Sampling stability above was tested with the declared `adapt_delta=.999`;
the public sampler default remains .95. This bounded result is not a guarantee
of convergence for every dataset or control setting. No package release, push,
full banking rerun, bridge sampling or manuscript update was performed.
