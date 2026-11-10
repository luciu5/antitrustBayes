# History log-share variance time-trend test (2026-11-09 UTC)

The authorized history-only test completed. Adding one regularized variance trend
improves the PNB 1960 baseline and part of the historical likelihood loss, with
stable sampling and unchanged margin fit. It does **not** fully restore historical
statistical fit. Production models, defaults and manuscript results are unchanged.

## Model and execution

```
sigma_share_i,t^2 = sigma_1960^2 * A_reference / institution_assets_i,t
                    * exp(gamma * (year - 1960) / 5)
gamma ~ Normal(0, 1)
```

The existing inverse-assets rule, global asset scaling and all original priors
remain fixed; the base residual-scale prior is explicitly at 1960 and reference
assets. Time runs from 0 to 1 over 1960–65. Gamma is the full-period log-VARIANCE
change, so the SD multiplier is `exp(.5*gamma*time)`. This adds exactly one HMC
coordinate and one normalized prior; the original total/fraction Jacobian is
unchanged. Quality stays analytically collapsed. There is no demand-mean trend,
new demand shock or extra measurement component. The complete vector SD enters
the contrast likelihood, share–margin covariance and conditional quality recovery.
Refitting can change 1960 parameters and quality even though time in 1960 is zero.

One new Bertrand history fit used four chains, 500 warmup + 500 retained draws per chain,
at most four cores, one thread per chain/numerical library, adapt_delta .999 and
max_treedepth 13. The 90-minute execution clock began 05:56:38 UTC; the controller
completed 06:00:37 UTC (3m58s), including preflight/compilation. There were no fit
retries, bridges, model averaging, cluster jobs or dependency refreshes.

Both constant-variance and inverse-assets history references passed exact data,
attributes, priors, source/kernel and original RDS/CSV integrity checks. Their
100 policy-draw realizations were copied byte-for-byte; snapshots were not refitted
or replayed. The new candidate uses the same draw-index seed and 100 policy draws.
Every comparison cell and failure/invalid-draw field is retained.

Development source is antitrustBayes 0.3.0.9001 at base 8afbd535d891a4385e69a30d6482bd15528378db
plus preserved local changes, not a clean Git build. The manifest records 16 actual
source/dependency/input hashes and package paths/versions. The existing library
provides antitrust 0.99.43.9000, coordination 0.1.0, CmdStanR 0.9.0.9002 and CmdStan 2.39.0.
No production R/Stan, specification 14, model weights or paper estimates changed.

## Results

Solved baseline mean shares (%), fixed PNB 1960 target:

| Bank | Observed | Constant variance | Inverse assets | Inverse assets + time |
| --- | ---: | ---: | ---: | ---: |
| First Pennsylvania |22.96|22.02|21.66|22.70|
| PNB |21.05|19.10|18.76|19.82|
| Girard |14.23|14.11|13.90|14.17|
| Three banks combined |58.24|55.23|54.32|56.69|

Candidate individual mean errors are −0.25, −1.23 and −0.06 percentage points, within the
requested 2–3 pp tolerance. Combined error improves from −3.92 to −1.55 pp versus the
inverse-assets reference. At observed rates the candidate means are 22.71/19.61/14.20,
combined 56.53, so the improvement already exists before solving equilibrium.

| Metric | Constant variance | Inverse assets | Inverse assets + time |
| --- | ---: | ---: | ---: |
| All-bank margin point RMSE |.42031|.42153|.42039|
| Margin predictive RMSE |.61921|.62174|.61787|
| Solved-share draw RMSE (pp) |.5212|.5541|.2814|
| Slope per original rate unit |.73007|.71505|.72824|
| Whole-market Paasche mean (%) |-14.35|-13.83|-15.24|
| Merging-party Paasche mean (%) |-50.33|-49.97|-52.04|
| Mean collapsed share log likelihood |-56.78|-178.54|-150.32|
| Mean collapsed joint log likelihood |-183.86|-306.63|-275.83|

Draw RMSE includes posterior dispersion; it is not RMSE of posterior means. No
small-bank accuracy gate was imposed. Margin fit and the slope are essentially
unchanged from the constant-variance reference. All likelihood entries are
normalized in-sample posterior-average diagnostics, **not** marginal likelihoods,
Bayes factors or model probabilities. The trend recovers 30.8 of the 122.8 joint
log-likelihood loss introduced by fixed inverse assets, leaving a 92.0 gap versus
constant variance. Independent dense calculations confirm that this is real.

Gamma posterior mean is 2.696 (95% interval 2.084–3.352). The posterior mean variance
time multiplier `exp(gamma)` is 15.69 (8.04–28.55), and mean SD multiplier is 3.90
(2.84–5.34). These hold institution assets fixed. An institution's actual variance
ratio additionally multiplies by its 1960/1965 asset ratio. These are transformed
draw summaries, not exponentials of the posterior mean. The strong trend exceeds
the central 95% range of its prior and should not be described as a mild adjustment.

The saved-state year-by-year residual check still shows excess dispersion near
both endpoints and lower dispersion mid-window. A monotonic trend therefore
does not absorb the full pattern, and cannot change relative size precision
within a year. The result supports keeping this sensitivity for further validation;
it does not establish that the fixed inverse-assets-plus-trend model is an adequate
common default. No extra specifications were fitted after inspecting the result.

## Verification and artifacts

Candidate maximum R-hat 1.00378, minimum bulk ESS 884 and tail ESS 1038; zero
post-warmup divergences or depth-limit hits. Three rejected warmup proposals
reported zero Gaussian scale; logs are preserved. Constant generated quantities
have undefined R-hat and are counted separately. All 300 comparison policy states
and both indices are valid; all no-op effects are zero. Candidate maximum independent
FOC residual is 1.4271e-8 (existing tolerance 1e-6). Every Paasche denominator uses
that draw's solved pre-merger rates. Cost/value-of-funds conventions and signed
rate handling are unchanged.

Eleven deterministic expectations passed with no failures/warnings/skips, covering
the log-variance/SD factor, 1960 anchoring, unchanged prior parameter block except
one coefficient, and dense conditional quality with time/size covariance. Stanc
syntax passed. A one-off exact `identical` check of the numeric time grid failed
on a 1.11e-16 representation difference from R's `seq`; a numerical check passed.
No saved-fit data/source compatibility tolerance was relaxed. The ordinary package
suite was not repeated for this isolated research-only change.

Independent review verified all 100 new states: SD identity error<=8.34e-15,
exact seeded dense quality replay error<=5.66e-13, direct-share error<=5.56e-17.
Eight independent pre/post roots reproduce Paasche ratios within1.91e-8.
Three posterior rows per each history fit match dense normalized Gaussian GQ
likelihoods within 5.12e-12. All 16 source guards, two reference hashes, two copied
draw files and 18 preserved snapshot artifact seals pass. Audit scripts/receipts
are in the output's `independent_audit/` directory.

```sh
python3 doc/validation/run_pnb_variance_time_trend_diagnostic.py --plan
R_LIBS_USER=../mergerBayes/.scratch/paper-library Rscript --vanilla doc/validation/variance_time_trend_prototype_checks.R
R_LIBS_USER=../mergerBayes/.scratch/paper-library Rscript --vanilla doc/validation/report_pnb_variance_time_trend_diagnostic.R .scratch/variance_time_trend_diagnostic_20261109
```

The module/config/runner/collector live alongside this report. Raw outputs are in
`.scratch/variance_time_trend_diagnostic_20261109/`: clock/completion receipt,
source/data manifest, generated kernel/evaluator, one native fit/CSV set, original
reference pointers/copied draw plans, complete comparison and parameter summaries,
time multipliers, year dispersion, likelihood/fit/validity diagnostics and seals.
All preceding fits and reports remain preserved. No changes were pushed.
