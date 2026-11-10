# Hierarchical year log-share variance test (2026-11-09 UTC)

The requested history-only test completed. Hierarchical year variances improve
historical statistical fit more than the linear trend, with no divergence
increase. The 1960 baseline meets the stated large-bank tolerance; PNB remains
about 1.82 percentage points low. The specification is a research sensitivity,
not an adopted production model or manuscript result.

## Model and prior

```
sigma_i,t^2 = sigma_1960^2 * A_reference / institution_assets_i,t
              * exp(u_t - u_1960)
u = tau_year * H * z
z ~ Normal(0, I_5)
tau_year ~ HalfNormal(scale = 1/sqrt(2))
```

H is the normalized 6-by-5 Helmert basis: H'H=I and HH'=I-11'/6. Five noncentered
contrasts remove the redundant common year direction. With tau_year fixed, the
prior variance of a year difference is 2*tau_year^2; the declared half-normal
scale therefore makes the prior second moment of the 1960–1965 log-variance
change equal to one, matching the preceding Normal(0,1) linear-trend coefficient.
This matches that moment only, not the entire prior distribution.

The SD multiplier is exp(.5*(u_t-u_1960)), so the existing base residual-scale
prior remains at 1960 and reference assets. The six hierarchy sampler parameters
(five contrasts plus one positive SD) replace the single linear trend coefficient.
The half-normal density is normalized, and no extra tau Jacobian or density is
attached to the derived u. Original quality, base-scale, margin, outside-share and
rho priors, total/fraction Jacobian, collapsed quality and demand mean are unchanged.
The full observation-SD vector enters the contrast likelihood, share–margin
covariance and conditional quality consistently. There are no new bank-year
demand shocks, measurement components or observed-share anchoring.

**Prior-predictive caveat:** this is a proper prior, but at the chosen scale its
arithmetic prior mean of exp(u_t-u_1960) is infinite. Conditional on tau that mean
is exp(tau^2), whose integral is at the half-normal moment boundary. The independent
review identified this after execution began; the declared scientific source was
kept frozen, and the caveat is preserved. This does not imply an improper posterior.
Actual posterior draws and diagnostics below are finite and well behaved. A tighter
hierarchy-prior sensitivity is warranted before adopting the specification; this
test is not a calibration of the hierarchy scale.

## Execution and validity

One new PNB Bertrand history fit used four chains, 500 warmup + 500 retained draws
per chain, at most four cores, one thread per chain/numerical library, diagonal
metric, adapt_delta .999 and max_treedepth 13. The fixed 90-minute clock began
11:22:29 UTC and the controller completed 11:27:45 UTC (5m16s), including preflight
and compilation. No retries, bridges, model averaging, cluster jobs or package
refreshes ran. Snapshots were neither refitted nor replayed.

Three saved comparisons—constant variance, inverse assets, linear trend—passed
exact data/attribute, original source/kernel and native RDS/CSV seal checks. Their
100-draw policy files are copied byte-identically. The new candidate uses the same
draw-index seed with 100 conditional quality/parameter states, each reused for its
baseline and merger. Every planned cell and invalid/error field remains explicit.

Candidate maximum R-hat is 1.00655; minimum bulk ESS 832 and tail ESS 1039; zero
post-warmup divergences and depth-limit hits. Three rejected warmup proposals
reported zero Gaussian scale; logs are preserved. Constant generated quantities
have undefined R-hat and are counted separately. All 400 comparison policy states
and both Paasche indices are valid, with zero no-op effects. Candidate maximum
independent FOC residual is 1.4600e-8, below the existing 1e-6 tolerance. Cost/value
recovery, signed rates and outside normalization remain unchanged. Each index uses
its draw's solved pre-merger denominator.

Source is antitrustBayes 0.3.0.9001 at base revision
8afbd535d891a4385e69a30d6482bd15528378db plus preserved dirty local changes, not a
clean Git build. The manifest records actual source/dependency/input hashes and
paths. The existing library supplies antitrust 0.99.43.9000, coordination 0.1.0,
CmdStanR 0.9.0.9002 and CmdStan 2.39.0. Production R/Stan, specification 14, public
defaults, model-weight calculations and manuscript results were not changed.

## Results

Solved posterior mean baseline shares (%), PNB 1960:

| Bank | Observed | Constant variance | Inverse assets | Linear trend | Year hierarchy |
| --- | ---: | ---: | ---: | ---: | ---: |
| First Pennsylvania |22.96|22.02|21.66|22.70|23.10|
| PNB |21.05|19.10|18.76|19.82|19.24|
| Girard |14.23|14.11|13.90|14.17|14.21|
| Three banks combined |58.24|55.23|54.32|56.69|56.54|

Hierarchy mean errors are +.14/−1.82/−.02 percentage points, all within the stated
2–3 pp tolerance. Combined error is −1.70 pp, versus −1.55 under the linear trend.
At observed rates the means are 23.22/18.92/14.25, combined 56.39. Thus PNB's
pre-solve mean error is −2.14 pp and its solved error is −1.82 pp. The hierarchy
improves First Pennsylvania/Girard accuracy while PNB is closer under the linear
trend. This is a tradeoff, not a uniform baseline improvement.

| Metric | Constant | Inverse | Linear trend | Year hierarchy |
| --- | ---: | ---: | ---: | ---: |
| Margin point RMSE |.42031|.42153|.42039|.42123|
| Margin predictive RMSE |.61921|.62174|.61787|.62370|
| Solved-share draw RMSE (pp) |.5212|.5541|.2814|.3686|
| Slope per original rate unit |.73007|.71505|.72824|.71923|
| Market Paasche mean (%) |-14.35|-13.83|-15.24|-14.89|
| Party Paasche mean (%) |-50.33|-49.97|-52.04|-50.84|
| Mean collapsed share log likelihood |-56.78|-178.54|-150.32|-111.06|
| Mean collapsed joint log likelihood |-183.86|-306.63|-275.83|-235.19|

Draw RMSE includes posterior dispersion, rather than being RMSE of posterior mean
shares. No small-bank accuracy gate was introduced. Margin point RMSE rises .2%
and predictive RMSE .7% relative to constant variance; the slope is about 1.5%
lower. Policy means remain close to the linear-trend and original-history results.

Likelihood entries are normalized in-sample posterior averages, **not** marginal
likelihoods, Bayes factors or model probabilities. The hierarchy gains 40.64 joint
log-likelihood units over the linear trend. A 51.33-unit deficit remains versus
constant variance, so improved year flexibility has not resolved all historical
fit loss from imposing inverse assets. Updated model probabilities were not
calculated or inferred from these diagnostic scores.

Posterior variance multipliers relative to 1960, holding institution assets fixed:

| Year | Mean | 95% interval |
| --- | ---: | --- |
|1960|1.00|[1.00,1.00]|
|1961|.278|[.127,.527]|
|1962|.085|[.011,.261]|
|1963|.484|[.204,1.036]|
|1964|3.20|[1.48,6.13]|
|1965|12.92|[5.88,24.17]|

These are summaries of transformed posterior draws, not exponentials of mean
effects. Actual institution ratios additionally multiply by asset growth ratios.
The underlying year-effect SD has mean 1.506, interval [1.005,2.218]; centered u
has marginal prior SD sqrt(5/6)*tau, whereas relative-year differences have SD
sqrt(2)*tau. Only six years inform this hierarchy. The fit captures the nonmonotonic
annual pattern without divergence problems in this pilot, but is not established
as the primary model or as universally stable across conduct families.

The year-wise size/margin-adjusted residual RMS is now .95–1.06 across all six
years, versus .64–1.38 with the linear trend. This descriptive check indicates
that the hierarchy captures the annual dispersion pattern; it is conditional
on the fitted states and is not an out-of-sample validation or a resolution of
the remaining within-year precision assumptions.

## Verification and replay

Fourteen deterministic expectations passed after fixing one test-only dimname
mismatch between a labeled Helmert covariance and an unlabeled dense oracle.
The final checks have zero failures, warnings or skips. They cover contrast prior
covariance/moment matching, exact 1960 anchoring, unchanged base parameter block,
proper half-normal normalization and dense conditional quality. Stanc syntax passed.
The ordinary package suite was not repeated for this isolated research change.

Independent audits reproduce all 100 candidate states: SD relative error <=9.10e-15,
seeded dense quality replay error <=2.17e-12, direct share error <=5.56e-17. Eight
independent pre/post roots reproduce market/party Paasche ratios within 1.11e-8.
Three posterior rows per all four models match independently normalized dense
share/joint log likelihoods within 1.41e-11. Scripts, receipts and the analytic
prior-tail note are preserved in `independent_audit/`.

All 17 source guards, three reference fit hashes, three copied draw files and
18 preserved snapshot artifact seals independently pass. Reconstructed year
effects and reference offsets agree with the saved fields within 2.67e-15.

```sh
python3 doc/validation/run_pnb_year_variance_hierarchy_diagnostic.py --plan
R_LIBS_USER=../mergerBayes/.scratch/paper-library Rscript --vanilla doc/validation/year_variance_hierarchy_prototype_checks.R
R_LIBS_USER=../mergerBayes/.scratch/paper-library Rscript --vanilla doc/validation/report_pnb_year_variance_hierarchy_diagnostic.R .scratch/year_variance_hierarchy_diagnostic_20261109
```

All new configuration, adapter, runner, checks and collector files are beside this
report. Raw outputs are in `.scratch/year_variance_hierarchy_diagnostic_20261109/`:
clock/completion receipt, manifest and generated kernel/evaluator, one native
fit/CSV set, original reference pointers/copied draw plans, full comparison with
uncertainty intervals, yearly variance multipliers, residual dispersion, parameter,
likelihood, fit and validity diagnostics, audit receipts and artifact seals.
All preceding results and uncommitted changes are preserved. Nothing was pushed.
