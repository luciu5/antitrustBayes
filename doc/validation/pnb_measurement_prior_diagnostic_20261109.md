# PNB jointly estimated measurement precision: bounded diagnostic

Completed 2026-11-09. The requested measurement-prior test passed its bounded
sampling and policy checks. Two new Bertrand fits were run sequentially, with
four chains, 500 warmup and 500 retained draws per chain; no divergences or
depth-limit hits occurred. Large-bank mean share errors are at most 0.30
percentage points in the snapshot and 0.59 in history. This is a research
prototype, not an adopted package default or a manuscript result.

## What changed

The preceding [fixed-.10 experiment](pnb_demand_shock_diagnostic_20261109.md)
used persistent bank quality plus a period demand shock and independent
log-share measurement error. This experiment replaces its fixed measurement SD
with one jointly estimated scalar:

\[\log m\sim N(\log 0.10,0.5^2).\]

The median is .10; its prior mean is .1133 and central 95% interval
[.03753,.26644]. The existing quality prior, period-shock SD prior (the former
share-residual prior), shock-margin correlation prior and total/fraction
Jacobian remain as in the fixed-.10 prototype. No inverse-size weights,
linear trend or hierarchical year effects enter these two candidate fits.

In each within-year relative-share contrast the statistical model is
persistent quality + period demand + independent measurement. If period-demand
SD is omega, margin-error SD is sigma_m and their correlation is rho_d, the
observable residual covariance is

\[\begin{pmatrix}\omega^2+m^2 & \rho_d\omega\sigma_m\\
\rho_d\omega\sigma_m & \sigma_m^2\end{pmatrix}.\]

The primitive sampled coordinate is log(m); its normalized Normal prior needs
no additional exponential Jacobian. Both quality and demand shocks stay
analytically integrated in the HMC likelihood. Joint conditional recovery draws
persistent quality first and then period demand from the same relative-share
contrasts and observed margins, using **that posterior row's m**. Only quality
and recovered demand enter the economic state. Measurement noise is excluded.
The identical realized state enters both equilibria; there is no requirement
for exact share reconstruction and no observed-share fallback. Original
observed-share cost/value inversion conventions are preserved.

The snapshot has 39 banks in 1960. History retains 204 distinct bank-year rows
in 1960–1965 and 41 persistent bank-quality states; deposits are not pooled
across years. Target prices, deposits, margins, merger parties and signed-rate
index rules are unchanged.

## Sampling and measurement posterior

| Information | Chains | Retained total | Max R-hat | Min bulk ESS | Min tail ESS | Divergences | Depth hits | Min E-BFMI |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| snapshot | 4 | 2000 | 1.00612 | 834 | 521 | 0 | 0 | 0.865 |
| history | 4 | 2000 | 1.00384 | 941 | 957 | 0 | 0 | 0.927 |

The new log-measurement coordinate itself has adequate mixing. Constant derived
zeros have undefined R-hat and are counted separately in the raw diagnostics;
they are not treated as failed varying parameters. Logs preserve three snapshot
and four history zero-scale rejected proposals at the first reported warmup
iteration; no such messages occur during retained sampling. These are not
post-warmup divergences. No rescue fits or sampler retuning were used.

| Distribution | Mean SD | Median SD | 95% interval | Pr(m<.10) | Pr(m>.20) |
| --- | --- | --- | --- | --- | --- |
| prior | 0.113 | 0.100 | [0.038, 0.266] | 0.500 | 0.083 |
| snapshot_estimated | 0.116 | 0.102 | [0.040, 0.274] | 0.483 | 0.093 |
| history_estimated | 0.081 | 0.079 | [0.033, 0.142] | 0.736 | 0.000 |

Snapshot barely updates the prior. History shifts toward a smaller measurement
SD, but its measurement and demand-shock log SDs have posterior correlation
-0.702. The fraction of draws with |rho_d|>.90 is .032 in snapshot and zero in
history. These diagnostics do **not** establish separate measurement-error
identification: the likelihood depends on omega^2+m^2 and rho_d*omega, so with
free rho_d a continuum of allocations has the same likelihood. Updating is
conditional on both variance priors and correlation support. The prior on m is
therefore an explicit modeling assumption, not an empirically established
calibration.

## Requested share comparison

Posterior percentages, mean [2.5%,97.5%], from 100 fixed policy states per cell.
These are in-sample conditional-state reconstructions, not held-out predictions.
The observed-price column uses the recovered state at observed rates; the
solved column uses that draw's pre-merger equilibrium. Error is the solved
posterior mean minus the observation, in percentage points.

| Information | Bank | Observed % | At observed prices | Solved baseline | Mean error pp |
| --- | --- | --- | --- | --- | --- |
| snapshot | First Pennsylvania | 22.96 | 23.00 [18.03, 27.55] | 22.97 [18.90, 26.54] | 0.01 |
| snapshot | PNB | 21.05 | 21.02 [17.34, 25.41] | 21.00 [17.87, 24.41] | -0.05 |
| snapshot | Girard | 14.23 | 13.90 [10.72, 16.64] | 13.93 [11.10, 16.29] | -0.30 |
| snapshot | Three banks combined | 58.24 | 57.91 [52.75, 61.68] | 57.90 [53.27, 61.17] | -0.34 |
| history | First Pennsylvania | 22.96 | 22.78 [20.13, 25.00] | 22.79 [20.63, 24.64] | -0.17 |
| history | PNB | 21.05 | 20.36 [18.33, 22.75] | 20.46 [18.77, 22.43] | -0.59 |
| history | Girard | 14.23 | 14.20 [12.78, 15.92] | 14.19 [12.93, 15.67] | -0.04 |
| history | Three banks combined | 58.24 | 57.35 [55.24, 59.72] | 57.44 [55.56, 59.54] | -0.80 |

Both candidates meet the requested 2–3 percentage-point large-bank tolerance
by their posterior means. No small-bank accuracy gate was imposed.
The three-bank combined error is -0.34 pp in snapshot and -0.80 pp in history.
Broad intervals alone were not the acceptance criterion.

The explicit benchmarks remain available:

| Information | Arm | First Pennsylvania % | PNB % | Girard % | Combined % |
| --- | --- | --- | --- | --- | --- |
| snapshot | reference | 18.67 | 17.37 | 9.12 | 45.17 |
| snapshot | fixed | 23.06 | 20.96 | 14.10 | 58.12 |
| snapshot | observed | 22.96 | 21.05 | 14.23 | 58.24 |
| snapshot | estimated | 22.97 | 21.00 | 13.93 | 57.90 |
| history | reference | 22.02 | 19.10 | 14.11 | 55.23 |
| history | fixed | 22.64 | 20.25 | 14.32 | 57.21 |
| history | observed | 22.96 | 21.05 | 14.23 | 58.24 |
| history | estimated | 22.79 | 20.46 | 14.19 | 57.44 |

`reference` is the original constant-variance model with quality-only policy
state, using stronger compatible four-chain reference fits. `fixed` is the
previous demand-shock prototype with measurement SD .10. `observed` reuses the
original reference fit and retains the full residual (zero independent
measurement) as an explicit calibration benchmark; it reproduces observed
inside shares and baseline rates. `estimated` is the new model in this report.
The earlier variance experiments, including their already-corrected
counterfactuals, are preserved in the
[experiment record](share_residual_experiment_record_20261109.md).

## Margins and merger policies

| Cell | Margin point RMSE | Margin predictive RMSE | All-bank solved share RMSE pp |
| --- | --- | --- | --- |
| snapshot_reference | 0.42543 | 0.59838 | 2.987 |
| snapshot_fixed | 0.42528 | 0.60447 | 0.429 |
| snapshot_observed | 0.42543 | 0.59838 | 0.000 |
| snapshot_estimated | 0.42531 | 0.60700 | 0.586 |
| history_reference | 0.42031 | 0.61921 | 0.521 |
| history_fixed | 0.41966 | 0.61797 | 0.389 |
| history_observed | 0.42031 | 0.61921 | 0.000 |
| history_estimated | 0.42000 | 0.61876 | 0.322 |

Margin errors use observed target-year log inverse margins and the same fixed
policy rows; predictive RMSE includes the fitted margin residual variance.
There is no material deterioration in these bounded comparisons. The share
RMSE averages all bank/draw errors, rather than just errors of posterior means.

| Information | Arm | Whole-market Paasche % | Merging-party Paasche % | Market valid | Party valid |
| --- | --- | --- | --- | --- | --- |
| snapshot | reference | -8.18 [-21.63, -0.10] | -37.60 [-111.95, -3.10] | 100/100 | 100/100 |
| snapshot | fixed | -16.66 [-22.20, -11.89] | -55.70 [-72.93, -39.73] | 100/100 | 100/100 |
| snapshot | observed | -16.29 [-20.24, -10.92] | -54.26 [-67.28, -36.25] | 100/100 | 100/100 |
| snapshot | estimated | -16.15 [-21.63, -10.61] | -54.41 [-74.63, -35.32] | 100/100 | 100/100 |
| history | reference | -14.35 [-17.01, -11.18] | -50.33 [-60.21, -40.07] | 100/100 | 100/100 |
| history | fixed | -16.26 [-20.54, -12.14] | -55.02 [-69.21, -41.23] | 100/100 | 100/100 |
| history | observed | -16.82 [-19.54, -13.07] | -56.02 [-64.96, -43.57] | 100/100 | 100/100 |
| history | estimated | -16.15 [-19.89, -12.25] | -54.39 [-66.72, -41.85] | 100/100 | 100/100 |

The input-market signs indicate reductions in deposit rates under the stated
Bertrand model. Every Paasche basket uses post-policy quantities and **that
draw's solved pre-merger rates** in its denominator. Policy intervals and
means above are single-model results; no updated model probabilities or
model-averaged estimates were computed or inferred from earlier fits.
All 800 requested benchmark/candidate policy draws are structurally and
index-valid. All 200 new-model policies are valid, with no draws discarded or
probabilities renormalized.

## Verification, provenance and limits

The extended deterministic checks passed 25 expectations for conditional
Gaussian recovery at m=0,.04,.20, the exact variance-allocation likelihood ridge
and a proper lognormal measure; a further three reporting checks passed for the
upper-tail probability. Stan syntax validation passed. Independent dense
Gaussian integration checks 64 native posterior rows across all eight cells;
maximum normalized likelihood discrepancy is 3.1e-12. State-implied shares
agree within 1.2e-16. The same 100-draw plans replay exactly for both new fits,
and no-op effects are zero. Maximum saved independent FOC residual is below
1.5e-8 across all cells. A separate reviewer solved ten mergers per new model
with an independent Broyden solver; index agreement was within 1.83e-6
percentage points, using solved-pre denominators.
The [separate policy audit](measurement_prior_independent_policy_check.R)
preserves that executable check and its ten posterior-row indices per model.

The original first test invocation failed because testthat changed the working
directory; source paths were repaired. A prior-quantile test literal was also
corrected before fitting. One report-only prior-tail complement typo was found
after fitting: the original prior Pr(m>.20) column incorrectly showed .91717
rather than .08283. Initial summaries, executed source and seals remain
preserved. `prior_posterior_summary_corrected.csv` is authoritative; the exact
one-line repair and old/new hashes are recorded in
`.scratch/measurement_prior_diagnostic_20261109/reporting_correction/receipt.json`.
The posterior, likelihood, states and policies are unaffected. Source checking
accepts only that exact reporting repair against the original archived hash;
no kernel, prior/data or inference seal is relaxed or rewritten.
The reviewer also corrected a fraction/percentage unit mismatch in an initial
audit assertion before the successful independent solver check; model outputs
were unaffected.

Execution outputs are in `.scratch/measurement_prior_diagnostic_20261109`.
The original controller ran for about 3.5 minutes, including preparation,
compilation and both fits, under its absolute 90-minute cap. One numerical
thread per library and one thread per chain were enforced. Four sealed
reference/fixed fits and their byte-identical policy plans were reused; only
two new HMC fits were made. The original snapshot reference now has four chains
and 500 retained draws per chain, rather than the short two-chain fit in the
first demand-shock experiment. The reference policy seed is 10902; fixed-.10
policies retain their original 10702 plan; new-model seed is 10902. Comparisons
across different fitted models do not represent identical parameter draws.

The manifest records current Git revision plus dirty source hashes, native
kernel/data SHA256 values, R package/library paths and CmdStan provenance.
Inference is from a separately generated research kernel, not an unmodified
specification-14 BayesFit. It is ineligible for production evidence/BMA.
The compatible existing library is
`/mnt/d/Projects/mergerBayes/.scratch/paper-library`, with development
antitrustBayes 0.3.0.9001, antitrust 0.99.43.9000, coordination 0.1.0,
cmdstanr 0.9.0.9002 and CmdStan 2.39.0. No new dependencies were installed.

Production R/Stan models, inference specification, defaults and the manuscript
were not changed by this diagnostic. The ordinary package suite was not rerun
for standalone research scripts; no full seven-model rerun, bridge sampling,
model averaging, cluster job, commit or push occurred. No candidate was
silently skipped; both new fits completed. Sampling failure did not trigger the
inverse-size/year-effects fallback. Its finite-tail hierarchy-prior revision
and a production API remain separate work.

The next decision is comparative: estimated measurement keeps good baseline
shares without a divergence penalty in this bounded case, while inverse-size
variance plus year effects is a distinct successful route already propagated
through its pilot counterfactuals. This pilot alone does not establish a
universal default or held-out/known-truth performance.

## Reproduction

From the antitrustBayes root, plan mode performs no fits or writes:

```sh
python3 doc/validation/run_pnb_measurement_prior_diagnostic.py --plan
```

A separately authorized repeat needs a distinct output directory:

```sh
python3 doc/validation/run_pnb_measurement_prior_diagnostic.py --execute --output=.scratch/measurement_prior_diagnostic_repeat
```

Collection uses saved draws only and never fills missing fits:

```sh
Rscript --vanilla doc/validation/report_pnb_measurement_prior_diagnostic.R .scratch/measurement_prior_diagnostic_20261109
Rscript --vanilla doc/validation/measurement_prior_independent_policy_check.R .scratch/measurement_prior_diagnostic_20261109
```

Deterministic tests use the existing compatible library:

```sh
Rscript --vanilla -e '.libPaths(c("/mnt/d/Projects/mergerBayes/.scratch/paper-library", .libPaths())); for (p in c("doc/validation/measurement_prior_prototype_checks.R", "doc/validation/measurement_prior_reporting_checks.R")) testthat::test_file(p, reporter="summary", stop_on_failure=TRUE)'
```
