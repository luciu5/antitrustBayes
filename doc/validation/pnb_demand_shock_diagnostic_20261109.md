# Bounded PNB period-demand-shock diagnostic (2026-11-09 UTC)

This is a separate Bertrand-only prototype, not the old residual-prior pilot and not a production-model or manuscript update. The naming/alias repair is independent of baseline shrinkage.

## Alias verification without estimation

A sealed, saved PNB Bertrand `total_fraction` fit passed ten policy draws. The selective economic quality SD equals the actual `sigma_b_strat_economic` values; the equivalent transformation from rounded historical coordinates agrees to relative error below `6.70e-8` (eight-significant-digit CSV precision). Reusing one draw plan reproduces both equilibria exactly, unchanged ownership gives exactly zero price effects, and the maximum independent input-Bertrand FOC residual is `6.676058e-9`, below the existing `1e-6` tolerance. The saved fit checksum is unchanged. This is no evidence that the alias defect caused share shortfalls.

## Declared prototype

Persistent quality retains its existing prior and centering. A new centered period demand shock has the former share-residual SD prior, `log(omega) ~ Normal(log(.5), .5)`. Measurement error is independent, constant and fixed at log-SD `.10`. There are no share-dependent measurement weights. The existing rho prior applies to the period-shock/margin correlation; the likelihood uses `S=sqrt(omega^2+.10^2)` and `rho_eff=rho_shock*omega/S`. The original quality/shock scale transformation and Jacobian are retained. The production Stan kernels, priors, public API and model specification 14 were not edited.

Quality is recovered conditional on shares and observed margins; each period shock is then drawn conditional on that same quality and parameter realization using the full relative-share contrasts. The target-period quality plus shock stays fixed through baseline and merger equilibria. Measurement noise is excluded from the structural demand state. Original observed-share-inverted value primitives are preserved, with signed input rates allowed. Paasche denominators use each draw's solved pre-merger rates, never observed rates.

## Execution and convergence

Four comparison cells completed: one compatible saved snapshot reference (2 chains, 250 retained per chain), plus **three new fits**, each with four parallel chains, 500 warmup and 500 retained draws per chain. Fits ran sequentially, with one thread per chain/numerical library and 100 policy draws per cell. There was no bridge sampling, model averaging, cluster job or package/library refresh. PNB history uses the existing 1960–1965 data and pooling conventions: 204 bank-year rows; snapshot has 39 banks. The fixed target remains PNB 1960.

The absolute cap ran from `02:44:57` to `04:14:57 UTC`, including preparation and compilation. The final sampling/policy worker finished at `03:15:23 UTC` (30m26s after start). The initial candidate compile failed **before sampling** because CmdStan's precompiled-header directory was read-only in the sandbox. The approved retry rebuilt the cache and used the same source, data and deadline. Original failed logs, controller manifest and first-attempt statuses are preserved. Successful-retry collection is separately named; no old inference or seals were rewritten.

- Snapshot reference: R-hat 1.010523, minimum bulk ESS 283; **convergence-inconclusive**, shorter reused fit, not a fresh four-chain reference.
- Snapshot candidate: R-hat 1.005174, minimum bulk ESS 956.
- History reference: R-hat 1.005571, minimum bulk ESS 1062.
- History candidate: R-hat 1.003282, minimum bulk ESS 1065.
- All cells: zero divergences and depth-limit hits. History logs retain sporadic rejected warmup proposals with scale-zero informational messages. All 400 baseline/policy draws and all 400 whole-market and party indices are valid; no draws were dropped. Maximum production pre/post FOC residual is `1.484e-8`; all no-op effects are zero.

## One comparison table

Shares and Paasche effects are percentages; brackets are 95% posterior intervals based on the same 100 selected draw states. Margin errors are on the existing log inverse-markdown response scale. `Margin RMSE` compares structural predictions with observations; `Margin PPC` includes marginal predictive residual variance. The combined row sums the three bank shares; its margin metric is a descriptive deposit-weighted RMS across those three banks, **not a likelihood weight or a small-bank accuracy gate**. The full-precision table is saved as `comparison_completed.csv`.

| Information / arm | Bank | Observed share | At observed rates | Solved baseline | Margin RMSE / PPC | Market Paasche | Merging-party Paasche |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Snapshot / reference | First Pennsylvania | 22.96 | 16.91 [1.76, 42.71] | 17.27 [2.20, 37.80] | 0.160 / 0.456 | -9.18 [-28.62, -0.11] | -43.66 [-143.69, -3.17] |
| Snapshot / reference | PNB | 21.05 | 19.07 [1.75, 37.92] | 18.74 [2.13, 33.76] | 0.108 / 0.440 | -9.18 [-28.62, -0.11] | -43.66 [-143.69, -3.17] |
| Snapshot / reference | Girard | 14.23 | 9.03 [1.24, 20.98] | 9.29 [1.39, 20.17] | 0.424 / 0.602 | -9.18 [-28.62, -0.11] | -43.66 [-143.69, -3.17] |
| Snapshot / reference | Combined | 58.24 | 45.02 [6.41, 68.88] | 45.30 [7.53, 66.25] | 0.241 / 0.490 | -9.18 [-28.62, -0.11] | -43.66 [-143.69, -3.17] |
| Snapshot / candidate | First Pennsylvania | 22.96 | 23.11 [19.65, 26.52] | 23.06 [20.23, 25.73] | 0.153 / 0.456 | -16.66 [-22.20, -11.89] | -55.70 [-72.93, -39.73] |
| Snapshot / candidate | PNB | 21.05 | 20.94 [18.36, 24.00] | 20.96 [18.87, 23.53] | 0.105 / 0.442 | -16.66 [-22.20, -11.89] | -55.70 [-72.93, -39.73] |
| Snapshot / candidate | Girard | 14.23 | 14.08 [11.69, 16.68] | 14.10 [11.98, 16.26] | 0.413 / 0.596 | -16.66 [-22.20, -11.89] | -55.70 [-72.93, -39.73] |
| Snapshot / candidate | Combined | 58.24 | 58.13 [55.22, 61.09] | 58.12 [55.56, 60.74] | 0.234 / 0.489 | -16.66 [-22.20, -11.89] | -55.70 [-72.93, -39.73] |
| History / reference | First Pennsylvania | 22.96 | 21.79 [19.32, 24.26] | 21.91 [19.89, 23.92] | 0.133 / 0.477 | -14.46 [-17.40, -11.62] | -50.74 [-60.27, -41.38] |
| History / reference | PNB | 21.05 | 18.45 [16.73, 20.23] | 18.81 [17.39, 20.29] | 0.076 / 0.464 | -14.46 [-17.40, -11.62] | -50.74 [-60.27, -41.38] |
| History / reference | Girard | 14.23 | 14.41 [12.85, 16.18] | 14.30 [12.96, 15.80] | 0.405 / 0.611 | -14.46 [-17.40, -11.62] | -50.74 [-60.27, -41.38] |
| History / reference | Combined | 58.24 | 54.66 [52.74, 57.04] | 55.03 [53.36, 57.12] | 0.222 / 0.509 | -14.46 [-17.40, -11.62] | -50.74 [-60.27, -41.38] |
| History / candidate | First Pennsylvania | 22.96 | 22.60 [20.60, 25.23] | 22.64 [21.02, 24.69] | 0.134 / 0.473 | -16.26 [-20.54, -12.14] | -55.02 [-69.21, -41.23] |
| History / candidate | PNB | 21.05 | 20.12 [18.09, 23.49] | 20.25 [18.61, 22.96] | 0.077 / 0.460 | -16.26 [-20.54, -12.14] | -55.02 [-69.21, -41.23] |
| History / candidate | Girard | 14.23 | 14.37 [12.27, 16.41] | 14.32 [12.49, 16.08] | 0.407 / 0.609 | -16.26 [-20.54, -12.14] | -55.02 [-69.21, -41.23] |
| History / candidate | Combined | 58.24 | 57.09 [54.36, 59.75] | 57.21 [54.88, 59.53] | 0.223 / 0.505 | -16.26 [-20.54, -12.14] | -55.02 [-69.21, -41.23] |

## Interpretation and limits

Under the declared `.10` measurement precision, all candidate large-bank mean errors fall comfortably inside the requested 2–3 percentage-point tolerance. Snapshot solved errors are +.10/−.09/−.13 pp; history errors are −.32/−.80/+.09 pp for First Pennsylvania/PNB/Girard. Combined errors shrink from −12.94 to −.12 pp (snapshot) and −3.21 to −1.03 pp (history). This assessment uses mean errors, not merely interval coverage.

Across all 39 target banks, unweighted margin point-RMSE changes only `.42381 -> .42528` (snapshot) and `.41986 -> .41966` (history). Marginal predictive RMSE changes `.60145 -> .60447` and `.62118 -> .61797`. There is no material margin deterioration in this bounded comparison. Structural slope means are nearly unchanged: `.73963 -> .73607` in snapshot and `.73007 -> .72946` in history (per original rate unit). There are no conduct probabilities because conduct is fixed to Bertrand.

The candidate policy intervals are finite, all index draws are valid, and independent equilibrium calculations agree. History party uncertainty is somewhat wider with period shocks; it should not be described as unchanged. The snapshot reference has a borderline convergence diagnostic and only 500 total retained draws, so its paired comparison remains exploratory/inconclusive for adoption. More fundamentally, a snapshot does not identify persistent quality separately from a transitory shock. Fixing measurement precision makes this a transparent diagnostic assumption, not an empirically established calibration. These results support keeping the prototype for further validation, not changing package defaults, conducting an unauthorized full rerun, or replacing manuscript numbers.

## Verification and reproducibility

The 112 focused package extraction/input-policy assertions passed with zero failures, warnings or skips. Nine deterministic prototype assertions passed, including an independent dense joint-Gaussian conditional-mean identity with partial margin observations and Stan one-dimensional-array inputs. Independent extended algebra checks compare full joint quality/shock covariance with a dense Gaussian oracle. The saved-results auditor recomputed all four cells' share summaries from 400 saved states and solved 20 independent pre/post FOC paths; maximum share-summary discrepancy is `5.69e-14`, independent FOC residual `6.31e-12`, and Paasche discrepancy `7.07e-8` in ratio units. Source/data hashes and reference seal verify.

Runner/config/module: `run_pnb_demand_shock_diagnostic.py`, `run_pnb_demand_shock_diagnostic.R`, `pnb_demand_shock_config.R`, and `demand_shock_prototype.R` in this directory. The Python entrypoint enforces an outer fixed-deadline timeout around the entire controller. The original executed driver is retained in the output as `executed_driver.R`; after completion, only CLI failure-exit handling and mutually exclusive flags were tightened in the maintained driver. Numerical source/model/data were unchanged.

```sh
# No inference or writes:
python3 doc/validation/run_pnb_demand_shock_diagnostic.py --plan
# A separately authorized repetition must use a distinct output:
python3 doc/validation/run_pnb_demand_shock_diagnostic.py --execute --output=.scratch/pnb_demand_shock_followup
# No-inference collection of this completed run:
Rscript doc/validation/report_pnb_demand_shock_diagnostic.R .scratch/demand_shock_diagnostic_20261109
# Deterministic checks with a compatible existing library:
Rscript doc/validation/demand_shock_prototype_checks.R
```

Durable raw artifacts are under `.scratch/demand_shock_diagnostic_20261109/`: `alias_check.*`, `clock.json`, `manifest.rds`, all three new `fit.rds`/CSV sets, four `policy_draws.rds` and validity CSVs, `comparison_completed.csv`, `execution_status_completed.csv`, `margin_fit.csv`, `parameter_summary.csv`, `completed_artifact_seals.csv`, `completed_run.json`, failed/retry logs, and `independent_audit/`. The original `execution_status.csv` describes the **first, compilation-failed controller attempt**; use the separately collected completed status for final results. These artifacts are diagnostic, not primary banking results or production-fit migrations.
