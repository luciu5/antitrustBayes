# Native profile v2 (baseline compiled kernels)

The profile used only the installed `antitrustBayes` objects from the baseline
library and `rstan::sampling(..., algorithm = "Fixed_param")`,
`log_prob`, `grad_log_prob`, and standalone `rstan::gqs`. It did not compile or
run HMC. The v1 invalid Gaussian-share/Student-t-margin sampled probes remain
in `native_profile_v1_failed_hybrid.csv`; v2 uses valid Gaussian/Gaussian
fixtures for both sampled and collapsed scaling.

## Native results

All 16 rows in `native_profile_v2_results.csv` had finite setup, target,
gradient, and GQ results. The archive rows used the four valid saved fits:
logit simultaneous, logit BF, logit BS, and CES simultaneous. Each standalone
GQ call used the five-column raw parameter schema selected by
`.bayes_bridge_raw_parameter_names`; the eight compared rows matched every
common output column exactly (`max_abs = max_rel = 0`, all finite). The common
output counts were 149, 279, 275, and 149 in that archive order.
`raw_schema_v2.csv` records five selected raw flat columns for all archives and
collapsed scaling, and 104/1,004 for sampled scaling at G=100/1,000; transformed
parameters, generated quantities, and `lp__` are excluded.

| fixture | N/G | setup s | log_prob s | grad s | GQ rows / s |
|---|---:|---:|---:|---:|---:|
| archive logit simultaneous | 24/24 | 0.027 | 0.00005 | 0.00010 | 8 / 0.039 |
| archive logit BF | 24/24 | 0.207 | 0.00010 | 0.00020 | 8 / 0.107 |
| archive logit BS | 24/24 | 0.029 | 0.00010 | 0.00015 | 8 / 0.025 |
| archive CES simultaneous | 24/24 | 0.034 | 0.00010 | 0.00010 | 8 / 0.021 |

Synthetic simultaneous scaling was valid Gaussian share plus Gaussian margin
for both quality modes. At N=100 and 1,000, G=N; at N=10,000, G=1,000.
The target/gradient/GQ times (seconds) were:

| demand / quality | N=100 | N=1,000 | N=10,000, G=1,000 |
|---|---:|---:|---:|
| logit collapsed | 0.00010 / 0.00030 / 0.026 | 0.00110 / 0.00180 / 0.060 | 0.02667 / 0.02000 / 0.550 |
| logit sampled | 0.00020 / 0.00040 / 0.024 | 0.00150 / 0.00270 / 0.087 | 0.01067 / 0.01600 / 0.378 |
| CES collapsed | 0.00010 / 0.00020 / 0.022 | 0.00120 / 0.00250 / 0.066 | 0.01067 / 0.02300 / 0.439 |
| CES sampled | 0.00010 / 0.00020 / 0.024 | 0.00160 / 0.00300 / 0.082 | 0.00800 / 0.02267 / 0.595 |

Entries are `log_prob / grad_log_prob / one-row GQ`; setup times and all
statuses are in the CSV. The complete v2 process took 44.28 s wall time and
reached 636,592 kB (621.7 MiB) maximum RSS under `/usr/bin/time -v`.

The GQ amortization fit measures repeated standalone calls at 1, 8, 50, and
200 rows (three repeats per count). The intercept includes per-call setup;
the slope is the fitted incremental row cost:

| fixture | intercept s | incremental s/row | R² | source |
|---|---:|---:|---:|---|
| logit BF, N=24/G=24 | 0.02962 | 0.00000838 | 0.004 | archived posterior rows |
| CES simultaneous, N=24/G=24 | 0.01669 | 0.00001122 | 0.576 | archived posterior rows |
| logit sampled, N=1000/G=1000 | 0.06110 | 0.00073971 | 0.971 | one Fixed_param row repeated |

The small archive slopes are near the timing noise floor; the large sampled
fixture shows a clear row-dependent cost. See
`gq_amortization_v2_observations.csv` for every call and
`gq_amortization_v2_fit.csv` for the fits.

## Helmert reference check

The independent dense matrix and the exact O(G) suffix recurrence agreed to
max absolute output errors of 0, 2.78e-17, 1.39e-16, 5.55e-16, and 3.33e-15
at G=1, 3, 10, 100, and 1,000. Directional output errors were 0, 0,
4.44e-16, 6.66e-16, and 4.66e-15. Centered sums stayed at or below
3.33e-15 for the recurrence. Central-difference checks (h=1e-6) were at most
1.10e-9 for dense and 1.82e-10 for recurrence. Dense/reference timing ratios
were 1x, 1x, 1x, 15.3x, and 161.7x over the same G sequence; G=1, 3, and 10
are below the timer resolution. No G=10,000 matrix was allocated. Its Q-only
bound is 799,920,000 bytes (762.86 MiB), recorded in
`helmert_storage_estimate_v2.csv` and `storage_estimates_v2.csv`.

## Thread readiness

For all four installed kernels, the source scan found no `reduce_sum`,
`map_rect`, OpenMP primitive, or `STAN_THREADS` use. The RStan
`threads_per_chain` option round-tripped successfully (default before probe:
1), and `RcppParallel::defaultNumThreads()` was 12. Local
`StanHeaders:::CxxFlags()` reported the existing `-D_REENTRANT -DSTAN_THREADS`
flags. The API and flags are present, but the kernels have no threaded
likelihood primitive, so this probe does not establish runtime parallelism.
