# Nutpie backend validation

> Closeout update, 2026-10-20: the historical results below are preserved.
> Version `0.2.0.9006` enables normalized Nutpie and mixed-backend evidence;
> sampler diagnostics remain reported rather than acting as permission gates.
> See the [current closeout results](nutpie_closeout_results.md) for executed
> public-API checks and their remaining scope. CmdStanR remains the default.

## Executed benchmark: 2026-10-20

The fixed 34-task protocol has completed on `nutpie` at
`6f7e7b7b0449e6c751ccfdd1869ff7a6b80b1aad`, package `0.2.0.9005`.
The implementation baseline is `fff419be8f4d9096eeff7b7f39dcc5e4481ef8ab`.
All 34 posterior fits completed and passed their save/reload draw-identity
checks. The resumed invocation reused 16 completed tasks and finished the
remaining 18, preserving the original fixture manifest, priors, sampler
budgets and raw warnings. No sampling failure or target-identity mismatch
occurred.

The original results contain **1,210 parameter and 164 policy comparisons
within the predeclared four-combined-MCSE tolerance**, with none outside.
Another 870 parameter and 144 policy checks remain unscored. All 870
parameter entries are fixed-k component/leader probabilities with exactly zero
backend differences. Unscored policy entries include constants, cost-support
attrition, and the separately identified CES Stackelberg validation-runner role
error. The largest absolute scored difference is 3.380 combined MCSEs. Execution completion therefore
does not mean every policy calculation or diagnostic passed.

All 34 saved fits were independently rehashed against the current unchanged
Stan sources and their declared Stan inputs. Within all 12 fixtures, source,
data, full prior metadata, observed rows and complete target identities agree.
The [provenance table](nutpie_backend_provenance.csv) records each hash,
backend/version, model-specification version, controls, native artifact seals
and output precision. CmdStan uses its unchanged eight-significant-digit CSV
default in this frozen benchmark; Nutpie retains native double arrays. The
resulting serialization difference is much smaller than the reported MCSEs.

The complete 15-task evidence protocol has executed, including all 18 bridge
repetitions. Four native kernel checks passed all 92 states (maximum density
difference `1.101341e-13`, zero Jacobian discrepancy). B/MC/C log-evidence
differences are 0.0241, -0.0428 and -0.0207 (Nutpie minus CmdStan), all within
the fixed uncertainty bounds; maximum model-probability difference is 0.0100.
Certification remains inconclusive because all six evidence fits have
2–11 divergences and the fixed protocol requires zero. No mathematical/API
capability is missing from the normalized evaluator; adequately certified
independent posterior samples remain the empirical blocker. Public Nutpie
bridge/BMA and mixed-backend evidence support therefore remain disabled.
Both corrected CES Stackelberg policy sets have also completed from exactly
the original fits/plans after every sealed evidence worker exited. The
[separate corrected audit](#corrected-ces-stackelberg-policy-audit) retains
all original failures and adds 11 scored policy comparisons within tolerance.
Across original and corrected measurements, all 1,385 scored checks agree
within the fixed tolerances. Independent final review found no substantive
correctness or regression issue. Existing opt-in CmdStan regressions exposed a pre-existing
CSV-precision mismatch in three strict algebra fixtures. Explicit 18-digit
output at those three call sites preserves all assertions, scientific controls
and production defaults. The resulting run passed 590 assertions across all
28 blocks, including all 12 opt-in sampling blocks, with zero failures/errors/
skips and nine existing simulation warnings. Both earlier 559-pass/31-failure
runs and native CSVs remain in the durable ignored regression directories.
Independent raw-CSV checks showed that all 1,560 noncentered and 240 share
identities were within the propagated eight-digit rounding bounds.

Recommendation: **A — CmdStan remains default; Nutpie experimental**.
Diagonal Nutpie reduced divergences in each pooled fixture, but structural ESS
per wall-second did not improve uniformly. Low-rank was faster in these pooled
fits but had less favorable R-hat/divergence behavior than diagonal under the
fixed budget. These bounded fixtures establish neither universal efficiency
nor broad statistical recovery. Individual diagnostic values appear below.

Cost-support attrition is retained: C/diagonal lost one draw; CF lost 10
CmdStan, eight diagonal and 18 low-rank draws; BS/CmdStan lost two draws.
Every inspected baseline failure is `bayes_finite_nonpositive_cost`.
Affected policy summaries describe the conditional paired-valid distribution,
with intact-chain MCSEs withheld. No failed state was converted into an effect
of zero. All nine pooled policy sets retained all 4,000 realizations.

Both original CES Stackelberg policy sets have 0/4,000 valid realizations per
regime because the runner supplied `corePost` to a `StackelbergCES` game.
All 16,000 failed records retain the exact message
`corePost is only meaningful for CoreFringe games`. This is a validation
adapter error; parameter sampling and the package economics layer succeeded.
The corrected audit below fixes the runner class dispatch in `c57aebb`
without overwriting these artifacts or refitting.

Full tables: [task/timing status](nutpie_backend_status.csv),
[chain diagnostics](nutpie_backend_diagnostics.csv),
[parameter summaries](nutpie_backend_parameters.csv),
[policy summaries](nutpie_backend_policy.csv),
[paired differences/MCSEs](nutpie_backend_comparisons.csv), and
[native evidence results](nutpie_evidence_results.md).
Exact retained trajectory-step counts and BFMI/step-size values remain in the
chain diagnostic CSV. Total gradient-evaluation counts, phase-separated
Nutpie warmup/sampling timing and peak memory were not available/measured;
no ESS/gradient or memory advantage is claimed. Compile/cache-preparation
values below have different cache histories; the separately executed cold
Nutpie compilation table is preserved in the earlier handoff.

## Executed benchmark table

All fits use four chains, 1,000 warmup and 1,000 retained draws per chain, two cores, target acceptance 0.95 and maximum depth 13. ESS/sec uses the common cached-model `bayes_est()` wall duration, including data preparation and backend extraction; it excludes separately measured model preparation and downstream policies. Native total includes warmup and retained sampling across all chains. CmdStan native time includes process execution; Nutpie includes callbacks/native extraction. Their separate phase fields differ: CmdStan warmup/sampling times are individual-chain maxima; Nutpie does not expose separate phases. Exact phase fields remain in the timing CSV. Cross-backend ESS/sec uses the common fit wall duration.

| Case | Backend/adaptation | Compile/cache prep s | Fit wall s | Native total incl. warmup s | Divergences | Depth hits | Max R-hat | Min bulk ESS | Min tail ESS | mu_log_a ESS/fit-s | Demand coefficient[1] ESS/fit-s |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| MC | cmdstanr/diag_e | 72.486 | 4.300 | 1.963 | 1 | 0 | 1.0033 | 1338.4 | 1041.5 | 311.2 | 311.2 |
| MC | nutpie_diag/diag | 0.924 | 2.154 | 0.989 | 0 | 0 | 1.0032 | 913.2 | 800.1 | 424.0 | 424.0 |
| MC | nutpie_low_rank/low_rank | 0.304 | 1.823 | 0.945 | 2 | 0 | 1.0054 | 734.5 | 691.4 | 578.6 | 578.6 |
| B | cmdstanr/diag_e | 0.072 | 4.166 | 2.095 | 9 | 0 | 1.0049 | 1150.9 | 1335.0 | 348.6 | 348.6 |
| B | nutpie_diag/diag | 0.309 | 1.957 | 0.947 | 11 | 0 | 1.0050 | 757.4 | 517.9 | 447.5 | 447.5 |
| B | nutpie_low_rank/low_rank | 0.339 | 2.083 | 1.065 | 1 | 0 | 1.0025 | 780.8 | 575.7 | 374.8 | 374.8 |
| C | cmdstanr/diag_e | 0.072 | 4.504 | 2.482 | 13 | 0 | 1.0032 | 1178.5 | 805.7 | 424.1 | 424.1 |
| C | nutpie_diag/diag | 0.320 | 2.103 | 1.202 | 8 | 0 | 1.0059 | 795.9 | 733.9 | 564.4 | 564.4 |
| C | nutpie_low_rank/low_rank | 0.311 | 1.915 | 1.017 | 3 | 0 | 1.0082 | 905.3 | 745.7 | 694.8 | 694.8 |
| BF | cmdstanr/diag_e | 71.084 | 4.354 | 2.049 | 8 | 0 | 1.0044 | 992.2 | 683.7 | 236.9 | 236.9 |
| BF | nutpie_diag/diag | 0.642 | 1.905 | 0.942 | 6 | 0 | 1.0096 | 548.0 | 360.5 | 287.7 | 287.7 |
| BF | nutpie_low_rank/low_rank | 0.340 | 1.924 | 0.999 | 3 | 0 | 1.0049 | 845.7 | 722.2 | 506.9 | 506.9 |
| CF | cmdstanr/diag_e | 0.070 | 4.663 | 2.394 | 22 | 0 | 1.0070 | 785.1 | 850.2 | 467.9 | 467.9 |
| CF | nutpie_diag/diag | 0.313 | 2.014 | 1.045 | 12 | 0 | 1.0040 | 534.9 | 618.1 | 736.0 | 736.0 |
| CF | nutpie_low_rank/low_rank | 0.321 | 1.902 | 0.982 | 1 | 0 | 1.0041 | 605.5 | 546.7 | 392.4 | 392.4 |
| BS | cmdstanr/diag_e | 0.070 | 4.237 | 2.047 | 6 | 0 | 1.0074 | 666.3 | 442.0 | 184.5 | 184.5 |
| BS | nutpie_diag/diag | 0.672 | 2.860 | 1.134 | 4 | 0 | 1.0092 | 665.9 | 310.0 | 256.3 | 256.3 |
| BS | nutpie_low_rank/low_rank | 0.619 | 2.173 | 1.092 | 1 | 0 | 1.0040 | 824.6 | 587.2 | 394.5 | 394.5 |
| CS | cmdstanr/diag_e | 0.209 | 5.222 | 2.794 | 19 | 0 | 1.0025 | 1138.2 | 1520.2 | 433.6 | 433.6 |
| CS | nutpie_diag/diag | 0.331 | 2.094 | 1.179 | 7 | 0 | 1.0087 | 469.3 | 455.5 | 411.1 | 411.1 |
| CS | nutpie_low_rank/low_rank | 0.335 | 2.056 | 0.986 | 0 | 0 | 1.0059 | 846.3 | 538.0 | 494.8 | 494.8 |
| C_pooled | cmdstanr/dense_e | 0.356 | 45.556 | 42.750 | 23 | 0 | 1.0091 | 554.3 | 679.8 | 34.2 | 43.9 |
| C_pooled | nutpie_diag/diag | 0.436 | 40.487 | 39.392 | 1 | 0 | 1.0083 | 641.0 | 592.7 | 15.8 | 84.3 |
| C_pooled | nutpie_low_rank/low_rank | 0.282 | 35.466 | 34.449 | 3 | 0 | 1.0112 | 524.7 | 633.6 | 22.0 | 73.5 |
| CF_pooled | cmdstanr/dense_e | 0.064 | 38.532 | 35.486 | 26 | 0 | 1.0120 | 430.9 | 238.5 | 14.1 | 84.2 |
| CF_pooled | nutpie_diag/diag | 0.308 | 25.216 | 24.007 | 1 | 0 | 1.0081 | 366.1 | 291.1 | 14.5 | 166.8 |
| CF_pooled | nutpie_low_rank/low_rank | 0.321 | 25.092 | 23.852 | 13 | 0 | 1.0144 | 288.7 | 299.2 | 14.3 | 153.9 |
| CS_pooled | cmdstanr/dense_e | 0.067 | 37.323 | 34.032 | 37 | 0 | 1.0060 | 724.4 | 461.2 | 28.3 | 105.9 |
| CS_pooled | nutpie_diag/diag | 0.302 | 27.824 | 26.573 | 2 | 0 | 1.0096 | 412.7 | 472.4 | 14.8 | 137.4 |
| CS_pooled | nutpie_low_rank/low_rank | 0.314 | 24.637 | 23.381 | 9 | 0 | 1.0126 | 361.5 | 678.4 | 24.0 | 133.2 |
| CES_B | cmdstanr/diag_e | 57.257 | 6.418 | 2.985 | 23 | 0 | 1.0038 | 1048.2 | 854.7 | 163.4 | 163.4 |
| CES_B | nutpie_diag/diag | 0.421 | 1.899 | 1.083 | 15 | 0 | 1.0109 | 513.8 | 255.0 | 291.0 | 291.0 |
| CES_CS | cmdstanr/diag_e | 69.501 | 4.998 | 2.756 | 35 | 0 | 1.0048 | 860.9 | 888.3 | 213.3 | 213.3 |
| CES_CS | nutpie_diag/diag | 0.445 | 2.199 | 1.316 | 18 | 0 | 1.0162 | 459.3 | 285.1 | 241.8 | 241.8 |

## Policy tables

Price effects are percent changes. Summary distributions use draws valid for both policies, with failures retained in the original results. An asterisk marks attrition or the identified CES runner failure, with intact-chain MCSE unavailable. The complete CSV records medians, intervals, threshold probabilities and all MCSEs.

| Case | Backend | Paired valid/planned | Merger mean (SD), 95% interval | Tariff mean (SD), 95% interval |
| --- | --- | ---: | --- | --- |
| MC | cmdstanr | 4000/4000 | 0.000 (0.000), [0.000, 0.000] | 6.760 (0.164), [6.440, 7.073] |
| MC | nutpie_diag | 4000/4000 | 0.000 (0.000), [0.000, 0.000] | 6.772 (0.154), [6.461, 7.114] |
| MC | nutpie_low_rank | 4000/4000 | 0.000 (0.000), [0.000, 0.000] | 6.766 (0.139), [6.485, 7.057] |
| B | cmdstanr | 4000/4000 | 7.991 (0.987), [6.097, 10.195] | 6.634 (0.204), [6.248, 7.070] |
| B | nutpie_diag | 4000/4000 | 8.002 (1.017), [6.065, 10.245] | 6.643 (0.211), [6.257, 7.124] |
| B | nutpie_low_rank | 4000/4000 | 8.006 (1.028), [6.062, 10.219] | 6.637 (0.219), [6.230, 7.110] |
| C | cmdstanr | 4000/4000 | 6.659 (0.613), [5.531, 8.017] | 2.868 (0.283), [2.320, 3.472] |
| C | nutpie_diag | 3999/4000* | 6.659 (0.634), [5.496, 8.103] | 2.870 (0.291), [2.297, 3.477] |
| C | nutpie_low_rank | 4000/4000 | 6.665 (0.600), [5.655, 8.032] | 2.874 (0.274), [2.330, 3.453] |
| BF | cmdstanr | 4000/4000 | 7.252 (1.186), [5.208, 9.827] | 7.521 (0.238), [7.059, 8.039] |
| BF | nutpie_diag | 4000/4000 | 7.219 (1.145), [5.173, 9.828] | 7.518 (0.247), [7.052, 8.105] |
| BF | nutpie_low_rank | 4000/4000 | 7.313 (1.166), [5.397, 9.989] | 7.497 (0.222), [7.050, 7.901] |
| CF | cmdstanr | 3990/4000* | 7.181 (0.801), [5.886, 9.043] | 5.074 (0.163), [4.707, 5.378] |
| CF | nutpie_diag | 3992/4000* | 7.133 (0.838), [5.633, 8.978] | 5.091 (0.186), [4.727, 5.453] |
| CF | nutpie_low_rank | 3982/4000* | 7.179 (0.794), [5.847, 9.218] | 5.069 (0.165), [4.709, 5.372] |
| BS | cmdstanr | 3998/4000* | 8.811 (1.144), [6.512, 11.214] | 6.397 (0.267), [5.978, 6.962] |
| BS | nutpie_diag | 4000/4000 | 8.814 (1.091), [6.672, 11.196] | 6.401 (0.227), [5.994, 6.930] |
| BS | nutpie_low_rank | 4000/4000 | 8.813 (1.080), [6.682, 11.091] | 6.396 (0.209), [5.997, 6.863] |
| CS | cmdstanr | 4000/4000 | 5.250 (0.384), [4.520, 6.096] | 4.874 (0.324), [4.210, 5.536] |
| CS | nutpie_diag | 4000/4000 | 5.260 (0.451), [4.439, 6.272] | 4.878 (0.375), [4.132, 5.583] |
| CS | nutpie_low_rank | 4000/4000 | 5.273 (0.387), [4.528, 6.148] | 4.851 (0.350), [4.075, 5.484] |
| C_pooled | cmdstanr | 4000/4000 | 6.674 (0.091), [6.484, 6.837] | 2.895 (0.076), [2.740, 3.040] |
| C_pooled | nutpie_diag | 4000/4000 | 6.676 (0.088), [6.488, 6.842] | 2.898 (0.073), [2.750, 3.036] |
| C_pooled | nutpie_low_rank | 4000/4000 | 6.677 (0.088), [6.494, 6.839] | 2.897 (0.072), [2.754, 3.035] |
| CF_pooled | cmdstanr | 4000/4000 | 7.203 (0.070), [7.066, 7.342] | 5.123 (0.023), [5.077, 5.168] |
| CF_pooled | nutpie_diag | 4000/4000 | 7.205 (0.070), [7.068, 7.348] | 5.123 (0.023), [5.077, 5.166] |
| CF_pooled | nutpie_low_rank | 4000/4000 | 7.206 (0.070), [7.068, 7.341] | 5.122 (0.023), [5.075, 5.168] |
| CS_pooled | cmdstanr | 4000/4000 | 5.178 (0.055), [5.068, 5.287] | 4.971 (0.082), [4.818, 5.141] |
| CS_pooled | nutpie_diag | 4000/4000 | 5.178 (0.053), [5.075, 5.285] | 4.971 (0.081), [4.812, 5.134] |
| CS_pooled | nutpie_low_rank | 4000/4000 | 5.177 (0.054), [5.072, 5.283] | 4.971 (0.082), [4.814, 5.137] |
| CES_B | cmdstanr | 4000/4000 | 5.194 (0.559), [4.049, 6.426] | 9.526 (0.086), [9.338, 9.696] |
| CES_B | nutpie_diag | 4000/4000 | 5.183 (0.557), [3.936, 6.330] | 9.526 (0.088), [9.336, 9.695] |
| CES_CS | cmdstanr | 0/4000* | — (—), [—, —] | — (—), [—, —] |
| CES_CS | nutpie_diag | 0/4000* | — (—), [—, —] | — (—), [—, —] |

## Paired policy differences

Cells show CmdStan minus Nutpie in percentage points, followed by the difference divided by combined MCSE in parentheses. A dash for the scaled value means the comparison remains unscored. Full quantile, median, threshold and SD comparisons are in the comparison CSV.

| Case | Nutpie adaptation | Merger mean difference (MCSEs) | Merger SD difference (MCSEs) | Tariff mean difference (MCSEs) | Tariff SD difference (MCSEs) |
| --- | --- | ---: | ---: | ---: | ---: |
| MC | nutpie_diag | 0.000 (—) | 0.000 (—) | -0.012 (-1.40) | 0.011 (0.49) |
| MC | nutpie_low_rank | 0.000 (—) | 0.000 (—) | -0.005 (-0.67) | 0.025 (1.08) |
| B | nutpie_diag | -0.011 (-0.33) | -0.029 (-0.97) | -0.009 (-0.88) | -0.007 (-0.39) |
| B | nutpie_low_rank | -0.015 (-0.48) | -0.040 (-1.36) | -0.003 (-0.31) | -0.015 (-0.89) |
| C | nutpie_diag | -0.000 (—) | -0.021 (—) | -0.002 (—) | -0.007 (—) |
| C | nutpie_low_rank | -0.006 (-0.30) | 0.013 (0.43) | -0.006 (-0.56) | 0.010 (1.04) |
| BF | nutpie_diag | 0.033 (0.88) | 0.041 (1.19) | 0.003 (0.23) | -0.010 (-0.49) |
| BF | nutpie_low_rank | -0.061 (-1.63) | 0.020 (0.48) | 0.025 (2.09) | 0.016 (0.66) |
| CF | nutpie_diag | 0.048 (—) | -0.038 (—) | -0.018 (—) | -0.023 (—) |
| CF | nutpie_low_rank | 0.002 (—) | 0.006 (—) | 0.005 (—) | -0.002 (—) |
| BS | nutpie_diag | -0.003 (—) | 0.052 (—) | -0.004 (—) | 0.039 (—) |
| BS | nutpie_low_rank | -0.002 (—) | 0.063 (—) | 0.000 (—) | 0.057 (—) |
| CS | nutpie_diag | -0.010 (-0.47) | -0.067 (-1.94) | -0.004 (-0.22) | -0.051 (-1.50) |
| CS | nutpie_low_rank | -0.023 (-1.56) | -0.003 (-0.17) | 0.023 (1.68) | -0.026 (-1.46) |
| C_pooled | nutpie_diag | -0.001 (-0.51) | 0.003 (1.12) | -0.003 (-1.43) | 0.003 (1.36) |
| C_pooled | nutpie_low_rank | -0.002 (-0.97) | 0.003 (1.20) | -0.002 (-0.94) | 0.004 (1.82) |
| CF_pooled | nutpie_diag | -0.002 (-1.21) | -0.000 (-0.22) | 0.000 (0.46) | 0.001 (1.07) |
| CF_pooled | nutpie_low_rank | -0.003 (-1.70) | -0.001 (-0.46) | 0.001 (1.19) | -0.000 (-0.40) |
| CS_pooled | nutpie_diag | 0.001 (0.42) | 0.002 (1.48) | -0.000 (-0.11) | 0.001 (0.64) |
| CS_pooled | nutpie_low_rank | 0.001 (0.94) | 0.001 (0.62) | 0.000 (0.04) | -0.000 (-0.10) |
| CES_B | nutpie_diag | 0.010 (0.59) | 0.002 (0.10) | 0.001 (0.25) | -0.002 (-0.52) |
| CES_CS | nutpie_diag | — (—) | — (—) | — (—) | — (—) |

## Key hierarchy scale posteriors

Each cell gives mean and 95% interval. These are posterior scale parameters, not sampler mass-matrix settings.

| Pooled case | Backend | Alpha-market SD | Outside-share SD | Market-intercept SD |
| --- | --- | --- | --- | --- |
| C_pooled | cmdstanr | 0.223 [0.062, 0.615] | 0.156 [0.019, 0.513] | 0.383 [0.154, 0.838] |
| C_pooled | nutpie_diag | 0.242 [0.068, 0.664] | 0.155 [0.017, 0.503] | 0.379 [0.153, 0.838] |
| C_pooled | nutpie_low_rank | 0.230 [0.067, 0.585] | 0.158 [0.018, 0.529] | 0.384 [0.152, 0.849] |
| CF_pooled | cmdstanr | 0.247 [0.083, 0.648] | 0.098 [0.008, 0.403] | 0.316 [0.122, 0.732] |
| CF_pooled | nutpie_diag | 0.248 [0.085, 0.670] | 0.100 [0.008, 0.387] | 0.320 [0.123, 0.723] |
| CF_pooled | nutpie_low_rank | 0.249 [0.079, 0.611] | 0.095 [0.008, 0.387] | 0.326 [0.121, 0.792] |
| CS_pooled | cmdstanr | 0.240 [0.085, 0.566] | 0.145 [0.013, 0.488] | 0.287 [0.094, 0.705] |
| CS_pooled | nutpie_diag | 0.245 [0.083, 0.603] | 0.150 [0.017, 0.500] | 0.286 [0.093, 0.680] |
| CS_pooled | nutpie_low_rank | 0.241 [0.080, 0.576] | 0.144 [0.016, 0.462] | 0.291 [0.093, 0.729] |

## Comparison counts

- parameter / unavailable_mcse: 870
- parameter / within_mcse_tolerance: 1210
- policy_conditional_paired_valid / unavailable_mcse: 144
- policy_conditional_paired_valid / within_mcse_tolerance: 164
- policy_unconditional_availability / different_availability_count: 10
- policy_unconditional_availability / same_availability_count: 34

## Corrected CES Stackelberg policy audit

The one-line runner fix recognizes `StackelbergCES` alongside
`StackelbergLogit` and forwards inherited `leadersPost`. The actual economic
reconstruction and trade solvers are unchanged. The focused no-HMC runner
suite passed **61 assertions in 10 blocks**, including 11 checks of actual
CES merger/tariff solves and inherited leaders, with zero failures/errors/
warnings/skips.

The policy-only replay used the original posterior states, seed 491871,
complete plans, indices, latent quality, k states and fingerprints. Exact
R/platform/package/source guards passed. All 60 original manifest, protocol,
status, fit, argument, task and result byte hashes remained unchanged. The
only changed source seal is the reviewed runner condition; all 32 production
R source files are identical to the original run. Both backends retained
4,000/4,000 valid draws in each regime, with zero baseline or policy failures
and physical-cost identity error exactly zero. No compilation or sampling
was needed.

| Backend | Policy | Mean | SD | Median | 2.5% | 97.5% | P(effect > 5%) |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| cmdstanr | merger | 4.6008 | 0.5298 | 4.5834 | 3.6465 | 5.7090 | 0.1610 |
| cmdstanr | tariff25 | 8.8077 | 0.2588 | 8.8462 | 8.2122 | 9.2121 | 1.0000 |
| nutpie_diag | merger | 4.6143 | 0.6622 | 4.5743 | 3.6049 | 5.8827 | 0.1583 |
| nutpie_diag | tariff25 | 8.8063 | 0.2591 | 8.8456 | 8.1872 | 9.2107 | 1.0000 |

All 11 scored policy comparisons are within four combined MCSEs; three
constant-indicator MCSEs remain unavailable. Maximum absolute scaled
difference is 0.972. Means agree closely, but the merger SD difference is
0.1324 with combined MCSE about 0.136. Nutpie's larger SD, maximum R-hat 1.0162
and 18 divergences leave tail precision uncertain. These results support
numerical compatibility under this budget, rather than tight distributional
certification.

| Policy | CmdStan minus Nutpie mean | Combined MCSE | Scaled mean difference | CmdStan minus Nutpie SD | Combined SD MCSE |
| --- | ---: | ---: | ---: | ---: | ---: |
| merger | -0.013501 | 0.042067 | -0.321 | -0.132433 | 0.136305 |
| tariff25 | 0.001409 | 0.007869 | 0.179 | -0.000333 | 0.005941 |

Tracked evidence: [original/corrected summaries](nutpie_ces_policy_replay.csv),
[comparisons](nutpie_ces_policy_replay_comparisons.csv),
[plan/fit identities](nutpie_ces_policy_replay_identity.csv), and
[source seals](nutpie_ces_policy_replay_sources.csv). Raw audit and executable
replay are retained in `.scratch/nutpie-validation/ces-role-audit/`, with
completed output `replay-25a3df78f/`. The first replay guard correctly rejected
Git's CRLF conversion before any policy work; materializing the exact reviewed
LF snapshot restored the required source bytes. Original inputs and seals were
preserved throughout.

The completed original backend/evidence manifests intentionally seal the old
runner source. Running their resumable entrypoints after this fix correctly
rejects that source mismatch. Preserve those roots; use a fresh
`AB_NUTPIE_OUTPUT` directory for a new full run. The separately sealed audit
is the supported reuse route for these completed CES fits.

The ordinary development-source suite previously passed 4,905 expectations,
zero failures/errors, 12 opt-in HMC skips and 21 existing warnings. The
subsequent opt-in run above executes all 12 skipped blocks and passes 590
assertions. The full package build/check completed with zero errors/warnings
and one pre-existing `tail()` namespace NOTE. Production defaults, the four
maintained Stan kernels and scientific specifications are unchanged; no
`mergerBayes` production analysis, DGP or sampler setting was modified.

## Earlier interrupted-run checkpoint






## Partial execution update: 2026-10-20

The user's local run reached the shell tool's 3,600-second timeout. Completed
artifacts survived: **15 of 34 tasks** completed, covering MC, B, C, BF and CF
with CmdStanR, Nutpie diagonal and Nutpie low-rank sampling. Each fit retained
four chains of 1,000 draws after 1,000 warmup iterations. An additional BS
CmdStan fit and its preliminary validation checkpoint were saved; its policy
work was incomplete at this inspection. The other 19 tasks lack complete
collected results. No native evidence-run artifacts exist yet.

For the completed pairs, there are **400 parameter and 65 policy comparisons
within the predeclared four-combined-MCSE tolerance**, with no scored comparison
outside tolerance or recorded target-identity mismatch. Another 120 parameter
and 75 policy comparisons lack MCSEs and remain unscored; constant quantities
and policy attrition must not be counted as passed precision checks.

Sampler divergences, summed over the four retained chains, are:

| Fixture | CmdStanR | Nutpie diagonal | Nutpie low-rank |
| --- | ---: | ---: | ---: |
| MC | 1 | 0 | 2 |
| B | 9 | 11 | 1 |
| C | 13 | 8 | 3 |
| BF | 8 | 6 | 3 |
| CF | 22 | 12 | 1 |

No completed chain hit its maximum depth. Maximum recorded R-hat across these
fits is 1.00963; these diagnostics do not remove the need to investigate
divergences. C/Nutpie diagonal lost one of 4,000 policy realizations. CF lost
10 CmdStanR, eight Nutpie-diagonal and 18 Nutpie-low-rank realizations. The
runner retains those failures and withholds intact-chain policy MCSEs when
attrition occurs. Six policy-availability comparisons therefore differ.

As an illustrative completed control, B merger mean price effects are 7.9911%,
8.0018% and 8.0061%, and tariff means are 6.6339%, 6.6430% and 6.6371%, in
CmdStanR/diagonal/low-rank order. These are bounded fixture results, not broad
posterior or policy certification. Full numbers remain in the linked CSVs.

The implementation remains suitable for integration as an **optional
experimental backend**, with CmdStanR the default and Nutpie evidence/BMA
disabled. Unfinished Stackelberg, hierarchical and CES comparisons, divergences,
policy attrition and missing native evidence validation prevent a production
recommendation or a demonstrated difficult-model performance claim.

The following 2026-10-19 handoff is retained as historical preparation and
verification evidence. Its statements that HMC had not run describe that date;
the partial execution update above owns the current sampling status.

## Original execution handoff: 2026-10-19

Status: 2026-10-19. Implementation and deterministic checks are available;
posterior, policy, difficult-model, and native evidence validation await the
user's local execution. **No HMC fit or paired sampler result is claimed here.**
The user explicitly elected to run the long jobs locally after recovery from a
computer reset. This report distinguishes recovered preparation, executed
compilation, and planned statistical checks.

The checkout started at `fff419be8f4d9096eeff7b7f39dcc5e4481ef8ab`, package
`0.2.0.9005`, with an uncommitted Nutpie implementation. The recovered manifest
`20261019T200221-0400` contained 12 fixtures and 34 backend tasks, no started
execution clock, and no completed or failed fits. Reusing a saved fit was
therefore impossible. The original fixture objects and their fingerprint remain
in `.scratch/nutpie-validation/manifest.rds`; preparation revisions are retained
under `preparation-history/`.

## Local execution

From the repository root:

```sh
bash doc/validation/run_nutpie_validation.sh --prepare
bash doc/validation/run_nutpie_validation.sh --run
bash doc/validation/run_nutpie_validation.sh --evidence
```

`--prepare` constructs/collects metadata without compilation or inference.
`--compile` compiles the four unchanged kernels without inference. `--run`
executes posterior and policy tasks, and `--evidence` performs the separately
gated target-density/evidence validation. The latter can create the explicitly
specified additional same-data comparison fits; see the
[evidence design](nutpie_evidence_design.md) and
[evidence runner](nutpie_evidence_validation.R).

The wrapper selects `.scratch/nutpie-library` and the existing compatible
`../mergerBayes/.scratch/paper-library` for **read-only dependency loading**.
Override `AB_NUTPIE_R_LIBS` with a colon-separated compatible library list on
another machine. No downstream production file or installed package is changed.
Each invocation writes a dated `invocation-*.log` under
`.scratch/nutpie-validation`; each task keeps arguments, native logs, status,
fit, and validation checkpoints in its own directory. `flock` prevents two
wrapper invocations from running against the same output root simultaneously.

The default budget is 14,400 seconds **per invocation**. For example:

```sh
AB_NUTPIE_MAX_SECONDS=28800 bash doc/validation/run_nutpie_validation.sh --run
```

Rerun the same command to continue tasks not yet started. Completed tasks are
reused; unfinished postfit work resumes from the saved fit without resampling.
Failed inference is preserved and requires explicit
`AB_NUTPIE_RETRY_FAILED=true` to retry; previous task/status/log artifacts are
archived. R source, validation-runner source, dependency locations/versions/
recorded revisions, and R/platform identity are sealed once execution starts.
A changed seal requires a fresh `AB_NUTPIE_OUTPUT` directory. This prevents
mixing stale economic-policy calculations with changed implementation code.
The source/data/model/prior identity check is additional to that runtime seal.

## Recovered native environment and compilation

The default R library was incompatible. The selected stack is R 4.6.1,
antitrust 0.99.42, coordination 0.0.0.9006, CmdStanR 0.9.0.9002,
CmdStan 2.39.0, posterior 1.7.0, trade 0.8.9.9002, nutpieR 1.9.1 at
`8e88a975d5e0b8ed3c76de84d22adb33de601a1e`, and bridgestan R 2.8.0.
The dependency CSV records the currently resolved optional prerequisites.

The original native compile attempt failed when nutpieR tried to download the
235 MB BridgeStan source release. Both sandboxed and escalated retries reported
`failed to lookup address information: Temporary failure in name resolution`.
This was a download/environment failure, not a demonstrated sampler defect.

Recovery used the approved read-only GitHub connector to retrieve the 11
official build/source/license files from BridgeStan 2.8.0 commit
`d2c499a7a85c9069f5b73b05d5ee2e9b8661a3e3`. Every file's exact Git blob
SHA-1 was verified. Sources and the verification manifest are retained in
`.scratch/nutpie-toolchain/bridgestan-source/bridgestan-2.8.0/`.
The official Makefile and all C/C++ files remain unchanged.

BridgeStan [documents reuse of an existing Stan installation](https://roualdes.us/bridgestan/v2.8.0/getting-started.html#using-pre-existing-stan-installations).
Read-only symlinks reuse the existing Stan 2.39.0 tree and `stanc` binary.
An approved symlink at `~/.bridgestan/bridgestan-2.8.0` provides the native
nutpieR discovery location. A make dry run confirmed that the existing
TBB/Sundials libraries were already built and compilation writes only workspace
artifacts. No `HOME` override, package-library write, upstream patch, or Python
translation was used. On an ordinary connected installation, nutpieR's normal
one-time BridgeStan download supplies these sources.

Compile-only measurements are retained in
`.scratch/nutpie-validation/compile-recovery.log` and `precompiled.rds`.
These are cold native model preparation times, including small adapter/cache
overhead, from a single execution; they are not sampler performance results.

| Kernel | Source SHA-256 | Native preparation seconds |
| --- | --- | ---: |
| logit_supply | `7456926489bb1ecae6cb0e1824a9fae3615cf3378113bd76853ca3b177a0d62d` | 100.524 |
| logit_stackelberg | `ac4a2f4a1f1192ebf84047f15e5b625271c1b3e9bc91746f32916f93f14bd229` | 103.493 |
| ces_supply | `ca3dbf4260de1927363af9febf87fc2188a49963b39cb94fdb97333ae9a72962` | 86.008 |
| ces_stackelberg | `d78577973bad96304ff73593700b267fc689c768af56b24bbb0f7bf73abaa37e` | 99.031 |

## Fixed statistical design

The [runner](nutpie_backend_validation.R) prepares seven logit single-market
conduct cases MC, B, C, BF, CF, BS, CS; three hierarchical cases C_pooled,
CF_pooled, CS_pooled; and two CES cases CES_B and CES_CS. All fringe/Stackelberg
components use fixed `k=2`, preserving the existing conditional-conduct model.
The hierarchical cases have three markets and two periods, alpha-market
pooling, outside-share pooling, market intercepts, active rho, and collapsed
market quality. Fixture equilibrium checks are in
[the fixture table](nutpie_backend_fixtures.csv).

Logit fixtures have CmdStanR, Nutpie diagonal, and Nutpie low-rank variants;
CES fixtures have CmdStanR and Nutpie diagonal variants. Each task uses four
chains, 1,000 warmup and 1,000 retained draws per chain, at most two concurrent
chains, seed 419771, target acceptance 0.95, maximum depth 13, and initialization
0. CmdStan uses its existing diagonal default for single-market cases and
`dense_e` for pooled cases; Nutpie variants explicitly select `diag` or
`low_rank`. These adaptation strategies are recorded separately, with no
claim that dense and low-rank matrices are identical.

This is a bounded paired fixture design, not a production grid. B serves as a
control, pooled C/CF/CS exercise the difficult hierarchical families, and
fixed-k CF supplies a single-market conditional component. A pooled historical
bank benchmark is omitted: the package's `pnb` data supply a single market of
39 banks and do not establish a reproducible bank panel. No pseudo-panel was
created and no mergerBayes production data/analysis was imported.

Every pair must have identical Stan source/data hashes, prior metadata,
model-specification identity, observed data, and required posterior variables.
The comparison records means, medians, 2.5%/97.5% quantiles and SDs; paired policy
summaries additionally record probabilities of price effects exceeding 0 and 5
percentage points. For a summary statistic `s`, the predeclared numerical
comparison is

```
abs(s_Nutpie - s_CmdStan) <= max(4 * sqrt(MCSE_Nutpie^2 + MCSE_CmdStan^2), 1e-10).
```

MCSEs use intact chain axes and the corresponding posterior-package mean, SD,
or quantile estimator. Four combined MCSEs are a per-statistic discrepancy
screen, not a multiplicity-adjusted certification or substitute for adequate
mixing. Identity mismatches, unavailable MCSEs, and outside-tolerance results
remain explicit rows. Signed CmdStan-minus-Nutpie differences and differences
scaled by combined MCSE are saved. Sampling does not require identical draws.

## Policy and diagnostics contracts

All posterior draws enter a fixed posterior-policy draw plan. The merger joins
the two largest firms. The tariff is 25% on the other firms, represented by
the established tax-inclusive wedge `0.25/1.25`. Both regimes share the same
latent-quality state, draw index, and plan fingerprint. The runner verifies
that the physical marginal-cost identity is preserved through each solve.
Definitions and posterior targets do not vary by sampler.

Mean, SD, median, interval, and threshold-probability comparisons are produced
for the quantity-weighted price percentage change in market m1/period t1.
All realization and policy failures are retained, including per-backend valid
fractions. When every planned draw survives, chain order is restored before
computing MCSE. With attrition, the paired-valid distribution is explicitly
conditional/descriptive and MCSE comparison is withheld; failure counts are
reported separately. Failed states are never silently reweighted away.

Diagnostics use the common backend adapter, including authoritative native
maximum-depth flags when available, chain divergence counts, BFMI, step sizes,
R-hat, and bulk/tail ESS. BFMI matches the installed CmdStanR definition
`(sum(diff(energy)^2)/length(energy))/var(energy)` per retained chain.
Native diagnostics remain in checkpoints. Missing metrics stay unavailable.
CmdStan reports warmup/sampling times separately; the current nutpieR route
exposes a combined native duration, so warmup-only and sampling-only Nutpie
efficiency columns remain unavailable. Per-parameter ESS divided by full
estimation time offers a common measured denominator; it includes backend
overheads and is labeled accordingly. Cold compilation is separate.

## Current quantitative results and decision

| Evidence | Executed result |
| --- | --- |
| Fixture fresh-equilibrium checks | All 12 prepared cases pass; numeric errors in fixture CSV |
| Native unchanged-kernel compilation | All four kernels compile; measured times above; no inference |
| Paired posterior agreement | Not run; all 34 task statuses remain explicit |
| Difficult-model ESS/speed comparison | Not run |
| Merger/tariff equivalence | Not run |
| Native target-density/evidence comparison | Not run; analytic bridge reference is separate evidence |

## Implementation and regression verification

The final development-source ordinary suite completed with **4,905 passing
expectations**, zero failures/errors, 12 explicit sampling skips and 21 solver/
simulation warnings across 352 test blocks. No HMC was enabled. The warnings
retain their original messages in
`.scratch/nutpie-validation/resumed-final-suite.log`; they include equilibrium
convergence notices, unchanged ownership and unused-share notices. They were
not suppressed or converted into successful sampling evidence.

A separate CmdStan-only library run, with both `nutpieR` and `bridgestan`
absent, completed 4,881 expectations without failures/errors, 13 skips and
the same 21 warnings. Its additional skip was the optional Nutpie compilation
test. Saved draw extraction, summaries, diagnostics and policy use remained
available without either optional package. The later focused adapter checks
passed 196 expectations with optional packages present. Evidence/conditional-k
focused checks passed 319 expectations; runner checks passed 50.

The new identical-state regression exercises all seven conduct families,
observed/latent-quality anchoring, role inheritance, structural reconstruction
and merger simulation using differently ordered native columns. It also tests
conditional-k backend forwarding and preservation of failed support states.
These deterministic checks establish the adapter contract, not agreement of
independent posterior samples or tariff distributions.

The built package's installed test suite passed. The initial check of a build
that deliberately omitted vignette outputs reported two missing-`inst/doc`
warnings. A full vignette build resolved both; the final packaging check
(reusing the already passed tests) finished with **zero errors, zero warnings,
one pre-existing NOTE** for unqualified `tail()` in unchanged
`R/cmdstan_backend.R`. Logs are retained under
`.scratch/nutpie-validation/resumed-package-check*.log`. Manuals were omitted;
vignette examples did not run HMC.

Fresh independent review found and resolved stale resume provenance, a missing
CmdStan comparison metadata reference, and absent paired SD comparisons. The
reviewer rechecked the fixes. Implementation commits are `c511251` (backend
seam), `81c5951` (Nutpie integration) and `8cdf11d` (guarded evidence); the
validation/handoff commit follows them on `nutpie`.

The generated status, parameter, diagnostic, policy, and comparison CSVs are
the durable numeric outputs. Empty result tables explicitly state `not_run`;
they must not be read as passing comparisons. After local execution, update
this report from those outputs, retain actual warnings and incomplete cases,
and assess MCSE adequacy before recommending a sampler.

The current recommendation is **A: CmdStan remains the default; Nutpie is
experimental**. Successful compilation and deterministic adapter tests do not
establish posterior or policy equivalence, a speed advantage, or evidence/BMA
comparability. Native evidence stays fail-closed pending its separate
equivalence gate. No mergerBayes production analysis, DGP, sampler setting,
maintained Stan kernel, scientific specification, prior, or policy definition
was changed for this validation work.
