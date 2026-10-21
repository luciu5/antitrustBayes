# Nutpie native-target and evidence validation

> Closeout update, 2026-10-20: the historical evidence and verdict below are
> unchanged. Version `0.2.0.9006` removes the global certification gate and
> permits mixed-backend evidence while retaining target/provenance checks.
> See the [current closeout results](nutpie_closeout_results.md).

## Execution status, 2026-10-20

The fixed protocol completed against Git `6f7e7b7b0449e6c751ccfdd1869ff7a6b80b1aad`,
package `0.2.0.9005`. All 15 tasks completed successfully: four native-target
checks, one analytic reference, four additional same-data fits and six
economic-model bridge calculations with three repetitions each. The runner
exited with status zero and the final collector wrote its protocol at
08:04:52 EDT on October 20.

All numerical target, evidence-precision and model-probability comparisons
passed their fixed checks. The overall certification status is **inconclusive**:
all six evidence fits have divergences, violating the prespecified zero-
divergence diagnostic requirement. No fits, priors or sampler controls were
changed to improve this outcome. The public Nutpie bridge/BMA gate remains
closed; these are bounded experimental validation results.

The scientific design and source hashes are frozen in
`.scratch/nutpie-validation/evidence/manifest.rds`. Fit seals, exact coordinates,
native model information and complete task results are retained beneath that
directory. Tracked compact results are
[`nutpie_evidence_status.csv`](nutpie_evidence_status.csv) and
[`nutpie_evidence_targets.csv`](nutpie_evidence_targets.csv).
The [design note](nutpie_evidence_design.md) specifies the coordinates,
normalization, Jacobian, uncertainty and complete-support contracts.

## Native-target comparison

The B fits use the identical installed `logit_supply.stan` source and Stan data,
model specification 11, baseline priors, Gaussian share/margin likelihoods,
collapsed market quality and active rho. Both contain four chains and 1000
retained draws per chain. The evaluator reconstructs all 4000 constrained draws
in the parameter block's declared order, verifies each constraint round trip,
and checks that CmdStan and BridgeStan use the same unconstrained names/order.

Twenty retained states and three deterministic interior states were evaluated
with both Jacobian settings and both BridgeStan proportionality settings.
The six unconstrained coordinates are `sigma_b_strat`, `mu_log_a_free.1`,
`s0_logit.1`, `rho_raw.1`, `log_sigma_share` and `log_sigma_margin`.

| Kernel | Points | Maximum normalized-target difference | Maximum constants-on/off difference | Maximum Jacobian difference | Result |
| --- | ---: | ---: | ---: | ---: | --- |
| logit_supply | 23 | 3.730349e-14 | 3.730349e-14 | 0 | Passed |
| logit_stackelberg | 23 | 2.842171e-14 | 2.842171e-14 | 0 | Passed |
| ces_supply | 23 | 1.101341e-13 | 1.101341e-13 | 0 | Passed |
| ces_stackelberg | 23 | 2.131628e-14 | 2.131628e-14 | 0 | Passed |

Differences are absolute log-density differences. The fixed point tolerance is
`1e-7 + 1e-9 * max(1, abs(targets))`. The normalized evaluation explicitly uses
`propto = FALSE, jacobian = TRUE`; the table does not infer normalization from
Nutpie's sampler log density. Native model information identifies BridgeStan
2.8.0, Stan 2.39.0, stanc3 v2.39.0 and `STAN_THREADS=true`. The R BridgeStan
client is the matching version 2.8.0; CmdStanR is 0.9.0.9002, nutpieR 1.9.1 and
bridgesampling 1.2.1.

The B pair's source SHA-256 is
`7456926489bb1ecae6cb0e1824a9fae3615cf3378113bd76853ca3b177a0d62d`;
its Stan-data SHA-256 is
`8e555ba90540b41bb2e46c237fee7c591e38966c800c28b9ec7510a67920d862`.
Source, saved source text, data, serialized JSON and native-library hashes are
checked before evaluator restoration. These numerical results establish this
fixture's native density agreement, not equivalence of every configuration.

The CS comparison uses fixed k = 2, the same six continuous coordinates and
the same Gaussian/collapsed-quality/active-rho specification. Its source
SHA-256 is
`ac4a2f4a1f1192ebf84047f15e5b625271c1b3e9bc91746f32916f93f14bd229`;
its Stan-data SHA-256 is
`806fe6e4925495f995986ed65a65f73c504a26a3f8de9e049213c28300ace23e`.
All four kernels therefore pass their 23-point checks independently, including
the bounded/positive parameter transformations and their Jacobians.

The CES supply source/data SHA-256 values are
`ca3dbf4260de1927363af9febf87fc2188a49963b39cb94fdb97333ae9a72962` /
`6d9f12daf4e59d4d04d50c06a69ae9084b35fef42d43dd1a1f0bea46e4cb40c9`.
The CES Stackelberg source/data values are
`d78577973bad96304ff73593700b267fc689c768af56b24bbb0f7bf73abaa37e` /
`97c146f5b25c347b20c42667211c5bae1a70c19f26b375950551585a074ac70c`.

## Analytic integration reference

A normalized Exp(1) times Uniform(0,1) density was integrated in unconstrained
coordinates, including both transformation Jacobians, using 4000 exact iid
draws and three normal bridge repetitions. This reference uses no Stan model
or HMC fit.

| Exact log evidence | Median estimate | Repetition SD | Maximum within-repetition MCSE | Fixed four-MCSE tolerance | Result |
| ---: | ---: | ---: | ---: | ---: | --- |
| 0 | 0.005866432 | 0.002834759 | 0.006304421 | 0.02521768 | Passed |

The three estimates were `0.005866432`, `0.006985260` and `0.001612453`.
This checks the integration API and analytic Jacobians; it does not validate
economic-model marginal likelihoods.

## Marginal likelihoods on common observations

The complete atomic B/MC/C universe uses B's unchanged observations and priors.
Four additional MC/C fits were run using four chains, 1000 warmup and 1000
retained draws per chain, two cores, seed 419771, target acceptance .95 and
depth 13. CmdStan uses its diagonal adaptation and Nutpie uses `diag`.
These fits differ from the backend benchmark's independently generated MC/C
fixtures, which cannot be compared as one observed-data universe.

Each cell below uses three normal bridge repetitions. The estimate is the
median; numerical uncertainty is the larger of repetition SD and the largest
within-repetition log-evidence MCSE. All 18 repetition estimates and their
MCSEs are retained in
[`nutpie_evidence_repetitions.csv`](nutpie_evidence_repetitions.csv).

| Model | Backend | Log evidence | Repetition SD | Maximum MCSE | Numerical uncertainty |
| --- | --- | ---: | ---: | ---: | ---: |
| B | CmdStanR | -3.390662 | 0.028564 | 0.025492 | 0.028564 |
| B | Nutpie | -3.366565 | 0.042761 | 0.033272 | 0.042761 |
| MC | CmdStanR | -5.571402 | 0.032891 | 0.024482 | 0.032891 |
| MC | Nutpie | -5.614153 | 0.022336 | 0.026759 | 0.026759 |
| C | CmdStanR | -4.748128 | 0.018625 | 0.021797 | 0.021797 |
| C | Nutpie | -4.768809 | 0.020417 | 0.020238 | 0.020417 |

The paired difference is Nutpie minus CmdStan. Its fixed tolerance is four
times the root-sum-square numerical uncertainty, and that tolerance must not
exceed .25 log units. Every model passes both the agreement and precision
checks; these numerical checks do not remove posterior-sampling bias.
Full values are in
[`nutpie_evidence_comparison.csv`](nutpie_evidence_comparison.csv).

| Model | Difference | Combined uncertainty | Difference / uncertainty | Tolerance | Numerical result |
| --- | ---: | ---: | ---: | ---: | --- |
| B | 0.024097 | 0.051424 | 0.4686 | 0.205695 | Passed |
| MC | -0.042752 | 0.042401 | -1.0083 | 0.169604 | Passed |
| C | -0.020681 | 0.029866 | -0.6925 | 0.119463 | Passed |

## Complete-universe probability calculation

The validation universe assigns prior probability 1/3 to each of B, MC and C.
All three candidates and every bridge repetition are retained on both
backends; no failed or unavailable positive-prior model was dropped. The
backend-neutral aggregation calculation gives:

| Model | CmdStan probability | Nutpie probability | Nutpie minus CmdStan |
| --- | ---: | ---: | ---: |
| B | 0.729783 | 0.739810 | 0.010027 |
| MC | 0.082435 | 0.078164 | -0.004271 |
| C | 0.187782 | 0.182026 | -0.005756 |

The maximum absolute difference is 0.010027, within the fixed .05 probability
tolerance. These values are validation calculations from the retained
uncertified fits; they do not enable public Nutpie or mixed-backend BMA.
The complete table and scientific prior are retained in
[`nutpie_evidence_universe.csv`](nutpie_evidence_universe.csv).

## Diagnostic limitation and public gate

All fits have four complete chains and zero depth hits. R-hat, bulk/tail ESS
and BFMI satisfy this protocol's fixed limits, but none satisfies its
zero-divergence requirement. Full chain values and available native diagnostic
fields are retained in
[`nutpie_evidence_diagnostics.csv`](nutpie_evidence_diagnostics.csv).

| Model | Backend | Divergences | Max R-hat | Min bulk ESS | Min tail ESS | BFMI range |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| B | CmdStanR | 9 | 1.004864 | 1150.90 | 1335.01 | 0.854–1.025 |
| B | Nutpie | 11 | 1.005001 | 757.37 | 517.91 | 0.843–0.899 |
| MC | CmdStanR | 2 | 1.004495 | 1021.77 | 794.69 | 0.836–1.012 |
| MC | Nutpie | 2 | 1.004653 | 700.11 | 416.62 | 0.788–0.949 |
| C | CmdStanR | 4 | 1.006130 | 1104.78 | 1123.07 | 0.863–1.090 |
| C | Nutpie | 3 | 1.003492 | 1161.39 | 1084.94 | 0.886–0.971 |

The missing certification is therefore posterior-quality evidence, rather
than an unavailable normalized-density API or a detected Jacobian/constant
discrepancy. The implemented evaluator has passed the native density checks
for all four kernels in these configurations; the bridge and aggregation
calculations agree numerically for the complete three-model universe.
An independently reviewed validation using adequate posterior sampling is
still needed before enabling public Nutpie evidence. No metadata flag or
validation-runner option changes the production gate.

No Stan kernels, priors, likelihoods, conduct equations or downstream production
analyses have been modified for this validation.
