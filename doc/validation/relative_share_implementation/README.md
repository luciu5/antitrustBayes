# Relative-share implementation validation

Completed bounded implementation validation, 2026-10-23. The relative-share
likelihood, posterior reconstruction, normalized evidence and conditional-k
paths pass the checks below. The original noiseless CES M4T4 stress attempt
remains incomplete; its partial artifacts are preserved. This is not a claim
of broad repeated-sample recovery or that all sampler diagnostics are good.

The implementation starts from `01e5a2250e5bc1205dcf20774d12010db4e8c508`,
package 0.2.0.9008/specification 13. The candidate is version 0.2.0.9009/specification 14.
Its default is `share_basis = "relative"`; explicit `"absolute"` retains the
previous free-intercept likelihood. CmdStanR remains the sampler default.
The validated implementation is local commit
`f23e39617ef730370ff5056eaa683ed0d1dabf91`; this report and the research runners
are recorded separately. No remote publication is claimed here.
[Kernel hashes](kernel_identity.csv), [R/helper identities](relative_share_source_identity.csv)
and [checked package files](checked_package_source.csv) identify the executed
source. All R, Stan, help and test files match the checked snapshot byte for
byte. Its 191 staged files include one later README example edit, which is
excluded from the package build; all inputs included in the archive match.

## Public behavior and scope

```r
bayes_est(..., share_basis = "relative")  # default, Gaussian shares/margins
bayes_est(..., share_basis = "absolute")  # prior free-intercept specification
```

Logit and CES use J-1 orthonormal log-share coordinates per observed market-year.
The projected joint Gaussian includes active rho and sampled or analytically
collapsed persistent centered quality. Every realized quality/conduct draw has
its own derived market-year intercept calibrated to its structural outside
share. No sampled intercept or intercept prior contributes to the relative
target. Supply equations, conduct roles, prior defaults, outside-share pooling
and counterfactual solvers are unchanged. Relative shares currently require
Gaussian shares and margins; robust likelihoods and sampled/collapsed
intercept coordinates require explicit absolute shares.

Within either basis, mixed CmdStan/Nutpie evidence and mixed-backend fixed-k
components are supported. Absolute and relative evidence cannot be combined:
they are densities for different observed coordinates. This is a scientific
compatibility rule, not a sampler diagnostic gate. No raw A-versus-B Bayes
factor is reported. Saved fits from earlier specifications require refitting.

## Executed deterministic and package checks

- Ordinary source suite: 443 test cases, 7,435 passing expectations, 0 failures, 0 errors;
  22 warnings and 15 opt-in skips retained.
- Full opt-in source suite: 443 test cases, 7,619 passing expectations,
  0 failures, 0 errors, 30 warnings and no skips. These totals exclude warning
  and skip records; counting all records gives 7,472 and 7,649, respectively.
- Final native targets: 5,903 checks pass, including 319 baseline configurations,
  heterogeneous market scales and total/fraction coordinates, exact 17-digit
  thin-noise invariants, CES price-unit identities, and 20 old/new absolute Logit
  target evaluations. The latter are bit-identical.
- All four kernels compile through CmdStan and Nutpie/BridgeStan.
- Package check from a clean source snapshot:0 errors, 0 warnings, 1 pre-existing
  NOTE about unqualified `tail()` in `R/cmdstan_backend.R`.
- Independent reviews cover projected Gaussian mathematics, conditional
  quality reconstruction, rotated cell arithmetic, conditional-k family
  contracts, and the saved-artifact collector. A fresh final review checked
  this report against raw artifacts, including unavailable-slot denominators,
  passing-expectation counts, partial CSV rows and cache restoration, and
  approved the closeout without a remaining substantive finding.
- Native M4T4 posterior/policy checks: 546 Logit and 546 finite-noise CES checks
  pass, across both backends, all 16 cells and four posterior states. They check
  same-quality outside-share calibration, reconstructed baseline demand,
  physical costs and supply FOCs, plus merger/tariff construction in selected
  cells. These are structural identities, not a recovery experiment.

See [native summary](native_summary.csv), [invariants](native_invariants.csv)
and [absolute-target preservation](absolute_target_preservation.csv),
[test totals](verification_summary.csv), [retained test warnings](verification_warnings.csv),
[package check](package_check.csv), and the full
[Logit](panel_logit_policy_checks.csv)/[CES](panel_ces_policy_checks.csv) panel checks.

## Frozen paired protocol

[The runner](../relative_share_backend_validation.R) uses four chains,
500 warmup and 500 retained draws per chain, target acceptance 0.99 and maximum
depth 13, with at most four sampling cores. Both samplers use diagonal/default
adaptation. These are validation controls; package defaults are unchanged.
The total native execution budget is four hours, excluding the overnight
interruption. Source/data/prior identities are checked within every pair.

The seven-model Logit exercise is a **fixed-k validation universe**: MC, B, C,
BF,CF,BS,CS, with rich models fixed at k=2 and explicit equal model weights.
The separate CF conditional-k family has original support 1:3 and prior
(0.2,0.5,0.3). It is not an extra candidate in that seven-model universe.
The runner also covers CES C/CS and Logit/CES M4T4. These small fixtures are
not the paper production experiment or a repeated-sample recovery study.

Logit fixtures retain small prescribed observation errors. The original CES
fixtures are noiseless structural stress cases. A separate finite-error CES
panel experiment is described below. Sampler warnings are reported;
finite evidence is not suppressed by divergence, R-hat or ESS thresholds.

[The collector](../relative_share_validation_collect.R) retains every requested
bridge attempt and all missing/failed candidates. Its evidence tolerance uses
four times combined within-bridge MCSE and repetition SD, without dividing
repetitions that share posterior draws by sqrt(number of repeats). Posterior
summaries use four combined MCSE. Policy results retain paired-validity mass;
conditional mean/probability/SD MCSE use ratio influences in the original
chain order. Quantile uncertainty remains unavailable under attrition.

Runtime: R 4.6.1, CmdStan 2.39.0, CmdStanR 0.9.0.9002, nutpieR 1.9.1,
BridgeStan 2.8.0, posterior 1.7.0, antitrust 0.99.42, coordination 0.0.0.9006,
and trade 0.8.9.9002. The selected compatible libraries were read only;
development antitrustBayes source was loaded explicitly. See the manifest
provenance and [source identities](relative_share_input_identity.csv).

## Paired posterior and evidence results

The original protocol has **24 complete fits out of 26 planned**, forming
12 complete pairs. Of its 78 planned bridge slots, 72 yield finite evidence
from the 24 completed fits; six remain unavailable with the missing pair.
All 72 actually executed bridge calculations are finite. The missing pair
is the original noiseless CES panel; the
[original status table](relative_share_status.csv) deliberately remains
`incomplete_saved_validation`. Nothing is normalized over survivors of a
positive-prior model universe.

Among the complete original pairs, 575 posterior-summary comparisons lie
within four combined MCSE and 210 match deterministic saved fields. All 12
paired evidence comparisons lie within four combined bridge uncertainties;
the largest absolute difference is 0.09718 log units (BF), or 1.078 combined
standard errors. These uncertainties include both within-repetition estimates
and variation across repetitions. See
[posterior comparisons](relative_share_posterior_comparisons.csv),
[evidence attempts](relative_share_evidence_attempts.csv) and
[evidence comparisons](relative_share_evidence_comparisons.csv).

| Pair | Nutpie minus CmdStan log evidence | Combined uncertainty | Difference / uncertainty |
| --- | ---: | ---: | ---: |
| Logit B | 0.08964 | 0.09271 | 0.967 |
| Logit CF, k=2 | 0.02798 | 0.05135 | 0.545 |
| Logit CS, k=2 | -0.00298 | 0.05104 | -0.058 |
| CES C | -0.03925 | 0.09088 | -0.432 |
| CES CS, k=2 | 0.01501 | 0.06653 | 0.226 |
| Logit C, M4T4 | -0.04728 | 1.03522 | -0.046 |
| CES C, M4T4, finite-error experiment | 0.21679 | 0.24135 | 0.898 |

The seven-model fixed-k universe is complete under each backend and the mixed
assignment. With explicit equal priors 1/7, its largest CmdStan/Nutpie posterior
model-probability difference is **0.004484** (0.448 percentage points).
Cournot receives 0.81201/0.81288. These fixture probabilities are not PNB or
paper estimates. A separate call through public `bayes_compare(..., universe=)`
using alternating backends and one fresh bridge repetition also returns seven
finite evidences and probabilities summing to one. See
[model probabilities](relative_share_model_probabilities.csv),
[differences](relative_share_model_probability_comparisons.csv) and
[public mixed comparison](public_mixed_comparison.csv).

The complete CF conditional-k family restores prior (0.2,0.5,0.3) exactly once:

| Assignment | Family log evidence | Posterior k=1 | Posterior k=2 | Posterior k=3 |
| --- | ---: | ---: | ---: | ---: |
| CmdStan | -12.44764 | 0.05138 | 0.32182 | 0.62680 |
| Nutpie | -12.44789 | 0.05216 | 0.33103 | 0.61680 |
| Mixed: Nutpie / CmdStan / Nutpie | -12.45707 | 0.05264 | 0.32487 | 0.62249 |

The homogeneous family evidence difference is -0.000251 log units (combined
uncertainty 0.05645); mixed minus CmdStan is -0.009428 (0.05672). Maximum
posterior-k difference is 0.009992 for Nutpie and 0.004305 for mixed sampling.
The identity `Z = sum(prior_k * Z_k)` and complete original support are checked;
no failed k is discarded. See [family results](relative_share_conditional_k.csv),
[weights](relative_share_conditional_k_weights.csv) and
[comparisons](relative_share_conditional_k_comparisons.csv).

## Policy distributions and attrition

All 52 policy-summary comparisons with estimable uncertainty lie within four
combined MCSE. Sixteen pairs have matching constant threshold indicators,
which do not establish population tail equivalence. Thirty median/quantile
comparisons lack MCSE because of attrition and remain descriptive, not certified
equivalent. Every planned draw, paired-validity mask and failure reason is
retained. For example:

| Case / policy | Backend | Valid / planned | Mean | SD | 2.5% | 97.5% |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Logit C merger | CmdStan | 1997/2000 | 6.2347 | 0.8487 | 4.3413 | 7.6996 |
| Logit C merger | Nutpie | 2000/2000 | 6.2562 | 0.8668 | 4.3564 | 7.8286 |
| Logit CF merger | CmdStan | 1937/2000 | 6.1044 | 2.3184 | 1.0072 | 10.1062 |
| Logit CF merger | Nutpie | 1922/2000 | 6.4433 | 2.0698 | 0.9952 | 9.7674 |
| Logit CS tariff 25% | CmdStan | 1997/2000 | 3.9402 | 1.1504 | 2.0321 | 7.0581 |
| Logit CS tariff 25% | Nutpie | 1998/2000 | 4.0128 | 1.1481 | 2.0882 | 7.0515 |
| CES C merger | CmdStan | 2000/2000 | 9.0001 | 1.5246 | 6.0803 | 12.3030 |
| CES C merger | Nutpie | 2000/2000 | 8.8429 | 1.6098 | 5.4076 | 12.0066 |

Values are percent price changes for these fixed validation policies. Rows with
attrition describe the paired-valid conditional distribution. The largest
scaled mean difference shown is CF merger, 2.046 combined MCSE; the CF merger
probability of exceeding 5% differs by 2.528 MCSE. These are not exact-equality
tests. See [all policy summaries](relative_share_policy.csv),
[comparisons](relative_share_policy_comparisons.csv),
[conditional MCSE](relative_share_policy_conditional_mcse.csv) and
[failure records](relative_share_policy_failures.csv).

The posterior-k empirical policy mixture is also constructed. CmdStan/Nutpie
merger means are 6.3851/6.4309 and tariff means 5.0029/5.0300 (see the exact
linked values); retained paired-valid mass is 0.98080/0.97962. Mixture summaries
are descriptive: no joint-chain or bridge-weight uncertainty is imputed.
See [conditional policies](relative_share_conditional_policy.csv),
[attrition](relative_share_conditional_policy_attrition.csv) and
[differences](relative_share_conditional_policy_comparisons.csv).

## M4T4 diagnostics and the separate CES experiment

| Case | Backend | Divergences | Depth hits | Max R-hat | Min bulk ESS | Min tail ESS |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Logit C M4T4 | CmdStan | 0 | 4 | 1.0071 | 541.3 | 657.3 |
| Logit C M4T4 | Nutpie | 0 | 0 | 1.0451 | 113.0 | 151.9 |
| CES C M4T4, finite errors | CmdStan | 0 | 0 | 1.0132 | 441.9 | 460.9 |
| CES C M4T4, finite errors | Nutpie | 1 | 0 | 1.0163 | 332.3 | 596.7 |
| CES C, noiseless single cell | CmdStan | 37 | 0 | 1.0249 | 188.3 | 73.8 |
| CES C, noiseless single cell | Nutpie | 3 | 0 | 1.0285 | 186.6 | 101.6 |
| CES CS, noiseless single cell | CmdStan | 4 | 0 | 1.0134 | 450.1 | 317.6 |
| CES CS, noiseless single cell | Nutpie | 2 | 0 | 1.0238 | 148.5 | 75.2 |

All diagnostics remain visible and do not suppress a finite evidence estimate.
[The full table](diagnostics_summary.csv) includes each complete case/backend
and fit wall time (warmup included; cache lookup/first-use work may also be
included). [Compile records](compile_records.csv) are separate; resumed cache
hits are explicitly identified. These timings do not establish a universal
performance ranking.

The original exactly noiseless CES M4T4 attempt stalled as inferred residual
scales approached numerical precision. It was stopped within the execution
budget. Saved CmdStan CSVs contain **6/0/0/0 numeric retained rows**, with no
complete retained fit; Nutpie was not attempted for that original panel.
Partial rows, warnings and manifest snapshots remain in the durable workspace.

[The separate finite-error runner](../relative_share_noisy_ces_panel.R) keeps
the structural truth, prices, priors, pooling, kernels and frozen controls.
Its prescribed perturbations are `0.015*sin(j+2*t+m)` on log quantities,
renormalized to each cell's original total, and `0.02*cos(2*j+t+m)` on log
margins. This changes observations and is a distinct experiment. Both backends
receive identical perturbed data. Both complete, all 140 posterior comparisons
lie within four combined MCSE, and all six bridge repetitions are finite.
Together the two experiments have 26 complete fits out of 28 scheduled and
78 finite bridge calculations out of 84 planned slots: all 78 executed
calculations are finite, and six original slots remain unavailable. The
failed original pair is not converted into a success. The additional public
mixed-comparison call above uses seven fresh finite calculations and is
reported separately. See [attempts](noisy_ces_panel_attempts.csv),
[identity](noisy_ces_panel_paired_identity.csv),
[posterior comparisons](noisy_ces_panel_parameter_comparisons.csv) and
[evidence](noisy_ces_panel_evidence_comparison.csv).

## Saved-fit portability

Actual relative-share CF Nutpie save/reload preserves posterior arrays and a
seeded policy draw plan exactly. Moving its original shared library out of the
cache triggers reconstruction from the matching maintained source and saved
data. Eight original/reconstructed normalized log densities agree exactly
(maximum error zero). The original library is restored and its saved byte hash
rechecked; the fit file and scientific provenance are unchanged. See
[checks](saved_fit_checks.csv) and [runner](../relative_share_public_checks.R).

## Preserved failures and interruptions

Initial native stress checks exposed QR roundoff at sigma=2^-34. Orthogonal
cell-component preconditioning fixed it without changing the density, adding
jitter or relaxing tolerances. A separate 15-digit JSON fixture-rounding issue
was isolated; only the analytic stress fixture uses verified 17-digit
serialization. Production serialization is unchanged.

The first pre-correction sampling attempt and a later interrupted session are
preserved under `.scratch/option-b-implementation-20261022`. Source-matched fits
are reused only after hash verification. Initial rich-model validation calls
incorrectly supplied `use_cutoff=FALSE` together with explicit rich structures;
the API rejected them before sampling. Those failures are retained and the
corrected calls omit the redundant flag. No scientific target was changed.

The final check/opt-in logs and all raw fits, partial fits, posterior draw plans,
manifest snapshots and source snapshots are retained at
`.scratch/option-b-implementation-20261022/final-validation`. Earlier numerical
and interrupted attempts remain beside it. Compact evidence is tracked here;
raw fits, native libraries and caches are not publication artifacts.

## Downstream integration

No mergerBayes production code, data, priors, sampler profiles, results or Slurm
jobs were changed. Its inspected banking/synthetic callers explicitly request
sampled/collapsed market-intercept coordinates. They therefore need an explicit
basis choice: `share_basis = "absolute"` preserves that specification; adopting
relative shares requires removing those intercept-only controls and any joint
intercept/s0 rotation. Current-fit specification/provenance must be refreshed.
Runner guards and cache fingerprints should record `share_basis`. The mismatch
between market-pooled outside shares and DGP outside shares that vary after
period shocks remains a separate modeling decision; this implementation does
not change time pooling. No paper estimates are regenerated here.
