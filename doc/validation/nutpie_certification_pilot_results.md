# Nutpie qualification pilot results, 2026-10-20

The separately frozen nine-fit pilot completed from `nutpie` commit
`8d6daa115bf2fb7888a26d1c3ca4ba6666355c64`, package `0.2.0.9005`.
No backend/adaptation passed all three cases. General qualification and the
independent evidence/policy confirmation stage remain blocked; Nutpie stays
experimental and public evidence/BMA remains disabled. CmdStanR is still default.

This is useful progress, with a narrower conclusion than promotion: the new
B and pooled CF fits had no divergences, but pooled Nutpie mixing and CES
Stackelberg geometry remain unresolved.
Low-rank adaptation did not consistently improve the diagnostic outcome.

## Fixed experiment and target identity

The [protocol](nutpie_certification_protocol.md) fixed four chains, two concurrent
cores, 2,000 warmup plus 1,000 retained draws per chain, seed 419977, target
acceptance .99, depth 13 and initialization 0. There were exactly nine fits,
no replacement seeds, no retries and no automatic control changes. CmdStan used
its diagonal metric for single-market fits and dense metric for pooled CF;
Nutpie diagonal and low-rank were recorded as distinct adaptations.

All nine fits matched their sealed original reference in source/data hashes,
observations, specification, priors, normalization and k support. Each retained
all 4,000 states, with exact posterior equality after serialization/reload.
All four maintained Stan source hashes remained unchanged. Compilation used
the existing external caches; the installed package was not made mutable.
The [parameter/provenance CSV](nutpie_certification_pilot_parameters.csv)
records each kernel, source/data hash, selected posterior summaries and MCSEs.

The frozen executable SHA-256 is
`e1de5bb93d1ba0295db5dbf373c8c472393c13219dac02968d8a5642e74795a7`;
its exact snapshot is retained in
`.scratch/nutpie-certification/pilot-accept99/frozen-runner.R`.
The manifest was sealed at 08:42:11 EDT. All fits/results, arguments, logs,
posterior hashes and native provenance remain under that output root.

## Actual diagnostics

Whole-fit seconds cover cached-model `bayes_est()` elapsed time, including
sampler setup and posterior extraction. Cache lookup/preparation is recorded
separately in the [original diagnostic CSV](nutpie_certification_pilot_diagnostics.csv);
neither measure is a new cold-compilation benchmark. Warmup and sampling phases
are not presented as equivalent backend-native timing fields.

| Case | Backend/adaptation | Fit seconds | Divergences | Depth hits | Min BFMI | Max R-hat | Min bulk ESS | Min tail ESS | Diagnostic assessment |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| B | CmdStan diagonal | 6.974 | 0 | 0 | .843 | 1.0061 | 974 | 910 | Pass |
| B | Nutpie diagonal | 3.299 | 0 | 0 | .808 | 1.0032 | 720 | 598 | Pass |
| B | Nutpie low-rank | 3.904 | 0 | 0 | .713 | 1.0092 | 562 | 267 | Tail ESS |
| CF pooled | CmdStan dense | 88.645 | 0 | 0 | .729 | 1.0098 | 603 | 608 | Pass |
| CF pooled | Nutpie diagonal | 69.553 | 0 | 0 | .730 | 1.0197 | 407 | 640 | R-hat |
| CF pooled | Nutpie low-rank | 79.782 | 0 | 0 | .699 | 1.0126 | 319 | 679 | R-hat, bulk ESS |
| CES CS | CmdStan diagonal | 9.109 | 69 | 0 | .751 | 1.0185 | 218 | 73 | Divergences, R-hat, ESS |
| CES CS | Nutpie diagonal | 4.699 | 4 | 0 | .756 | 1.0048 | 615 | 543 | Divergences |
| CES CS | Nutpie low-rank | 4.699 | 3 | 0 | .755 | 1.0076 | 701 | 395 | Divergences, tail ESS |

The assessment uses the declared zero-divergence/depth-hit, BFMI .3, R-hat 1.01
and bulk/tail ESS 400 requirements, separately for every fit. Metrics summarize
all retained posterior columns; constants are handled explicitly below.
Actual values remain in the report even when a criterion fails.

Pooled Nutpie diagonal's largest R-hat is in `mu_log_a` (and its free coordinate),
with market-alpha raw coordinates also above 1.01. Low-rank's largest R-hat and
smallest bulk ESS are in the alpha hierarchy SD. CES CmdStan mixing difficulties
include the margin residual scale, rho and collapsed likelihood terms. This
localizes follow-up; it does not establish a model error or an upstream sampler
correctness bug. CES divergences occurred on both samplers, and the independent
seed/control change cannot be interpreted as a controlled causal comparison
against the original study.

Selected pooled CF hierarchy posterior means and 95% intervals are below.
These describe the pilot and do not override its diagnostic failures.

| Quantity | CmdStan dense | Nutpie diagonal | Nutpie low-rank |
| --- | --- | --- | --- |
| Alpha hierarchy SD | .2353 [.0790, .5723] | .2408 [.0829, .6225] | .2476 [.0823, .6030] |
| Outside-share hierarchy SD | .1088 [.0080, .4569] | .1042 [.0081, .4363] | .1011 [.0094, .4024] |
| Market-intercept hierarchy SD | .3141 [.1232, .7278] | .3196 [.1239, .7630] | .3208 [.1219, .7525] |

## Constant-reference audit and preservation

Independent review identified a validation-only precision issue after the
fixed pilot had started. Original CmdStan CSVs retain eight significant digits,
while saved native draws retain full precision. Exact constant-value checks
therefore falsely flagged nine pooled CF and two CES CS generated-quantity
columns per native fit. No posterior state, sampler control, target or numerical
diagnostic was affected.

The [no-HMC audit](nutpie_pilot_diagnostic_audit.R) sourced the exact frozen
runner and reproduced every original availability assessment. It then used
sealed original native-diagonal reference constants for both native adaptations,
retaining CmdStan references for CmdStan. Target/data/prior identities matched
for every reference. Exact name/value checks, missing-diagnostic rejection and
all thresholds stayed fixed. Newly stuck parameters cannot gain a constant
exemption. The [audited chain table](nutpie_pilot_diagnostic_audit_chains.csv)
and [variable table](nutpie_pilot_diagnostic_audit_variables.csv) retain original
and audited results side by side.

All nine pass/fail decisions were unchanged: three cells pass and six fail;
no adaptation qualifies across all cases. All 67 sealed inputs remained
unchanged, including original fit/result files and the original tracked CSVs.
The audit did not invoke sampling or generate replacement states.

The current pilot runner fixes reference selection and seals the additional
native references before inference. New runs keep reports inside their own
output root, preventing a fresh preparation from replacing this study's CSVs.
Its changed source cannot resume the old manifest. A separate preparation
under `.scratch/nutpie-certification/preflight-native-constants` checks the fix
without HMC. To inspect the executed study, use the preserved artifacts and
the explicit no-HMC audit, rather than resealing its original manifest.

## What advanced and what remains

The [five-fit fresh-process audit](nutpie_serialization_results.md) passed
without Nutpie/BridgeStan available, including all four kernels, pooled CF,
exact extraction/diagnostics/quality/plan checks and actual merger/tariff solves.
Missing or changed native evidence artifacts fail explicitly. The
[dependency-stack CI changes](nutpie_ci_readiness.md) fix incompatible pins,
explicitly test an optional-free installation and retain the opt-in CmdStan tier.
A hosted CI run and a clean native installation tier remain unverified.

Source checks for backend/policy/runner behavior passed 562 assertions in 54
blocks. The same filter without optional native dependencies passed 555, with
one explicit missing-nutpieR skip and no failures, errors or warnings. The final
three new runner/audit test files passed 49 assertions in 11 blocks; corrected
manifest preparation performed no HMC. YAML/embedded-shell/R static checks and
`git diff --check` passed. Logs and test result objects are retained under
`.scratch/nutpie-certification/`. Production R and package dependencies were
unchanged in this follow-up, so the earlier full CmdStan regression evidence
was preserved rather than rerun for research-runner and CI edits.

No new policy distributions or bridge estimates were computed from this pilot.
Those were intentionally reserved for independent confirmation after diagnostic
qualification; the prerequisite failed. The earlier posterior/policy/evidence
results remain available with their original uncertainty and divergence caveats.

The next statistical work is a focused examination of the unchanged CES target's
divergent trajectories and pooled alpha hierarchy mixing, before another
separately reviewed sampler experiment. No Stan or economic specification change
is authorized by this result. Independent confirmation, precise CES policy tails,
public evidence dispatch, full-support conditional-k and mixed-backend evidence
remain required before broader promotion. CmdStanR remains default and Nutpie
remains experimental. No mergerBayes production analysis was changed.
