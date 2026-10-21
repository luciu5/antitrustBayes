# Nutpie promotion follow-up, 2026-10-20

Starting commit: `8d6daa1`. This is a new sampler-control experiment, separate
from the completed [paired study](nutpie_backend_results.md) and
[evidence study](nutpie_evidence_results.md). Those results and their manifests
remain unchanged. CmdStanR remains the default; Nutpie evidence/BMA stays disabled.

Execution update: the nine-fit pilot completed. No adaptation qualified across
all three cases; independent confirmation was not launched. See the
[actual diagnostics and preservation audit](nutpie_certification_pilot_results.md).

The previous study established target-density identity and compatible posterior,
policy and evidence estimates within Monte Carlo uncertainty. It did not certify
sampling quality: every evidence fit had divergences, and CES Stackelberg merger
SD/tail uncertainty was substantial. This follow-up changes sampler controls only.
All four Stan kernels, data, priors, transformations, conduct and policy definitions
remain fixed. No downstream production analysis is part of this study.

## Frozen diagnostic pilot

The [runner](nutpie_certification_pilot.R) takes the original sealed B, pooled CF
and CES CS fixture arguments, rather than generating new observations. It compares
CmdStanR, Nutpie diagonal and Nutpie low-rank adaptation: nine fits total.

| Control | Value |
| --- | --- |
| Seed | 419977, distinct from the original study |
| Chains / concurrent cores | 4 / 2 |
| Warmup / retained per chain | 2,000 / 1,000 |
| Target acceptance / maximum depth | .99 / 13 |
| Initial values | 0, as in the original study |
| CmdStan metric | diagonal for single-market; dense for pooled CF |
| Nutpie adaptation | diagonal and low-rank, recorded separately |
| Other adaptation controls | existing backend defaults; divergence threshold unchanged |
| Budget | one attempt per task; 30-minute worker timeout; no automatic retries |

The pilot collects posterior and sampler diagnostics, with no policy solves or
bridge estimates. Every fit is saved before diagnostic assessment. Reloaded draws
must match exactly; original reference identity must match in Stan source/data,
observations, model specification, priors, normalization and k support. The new
manifest seals the runner, production R source, dependency runtime, original
manifest, original reference fits and four Stan sources. Changes require a new
output root; the original validation root cannot be reused.

Each individual fit must have zero divergences and depth hits, BFMI at least .3,
maximum R-hat at most 1.01 and minimum bulk/tail ESS at least 400. Missing metrics
or incomplete chains fail this assessment. A backend/adaptation is eligible for
the proposed cross-case confirmation only if it passes all three cases. This
eligibility is an experiment outcome, not a backend permission check.
Exactly constant columns are exempt from undefined R-hat/ESS only when their
names and values match constants in a sealed original reference fit. The current
runner uses original native references for native fits, preserving full numerical
precision, and original CmdStan references for CmdStan. The executed first runner
used CmdStan constants for all fits; a separate no-HMC audit corrected its CSV
rounding artifact without changing any pass decision. A newly stuck parameter
cannot obtain an exemption merely because all its draws agree.

Failed fits remain in the report. There are no replacement seeds, automatic
target-acceptance increases or removal of difficult cases. If none qualify, report
the failure and diagnose it before designing another separately labeled experiment.

The native diagnostic audit identified nuts-rs 0.18.2 in the installed binary and
the pinned nutpieR lockfile. Its `max_energy_error = -Inf` on divergent trajectories
is an intentional sentinel; `diverging`, `divergence_energy_error` and
`divergence_message` remain authoritative. Recorded divergences contain genuine
energy-error failures. See the [upstream collector](https://github.com/pymc-devs/nuts-rs/blob/v0.18.2/src/stepsize/dual_avg.rs#L131)
and [exact dependency lockfile](https://github.com/andytimm/nutpieR/blob/8e88a975d5e0b8ed3c76de84d22adb33de601a1e/src/rust/Cargo.lock).

The frozen first execution used:

```sh
bash doc/validation/run_nutpie_certification.sh --prepare  # seals inputs; no HMC
bash doc/validation/run_nutpie_certification.sh --pilot    # explicit nine-fit experiment
bash doc/validation/run_nutpie_certification.sh --collect  # no HMC
```

Raw fits/logs/manifests persist in the ignored
`.scratch/nutpie-certification/pilot-accept99/`. The wrapper uses the existing
optional library and a compatible structural dependency library; the latter is
read only. Override `AB_NUTPIE_LIBRARY` and `AB_NUTPIE_STRUCTURAL_LIBRARY` when
using another compatible installation. Compiled artifacts are reused from the
original external cache. Cache-lookup time is reported separately; this is not a
new cold-compile benchmark. `estimation_seconds` is whole-fit elapsed time, with
the same scope for all backends. Native warmup and sampling phases are not falsely
split into comparable timings.

The current runner contains the reviewed reference-selection correction and
therefore rejects the first execution's immutable manifest. Do not reseal it or
refit these states to correct an assessment. The separate audit reads its frozen
runner and preserves all original artifacts:

```sh
AB_NUTPIE_DIAGNOSTIC_AUDIT_EXECUTE=true \
  Rscript --vanilla doc/validation/nutpie_pilot_diagnostic_audit.R
```

Future experiments require a new `AB_NUTPIE_PILOT_OUTPUT` root; their CSVs are
written under that root's `reports/` directory. Preparing the corrected runner
under `.scratch/nutpie-certification/preflight-native-constants/` performed no
HMC and did not overwrite any first-execution report.

## Independent confirmation after pilot assessment

The manifest records a proposed confirmation budget with independent seeds
420101 and 420211, four chains, 2,000 warmup and 2,000 retained draws per chain,
target acceptance .99 and depth 13. The pilot does not execute this stage or
choose an adaptation automatically.

The minimal confirmation set is the same-observation B/MC/C evidence universe
and CES CS policy fixture on both backends, for both seeds: 16 fits. All individual
fits must meet the diagnostic criteria; combining seeds cannot hide a failure.
Any claim favoring difficult hierarchical models additionally needs pooled
C/CF/CS confirmation, adding 12 fits. Neither stage changes the package default.

Evidence uses three bridge repetitions per model/backend/seed, retaining the
four-combined-uncertainty comparison, .25 log-unit precision cap and .05 maximum
model-probability difference. Public comparison dispatch, a complete conditional-k
family and any mixed-backend universe still need real integration validation.

Policy comparison must establish precision as well as agreement. Proposed
practical caps for four combined MCSEs are .25 price percentage points for means
and interval endpoints, 10% of reference posterior SD for SD comparisons, and .03
for threshold probabilities. Differences must also lie within four combined
MCSEs. Exact constants are checked as identities; missing MCSEs cannot pass.
Attrition is reported explicitly and prevents a claim about the unconditional
policy distribution. These caps are follow-up design choices, not retroactive
criteria for the original study. The confirmation runner and independent review
must freeze the final policy budget and selection before confirmation sampling.

## Remaining production qualification

In addition to clean posterior/evidence confirmation, qualification requires a
compatible dependency-stack CI route, saved-fit use without optional native
dependencies, public evidence dispatch, complete conditional-k evidence/policy
mixing, and target checks for pooled arrays and sampled-quality predictive paths.
Mixed-backend evidence is a separately validated capability. No fit metadata or
user option may unlock the evidence gate. A future package-owned certification
record must state its supported kernel/specification/evaluator scope.
