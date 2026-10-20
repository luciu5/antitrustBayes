# Nutpie backend validation: execution handoff

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
