# Performance audit results

> **Dated audit evidence: September 2026.** For current development, use the
> [documentation index](../README.md). Revision and session observations below
> describe the recorded audit; the retained measurements and limitations remain.

Execution resumed after an environment reset on September 17. Workspace
sources and the compact evidence below survive; the latest `/tmp`-only paired
vectorization HMC arrays do not. Their comparison is incomplete and establishes
no sampler speedup. Recovery and retained partial telemetry are documented in
`recovery_status.md`.

The subsequent September 17 recovery now builds and verifies the four Helmert
kernels. The main Stan/header files match the tested candidate and have been
applied after the ordinary installed suite passes (13 skips, 13 warnings;
the added skip is the optional antitrust lifecycle absent from recovery deps).
Seven final boundary fixtures have matching targets, full gradients, schemas,
outputs and nonvacuous quality-centering checks. The original central-difference
gate is 6/7: one cancellation failure is retained and resolved separately with
fourth-order differences agreeing with each kernel's analytic derivative within
8.87e-12, under the same 1e-10 threshold. See
`recovery/native_review/fd_resolution.md` and its original/resolved CSVs.
The older chronological candidate-hold statements below describe earlier stages.
This establishes numerical equivalence, not a new HMC ESS/second improvement.

At the September 17 checkpoint, the user had committed and pushed the four
reviewed source/report units; main and the locally recorded origin/main matched
ea29e96. The working tree was clean
before this next audit unit. Optional parallel bridge execution remains outside
production: eight archived benchmark jobs (48 bridge estimates, 16 separate
public comparisons) have exact paired evidence and iteration counts, but only
1.137–1.184 speed ratios with 1.783–1.784 peak PSS ratios. The candidate warning
loss was reproduced and repaired in isolation, with 13 passing assertions and
a fresh Luna read-only review. Current serial defaults and economic targets
remain unchanged. See recovery/bridge/archived_benchmark_results.md.

## Independent chain scheduling

The first accepted unit resolves the default process count after merging user
options: `min(chains, getOption("mc.cores", 1L))`. An explicit `cores` wins.
This does not enable within-chain threading or change any density or HMC setting.

The frozen comparison comprises twelve fits: simultaneous Bertrand, exact BF,
and exact BS; two sampler seeds; one or four processes. Each fit uses the same
N=24 noisy-share fixture, four chains, 1,000 iterations/500 warmup,
Gaussian collapsed shares, fixed-t5 margins, rho off, and existing priors.
Every pair has identical full unpermuted posterior hashes and leapfrog counts.
All twelve fits have zero divergences and zero treedepth hits. These facts
establish scheduling equivalence for the fixtures, not general sampler recovery.

Sampling-call wall speedups range from 1.561 to 2.169. Peak process-tree PSS
ratios range from 1.891 to 2.080. PSS is polled at .25-second intervals and may
miss brief peaks; summed RSS double-counts shared fork pages. Native chain
sampling clocks include generated-quantity execution. This experiment does
not separately measure GQ or establish large-N/many-market scaling.

The focused policy tests pass 41 expectations, including interception of the
public `bayes_est()` sampling call and restoration of the caller's options.
A fresh independent review identified a test-option leak; it was repaired and
retested before source integration. The installed ordinary suite also passes
with twelve opt-in sampling skips and thirteen expected warnings. Two initial
test-harness failures are retained in the ignored failure ledger. No Stan or
economic code changed.

Compact metrics are in `chain_results.csv` and the requested comparison table
is in `summary.csv`. Full fits, logs and per-process monitors remain in the
ignored statistical-toolkit run directory.

## Remaining measurements

Executed lifecycle counters confirm one nonlinear price solve during baseline
construction for each Logit/CES Bertrand/Cournot Bayes draw, followed by the
actual policy solve. Trade's supplied-parameter Bertrand constructor solves
both an unchanged pre and post state. `baseline_solver_counts.csv` retains
baseline measurements; a public observed-baseline mode is being implemented
and must pass state/counterfactual equality before acceptance.

Selective extraction reference measurements are in
`extraction_reference_results.csv`. On the two small archived fixtures,
required copied fields use 33,736 and 130,536 bytes versus unrestricted
3,195,776 and 5,273,120 bytes. Required values match. Timings below .04 seconds
are noisy; process peak RSS decreases only modestly because library/fit memory
remains resident. These are reference extraction results, not yet end-to-end
production startup savings.

Native profiling completed sixteen finite target/gradient/GQ probes. Eight
standalone GQ rows match every common archived output exactly for each of four
fixtures. The compact report is in `native_profile/native_profile_report_v2.md`.
The small-fixture incremental GQ cost lies near timing noise; it does not justify
prioritizing a GQ split ahead of extraction and construction. The independent
Helmert reference agrees within 3.4e-15 and supports a compiled production gate.

## Selective extraction integration

The selected-draw contract was independently reviewed and applied to the Bayes
working files after repairs in temporary commits `32ecc4a` and `3ac3d42`.
The review reproduced and verified fixes for legacy outside-share logits,
required/optional missing fields, retained-versus-declared Stan variables,
stale draw-count metadata, historical Stackelberg cutoff summaries, and numeric
quality-mode codes without descriptive metadata. The focused tests pass 32
checks; canonical full/selected values match on archived Logit, CES, BF and BS.
Version-1 draw-plan fingerprints retain their explicit full-extraction path;
new version-2 plans bind the canonical economic/conditioning field union.

An initial combined full-suite attempt used an incoherent mixed source/installed
namespace and failed. The failure is retained in the ignored run directory.
An installed R package containing both accepted units passes the ordinary suite
with the same twelve opt-in skips and thirteen expected warnings. Its native
library is the unchanged baseline library, copied explicitly for this R-only
integration check; this is not a new Stan compilation.

The first 200-fixed-draw small-fixture benchmark has extraction times .008 to
.003 seconds and copied arrays 3,195,776 to 33,736 bytes. Total simulation is
2.630 to 2.956 seconds in this single pair, with all 200 draws valid and identical
draw indices, prices, shares and failure records. Peak process RSS is essentially
unchanged, 576,564 versus 577,392 KiB. This establishes a copied-array reduction,
not an end-to-end speedup or a welfare comparison.

The subsequent independent benchmark completed 36 paired policy rows across
ordinary Logit, ordinary CES, exact BF and exact BS: 200 fixed draw positions,
three policies and three warmed repetitions per fixture. Prices, shares, CV,
producer surplus, total welfare, draw identities and failure records match
exactly; every draw is valid. End-to-end gains remain inconsistent. A separate
200-row collapsed-quality Fixed_param extraction stress test at N=10,000 has
.834 versus .025 second extraction, 112,186,080 versus 19,336 copied bytes,
and peak RSS 881,796 versus 772,632 KiB. These are serialization measurements,
not HMC evidence. Compact results and the retained runner failure are in
`selective_extraction/`.

Shared baseline construction, threading, counterfactual scheduling and optional
parallel bridge evidence remain active audit units. Their speedups and
equivalence are not established by the chain benchmark.

The sampled-quality Helmert transform currently constructs a dense G by G-1
matrix. A mathematically equivalent linear-time recurrence is specified in
`sampled_quality_contrast.md`; independent reference and native gates precede
any kernel change. This matrix is separate from economic ownership.

The candidate compiled recurrence passes twelve isolated-process small nonzero
quality checks across all four kernels, BF/CF/BS/CS and the singleton boundary.
Target, full gradient-vector and complete fixed-output comparisons are zero;
coordinate finite differences agree within 1e-9. The full comparisons are in
`native_profile/root_transform_fixed_gate.csv`; they include transformed
parameters and generated quantities. Panel/rho gates and a matched-compiler
timing comparison remain required.
An earlier candidate report described compiler flags as matched, but clarification
established that the trusted old binary's flags were not verified. Those timings
are retained as unmatched-build measurements and do not establish a causal
speed ratio.

A subsequent matched rebuild uses identical GCC 11.4/R 4.6.1 compiler commands
and flags. Twenty full-vector fixed gates have target/output differences zero
and maximum gradient difference 4.02e-13. At N=G=10,000, target evaluation is
.045 versus .007 seconds, gradient .085 versus .011 seconds, and process peak
RSS 1,332,104 versus 554,920 KiB; target and gradient norm match. This is one
fixed-point large probe, not an HMC ESS/second comparison. Raw-schema checks
match for six fixtures. Canonical migration metadata uses installed Stan
`model_code` digests; file-byte digests are only labels. Independent panel/rho
coverage and saved-fit restoration still gate integration. No candidate Stan
source has been applied to the main working files.

A fresh independent reviewer now verifies all eighteen saved panel probes:
normalized and raw target differences zero; maximum full gradient, raw-gradient
and TP/GQ differences 3.55e-15; raw schemas, coordinate vectors, names, shapes,
finite patterns and finite-difference vectors match. Coverage includes sampled
and collapsed Logit/CES, missing-margin rho, legacy/rotated pooling, centered
and noncentered hierarchies, `X_s0`, and BF/CF/BS/CS. Additional combinations,
large-G nonzero stability and saved-fit compatibility still gate integration.
See `native_profile/helmert_independent_panel_comparison.csv`.

Saved-fit compatibility is independently reviewed and its helper is applied
to the main working files, with 106 focused checks passing. It explicitly
preserves rotation-era private coordinate flags instead of forcing legacy
coordinates. Forty saved rows across four kernels and both flags agree within
1.07e-14; 100 rows across five original-kernel archived fits agree exactly.
One short normalization fixture retains a divergence, so this is target
restoration evidence rather than an HMC or bridge-convergence claim. Native
source integration remains held for the final boundary gates. See
`kernel_compatibility.md` and `kernel_migration_*.csv`.

## Public baseline consumer

The Bayes consumer now selects `baseline = "observed"` only when the public
`antitrust::specify()` exposes that option; an explicit user option wins and
older antitrust installations retain their existing path. The integrated
installed R-only Bayes suite passes against the reviewed lifecycle candidate,
with twelve opt-in skips and thirteen expected warnings. It uses the unchanged
trusted native Bayes library. Counterfactual/state checks and solve counters
were executed separately.

The paired runtime comparison now completes 27 policy rows per path: Logit
Bertrand, Logit Cournot and CES Bertrand; merger, +10% cost shock and both;
three warmed repetitions, each with 200 matched posterior draw positions.
Prices, shares, CV, producer surplus, total welfare, draw identities and failure
records match exactly, with all draws valid. Construction price/solver calls
fall from 200 to zero per row; actual policy solver calls, iterations and
function evaluations agree. The unchanged baseline solves have zero iterations
on these fixtures, so removing them does not save a difficult nonlinear solve.
Median whole-policy speed ratios range from .772 to 1.217, including regressions;
no general wall-time gain is established. Peak process RSS is 607,900 versus
609,044 KiB over each complete nine-policy-batch process, not per policy.
The architecture separates baseline realization from policy solving, but the
performance evidence does not justify a claim of large startup savings. Compact
paired counts and numerical gates are in `baseline_consumer_*.csv`, reproduced
by `run_baseline_consumer_posterior_benchmark.R` and
`compare_baseline_consumer.R`. Two initial runner setup failures remain recorded.

The current trade branch advanced to `6609e28` during the audit, adding fitted
model promotion and requiring antitrust >=0.99.40. The earlier lifecycle
candidate installs were validated on the previous trade base. Its patch must
be reconciled with the newer code and dependencies before source integration;
neither actual cross-package worktree has been changed by this audit.
The user reports another modifying session may exist. Cross-package reconciliation
is therefore tested in temporary clones; actual working branches remain untouched.
The temporary current-base reconciliation installs successfully and its ordinary
suites exit zero: antitrust (77 documented extended/nightly tier skips, 122
warnings), trade (no skips, 71 warnings, including promotion tests), and Bayes
(12 opt-in skips, 13 warnings). It includes the pre-existing external antitrust
promotion patch as a separate compatibility prerequisite. Lifecycle development
versions preserve that newer dependency at antitrust 0.99.40.9000 and trade
0.8.8.9000. `current_refactor_suite_gate.csv` records this scope; it is not a
fresh 200-draw timing on the combined current-base tree.

Full firm ownership is economically block sparse, but legacy `ownerToMatrix()`
materializes a dense product matrix. Firm/block calculations are candidates
for a separately measured optimization. General partial-ownership behavior
must remain available. No ownership representation change has been accepted.

The separate well-conditioned ownership reference now has three warmed
repetitions for N=100, 300 and 1,000, including three-product firms. At N=1,000,
known-block solves take .002--.007 seconds versus .184--.220 seconds for the
dense inverse, with margin differences below 2e-16. Direct RHS solving alone
has little benefit in this fixture. These are linear-algebra timings; neither
counterfactual equivalence nor singular/partial-ownership behavior is established.
The reference still materializes dense matrices. See `ownership_rhs_reference.csv`.
