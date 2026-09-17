# Performance audit results

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
not an end-to-end speedup or a welfare comparison. Warmed repetitions, further
conduct fixtures, welfare checks and larger extraction-only probes are pending.

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
speed ratio. No candidate Stan source has been applied to the main working files.

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
