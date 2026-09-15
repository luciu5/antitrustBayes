# Computational diagnostic rerun

Written after seven completed initial fits failed their prespecified gates,
before any PNB fit. Initial call/data checkpoints and traces remain intact.
The sample includes CES simultaneous Bertrand with 244 depth hits and
minimum BFMI .2688, and simultaneous Logit cases with depth/BFMI failures.
No posterior recovery or model odds are inferred from these fits.

A controlled sampling experiment may use four chains, 2,000 iterations,
1,000 warmup iterations, and `control$metric="dense_e"`. Adapt delta stays
.95, maximum depth stays 13, and the original acceptance gates stay unchanged:
zero divergences, zero depth-limit hits, BFMI > .3, and maximum Rhat <= 1.01.
The target, compiled model, priors, DGP, seeds, cutoff support, and role
mapping are unchanged. The [RStan control documentation](https://mc-stan.org/rstan/reference/stan.html)
supports this metric. It can address linear correlation; it is not a promised
remedy for nonlinear hierarchical geometry.

Use a separate output directory and fingerprint the full changed sampling
call. A `run_sampling_experiment.R --sampling-profile=dense --fit-only --limit=4` pilot covers all four
demand/conduct cases at the first seed and layout, each with all three
candidate structures. It is a computational pilot, not a recovery study.
Apply a successful computational configuration across all prespecified seeds
and layouts; do not select attractive replications. Keep failed reruns and
unknown diagnostics. The other seeds are mandatory before any recovery claim.

The initial short run may be paused or stopped if the diagnostic review
establishes that continuing it adds repeated failures without resolving the
computational gate. Record every completed, interrupted, and not-started
configuration separately. Do not classify interrupted or unstarted work as
successful or failed posterior recovery. Any full followup experiment keeps
all originally prespecified datasets and candidate economics.

The initial process was stopped with exit status 130 after ten completed
fits all failed diagnostics. Timestamped fit-start/fit-done records and
checkpoint inventory distinguish ten completed, six interrupted without a
checkpoint, and 128 not-started fits across the full 144-fit plan. All initial
checkpoints, logs, and before/after ledgers are retained. New experiments
persist a planned manifest and a separate per-job start/completion ledger.


## Equivalent quality coordinates: native specification 8

The longer Native7 dense pilot has eight completed simultaneous/fringe fits, all failing the original diagnostic gates; four Stackelberg fits are still running. Before another posterior experiment, the lead derived and independently numerically verified an equivalent centered orthonormal quality representation. See [quality_coordinates.md](quality_coordinates.md) and [native8_gate_summary.csv](native8_gate_summary.csv). Fresh independent implementation review passes. The native coordinate density and dimension differ, but the physical demand-quality prior, all economic likelihoods, normalized evidence, original priors, DGP and role rules are preserved exactly.

After full package gates pass, run a new isolated twelve-fit computational pilot using the same four-chain, 2000-iteration, 1000-warmup dense profile, .95 adaptation, depth13 and original acceptance gates. Use strict Native8 fit provenance, a new frozen runner directory, new output directory and installed-code fingerprint; do not reuse Native7 checkpoints or mutate its running files/libraries. Keep all original failed/interrupted/not-started rows. Performance improvement is conjectural until measured. A passing configuration must then be applied to all original seeds/layouts/candidates before any recovery claim or PNB fit. The full-rank 36-panel conditional-demand design diagnosis rules out a simple design-rank deficiency, but establishes neither conduct identification nor posterior recovery.
