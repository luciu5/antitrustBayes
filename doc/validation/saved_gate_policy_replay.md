# Sealed-fit and current-policy replay

The November 1 specification-14 gate has two completed fits and frozen
inference thresholds. Later package metadata may differ: the current package
is `0.3.0`, while those fits retain `0.2.0.9009` provenance. The original
`coordinated_release_gate.R` correctly rejects any full-source seal mismatch.
Leave that runner, its manifest and its dated report unchanged.

The sampling-free `saved_gate_policy_replay.R` verifies the committed artifact
seals, manifest call hash, complete two-fit support, saved fit digests,
specification, kernels, data and frozen sampler controls. Every sealed R/Stan
and original research-runner file must match byte for byte. It refuses changed
or added production files. A reviewed DESCRIPTION-only change requires the
explicit acknowledgment flag and appears in the source comparison table;
it is not hidden by updating a historical hash.

Run from this repository with an explicitly compatible structural library:

```sh
R_LIBS=/path/to/compatible/library \
  Rscript --vanilla doc/validation/saved_gate_replay_checks.R
R_LIBS=/path/to/compatible/library \
  Rscript --vanilla doc/validation/saved_gate_policy_replay.R \
  /path/to/read-only/saved-gate /path/to/fresh/replay-output \
  --allow-description-change
```

The output must be separate from the saved gate and empty. Existing output
cannot be overwritten. The script reads large fits in place and writes compact
policy residuals, saved-fit diagnostics, source comparison, warnings and a
provenance RDS. It checks all 32 baseline/merger FOCs at the original truth
states and fixed draw 83, with the original `1e-6` residual limit. Demand,
cost and ownership equations supply the oracle; package margin/share helpers
are not an oracle. It also checks source-object immutability and recomputes
the retained fits' original diagnostic bounds. Runtime package hashes are
recorded separately from the historical inference environment.

This is a policy replay using verified saved inference, not new HMC, broad
recovery, conduct-identification certification, or a replacement for the
incomplete noiseless CES panel attempt. Historical unavailable results remain
unavailable. A change to scientific source requires a separate reviewed
validation design rather than this replay route.
