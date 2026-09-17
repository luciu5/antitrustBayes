# Restart-safe benchmark jobs

The helper in `checkpoint.R` is a benchmark utility, not a package API or a
change to inference. One worker owns each job id. Use a distinct job for each
fit, target fixture, bridge repetition or counterfactual batch so a shutdown
does not lose an entire ladder.

Before starting, persist the complete plan: data or data digest, source and
installed-kernel digests, compiler/build flags, priors, likelihood/conduct and
pooling options, initial values, explicit seeds, sampler/solver settings and
posterior draw/quality/k plans where applicable. Pass these as `specification`.
Changes create a new fingerprint directory; they cannot reuse an old result.
The wrapper cannot infer a closure's captured inputs, so they must be included.

Use the ignored workspace artifact directory, not `/tmp`, for results,
libraries, source worktrees, compile objects and logs. Keep only source,
compact CSVs and reports in Git. Do not rely on a sequence of global RNG draws:
freeze each job's seed and initial values before the first fit.

Successful jobs save their payload before the completion record, with atomic
same-directory renames and a SHA-256 checksum. Matching completed jobs load
their existing result. Errors and warning text are retained; failed jobs are
not rerun by default. Explicit retries create new attempts and retain failures.
No job is declared successful merely because its process exited: the caller
still applies the unchanged scientific/diagnostic gates to the returned data.

An unfinished attempt is not automatically assumed dead. After a known host
restart, confirm no other writer owns that job and use
`resume_incomplete = TRUE`; its previous start record is retained and marked
`INTERRUPTED`. Active writers must not be displaced. The utility introduces no
OS, power-management, account or scheduler changes.

Compute timing excludes checkpoint writes and records UTC endpoints as well
as `system.time()` components. Host sleep or unusual pauses can invalidate a
wall-time comparison even if results complete. Report those cases and rerun
the identical prespecified job as a separate timing repetition; do not adjust
priors, tolerances or sampler settings to obtain a better result.

Fresh native comparisons are required after a source, compiler, dependency or
kernel change. Persistent object files are reusable only when the recorded
build specification matches; a surviving file is not proof of correctness.
