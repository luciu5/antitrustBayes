# Optional bridge execution specification

Serial bridge execution remains the default. A future accepted implementation
may expose a positive integer `bridge_cores` and an optional `bridge_seed`.
Each job is one model/repetition; a model context is restored and validated
once, and repetitions retain the current maxiter and evidence diagnostics.

A reproducible comparison must allocate explicit per-job random streams before
dispatch and use the same streams in serial and parallel execution. The
existing serial call without an explicit seed must retain its current random
number behavior. Restore the caller RNG state when an explicit seed is used.
Do not interpret failed worker outputs as successful bridge results, renormalize
away positive-prior failed hypotheses, or lower repetition safeguards.

Fork workers can share read-only fit/context pages on supported platforms,
but bridge working memory is additional. Keep the worker count explicit and
conservative, collect model/repetition error records, and benchmark process-tree
PSS as well as runtime. Unsupported platforms should use documented serial
execution rather than silently promise parallelism.

A source optimization changes source hashes even when its normalized target
is equivalent. Saved-fit restoration requires explicit validated source-hash
pairs and unchanged positive-width raw parameter schemas. A pure Helmert
algorithm change must preserve both legacy and rotated market-s0 flags; the
older original-to-rotation migration still fills missing zero-width fields
and requires the original legacy geometry. Unknown source hashes must fail.
# Saved-kernel compatibility after an exact transform

The linear-time Helmert implementation changes source hashes while preserving
the raw parameter schema and normalized target. Add its hashes to the explicit
allow-list only after the compiled target/gradient/output gate passes for all
four kernels. Validate both rotation-to-Helmert and the composed original-to-
Helmert migrations; arbitrary source changes must still fail closed.

An original pre-rotation source requires the private legacy rotation flag of
1, as before. A saved rotation-era source may legitimately have flag 1 or 2:
preserve that flag when migrating to the Helmert source. Do not reuse the old
blanket `migration => flag 1` rule for every approved performance migration,
because that would change a rotated fit's sampling coordinates and density.
Keep serialized fits immutable and retain all raw-schema/zero-width checks.
