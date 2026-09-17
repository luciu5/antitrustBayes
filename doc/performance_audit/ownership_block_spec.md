# Ownership blocks: candidate, not an accepted implementation

The legacy ownership slots store base R matrices, even for full firm ownership
whose off-firm blocks are zero. Consequently, economic block sparsity does not
currently imply sparse storage or a sparse margin solve. Changing the slot type
is a broader compatibility change and is not approved solely from this finding.

The generic Bertrand margin calculation forms

\[
A=(E')\odot\Omega,\qquad b=r\odot\operatorname{diag}(\Omega),
\]

then computes a signed `solve(A) %*% b / r`, with a global `MASS::ginv(A)`
fallback. If ownership has disjoint blocks, A has the same block structure and
the nonsingular calculation can be performed block by block. This is a linear
algebra optimization, independent of the demand or conduct equations.

First measure the smaller change `solve(A, b)` against the current
`solve(A) %*% b`. It avoids materializing a full inverse without changing
ownership representation or partition detection. Retain the same global
`MASS::ginv(A) %*% b` fallback. Roundoff can differ, so compare actual margin
vectors and counterfactual equilibria before acceptance. Commit this change
separately from any block implementation if it passes and has useful savings.

Any candidate must support arbitrary current weighted ownership through the
existing path. Do not assume every binary matrix is a valid full-firm partition,
or infer firm IDs from product labels. Keep the global generalized-inverse
fallback: separate block pseudoinverses can use different singular-value cutoffs
and change the original answer. Likewise, global `solve()` may reject a matrix
with large scale differences between individually well-conditioned blocks.
Require a conservative global-conditioning gate or preserve the original path
for such cases. A candidate that cannot reproduce these decisions stays a
research benchmark rather than a default.

Benchmark margin evaluation and full merger/cost-shock solves, with solver
counts/iterations, for flat and nested Logit/CES, input/output markets,
multiproduct firms, weighted/cross ownership, active subsets, singleton blocks,
and near-singular fixtures. Compare against the original prices/shares/welfare
at existing tolerances. Report dense ownership and elasticity allocations
separately: reducing factorization work alone does not remove their N-squared
storage. No production implementation is established by this specification.

`ownership_rhs_reference.R` and its CSV record three warmed alternating-order
repetitions at N=100, 300 and 1,000, with one or three products per firm and
known full-ownership partitions. At N=1,000, dense inverse evaluation took
.184--.220 seconds, a direct right-hand-side solve .180--.192 seconds, and the
known-block reference .002--.007 seconds. All tested margin differences are
below 2e-16. Small timings have millisecond resolution. This only establishes
linear-algebra agreement on well-conditioned output-Logit fixtures; it does
not test singular fallback decisions, partial ownership, full equilibrium,
or memory reduction. The reference still constructs dense matrices.
