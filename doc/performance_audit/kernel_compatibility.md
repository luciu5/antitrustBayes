# Saved-fit compatibility gate for native performance changes

This records a private candidate and its executed restoration checks. It does
not approve arbitrary Stan source changes, establish bridge convergence, or
change any statistical or economic specification.

## Coordinate policies

The saved source digest is computed from the installed `stanmodel@model_code`,
not the bytes of a source file. A digest pair is approved separately for each
registered native kernel. Unknown, reverse and cross-kernel pairs remain
errors, and the existing raw parameter-schema checks remain required.

There are two policies for the currently checked generations:

* Original kernels precede the market/outside-share coordinate option. A
  private reconstruction uses legacy coordinates, adding data flag 1 only
  when absent and rejecting a contradictory or malformed flag.
* Rotation-era kernels already identify the selected coordinates with flag 1
  or 2. A performance-only migration must preserve that saved flag. Replacing
  a rotated fit's flag with 1 would evaluate its draws in different coordinates.

The linear-time Helmert recurrence preserves the sampled contrast coordinates,
their prior, the density measure and all economic outputs. It needs no new
Jacobian. That algebra alone is insufficient to register a source migration;
native normalized targets, gradients, output schemas and archived restoration
are checked independently.

## Executed restoration checks

Eight small fixtures cover the four native kernels with both coordinate flags.
Five saved posterior rows per fixture are restored against the Helmert
candidate. All 40 stored `lp__` comparisons agree within 1.07e-14. The coordinate
flags, raw schema, original fit and caller RNG state are preserved.

Those very short chains are normalization fixtures, not inference benchmarks.
One CES Stackelberg rotated fixture has one divergence, and the short-chain
ESS warnings are retained. No HMC or marginal-likelihood success is inferred
from them.

Five proper archived fits cover Logit Bertrand, Logit Cournot, CES Bertrand,
exact Bertrand fringe and exact Bertrand Stackelberg. They use the original
pre-rotation kernels. Twenty stored rows per fit are restored through the
composed original-to-Helmert path. All 100 stored `lp__` comparisons are exactly
equal; the original fits and caller RNG state remain unchanged.

Fresh independent review approves the compatibility helper within this scope.
The helper, tests, compact results and executable runners are now in the main
working files. Its 70 migration checks and 36 existing context-restoration
checks also pass against those files. Native kernel integration remains a
separate gate; merely registering its exact digest does not select that kernel.
These checks verify target restoration at saved rows; repeated bridge runs and
their error/convergence diagnostics are a separate benchmark.

Any later vectorization or threading generation requires its own explicit
source-digest pairs and native gates. It must not be admitted by a generic
"performance change" exception.
