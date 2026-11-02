# Restoring a saved posterior's native density context

A serialized Stan fit retains draws and sampler diagnostics, but its live C++
model instance may be unavailable after restart. This is a runtime limitation,
not evidence that its distribution is improper. Preserve the original normal
`stanfit` object and evaluate saved draws through the separate `stanfit_model`
argument supported by bridgesampling.

Construct a context from the installed matching kernel and a private copy of
the original Stan data, using one `Fixed_param` iteration and no HMC. The context
initial state is not a posterior sample and never enters the evidence sample.
Reuse the saved posterior sample exactly; retain every likelihood constant,
prior constant, truncation normalizer, and native-coordinate Jacobian.

Kernel identity must come from the saved Stan model name: exact ordered-k
fringe uses the Stackelberg kernel as well. Inferring the kernel from a broad
fringe label can construct a different target and is prohibited.

Require identical model-code hashes, or the exact documented original-to-
rotation hash pair for one of the four approved kernels. The latter migration
adds only `market_s0_parameterization=1` to a private copy of old data, selecting
the algebraically verified legacy path. Other source mismatches must reject.
This adapter does not change posterior coordinates. Nonzero parameter names
and dimensions must match; additional zero-width declarations are harmless.

Before bridge evaluation, compare the restored normalized log target against
the stored `lp__` at fixed, evenly distributed saved rows from each chain.
All rows must be finite and agree within a small absolute tolerance plus
machine-precision accumulation error. Current kernels deliberately use
explicit normalized lpdf contributions, so `lp__` retains the required
constants. This check supplements the source/hash gate; it does not prove
equivalence of arbitrary unknown programs.

The resulting bridge uses the same pair

\[
\{\theta^{(r)}\},\qquad
\ell_{\rm restored}(\theta^{(r)})
=\ell_{\rm saved}(\theta^{(r)}).
\]

There is no refitting, posterior resampling, new prior, or new parameterization.
Preserve the caller's R random-number state while constructing and validating
the temporary context. A restoration or density-check failure remains an
explicit incomplete evidence calculation, never a successful comparison.

Status: root-approved implementation contract. Native reconstruction and bridge
checks remain executed validation gates; initial failed reload and wrong-kernel
reference attempts are retained locally.
