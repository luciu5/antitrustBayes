# Stan performance audit and measurement plan

This audit covers the four Stan kernels at the performance-audit baseline
(`3ef5ee6`). It is an architecture review and measurement plan; it contains
no sampler run, Stan compilation, or change to the likelihood.

## Findings

- The kernels do not use `reduce_sum`, `map_rect`, OpenMP, or a threaded Stan
  primitive. `RcppParallel` is linked by the package, but it is not called by
  the model code. The existing `threads_per_chain` and `grainsize` path is
  therefore not an active target-execution control.
- The ordinary sampled-quality path declares a dense Helmert matrix
  `Q[N_group, N_group - 1]` and evaluates `Q * b_group_centered`. This is
  quadratic in the number of groups in both storage and multiplication. At
  `N_group = 10,000`, the matrix alone is about 800 MB using ordinary doubles,
  before Stan's copies and autodiff allocations. A large sampled-quality
  probe must be resource bounded; collapsed-quality probes avoid this matrix.
- Simultaneous models have linear observation work after transformed data.
  Collapsed Gaussian shares use structured sufficient statistics and remain
  `O(N + G)` under their supported factorization. Student-t observation terms
  remain `O(N)`.
- Ordered-k models evaluate a complete market-panel mixture. Their work grows
  with the supported `k` states and market-panel observations. Stackelberg
  generated quantities allocate per-market `N_m x K_m` matrices and repeat
  response calculations, so generated-quantity cost should be measured apart
  from model target/gradient cost.
- RStan `cores` controls independent chain processes. The package chain policy
  now defaults to `min(chains, getOption("mc.cores", 1L))`, retains an explicit
  `stan_options$cores`, and does not infer a host CPU budget or alter within
  chain threading. The default remains one process when the option is absent.

## Measurement design

The accompanying `stan_performance_reference.R` creates a frozen manifest for
`N = 100, 1,000, 10,000`, sampled versus collapsed quality where supported,
the four kernels, and separate model versus generated-quantity measurements.
It only performs environment and compiled-model preflight. A future runner
must add valid fixed-parameter target/gradient calls before recording timings;
it must not substitute HMC timing for target/gradient timing.

The first reference pass should use `N_group = N` for the sampled-quality
stress cases and record memory failures as results. Run `N = 100` and `1,000`
for the dense sampled transform. For `N = 10,000`, run the collapsed case and
retain a bounded dense-case failure or memory estimate rather than risking a
machine-wide allocation. The same economic data and parameter values must be
used for model and generated-quantity measurements.

Before considering an `O(G)` Helmert transform, independently compare its
values and autodiff gradients with the dense reference at small and moderate
`G`. The candidate recurrence is an algorithmic optimization, not evidence
for a changed prior or coordinate system. Its target, gradient, and generated
quantities must pass the existing native-equivalence checks before any Stan
source change.

## Interpretation and dependencies

Independent-chain scheduling is measured first because it does not change
target execution. The paired measurements are recorded in
`doc/performance_audit/chain_results.csv`. Profile the model target/gradient and
generated quantities before modifying either native code path. Interpret
chain speedups together with process-tree PSS and convergence diagnostics. Automatic CPU detection is not justified until those measurements
and a memory budget exist. A future `reduce_sum` implementation also requires
an explicit factorization proof for each likelihood mode; collapsed Gaussian
quality cannot be split into rowwise terms merely because its final density is
fast to evaluate.
