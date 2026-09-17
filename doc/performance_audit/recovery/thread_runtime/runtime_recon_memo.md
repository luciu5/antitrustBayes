# Thread runtime reconnaissance

This is a bounded source audit for the threading implementation contract in
`doc/performance_audit/threaded_likelihood_spec.md`. It does not compile a
model or claim that current sampling uses within-chain threads. The only
runtime operation still needed is the separate validation script in
`validate_runtime.R`, run in a fresh R process with an already compiled model.

## Findings confirmed from the installed stack

The installed stack is R 4.6.1, rstan 2.32.7 (Stan 2.32.2), StanHeaders
2.32.10, and RcppParallel 5.1.11.2. The relevant installed headers are under
`/home/ctaragin/R/x86_64-pc-linux-gnu-library/4.6/`.

* `rstan::rstan_options("threads_per_chain")` is the supported RStan control.
  Its setter stores the value in RStan's private option environment and also
  writes `STAN_NUM_THREADS`; the base option
  `getOption("rstan.threads_per_chain")` is not the RStan control. The
  installed setter clamps to at least one integer.
* The installed `rstan:::stanc` body emits the threading prologue only when
  `rstan_options("threads_per_chain") > 1`:

  ```cpp
  #ifdef STAN_THREADS
  #ifndef RSTAN_THREADING
  #define RSTAN_THREADING
  #include <stan/math/prim/core/init_threadpool_tbb.hpp>
  auto tbb_init = stan::math::init_threadpool_tbb();
  #endif
  #endif
  ```

  `rstan:::stan_model` sets `STAN_NUM_THREADS` before invoking the Stan
  translator. `StanHeaders:::CxxFlags()` supplies `-DSTAN_THREADS` and the
  RcppParallel TBB include; `StanHeaders:::LdFlags()` supplies the TBB
  libraries. The local build recipe and these flags are documented in
  `StanHeaders/doc/stanmath.Rmd:43-47` and `:309-314`.
* The installed TBB is the legacy interface (`TBB_VERSION_MAJOR 2019` in
  `RcppParallel/include/tbb/tbb_stddef.h`), so the active branch of
  `StanHeaders/include/stan/math/prim/core/init_threadpool_tbb.hpp` is the
  `tbb::task_scheduler_init` implementation at lines 122-161. Its thread
  count comes from `STAN_NUM_THREADS` through `get_num_threads()` at lines
  30-73. The new-TBB branch at lines 77-120 has the same first-initialization
  caveat: if a scheduler has already been instantiated, later initialization
  is ignored (lines 79-83). The legacy branch carries that caveat at lines
  123-128.
* The installed `rstan::sampling,stanmodel` method has no
  `threads_per_chain` formal. When the option is greater than one it sets
  `STAN_NUM_THREADS` (the method body does this before sampler construction),
  then calls the compiled sampler. Its Unix multi-chain path uses
  `parallel::mclapply(..., mc.preschedule = FALSE,
  mc.cores = min(chains, cores))`; each child invokes `rstan::sampling` with
  one chain. This is the concrete fork interaction: child processes inherit
  the option/environment state established by the parent. The project default
  is `cores = 1`, and the bridge context explicitly uses `cores = 1`.
* The installed `libStanServices.a` contains the command-level
  `stan::math::init_threadpool_tbb` call in `stan_fit.o`; its DWARF maps the
  call site to `stan_fit.cpp` lines 397-401. This is binary/source-line
  evidence from the installed runtime, not a shipped textual `stan_fit.cpp`.
  Together with the header above, it establishes that the sampler command
  initializes TBB from the environment. It does not establish that changing
  the option after the first model/DSO initialization can resize an existing
  arena.

The source-level lifecycle is therefore: set the RStan option before model
translation/DSO loading; translation injects the TBB initialization only for a
threaded build; sampling passes the environment into the compiled command; the
first TBB scheduler initialization in a process controls that process. A
thread count must be validated in a fresh R process. Restoring the R option
alone is insufficient because its setter also changes `STAN_NUM_THREADS`.

## Current repository plumbing

`R/estimate.R:569-580` reads `default_stan$threads_per_chain`, falls back to
`getOption("rstan.threads_per_chain", 1L)`, uses it only to derive a default
`grainsize`, and then removes both `threads_per_chain` and `grainsize` from the
list passed to `rstan::sampling` at `R/estimate.R:699-703`. Thus the current
caller value is not an RStan runtime setting. It is also using the wrong option
namespace for the installed RStan control. No implementation change belongs
in this reconnaissance assignment.

`R/helpers.R:379` and `:782` accept `grainsize`; the data lists at
`R/helpers.R:644` and `:1076` carry it. None of the current Stan data blocks
declares or consumes `grainsize` (repository search found no `reduce_sum`,
`reduce_sum_static`, `map_rect`, `STAN_THREADS`, or `grainsize` in
`inst/stan`). For example, the sampled-quality scalar likelihood remains in
`inst/stan/logit_supply.stan:554-614`: the share loop is `:559-566` and the
margin loop is `:586-614`. The CES source has the analogous scalar paths.

The real data/dispatch additions, once the serial vectorized target is gated,
are consequently:

1. For simultaneous and smooth-fringe sampled-quality likelihoods, add a
   consumed observation index vector and a `grainsize` data integer to the
   Stan declaration, and pass both to a reducer. The reducer must receive the
   shared predictors, market-indexed scales, margin mask/index data, rho
   inputs, and degrees of freedom needed to reproduce the current scalar
   terms. A one-thread reducer call remains the numerical reference.
2. For exact ordered-k Stackelberg kernels, use a separate market-level
   grainsize and a market slice that owns all periods and all admissible k
   states. The existing market order/start/count arrays are declared at
   `inst/stan/logit_stackelberg.stan:297-302`; the market mixture ends in one
   `log_sum_exp` at `:716-740`. A task must own that complete mixture, with
   share outside it and evaluated once. Observation grainsize and market
   grainsize are different units.
3. Keep collapsed global and supported rho branches serial as required by
   `threaded_likelihood_spec.md:48-52`. A zero-length margin index must remain
   an exact zero contribution, and one-observation/one-market cases must not
   create an invalid reduction range.

`R/bma.R:473-485` creates the fresh bridge context with one chain and
`cores = 1`; restoration is at `R/bma.R:585-601`, and bridge sampling is at
`:825-861`. Preserve this serial bridge path while the sampling implementation
is being normalized. The fork evidence above is a reason to keep any future
parallel bridge experiment in its own fresh-process gate, not to change the
archived bridge path here.

## State restoration requirements

Any future wrapper that changes the RStan setting must snapshot both
`rstan::rstan_options("threads_per_chain")` and the exact prior state of
`Sys.getenv("STAN_NUM_THREADS", unset = NA_character_)`. An `on.exit` handler
must restore the option and then either unset the environment variable (if it
was absent) or restore its prior string, on both success and error. The TBB
scheduler itself is process-global/static, so this restores caller-visible R
state but cannot undo an already initialized arena; separate processes are
required to compare thread counts. Do not infer the number of TBB workers from
the option name or from `chains * cores`.

The installed `reduce_sum` implementation validates a positive grainsize at
`StanHeaders/include/stan/math/prim/functor/reduce_sum.hpp:198-218`. Its
automatic partitioning is not deterministic (`:123-136`); a static reducer is
needed where the exact ordered-k summation contract requires deterministic
partitioning. This is a later likelihood implementation decision, not a
runtime result from this audit.

## Validation status and next gate

No model was compiled and no benchmark was run in this audit. A parse-only
`rstan::stanc` probe confirms the expected distinction: with the installed
default option of one, a `reduce_sum` source is translated without the
`RSTAN_THREADING`/`init_threadpool_tbb` prologue; setting the RStan option above
one before translation is the compile-time gate. The actual worker count and
first-constructor behavior remain awaiting execution in a clean process.

Run `validate_runtime.R` once per requested thread count with an already
compiled threaded model and its data RDS. It records the RStan option and
environment before/after sampling, parse-only prologue flags, sampling success
or error, and exact state restoration. It deliberately uses one chain and one
core. An instrumented reducer or equivalent TBB probe must be supplied in the
later runtime gate to record actual arena concurrency; a successful RStan fit
alone is not evidence of a particular number of TBB workers. Run the success
and forced-error cases in separate R processes and repeat with a malformed
data fixture to exercise the sampler's own error path.

