# Performance archaeology (read-only, 2026-09-16)

Scope: immutable baseline `HEAD=12b26c0` plus a separate inspection of
`/tmp/rotation-fresh-build`. No Stan build, install, sampling, or source edit.

## Historical claims

- `a851493` parent (`a851493^`) contains `partial_sum_fn` and one
  `reduce_sum(partial_sum_fn, ..., grainsize, ...)` in both
  `inst/stan/logit_supply.stan` and `ces_supply.stan`. `a851493` replaces
  those with centered logistic-normal share blocks (`clr_shareIn`) and direct
  loops; its four source files contain neither `reduce_sum`, `partial_sum_fn`,
  nor `grainsize`.
- `962a130` starts from that centered-CLR implementation and restores
  `log_shareIn` data and direct per-observation share terms, adding the
  bivariate Student-t share/margin factorization and rho path. This is a
  likelihood/statistical revision, not a threading restoration.

## Current threading path

All four baseline kernels (`logit_supply`, `ces_supply`,
`logit_stackelberg`, `ces_stackelberg`) and all four separate rotation-source
kernels contain zero `reduce_sum`/`partial_sum` occurrences. The same is true
in `/tmp/rotation-fresh-build`. There are no `BUILD_STAN_THREADS` or
`STAN_THREADS` references in Git history, source/configuration, generated
exports, or `src/Makevars*`; Makevars only contains ordinary RcppParallel flags.

`R/estimate.R:504-510` defaults to `chains=4`, `cores=1`. At lines 527-538 it
reads `default_stan$threads_per_chain` (falling back to
`getOption("rstan.threads_per_chain", 1L)`), computes `target_chunks` and a
`grainsize`, then deletes both `threads_per_chain` and `grainsize` before the
`rstan::sampling` call at lines 658-661. `grainsize` is still passed into both
helpers, and helpers still document/store it, but no current Stan data block
has a `grainsize` declaration. Therefore neither option controls target
execution. A caller setting global `rstan_options("threads_per_chain" > 1)`
would set RStan's `STAN_NUM_THREADS`, but there is no threaded Stan primitive
for it to accelerate. `threads_per_chain` in `stan_options` is stripped before
RStan sees it.

The only active package parallel default is simulation `mc.cores`:
`R/simulate.R:85-93` uses `parallel::detectCores(logical=TRUE)`, caps at four,
and forces one on Windows. That is forked counterfactual-draw work and is
separate from RStan chain cores. `RcppParallel` is linked but no current
kernel invokes `map_rect` or `reduce_sum`.

## Factorization units

Ordinary sampled kernels add share and observed-margin terms in direct
observation loops (`logit_supply.stan:489-550`, with the same structure in
CES). Conditional on the sampled global quality vector, the product of those
row terms is the target factorization. The generated `log_lik_period`,
`log_lik_market`, and `log_lik_cluster` fields are accounting aggregates; they
do not change the target's factorization.

Collapsed Gaussian shares call `collapsed_share_gaussian_density` once over all
share rows and group labels (or its sufficient-statistic equivalent), so the
shared quality is integrated jointly over all observations sharing each global
group and the global centered contrast. Active Gaussian rho augments that
joint block with matched observed margins. This is a shared-group/global
multivariate density, not a product of row-level collapsed terms. Missing
margin rows remain in the share block and are omitted from the margin block.

Exact ordered-k kernels evaluate the common share likelihood once outside the
mixture, then for each market initialize each component score with its supplied
prior and add the observed-margin terms for that market before one
`log_sum_exp(score)` (`logit_stackelberg.stan:572-675`; CES is analogous).
Thus k is one complete market-panel state, not an independent per-row or
per-period mixture; exact data construction has one partition unit per market
and fixed ranking across the panel. In collapsed Gaussian rho exact-k, the
joint quality/margin score is evaluated per k only for one market and the
share base is subtracted once (`:635-651`). The transformed-data guard rejects
that combination whenever `N_market > 1` (`:338-341`; CES `:345-348`) because
under global centered quality the market mixtures cannot be multiplied.

## CPU/default conclusion

Current package default `cores=1` for RStan is the safe default. RStan's
`cores` forks separate chain processes; automatic `parallel::detectCores()`
would be a host-CPU answer, not a memory-safe chain budget (and can overreport
WSL/container capacity). This environment reports 12 logical CPUs but only
~7.7 GiB visible memory; no package dependency provides an empirical per-chain
peak-memory model or cgroup-aware available-core helper. A future opt-in
helper could cap by `min(chains, effective CPU quota, floor(memory budget /
measured per-chain peak))`, with one as fallback, but changing the default to
all detected CPUs is not justified by current evidence. Keep simulation's
existing cap of four separately; do not advertise `threads_per_chain` or
`grainsize` as active performance controls until a threaded Stan contract is
actually implemented and tested.

## Verification commands

`git show a851493^:inst/stan/{logit_supply,ces_supply}.stan`,
`git show a851493:...`, and corresponding `962a130` objects; `git grep` over
HEAD; source scans over all four kernels in HEAD and `/tmp/rotation-fresh-build`;
RStan 2.32.7 `sampling(stanmodel)` method inspection; package/config scans;
and read-only CPU/cgroup probes.
