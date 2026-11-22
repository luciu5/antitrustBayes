# Relative-share panel validation, 2026-11-22

Nine CmdStan 2.39.0 executables built one at a time on a WSL2 instance
with 7.7 GiB RAM and 2.0 GiB swap. The build used GNU Make with `-j1`,
thread support enabled and precompiled headers disabled. Peak memory below
is the maximum sampled sum of proportional set size (PSS) over the command
and its descendants. The [monitor](../performance_audit/process_memory.py)
polls every 0.25 seconds and can miss brief spikes.

| Stage | Peak process-tree PSS | Result |
| --- | ---: | --- |
| Four generic and five optimized kernel builds | 2.60–3.26 GiB each | All nine completed. The largest was the CES Stackelberg relative-panel kernel. |
| Native panel comparisons and bridge regressions | 0.45 GiB | 687 expectations passed, zero skips or warnings. |
| Ordinary development-source suite | 0.75 GiB | 8,502 expectations passed, 17 opt-in skips, six warnings, zero failures or errors. |
| Four sequential CmdStanR chains on one CES fixture, including fit serialization and diagnostics | 0.44 GiB | 2,000 post-warmup draws, zero divergences and zero treedepth hits. |

The native checks exercised both demand systems, three supply structures and
12 panel designs. At the same initialization they compared every generated
output, normalized log density and gradient against the generic kernel, with
absolute tolerance `1e-7` and explicit checks of nonfinite values. In the
four 600-point ordered-k comparisons the largest log-density difference was
`9.09e-13` and the largest gradient difference was `3.27e-10`.
The first run completed the ordered-k checks and the Logit comparisons, then
stopped preparing the first CES fixture. Supplying its required arguments
allowed the full rerun against the same executables.

The ordinary suite emitted five equilibrium-solver convergence warnings and
one marginal-cost direction warning. The source audit also reproduced all
45 manifest hashes, detected an edited Stan include, and flagged a missing
generated-kernel seal.

For HMC, the CES executable used the sealed Stan data from the
[November 9 demand-shock fixture](residual_model_integration_results.md)
with two markets and three periods. Each chain used the original seed setting,
500 warmup iterations, 500 sampling iterations, target acceptance 0.999 and
maximum treedepth 13. One chain ran at a time on one thread. Across the
13 active parameters, maximum R-hat was 1.00933 and minimum bulk ESS was 708.9.
All 2,000 draws of each active parameter exactly matched the archived fit in chain order.
The fixture has multiple products per quality group and therefore exercises
the general likelihood inside the new executable. Sampling memory for the
rectangular fast path and multi-fit collection remains unmeasured.

The durable local run directory is `doc/statistical_toolkit_runs/panel-20261122/`.
It contains the compiled sources and binaries, per-command memory measurements,
the initial failure log and the rerun logs. The measured peaks identify
compilation as a source of memory pressure. The cause of the earlier freeze
remains unconfirmed.

## Skipping counts with zero prior weight

On the four-market, four-period Logit Cournot Stackelberg fixture, all prior
mass was assigned to three leaders. Each timed invocation evaluated 2,000
identical parameter points using the existing CmdStan toolchain. The three
configurations ran in rotating order. Times are medians over six runs after
one warmup run per configuration.

| Configuration | Declared counts | Seconds |
| --- | --- | ---: |
| Original ordered-panel kernel | 1 through 6 | 0.370 |
| Kernel that skips zero-weight counts in the likelihood | 1 through 6 | 0.226 |
| Original ordered-panel kernel with singleton support | 3 | 0.222 |

The three runs returned identical normalized log densities and every gradient
coordinate. The six-count fit retains all six counts in generated quantities.
The table measures fixed-point evaluations. Local measurements and the
comparison script are in `doc/statistical_toolkit_runs/sparse-ordered-k-20261122/`.

### Four-chain sampling comparison

On the same 112 observations and six declared counts, the original kernel took
51.55 seconds and the kernel that skips zero-weight counts took 23.76 seconds
(2.17× as fast, 54% less time). These are summed wall times for four sequential,
single-threaded CmdStan chains, including warmup, sampling and writing generated
quantities. Both arms used precompiled binaries, the same data, seed 452891,
chain IDs 1–4, initialization at zero, 500 warmup and 500 sampling iterations
per chain, target acceptance 0.99 and maximum treedepth 13. Alternating which
kernel ran first within each chain pair gave the following times.

| Kernel | Chain 1 (s) | Chain 2 (s) | Chain 3 (s) | Chain 4 (s) |
| --- | ---: | ---: | ---: | ---: |
| Original | 13.63 | 14.47 | 12.51 | 10.95 |
| Skips zero-weight counts | 6.40 | 6.96 | 5.44 | 4.96 |

All 2,000 draws of each of the 16 active parameter coordinates and all retained
sampler diagnostics matched exactly between kernels. Both runs had zero
divergences and zero treedepth hits, but maximum R-hat was 1.536 and minimum
bulk ESS was 7.24. The timing measures identical sampler work on this fixture.
Time to reach a well-mixed posterior remains unmeasured. Inputs, per-chain
timings and diagnostics are recorded under `hmc-pair/` in the local run directory
above, with input and executable hashes in `run.log`.
