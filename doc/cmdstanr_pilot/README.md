# CmdStanR migration pilot

See [results.md](results.md) for the bounded run's findings and migration
decision.

This is an archived pre-migration compatibility probe. The runner requires an
installed pre-migration `antitrustBayes` package; it stops if the installed
package already uses CmdStanR, so it cannot silently compare CmdStanR against
itself. The frozen 36-row inputs are the Logit Bertrand (`B`) and CES Bertrand
(`CES_B`) synthetic fixtures from the market-quality benchmark, seed 271828.
Both backends use the current `bayes_est()` Stan data and the same four-chain
sampling budget: 300 warmup and 300 retained iterations, two concurrent chains,
`adapt_delta = 0.95`, maximum tree depth 13, diagonal metric, and `init = 0`.
They use the same seed, though the two Stan interfaces need not produce
identical random draws.

Run from the repository root with an installed copy of this package plus
CmdStanR, CmdStan, `posterior`, `loo`, and `bridgesampling`:

```sh
Rscript doc/cmdstanr_pilot/run.R logit /path/to/writable/logit-output
Rscript doc/cmdstanr_pilot/run.R ces /path/to/writable/ces-output
```

Each run writes sampling, compilation, LOO/WAIC, bridge, and fixed-coordinate
target checks to `pilot-results.csv`, with detailed objects and version/hash
provenance beside it. The CmdStan CSV files and fit objects stay in the output
directory. Read the compilation and sampling columns separately: RStan's
package-installed models are already compiled when the script begins, while
CmdStan compiles the copied kernel during the timed run. The RStan timing is
the public `bayes_est()` call and includes its data preparation. The
ESS-per-second columns are labeled with these different timing scopes and
should not be compared as an equal-scope throughput benchmark.

This pilot covers the two simultaneous Bertrand kernels with sampled quality.
It does not establish behavior for the exact-k kernels, archived RStan fits,
other parameterizations, or a package distribution strategy. Bridge estimates
from short pilot chains are a compatibility signal, not production evidence.
