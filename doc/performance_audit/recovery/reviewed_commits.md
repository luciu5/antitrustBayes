# Durable reviewed source snapshots

The main working files contain the reviewed changes; actual main's Git HEAD
has not been advanced by this recovery. No remote has been pushed.
Three separable commits are preserved in the ignored local repository
`doc/statistical_toolkit_runs/20260916/performance-audit/recovery/worktrees/bayes-reviewed-integration`,
based on main `e67cf31`:

| Commit | Unit |
|---|---|
| `ff926ad` | Explicit saved-fit kernel migration and coordinate preservation |
| `437e93c` | Consume the public observed-baseline lifecycle when available |
| `b4180fa` | Linear-time sampled-quality contrast transform in four kernels |

Their format patches are in the adjacent ignored `reviewed-patches/` directory.
These are local recovery snapshots, not additional GitHub branches. Tests,
source-generation checks and original/resolved numerical gates are documented
in `native-build-checkpoint.md` and `native_review/`. Heavy fit/library/log
artifacts remain ignored. The lifecycle remains conditional on the installed
antitrust public API; the actual antitrust and trade branches are unchanged.

Parallel bridge dispatch and serial vectorization remain candidates. Their
focused tests and interrupted-run telemetry do not establish a real bridge
speedup or vectorization HMC gain, and neither production dispatch is applied.
