# Bridge recovery checkpoint

This checkpoint contains the recovered parallel bridge candidate and its
regression tests. The isolated durable worktree is
`doc/statistical_toolkit_runs/20260916/performance-audit/recovery/worktrees/bridge`;
it is based on `e67cf31206e7662044f811e50ba4d26a0c54937a` and has repair commits
`c44dbc52084d48030f7ebf2a4734709279892c78` and
`1bf79e33c045f46616d590657780915b4dabef22`, with the focused parallel test
extension at `50ad896c6f481d69a21a1fa0bf054811c271a27`.

The concrete repair is at `R/bma.R:1123`: context storage uses
`bridge_contexts[i] <- list(bridge_context)`, so a failed `NULL` restoration
keeps its model slot and cannot shift later models. `R/model_universe.R` forwards
`bridge_cores` and `bridge_seed` into the ordinary comparison call. The public
regression covers first, middle, and last failed restorations, every repetition,
model labels, finite successful evidence, positive-prior missing evidence, and
zero-prior missing evidence.

Static R parsing passed for the two source files and the recovery test. The
focused recovery suite now passes with exit status 0 using the durable R-only
candidate package, the rebuilt dependencies, and the unchanged baseline native
DLL. Six displayed worker warnings are intentional synthetic failure warnings;
see `test-bridge-performance-recovery.log` and `bridge-focused-results.md`.

The complete source/test patch is `bridge-recovery.patch` and its SHA-256 is
`9cb9a78add349568fc4e43d560c2f072dd2199b3a3a2e4e9cf31f9beab1dcfe2`.
