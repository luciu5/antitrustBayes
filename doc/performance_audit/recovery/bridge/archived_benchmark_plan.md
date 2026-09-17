# Frozen archived bridge benchmark (2026-09-17)

No posterior fitting or native compilation. Use the unchanged, archived
`logit_BF_noisy_271828__D.rds` and `logit_BS_noisy_271828__D.rds`: N=24,
four chains, 500 saved posterior draws per chain, Gaussian collapsed shares
and fixed-df Student-t margins. Keep each recorded prior, data and k support.
They use the same baseline prior manifest. BF supports
0:5 and BS supports 1:5; this benchmark does not change either support.

**Plan correction before any evidence estimate:** the recorded BF and BS
DGP datasets differ. The original combined public comparison was correctly
rejected by the same-observations validation, before checkpoint or bridge
sampling. Its frozen runner and plan are preserved under the ignored
`bridge-first-plan/` directory and its failed command log is retained.
Do not change either dataset or relax comparison validation. Version 2 runs
each fit's public comparison separately and sums the timing; it cannot produce
a scientific BF-versus-BS model comparison. In the parallel mode, each fit's
three bridge repetitions are forked separately. Cross-model scheduling is
covered by the focused helper tests, not by this public archived timing.

Run the two complete public comparisons with 5,000 maximum bridge iterations
and three repetitions per model. Serial remains the production default.
Seeds are prespecified as 707101 and 707102. For each seed, run:

1. current main, serial, original caller RNG;
2. candidate serial, bridge_seed=NULL, same caller RNG;
3. candidate serial, explicit seed, preallocated per-job L'Ecuyer streams;
4. candidate two-worker fork, same explicit seed and per-job streams.

Use separate R processes and reverse mode execution order for the second
seed. Check public output/evidence/default caller RNG equality for 1 versus 2;
check every model/repetition logml and iteration count for 3 versus 4. Assign
the supplied seed to BF and supplied seed + 10,000 to BS for explicit streams,
so the separate model calls have independent stream plans.
Retain warnings, failed repetitions and missing evidence. Record process-tree
PSS/RSS with the existing memory monitor and end-to-end public comparison time.
Do not present the two repetitions of this performance design as independent
statistical recovery replications or new HMC evidence. Process-tree memory
includes library loading and all tasks; public elapsed time excludes startup
and source overlays but includes validation, restored contexts and assembly.

No production source changes are authorized by passing a timing benchmark.
Fresh independent review remains required before integrating the candidate.
Checkpoint each complete comparison with source, native DLL, runner, input
fit hashes and fixed settings. Artifacts are under the ignored durable
recovery directory, never /tmp. Do not edit the runner during an active run.

Setup failure retained: the first metadata-only R command omitted the standard
user library from R_LIBS_USER and could not find Rcpp. Explicitly adding the
existing user library restored package loading; no dependencies were changed.
