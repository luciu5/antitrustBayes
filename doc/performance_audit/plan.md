> **Historical checkpoint: September 2026 frozen audit plan.**
> For current development and the consolidated agenda, use the [documentation index](../README.md).
> Baseline revisions, execution constraints, and implementation steps below describe the original audit. Retained result reports determine which work was executed; this plan does not select new work.

# Performance and structural-lifecycle audit

Baseline: antitrustBayes main `3ef5ee6`, antitrust refactor `684c009`, and
trade refactor `584f16a`. All three worktrees were clean at audit start.
This work changes execution and structural lifecycle, not economic equations,
likelihood definitions, priors, role supports, or solver/sampler tolerances.

## Ordered implementation units

1. Resolve default chain processes from `min(chains, getOption("mc.cores", 1L))`
   after user options are merged; explicit `stan_options$cores` takes precedence.
   Do not infer within-chain threading from process counts. Benchmark identical
   starts, seeds, and HMC budgets with one and four chain processes.
2. Centralize the posterior-variable contract for simulation, native role
   reconstruction, bargaining, and latent-quality planning. Benchmark extraction
   and complete 200-draw simulation startup before replacing unrestricted
   extraction. Retain draw ordering, exact-k market coupling, and plan integrity.
3. Measure redundant structural-construction solves. Specify a public antitrust
   baseline lifecycle used by antitrustBayes and trade. Parameter loading must
   retain the same demand and cost state; subsequent counterfactual results must
   match. Legacy calibration and policy simulation remain separate operations.
4. Profile target, transformed parameters, and generated quantities before
   selecting a minimal-output/standalone-GQ design. Selected R extraction alone
   does not eliminate native GQ computation or stored sampling outputs.
5. Vectorize valid sampled-quality likelihood sums, validate normalized native
   targets/gradients, then implement observation partial sums with an actual
   threaded build. Measure N=100/1000/10000 and several frozen grainsizes.
6. For exact-k likelihoods, a task owns a complete market/partition mixture.
   Conditional independent margin mixtures may be additive even when a global
   collapsed share likelihood is evaluated separately. Global rho conditioning
   is not split across markets without a valid factorization.
7. Profile the existing O(N+G) collapsed kernel; preserve its centering and
   global quality structure. No dense production covariance or naive
   observation-level parallel density is introduced.
8. Benchmark posterior scheduling across draws, units, and flattened tasks.
   Warm starts may only change initial guesses. Retry canonical starts after
   failures; retain solver counters, tolerance checks, and original row order.
9. Parallel bridge execution is optional and must preserve normalized targets,
   repeat diagnostics, seeds, failure reporting, and memory constraints.
10. Audit canonical model identities and public StructuralFit extension points.
    Do not create an antitrustBayes dependency on trade or private cross-package
    structural shortcuts.

## Frozen chain benchmark

Three previously fitted noisy-share fixtures are reused: simultaneous Logit
Bertrand, exact BF, and exact BS, each with N=24. All use collapsed Gaussian
shares, fixed-df5 Student-t margins, rho off, common residual scales, and the
existing priors. There are two sampler seeds, 271828 and 314159, and one/four
chain-process settings: twelve fits. Starts and data are identical within each
pair. Four chains, 1,000 iterations, 500 warmup, adapt_delta .95, maximum
treedepth 13, and diagonal metric are unchanged. Larger likelihood/threading
benchmarks are separate experiments; these small fixtures do not establish
large-N or many-market scaling.

Normalized target code and all density inputs remain identical. Compare full
unpermuted posterior arrays, sampler diagnostics, leapfrog counts, elapsed
sampling-call time, and summed native warmup/sampling times. Summed chain time
is distinct from wall time. Native sampling time includes GQ work until a
separate measurement establishes its contribution.

## Measurement and storage

Large fits and logs remain under ignored
`doc/statistical_toolkit_runs/20260916/performance-audit/`; compact sources,
manifests, tables, and reports are reviewable Git artifacts. A process-tree
monitor polls `/proc` every .25 seconds and records peak summed RSS and PSS.
PSS accounts for shared pages; summed RSS can double-count forked pages.
Polling can miss short-lived peaks, and memory scans are not simultaneous.
Record this limitation rather than presenting these as exact allocations.

For each unit retain old/new time and memory, fixed-target differences where
applicable, counterfactual price/share/welfare differences, solver work,
diagnostics and failures, and a recommendation. Failed runs are retained and
never relabeled as sampler successes. Settings are not tuned after results.

The session previously could not mutate the actual repository's read-only
Git directory. Implementations can be committed in isolated temporary
worktrees for review, then applied to the requested working branches. Do not
rewrite history or publish remote changes without specific authorization.
