# Nutpie closeout regression checks

Executed 2026-10-20 EDT / 2026-10-21 UTC against production candidate
`e8e16f31cfe2374bc9e8d7ceedf3642a9e6e3ff5`, version `0.2.0.9006`.
Later report-only changes do not alter the tested package source.

| Tier | Passed assertions | Test blocks | Failures/errors | Existing test warnings | Skips |
| --- | ---: | ---: | ---: | ---: | ---: |
| Full development-source suite, optional native dependencies installed | 5,139 | 388 | 0 / 0 | 21 | 12 opt-in HMC |
| Full development-source suite, Nutpie/BridgeStan unavailable | 5,132 | 388 | 0 / 0 | 21 | 13 |
| Existing opt-in CmdStan regression filter | 590 | 28 | 0 / 0 | 9 | 0 |
| Installed source-package tests during `R CMD check` | 4,938 | — | 0 / 0 | 21 | 23 |

The installed package excludes repository-only research runners, accounting for
its lower assertion count and extra skips. The CmdStan run executes the HMC
tests skipped by the ordinary suite. Existing warning details remain in the
raw logs, including deliberately short test-chain diagnostics and structural
fixture warnings; no tests or production sampler defaults were weakened.

The source package built successfully, including vignette creation.
`R CMD check --no-manual` completed with **zero errors, zero warnings and one
pre-existing NOTE**: `.cmdstan_assert_fit()` uses unqualified `tail`.
Examples, help cross-references, installation/loading, tests and vignette
rebuilding passed. Dependency-index probes reported unavailable network access;
the installed dependency checks completed successfully.

The successful build used a clean `git archive` snapshot of the committed
candidate so that large ignored native caches were excluded before directory
traversal. An earlier direct-workspace build was interrupted during traversal;
its log is retained separately and it supplied no validation result.

Durable ignored evidence:

- `.scratch/nutpie-closeout/ordinary-suite.{log,rds}`
- `.scratch/nutpie-closeout/optional-free-suite.{log,rds}`
- `.scratch/nutpie-closeout/cmdstan-regressions/` and
  `.scratch/nutpie-closeout/cmdstan-regressions.log`
- `.scratch/nutpie-closeout/package-check-staged.log` and
  `.scratch/nutpie-closeout/package-check-staged/antitrustBayes.Rcheck/`

The selected runtime uses R `4.6.1`, CmdStan `2.39.0`, `nutpieR` `1.9.1`,
BridgeStan R/client `2.8.0`, `bridgesampling` `1.2.1`, `antitrust` `0.99.42`,
`coordination` `0.0.0.9006` and `trade` `0.8.9.9002`. Exact package revisions,
library paths and all four unchanged Stan source hashes are sealed in the
closeout runtime manifest. Downstream package libraries were read as dependencies;
no `mergerBayes` analysis, fit, specification or production artifact was edited.

The independent review of the committed candidate found no substantive code
blocker. Native statistical/public-API validation is reported separately in
[closeout results](nutpie_closeout_results.md), under the
[frozen protocol](nutpie_closeout_protocol.md). Hosted CI and a clean native
toolchain installation in CI have not been executed locally.
