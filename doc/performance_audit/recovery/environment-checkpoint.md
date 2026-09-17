# Native recovery environment checkpoint

Recorded 2026-09-17 before isolated source reconstruction and native builds.

## Runtime

- R 4.6.1
- `rstan` 2.32.7
- `StanHeaders` 2.32.10
- `Rcpp` 1.1.1.1.1
- `RcppParallel` 5.1.11.2
- `testthat` 3.3.2
- `posterior` 1.7.0
- installed `antitrust` 0.99.34 (obsolete for this refactor)
- `coordination` not installed
- `Bayes` not installed

## Recovered source bases

- Bayes main: `e67cf31`, with the reviewed working-tree compatibility files
  `R/bma.R`, `R/params.R`, and the associated baseline-consumer tests/docs.
- antitrust refactor: `/mnt/d/Projects/antitrust`, HEAD `684c009`.
- coordination source: `/mnt/d/Projects/coordination`.
- trade source: `/mnt/d/Projects/trade`.
- preserved patches: `doc/statistical_toolkit_runs/20260916/performance-audit/review-patches/`.
  The Helmert candidate is `bayes-helmert-candidate.patch`; the public antitrust
  promotion base is `existing-antitrust-promotion-base.patch`.

The current Bayes working tree remains protected. Candidate and baseline copies
will be built under the ignored
`doc/statistical_toolkit_runs/20260916/performance-audit/recovery/worktrees/`
directory with libraries under its `Rlib/` sibling. Builds use
`MAKEFLAGS=-j1` and one compiler process at a time. Heavy build outputs and
RDS/log files stay under that ignored directory; compact scripts and evidence
are copied here as they become available.
