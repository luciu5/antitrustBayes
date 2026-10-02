# antitrustBayes test-suite audit (2026-10-01)

The before-state, block-level inventory and dispositions are in
[`test_audit_inventory.csv`](test_audit_inventory.csv). It records all 286
`test_that()` blocks in 49 files, branch `main`, descriptions, measured block
time from one local run, behavior category, potential duplication, migration
classification, and final disposition. Runtime figures below are wall times on
one local host and are approximate, not a benchmark guarantee.

## Finding

The ordinary suite had some development scaffolding, especially a 124-cell
research benchmark's own harness test (4.394 seconds of block time, including
3.686 seconds to smoke every synthetic cell). Most of the remaining large
suite protects distinct economic, posterior, public API, or downstream behavior.
A large percentage reduction would remove useful protection. The more serious
problem is that full statistical validation is not currently an executable
release gate: the old enterprise HMC runner still requires RStan objects while
`bayes_est()` now returns CmdStanR fits.

| Measure | Before | After |
|---|---:|---:|
| Ordinary `testthat` files | 49 | 48 |
| Ordinary `test_that()` blocks | 286 | 272 |
| Routine wall time | 33.507 s | 28.731 s and 32.193 s (two runs) |
| Explicit extended fixture blocks | 0 | 6 |
| Blocks permanently removed | — | 8 |
| Blocks rewritten for CmdStanR test options | — | 1 |

The baseline run had zero failed expectations, zero errors, 13 skips, and 17
warnings. Twelve skips were opt-in Stan sampling tests; the thirteenth was an
untracked saved-fit cache snapshot, now removed. The warnings were expected
messages from the structural-model paths (ignored Logit BLP shares, identical
ownership, and positive cost shocks). The first post-edit ordinary run had 273
blocks, zero failures/errors, 13 skips, and 17 warnings. The final run had
272 blocks, zero failures/errors, 12 opt-in skips, and 17 expected warnings.
The apparent runtime gain varies across single runs; the 124-cell smoke was
3.686 seconds of measured baseline block time. Six extended fixture checks
passed in 6.916 seconds.

## What the suite protects

| Cluster | Best permanent protection | Main weakness |
|---|---|---|
| Conduct and structural bridge | `test-conduct.R` independently checks five Logit and CES FOCs, ownership effects, and antitrust methods; `test-params.R` checks costs and cost shocks. | Estimated-fit full simulation is opt-in; most routine public simulation uses supplied draws. |
| Posterior reconstruction | `test-posterior-quality.R` and `test-market-intercept-collapse.R` compare conditional means/covariances with dense Gaussian algebra; `test-posterior-predictive.R` checks Gaussian and Student-t laws. | No current, automated frequentist coverage or bias experiment for recovered posterior quantities. |
| Exact-k and model evidence | `test-posterior-k.R`, `test-restriction-comparison.R`, `test-model-universe.R`, and `test-bma.R` test component mappings, support, prior/evidence identities, and incomplete evidence. | Full posterior model-selection calibration against known truth is not an active gate. |
| Public API and cross-package integration | `test-simulate-dispatch.R`, `test-stackelberg-r-contract.R`, `test-auction-policy-models.R`, `test-latent-quality-simulation.R`, and `test-real-rate-simulation.R` exercise antitrust and coordination construction and simulation. | Few end-to-end sampled-fit policy draws cross both package boundaries. |
| Stan source and backend | Source contract tests inspect all kernels, `test-market-quality-stan.R` parses Stan when CmdStan is available, and opt-in sampling checks generated quantities. | Source-string tests can pass even when the current sampled posterior is wrong; short HMC checks cannot establish recovery. |
| Benchmark harness | Six manifest/DGP checks now live at `doc/market_intercept_enterprise/fixture_checks.R`. | The historical 124-task runner is incompatible with the CmdStanR backend. |

The strongest tests are the independently calculated FOCs and Gaussian
posterior moments, exact-k finite-mixture evidence identities, and public
`BayesFit`/`BayesSim` integration with antitrust and coordination. They catch
wrong economics that object-construction or source-string tests would miss.

## Changes and deletion red-team

Every proposed deletion was challenged for a plausible material failure before
editing. These are the final dispositions; the CSV specifies each block.

| Before location | Disposition | Evidence and lost-bug challenge | Replacement or remaining guard |
|---|---|---|---|
| `tests/testthat/test-market-intercept-enterprise-harness.R` (7 blocks) | Move 6 manifest/DGP checks to extended; delete 1 worker-path check | A DGP error could invalidate statistical claims, so ordered-k, stress, intercept, and conduct checks remain executable. Because the full runner is archived, manifest uniqueness and 124-cell smoke are also retained in the extended checks. Only worker-path uniqueness is removed: paths are constructed from distinct task IDs and checked by the strict aggregator. | `doc/market_intercept_enterprise/fixture_checks.R`; README command. |
| `tests/testthat/test-classes.R` (first 2 blocks) | Delete as redundant | A broken S4 extension would matter, but public `bayes_specify()` and `simulate()` tests assert `StructuralFit`, `BayesFit`, and `BayesSim` and actually run the methods. | `test-simulate-dispatch.R`, `test-params.R`; stale-fit validation stays. |
| `tests/testthat/test-helpers.R` (2 randomized structure blocks) | Delete as redundant | Malformed Stan data matters; numerical quantity/revenue-share, market-year, prior, model-version, and one-market shape tests already guard the contract more directly. | Remaining helper/economic tests and opt-in Stan compilation. |
| `tests/testthat/test-rho-studentt.R` (`nu -> Inf`) | Delete as low value | The test algebraically compares a formula written in the test with its own limit; it never calls production code. | Joint bivariate density identity and Gaussian/Student-t posterior prediction. |
| `tests/testthat/test-stan-core-policy.R` (nonmutation block) | Delete as redundant | Global option mutation matters, but the explicit-core block already asserts the option is unchanged. | Remaining core-policy and `bayes_est()` forwarding checks. |
| `tests/testthat/test-posterior-predictive.R` (saved-fit byte snapshot) | Delete as low value | A private CmdStanR lazy-cache mutation is not a promised contract. This test always skipped without an untracked RDS fixture and serialized internal state byte for byte. | Adjacent public posterior-prediction test checks unchanged fit and outcome behavior. |
| `tests/testthat/test-s0-information-ladder.R` (opt-in sampling block) | Rewrite | The old `algorithm = "Fixed_param"` option is unsupported by the CmdStanR backend, so the opt-in test failed before sampling. | Supported short HMC options retain the one-market panel outside-share dimension and equality check over every posterior draw. The compiled Logit/CES check passed with a writable CmdStan cache (54.35 seconds after compilation; 6.787 seconds with the model cache warm). |

## Migration and oracle classification

- **A, permanent regression benchmark:** public default/legacy-mode contracts
  (`test-default-model.R`, `test-collapsed-gaussian-api.R`), old ownership input
  equivalent to `Counterfactual` (`test-simulate-dispatch.R`), saved fit scope
  interpretation (`test-market-quality-R.R`), observed-share default
  (`test-latent-quality-simulation.R`), one-year exact-k plan selectors
  (`test-plan-k-key-compatibility.R`), HMC-only bridge/LOO evidence
  (`test-bridge-performance-migration.R`), and dense Gaussian oracles
  (`test-posterior-quality.R`, `test-market-intercept-collapse.R`). These
  represent intentional API compatibility or independent mathematics.
- **B, temporary migration assurance:** the enterprise benchmark's 124-cell
  RStan-era runner is historical validation infrastructure, not a current
  test of the CmdStanR architecture. Its DGP invariants are retained as
  extended fixture checks, while the runner requires a reviewed port.
- **C, legacy-behavior trap:** none found among retained regression benchmarks.
  The removed fixed-parameter sampler option in the opt-in test was obsolete
  test code, not a behavior to freeze.

## Tiers and follow-up

1. **Fast/core:** ordinary `testthat::test_local()` with deterministic
   economic/Bayesian algebra, public specified-fit integration, data contracts,
   and Stan syntax checks. This is 272 blocks in 48 files.
2. **Stan sampling (scheduled/manual CI today):** set
   `ANTITRUSTBAYES_RUN_STAN_TESTS=true` and use a writable
   `options(antitrustBayes.cmdstan_cache = ...)`. Twelve additional compiled
   sampling blocks cover rho/Student-t generated quantities, outside-share
   pooling, and estimated-fit dispatch. The current workflow runs this job on a
   schedule or explicit dispatch, not on pull requests. Short chains are smoke
   checks only.
3. **Extended validation:** run
   `testthat::test_file("doc/market_intercept_enterprise/fixture_checks.R")`
   for manifest and DGP contracts, then restore and review the historical 124-task HMC
   enterprise protocol against CmdStanR before treating its results as current
   statistical recovery evidence. Independent density/gradient validation
   scripts under `doc/` are research checks, not ordinary `testthat` gates.

The enterprise runner currently fails a preflight by looking up
`get("stanmodels", asNamespace("antitrustBayes"))`, which does not exist in
the installed package. Later it calls `rstan::get_sampler_params()` and
`rstan::extract()` on `fit@stanfit`, while the current fit is a CmdStanR fit.
The frozen kernel hashes also refer to the old representation. A safe rebuild
needs a CmdStanR adapter, refreshed hashes, verification of posterior summary
and diagnostic extraction, and a small known-truth pilot before running the
124 tasks. This is an important under-tested area despite the suite's size.

Permanent-test policy: every new test should name the economic, numerical,
public API, or downstream failure it detects. Use independent identities and
small known-truth examples first. Deep-test distinct economic models, and
smoke mechanically similar variants. Label migration checks with an expiry
condition. Put costly recovery and Monte Carlo validation in a scheduled or
release tier with explicit acceptance criteria, not in routine `testthat`.
