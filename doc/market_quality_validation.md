# Market-quality validation contract

Status: prespecified gates; execution results will be appended without removing
failed or interrupted runs. Existing defaults remain unchanged until all gates
below pass.

## Mathematical and native gates

1. Independently calculate small-N dense covariance densities and posterior
   moments, including all Gaussian constants. Compare to a separate structured
   reference and the production target. Ordinary-scale density tolerance is
   2e-10 in the independent oracle; numerical stress cases are separately
   labeled and never substituted for the ordinary cases.
2. Check one-market global/market equality with identical basis ordering,
   normalized priors, economic parameters, and sampling coordinates. Compare
   native log density, gradients, predictions, bridge target, quality moments,
   and paired policy draws.
3. Check expanded/collapsed equality by Gaussian completion of the square and
   independent low-dimensional quadrature. Expanded and collapsed targets at
   arbitrary beta are not equal: equality requires integrating the expanded
   prior/likelihood over beta.
4. For multiple markets with unequal group counts, repeated group IDs,
   multiproduct/panel observations, and missing margins, compare the production
   rho-active market density against the dense joint covariance oracle.
5. For small M,K, enumerate the Cartesian product only in the oracle and
   compare to market-wise log_sum_exp for BF, CF, BS, CS. Verify each market's
   component probabilities, not just the total likelihood.
6. Compare native automatic gradients with central finite differences for
   demand, outside share, tau, both residual scales, rho, and intercepts. Use
   several step sizes to distinguish truncation error from a defect. Record
   absolute and relative error; ordinary cases require relative error below
   2e-5 with absolute allowance 2e-6 near zero.
7. Run focused tests and the full package suite. Inspect normalization and
   bridge compatibility independently. Run no HMC before substantive findings
   in these gates are resolved.

Edge cases include singleton/two-group blocks, uneven markets, entry/exit,
tau near zero, tiny share scale, rho zero/nonzero, Logit input/output, CES,
sampled/collapsed quality, and all exact conduct families. Production large-N
tests must never use a dense global covariance.

## Frozen HMC and default gates

Compare legacy Student-t/sampled/global, Gaussian/sampled/market, and
Gaussian/collapsed/market on identical data and original prior settings where
comparable. All primary comparisons use active rho. Case-specific data,
initial states, seeds, and sampler budgets are frozen before fitting and kept
in the ignored durable run directory. No fit-specific sampler/ prior tuning
is permitted in the primary comparison. Any extra experiment is separately
labeled sensitivity.

The first benchmark covers simultaneous Logit Bertrand and Cournot, exact BF
and BS, and one CES case, with multiple markets and periods. Report every run's
wall/native time, divergences, treedepth, BFMI, R-hat, bulk/tail ESS and ESS/sec.
Legacy failures remain evidence; they do not make the proposed model pass.

A usable proposed-default fit requires zero divergences and treedepth hits,
all-chain BFMI above .3, maximum R-hat at most 1.01, and minimum bulk/tail
ESS at least 100 across estimated model coordinates and reconstructed economic
parameters. Sampled quality coordinates and nuisance scales are included;
poor nuisance mixing is not hidden by a structural-only summary. A short pilot
that cannot meet these precision requirements cannot authorize a default
change.

Scaling target/gradient tests use M=10,100,500,2000 with fixed firms/periods,
active rho, and an exact-k family. HMC scale runs follow only usable small
benchmarks and are reported as feasible, failed, interrupted, or not run.

PNB first runs B, MC, C, A, then the non-overlapping rich universe with BF,
BS, CF, CS. Gaussian global/market scope is identical for PNB; Gaussian versus
Student-t need not be. Report bridge repetitions/stability, model and k
probabilities, BMA policy distribution, and variance decomposition. Failed
PNB fits or bridge estimates remain visible and prevent default promotion.

Only after density, gradient, package, usable HMC, and PNB gates pass may all
public defaults and help pages change together. Outstanding evidence must be
reported explicitly rather than characterized as validation success.

## Pre-change package baseline (2026-09-17)

The installed specification-9 native kernels passed 136 fixed-target,
normalization, FOC-prediction, and finite-difference checks. Their fixed states,
native targets, gradients, and predictions are retained as the native baseline.
The pre-change unit suite exposed eight stale source-text assertions in
`test-quality-coordinates.R`: they still require a dense Helmert matrix, although
the accepted production implementation already uses the validated suffix
recurrence. These failures are retained. Updating those assertions to the sparse
coordinate contract does not replace or relax the executed density/gradient
gates. The first CSV export of unit results failed on its list-valued result
column; the complete RDS/log and a separate atomic-column summary are retained.

## Source contract migration

The combined private candidate passed 187 source/orthonormal-coordinate checks
(7 tests, zero failures). These include parsing all four scoped Stan programs
and the independently computed Helmert covariance and normalized coordinate
densities. They are not compiled target/gradient or sampler evidence. The
obsolete dense-matrix source assertions were replaced by assertions for the
accepted sparse per-block suffix recurrence. The single-writer native build
subsequently completed successfully; the compiled checks are reported below.

## PNB input provenance

The pulled rich-main comparator tables/settings are available. Its required
`supreme_data.RData` is absent, and none of its ignored fitted objects is present
locally. The bundled PNB data match the older cached fits in group identities
and quantities, but prices and margins differ (maximum absolute differences
2.172417 and 4.901329). Consequently these inputs cannot be silently substituted
for a direct recent-run comparison. PNB preparation's failed original-input
attempt, before any sampling, is retained.

The user subsequently approved the bundled 39-bank data for a separately
labeled application diagnostic. Its prices have the same summary statistics
as the recent rich-main run, but its mean margin is 0.70985 versus 0.47369 in
rich-main; row-level input equality is not established. Bundled data and
application settings are frozen in `pnb-bundled`. This comparison cannot isolate
the Gaussian/Student-t change from the margin-input change. No PNB sampling has
yet been executed on the bundled alternative.

A read-only provenance audit subsequently recovered the tracked historical
input blob from mergerBayes commit `4be1c47`. Rich-main preparation on that
archived snapshot matches the recorded price/share/margin summaries, including
mean markdown .4736931142. The reconstructed data and original application
budget were frozen in a separate `pnb-archived` directory before fitting. This
is an archived-input reconstruction, not proof of the original process's
serialized input. The first Bertrand fit completed with two divergences, so
PNB usability has not passed.

## First compiled native run

The immutable first native run completed 1,199 checks: 1,195 passed and four
failed, all with sigma_share=1e-4. The failed ordinary-double dense comparisons
and short-step finite difference are preserved. A separate 50/80-digit Decimal
dense Cholesky oracle, using the same recorded residuals, covariance and
parameters, gives density error 8.48e-9 (stress tolerance 2e-5). Its tau finite
differences agree with native AD to relative errors 5.37e-10 and 1.88e-10 at the
same two prespecified steps (tolerance 2e-5). The production kernel, likelihood,
priors, parameters and tolerances were unchanged. This resolves roundoff in the
dense-double oracle and subtraction of large log-density carriers; it does not
turn the archived first run into a zero-failure run. Independent Luna review
confirmed the high-precision resolution and found no blocking native economics,
centering, normalization or factorization defects. The additional MC/auction
native cases passed 96/96 checks. An R model-universe fingerprint repair was
subsequently verified, and the independent full-suite
rerun passed 1,855 expectations with no failures or errors (13 skips, 13 warning
expectations). The native .so hash was unchanged. HMC comparisons remain required.

## Fixed-state policy regression

The corrected public-API fixture passed 144/144 checks: single-market realized
quality, pre/post prices and shares agree across scopes for Logit/CES and
Bertrand/Cournot, through mergers, cost shocks and combined policies. All
solves in that fixture succeeded. Two-market BF/CF/BS/CS plans center quality
locally and share one realized k across periods. These are fixed structural
states, not posterior distribution or HMC evidence. The first harness attempt
used unsupported scalar cost scenarios and omitted ranking metadata; its
failures remain retained. Corrected scenarios use the supported per-unit
product cost vectors and the stored predetermined ranking. No economic equations,
solver tolerances, parameters or prior settings were changed.

Retained two-market BF/CF/BS/CS simulations additionally passed 272 checks for
solver success and reuse of the complete quality state. Single-market consumer,
producer and total-welfare comparisons passed 96 checks with exact zero
difference. The first extended harness compared each complete diagnostic quality
vector against a local subset, producing 256 retained test failures; a corrected
replay checks the same saved simulations without rerunning or changing any solve.
The combined reproducible replay passed 368/368 checks.
