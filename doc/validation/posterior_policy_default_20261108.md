# Posterior-state counterfactual default (2026-11-08)

Version 0.3.0.9001 changes policy simulation, not estimation specification 14.
An estimated `BayesFit` now resolves omitted `demand_anchoring` to
`"latent_quality"` in `simulate()`, `bayes_sim()` and
`posterior_policy_models()`. Explicit `"observed"` retains conditional
observed-share calibration. `BayesSpecifiedFit` has supplied primitives and
retains that supplied-state/calibration default; no estimated quality posterior
is invented for it.

## State and economic contract

`posterior_policy_draws()` fixes the posterior parameter index, recovered
conditional quality realization and exact-k role state jointly. Collapsed
quality is drawn from its conditional posterior, including the supported
margin/role conditioning, rather than replaced by its mean or zero. A reused
plan carries the same state through pre-policy/post-policy equilibria and
sequential policies. Relative quality retains the existing outside-share
projection. The original observed-state cost/value-of-funds recovery is held
fixed when demand is replaced; costs are not inverted again to force observed
prices or shares. Model probabilities, likelihoods, priors, kernels and role
support are unchanged.

Simultaneous output Logit/CES and input Logit support Bertrand, Cournot,
monopolistic competition, second-score auction and bargaining. Native exact
core-fringe/Stackelberg paths preserve their supplied plan's role controls.
Logit auctions use their cost-based value index with the appropriate signed
coefficient, not posted-price Logit inversion. Bargaining retains the existing
firm-share equations. Reusable `AntitrustFit@parameters` matches its realized
native model after quality replacement.

The paired `antitrust >= 0.99.43.9000` development dependency repairs plain
Logit-auction sequential promotion: the solved value index, incremental cost
state and active-product map become the next baseline. Input-rate validation
uses that active map rather than treating exited products' NA prices as failed
active equilibria. Nonzero pending outside-cost-shock promotion is explicitly
unsupported; the Bayesian latent auction route has zero outside cost shift.

## Precise remaining limitation

Smooth-cutoff Logit auctions use a custom blended conduct model whose existing
posted-price demand representation cannot supply the native auction cost-based
latent value index. Estimated default and explicit latent requests therefore
fail clearly for this family. Its historical explicit `"observed"` route is
retained; it is never selected as a fallback. This latent family remains
incomplete and requires a separately correct structural representation.
`posterior_policy_models()` also retains its prior rejection of smooth-fringe
model realization. Input CES latent reconstruction remains unsupported.

## Local verification and downstream adoption

Focused regressions cover nonzero sampled/collapsed quality, explicit observed
compatibility, fixed plans and invalid plans, exact-k roles, input Logit,
auction/bargaining economic identities, independent FOCs, no-op effects and
sequential state handoff. The ordinary development-source suite completed with
**8,110 passed assertions, zero failures, zero errors, 22 warnings and 15
skipped tests**. Every skip is an explicitly opt-in Stan execution/sampling
test. The warnings comprise constructor notices (same ownership, ignored Logit
share inputs and cost-shock direction) and five existing solver convergence
notices; they are recorded rather than treated as passing sampling evidence.
The first broad run exposed two obsolete auction support-rejection expectations;
the corrected test retains rejection of unsupported input CES, while the new
family tests exercise supported input Logit and output CES auctions.

The suite ran with `ANTITRUSTBAYES_RUN_STAN_TESTS=false` and
`R_LIBS=/mnt/d/Projects/mergerBayes/.scratch/paper-library`, using
`testthat::test_local(reporter="summary")`. Focused structural auction and
nightly sequential checks also passed. An independent review using only the
installed packages (no source adapters) passed sampled/collapsed Logit/CES
default/observed comparisons, specified-fit defaults, fixed-k/offset plans,
tampering rejection and eight input/output auction sequence paths with
independent probability, value-index, markup, exit and advancement identities.

The existing mergerBayes library workflow refreshed the single shared library
to `antitrustBayes 0.3.0.9001`, `antitrust 0.99.43.9000`,
`coordination 0.1.0` and `trade 0.9.0`. CmdStanR is `0.9.0.9002`
(revision `55342f4f6c9300f9a47046f107b4db9b5e468aa6`), bridgesampling is
`1.2.1`, testthat is `3.3.2`, and execution used R `4.6.1`.
The antitrustBayes checkout remains based on
`8afbd535d891a4385e69a30d6482bd15528378db` and antitrust on
`ed8be018714e12ca8b3b7de1b7e86d213c72bc06`, both with the uncommitted changes
described here. mergerBayes `renv.lock` records those base revisions and the
actual dirty-source `PaperSourceHash` values; a base revision alone cannot
recreate these local changes. No inference specification or saved-fit seal
was rewritten, and no changes were pushed.
The local suite log, structured test results and library-refresh log are
preserved in `.scratch/posterior_default_validation_20261108/`.

The first actual bounded banking pilot exposed an additional extraction bug
that full-draw fixtures missed. Under `quality_scale_parameterization =
"total_fraction"`, selective extraction correctly obtained the economic quality
SD but normalization dropped its explicit `sigma_b_strat_economic` alias.
Nested conditional recovery then interpreted the historical canonical name as
a raw total scale and rejected its missing fraction parameter. Normalization
now retains the economic alias. A new selective-adapter regression checks
idempotence, raw-scale equivalence and independent conditional Gaussian
moments; it adds 14 assertions to the earlier 8,096-assertion run.
The latest suite evidence uses the `posterior-default-suite-corrected.*` files.
The actual sealed MC fit also passed development-source recovery and no-op
replay: three valid draws, nonzero quality (SD 1.553579), exactly zero price
effects and maximum independent markdown FOC residual `4.440892e-16`.
Its saved fit checksum did not change. No likelihood, kernel, prior, fitting
or bridge code changed for this correction. The installed library will be
refreshed after the active snapshot workers finish, followed by a distinct
no-inference policy replay; the original pilot records remain preserved.

The mergerBayes banking call explicitly selects latent quality and a fixed
model-stable draw plan. Whole-market and merging-party Paasche indices use
each draw's solved pre-policy rates in their denominator. The available banking
selection contains compact CSV summaries but lacks the original fit/evidence
artifacts and seals; the historical referenced external fit root is absent.
No banking policy replay or numeric manuscript adoption can be claimed from
those summaries alone. The existing observed-anchor results remain historical
while corrected primary results await restored compatible saved inference.

The separate residual-prior sensitivity is opt-in in mergerBayes; it does not
change the package's global prior. A bounded local pilot was subsequently
authorized, with no full production rerun or cluster jobs.
Its candidate is `log(sigma_share) ~ Normal(log(0.2), 1)` versus the current
`Normal(log(0.5), 0.5)` prior. Their probabilities of a residual scale below
0.04 are 5.376031% and 0.0000219226%, respectively. This prior-mass calculation
motivates a sensitivity experiment; it does not establish improved posterior
share or margin predictions and does not adopt a new default.
