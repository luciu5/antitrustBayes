# Plan: integrate Gaussian market intercepts

Status: implemented as an opt-in experiment; sampler promotion is still
pending the frozen HMC comparison. The earlier geometry work and its failed
runs are checkpointed in `b28384a`, following `64fc60c`. The implementation
starts from checkpoint `d75fc9c`. Package defaults remain unchanged.

Implemented scope: all four Logit/CES simultaneous and exact ordered-k Stan
kernels, the public R data/API contract, conditional intercept-and-quality
reconstruction, reusable policy draw plans, normalized market-level generated
likelihoods, bridge migration guards, and specification version 11. A
single-market request remains the historical zero-intercept model and is
recorded as a no-op. The new path requires multiple markets, Gaussian shares
and margins, collapsed market-scoped quality, the joint likelihood, and base
(`legacy`) s0 coordinates.

## Objective and ownership

Test whether analytically integrating the existing Gaussian market intercepts,
in addition to market-specific Gaussian qualities, improves sampling without
changing the statistical model or policy counterfactuals.

Sol xhigh owns the specification, mathematics, acceptance gates, integration,
and final assessment. Luna xhigh workers perform bounded implementation and
independent verification. The bounded implementation and an independent
read-only red-team review are complete. Sol retains responsibility for final
integration and the HMC decision.

The outside-good utility convention, centered quality prior, all scale priors,
outside-share pooling, demand and conduct equations, exact-k support, and policy
definitions remain unchanged. Reference-product normalization is a separate
question and is not part of this experiment.

## 1. Freeze the comparison and approve the mathematics

Record the current source/diff hashes, package binaries, fixtures, initial
economic states, seeds, settings, and output schemas before workers edit files.
Use isolated worktrees from the reviewed checkpoint, with disjoint file
ownership. Do not build workers from the older `64fc60c`: that would omit the
existing experimental changes. No concurrent edits to the same production
files or generated headers.

For market m, conditional on continuous parameters and a conduct state k, write
the share equation as mu + b_m 1 + Z q, plus Gaussian residual noise. Retain

    b_m ~ Normal(mu_b, sigma_market^2)
    q ~ Normal_on_centered_subspace(0, tau^2 P).

After conditioning on observed Gaussian margins, use their normalized marginal
density and a share residual z with the conditional margin adjustment and mu_b
removed. The integrated share covariance is

    V_m(k) = D_m + tau^2 Z_m P_m Z_m' + sigma_market^2 11'.

D contains sigma_share^2 for missing-margin rows and
sigma_share^2 (1-rho^2) for matched rows. Market-specific residual scales, when
selected, enter D without changing their priors. The margin mean depends on k;
the covariance need not depend on k for the present common residual model.
Here sigma_market denotes the random-intercept SD, `sigma_b_market`, not a
share or margin residual scale. Integrate b_m - mu_b, retaining mu_b in the
conditional mean and keeping its original prior.

Derive the normalized density, stable sufficient statistics, and conditional
posterior moments/draw algorithm for (b_m, q_m). Preserve approximately linear
observation/group cost per market-state; dense matrices are small-fixture
oracles only. No global covariance and no enumeration of K^M states.

The single-market package fixes b_market to zero. Keep that path unchanged;
do not introduce a random intercept or new hyperparameters there.

Existing reference evidence, rerun during planning:

* 972/972 dense-versus-structured likelihood/posterior-moment checks passed;
  maximum absolute log-density discrepancy 4.28e-9.
* 16/16 numerical integration checks passed; maximum discrepancy 1.50e-13.
* 16/16 numerical-gradient comparisons passed; maximum scaled discrepancy
  2.15e-6. These are not Stan automatic-gradient tests.

These establish feasibility, not a sampler improvement or complete package
equivalence. Extend them to nonzero intercept means, hyperpriors, s0/Cournot
coordinates, and exact-k factorization before production validation.

## 2. Explicit experimental API and compatibility

Proposed argument:

    market_intercept_parameterization = c("sampled", "collapsed")

Keep sampled as the default. Initially support the new mode for Gaussian
shares, Gaussian margins, collapsed market-scoped quality, and optional rho.
Retain all existing sampled-intercept specifications, including Student-t and
hybrid modes. Reject unsupported combinations for the new mode explicitly;
do not silently alter the likelihood or quality scope.

Use the unchanged base s0 hierarchy for the collapsed path. The old rotated and
joint coordinates contain b_market and cannot be retained as sampling
coordinates after it is integrated out. Require the legacy/base s0 coordinate
mode for this initial opt-in path and explain incompatible requests. This is
a coordinate choice, not a change in s0 pooling or prior. Verify the existing
Cournot alpha/kappa transformation and all its Jacobians separately.

Store requested/effective intercept mode, normalization/version information,
and source hashes in fit metadata, Stan-data attributes, and cache contracts.
Old serialized fits retain their original interpretation. Keep stanfit ordinary
and inspectable. Reconstructed draws and conditional expectations must have
distinct output names; neither may masquerade as a sampled intercept.

## 3. Worker allocation and implementation sequence

| Worker | Bounded responsibility | Owned files |
|---|---|---|
| Luna A | Independent math/oracle and fixed-state tests | New reference and test files only |
| Luna B | Production Stan helper and integration | Four Stan kernels and regenerated headers |
| Luna C | R API, posterior recovery, draw plans, metadata, docs | R helpers/estimate/reconstruction and corresponding help/tests |
| Luna D | Fresh independent red team | Review report and separate verification scripts |

Run at most two substantial worker tasks concurrently, with the lead reserving
review capacity. Only one Stan compilation or HMC job at a time locally. Worker
D must not have authored the implementation. Sol approves A's specification
before B starts; B and C work in parallel only after their data/output contract
is fixed. Integrate after equation-by-equation review.

Start production integration with simultaneous Logit. The same reviewed helper
then enters Logit exact-k and CES ordinary/exact-k kernels; do not copy economic
equations into new solver or likelihood implementations. Complete this plumbing
and its small deterministic tests before expensive sampler comparisons.
The Logit exact-k target lives in logit_stackelberg.stan: update both its target
component scores and its generated posterior component probabilities before
testing R exact-k recovery. Updating probabilities alone is insufficient.

Primary files include R/estimate.R, R/helpers.R, R/posterior_draw_contract.R,
R/posterior_quality.R, R/latent_quality_simulation.R, R/bma.R, and the four
inst/stan supply/stackelberg programs. Keep coordination, trade, and mergerBayes
production code outside this experiment.

## 4. Exact-k, output, and policy recovery contract

Integrate b and q inside each full market/state contribution:

    log L = sum_m log_sum_exp_k(log pi_mk + log L_mk).

Preserve one role state per market across its periods. Common hyperparameters
do not obstruct factorization conditional on their sampled values. Do not
extend this factorization to legacy globally shared qualities by assumption.

For each selected structural posterior draw:

1. Draw k from that market's fully integrated component probabilities.
2. Draw b_m conditional on the market data, parameters, and that k.
3. Draw q_m conditional on the same state and realized b_m.

This sequence generates the joint posterior. Never draw independent marginal
intercepts and qualities, plug in only conditional means, or redraw k between
recovery and simulation. Report quality moments with their conditioning stated;
moments claimed to integrate b must include its variance contribution.

Extend posterior_policy_draws() to retain this single coherent realization and
reuse it for merger, cost/tariff, and combined scenarios, pre and post. Preserve
selective extraction. Audit share_mean, predictions, and posterior_quality():
the former explicit-intercept mean cannot silently acquire a different meaning.
Test contracts across market-years and multiproduct ownership.
Store the reconstructed per-market intercept matrix in the draw plan, validate
its dimensions and market labels, and include it in the plan fingerprint along
with quality and k. A modified intercept must invalidate a cached policy state.

Do not expose the prior center mu_b as if it were a posterior b_market draw.
Use explicit integrated-mean and reconstructed-intercept output names. Define
pred_shareIn, pred_shareUncond, gamma_market and the other derived quantities
through the same coherent recovered state or clearly labeled integrations.
Keep unsupported pointwise log-likelihood/LOO outputs disabled; the normalized
integrated likelihood is a market-level object.

## 5. Mandatory validation before HMC

* Compare the production normalized marginal target with the integral of the
  old target, including all priors and Jacobians. Targets have different
  dimensions: pointwise equality to one old intercept value is not the test.
* Native Stan values against independent dense and quadrature oracles; check
  automatic gradients against finite differences, including hyperparameters,
  alpha, s0, tau, both residual scales, rho, and the intercept scale/mean.
* Test small scales, large positive/negative rho, missing/all-missing margins,
  unequal group counts, repeated periods, entry/exit, and one-group markets.
* Verify the single-market no-op and market-year/market/unpooled s0 structures,
  X_s0, centered/noncentered base hierarchies, and Cournot's coordinates.
* Verify BF/CF/BS/CS market component scores and posterior k against a tiny
  brute-force enumeration oracle. Retain normalized bridge targets.
* Verify joint posterior recovery means, variances, cross-covariances, and
  predictions. Validate policy results at matching reconstructed economic states.
* Run targeted and full package tests. A fresh reviewer checks evidence and
  production diffs for omitted constants, duplicated priors, output-semantic
  changes, missing uncertainty, and unintended economic changes.

Freeze numerical tolerances before seeing candidate results, based on the
existing oracles and conditioning of each test. Resolve substantive findings
before sampling. Save failures rather than relaxing their gates.

## 6. Frozen HMC experiment and decision gate

First compare the current best sampled-intercept implementation (joint s0
coordinates, total/fraction quality scales, margin-first likelihood) with the
collapsed-intercept mode on the same six-market, 176-observation SC Bertrand
cross-section. Retain the same economic priors, data, rho, quality-scale mode,
seeds 19091/29091, and 4 chains x 2000 iterations with 1000 warmup, adapt_delta
.99, and max_treedepth 13. The collapsed path necessarily uses base s0
coordinates; map the common economic initial values and do not invent an
initial value for an integrated variable.

Use an additional sampled-intercept/base-s0 control if needed to distinguish
effects of eliminating b from effects of changing s0 coordinates. Prespecify
that control before candidate outcomes. Archived results supply context;
runtime claims require matched binaries, hardware/load, and comparable runs.

Primary gate for each seed: zero divergences, zero depth hits, all-chain BFMI
at least .30, max rank Rhat at most 1.01, minimum bulk/tail ESS at least 100.
Report runtime and ESS/sec, including recovered economic quantities rather
than only remaining HMC coordinates. Retain failed runs and compare usable
posterior summaries using MCSE. Do not adjust priors or sampler settings to pass.

After this gate, extend to MC, Cournot, auction, representative exact BF/CF/BS/CS
and CES fixtures, and a small synthetic panel. Then try one real full-panel
model, followed by Cournot if warranted. Do not launch the full ladder first.
For usable paired fits, compare bridge evidence with repeated-bridge uncertainty
and policy distributions with Monte Carlo error. Extra conditional recovery
draws do not count as independent HMC information or artificial ESS gains.

Do cheap target/gradient scaling tests before large HMC jobs. Use the user's
Slurm cluster for expensive expansion. Preserve chain-level checkpoints and
resumable run ledgers against WSL interruption; keep large fits local/cluster
side and commit only compact evidence.

## 7. Deliverables and stopping rule

Deliver mathematics/architecture notes, independent oracles, native target and
gradient tables, joint-recovery tests, a two-seed HMC table, policy/evidence
comparisons where fits are usable, a red-team report, and a recommendation.
Separate commits by mathematics/tests, Stan/API, recovery, and benchmark results
when commit authorization is provided. Do not push without target authorization.

Keep this opt-in unless exactness, recovery, and HMC evidence all support wider
use. If performance does not improve materially, retain the evidence and avoid
promoting an implementation that adds complexity without benefit. Do not claim
that integrating intercepts must remove all variance-component or s0 geometry.

## Executed validation status

The independent dense oracle covers ordinary and boundary-scale covariance
states, missing margins, nonzero rho, unequal groups, the one-group case, and
explicit enumeration of all 27 states in a three-market exact-k fixture. Its
largest structured-versus-dense log-density error is
`7.105427e-15`; the exact-k evidence error is zero and the largest state
probability error is `3.330669e-16`.

A native Stan harness compiles the exact production likelihood helpers. Across
five ordinary and boundary cases its largest normalized target error is
`3.86e-10`. This target includes the normalized production hyperpriors for the
integrated intercept, quality scale, residual scales, and rho. Stan automatic
derivatives agree with an independent dense-oracle central difference in all
165 comparisons under the prespecified tolerance; the largest scaled
discrepancy is `5.07e-6` at a deliberately difficult
small-share-noise state.

All four full kernels pass `stanc`, and the complete package compiles and
installs. The intercept-collapse tests and the focused API, market-quality,
model-universe, policy-reconstruction, bridge, and Stackelberg contract tests
have been executed. The complete non-sampling package suite passes after
updating the policy-plan signature and adding an explicit forward-only bridge
migration stage for the new kernels. Historical kernel stages remain accepted
only through their reviewed coordinate policies; collapsed-intercept saved
sources cannot be reinterpreted as sampled-intercept fits. The frozen Supreme
Court HMC comparison has not yet been run, so no claim about improved sampling
or a new default is made.

The independent red team found no economic-target defect. It confirmed that
sampled and integrated coordinates remain one statistical hypothesis in the
model universe, while cache and policy-plan identities retain the
parameterization and reconstructed latent state. It also confirmed that exact
ordered-k reconstruction uses the same realized market state for the
intercept, centered qualities, and policy simulation.
