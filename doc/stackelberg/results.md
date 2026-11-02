# Diagnostic record

This file is updated as gates execute. A mathematical result, a numerical
check and an estimated posterior are different kinds of evidence.

## Analytically established

The audit derives all four firm-level Logit/CES price/quantity Stackelberg
FOCs and follower response matrices. All provisional Logit formulas survive
on the stated interior domain. Input Logit requires positive markdowns and
marginal values p+mu. CES uses unconditional revenue shares and fixed nominal
expenditure; common CES margins are Lerner ratios rather than dollars.

The naive smooth weighted-sector formulas are rejected. A public random
boundary over nested hard games supplies a coherent differentiable likelihood,
with a shared partition per market-year. The direct ordered-count prior also
marginalizes exact hard games, sharing one count across all years when
market-scoped; it does not introduce a sampled discrete parameter. The all-leader endpoint nests ordinary simultaneous conduct. The all-follower
endpoint is retained only as a deterministic nesting diagnostic; Bayesian
Stackelberg support is now restricted to positive k. The public-boundary
lottery is conditioned on a nonempty leader core with its parameter-dependent
normalizer inside the likelihood. Global uniqueness and
point identification of timing/cutoff are not established.

## Numerically verified before production

The independent oracle passes 42 main demand/conduct/role/sign configurations
plus 26 boundary cases, and direct leader and follower own-profit optimization.
Maximum main-case follower FOC error is 4.5e-16; reaction finite-difference
error 1.1e-8; implicit-response error 5.3e-11; reduced leader gradient 6.3e-9.
The largest direct-optimization action discrepancy is 9.6e-8. Own-profit
Hessians are strictly negative definite on these fixtures. A fresh independent
CES inversion/derivative review also passes. These are local numerical checks,
not global theorems.

Thirty-six independent synthetic panel DGP configurations pass normalized
demand round trips with fixed firm qualities, positive structural costs and
valid margins. Three seeds and both multiproduct role layouts are retained.
This validates the data generator; it does not establish posterior recovery.
An additional 72 configurations satisfy the inherited single-market Bayesian
quality centering/share normalization, including cutoff and fixed-count timing.
A fresh independent algebraic review reproduces Logit output/input and CES
share normalization within 1.8e-15 and confirms the CES price-unit caveat.
Their outside truths are derived before fitting; the general .30-outside DGP
is a separately labeled statistical misspecification stress case.

Legacy coordination tests pass against antitrust 0.99.39 at commit 3ab5f37.
Legacy Bayesian tests pass with ANTITRUSTBAYES_RUN_STAN_TESTS=true and
NOT_CRAN=true, including all pooling tests. The runner reports twelve existing
warnings and no skipped tests. Twenty-seven legacy specified-demand/conduct/
cutoff/sign configurations have saved baseline and three intervention outputs
for an executable after-change comparison. That comparison is still pending.

The new Logit implementation passes the lead's independent 28-case comparison
against the signed primitive oracle: equilibrium prices, quantities, costs,
margins, analytic and generic follower responses, finite differences of fully
solved followers, generic equilibrium prices and reduced leader gradients.
Maximum implicit-response discrepancy is 5.1e-11, follower-response finite
difference discrepancy 1.8e-11, generic price discrepancy 4.6e-12 and reduced
leader gradient 3.3e-11. The committed implementation also passes the complete
coordination suite: 294 assertions, zero failures/skips and the same 24 warnings
as the legacy baseline. The reviewed Logit change is integrated locally as
coordination commit 32c9652. CES followed its approved contract and the final combined deterministic
gate is now complete; Bayesian production work has started only afterward.

Twenty-four additional off-equilibrium tests check actual nonzero raw Nash and
reduced leader profit gradients at market size 1,000. They pass with maximum
absolute error 1.6e-9, verifying output/input signs, actual quantity derivative
units and explicit candidate-state diagnostics. The lead caught and required
correction of the input sign and extraneous Cournot market-size factor before
this check. Zero-at-equilibrium comparisons alone would not detect either.

Twenty-four additional Logit endpoint comparisons use the unchanged antitrust
Logit/LogitCournot classes at identical fixed demand and costs. All-leader and
all-follower timing partitions agree with ordinary simultaneous conduct for
single- and multiproduct ownership, output and input markets. Maximum margin
and cost differences are 4.5e-16; ordinary solver price differences are below
8.5e-9. This checks nesting through an independent production code path.

The independent ordered-boundary likelihood reference also verifies the
existing correlated Student-t conditional margin equation by dividing a
direct bivariate density by its share marginal. It covers negative/zero/
positive correlation, masked rows, no-margin prior-only cells, ties and extreme
cutoff tails. Component probabilities, whole-cell likelihoods and sequential
accounting factors agree to 3.6e-15. The subsequent native compiled gate is
recorded below.
The direct-k two-period reference also matches whole-market logsumexp and
posterior component probabilities, while differing from separately marginalized
periods. It establishes the algebraic unit of marginalization numerically.

## Final deterministic gate, before Bayesian production

The reviewed Logit/CES implementation is integrated at coordination commit
b21b1ec, API version 0.0.0.9003. Its complete suite passes 383 assertions,
zero failures/errors/skips, and the same 24 warnings as baseline. Source build,
installation and R CMD check complete with Status: OK. Roxygen regeneration
is byte-identical on its second run. PriceLeadership/BLP economics are unchanged.

The lead reran independent checks against the final immutable worker commit:
42 baseline derivative/equilibrium cases; 36 nonzero off-equilibrium raw-gradient
cases; 36 calibration cases with noisy/partial/single margins; 36 ordinary
simultaneous endpoint comparisons; and 72 ownership/cost/combined interventions.
Every existing tight gate passes. Final CSVs and complete logs use the
coordination_final_ prefix. The policy check independently solves the raw
follower Nash game and differentiates each leader's own profit, including two
simultaneous leaders and both multiproduct roles. Maximum normalized leader
finite-difference gradient is 2.8e-11 and reaction discrepancy is 1.2e-10.

Tightening CES follower solving and scalar inverse precision reduced its
baseline fully solved reaction finite-difference discrepancy from 1.15e-7 to
4.1e-11. Generic versus analytic equilibrium prices agree within 2.4e-12
across all 42 fixtures. These are local numerical verifications; global
existence/uniqueness and identification remain unproved. Calibrated constructors
now require successful perturbed-start baseline reproduction. Both raw and
scale-normalized diagnostic gradients are returned. Optional nested implicit
checks may be skipped explicitly during posterior draws, with not-run status;
mandatory analytic equilibrium and baseline reproduction gates remain active.

## Independent-review numerical repair

Fresh randomized review found a numerical domain defect in the initial Logit
price root: log(p/target(p)) rejected negative intermediate markup-implied
targets, although negative Logit cost primitives can support an interior
positive-price baseline. The corrected root uses 1-target(p)/p. This preserves
the exact FOC and the positive final price requirement. Bertrand followers now
solve their profit gradient divided algebraically by sigma*q_j, retaining raw
gradients as diagnostics; this prevents tiny shares from masking a false root
at very high follower prices. The final normalized follower residual gate is
2e-10. No economic formula was changed.

Coordination commit b21b1ec, API version 0.0.0.9003, includes both corrections
and strict logical output validation
and regressions. The full suite and all five independent primitive validation
scripts pass again; its built source package passes R CMD check with zero
errors, warnings and notes. Logs and derivative/calibration/policy tables use
the coordination_solver_repair_ prefix. Testing the installed source package
recovers all 120 randomized baseline fixtures at the unchanged 2e-10 normalized
equilibrium gate.

One of the moderate price-leader deviations remains a recorded failure: the
unconstrained follower optimum has price -0.09747337, with normalized FOC
2.64e-14. It is outside the implemented positive-price Bertrand action domain,
not a positive interior equilibrium to force the solver to recover. Corner
price solutions and negative-price Bertrand simulations remain unsupported.
The original failure is retained in the randomized CSV and full log.

## Native Bayesian density gate

Specification 7 retains every prior PDF constant and positive-support
normalizer in all four Stan programs. The user waived compatibility of old
targets; earlier fits require refitting. There is no external evidence offset.
The physical CES outside-revenue conversion is corrected across structures,
and inactive share caps are removed from simultaneous equations.

The rebuilt installed programs pass 180 full normalized-prior/likelihood/
constraint-Jacobian checks over hierarchical flags, plus 72 comparisons with
archived version-6 targets. Logit target shifts equal the derived prior
constants and its gradients remain unchanged. CES target shifts equal those
constants plus the separately measured outside-conversion likelihood change.
Errors are below 2.0e-13. These are fixed-density checks, not posterior fits.

The 96-case independent whole-unit Stackelberg oracle also passes, including
missing margins, share-only/margin-only modes, ties, pooled panels, extreme
cutoff tails and both prior options. Twenty automatic-gradient comparisons
with Richardson finite differences pass the existing 2e-7 relative gate
(maximum 4.1e-8). All likelihood/prediction/role/density discrepancies pass the
existing absolute 1e-9 gate. A raw-quantity audit input fixes a one-ulp revenue
normalization-path mismatch amplified at tau=1e-7; the tolerance was unchanged.
The independent reference still reconstructs unique firms from product rows.

A single-period compiled probe caught missing length-one array dimensions in
new Stan data. The production data repair is installed. Twenty-four independent
compiled endpoint comparisons now pass: all-leader Stackelberg matches the
ordinary simultaneous program and independently computed physical-share FOCs
for both conducts, Logit output/input, CES output, and both ownership layouts.
The fixtures include dominant firms and CES outside revenue shares below 3e-9.
Maximum inverse-margin error is 1.8e-15; proper-target nesting error is
1.14e-13. The fresh review also passes 6,028
native source/math assertions and 1,020 independent primitive assertions;
its corrected report excludes the earlier tautological endpoint comparisons.
The first updated Bayesian suite passes 885 assertions, zero failures/errors/
skips and twelve existing warnings, after its eighteen old outside-conversion
expectations were replaced by the independently derived expenditure-accounting
ratio at the unchanged 1e-10 gate. The failed initial log is retained. Thirty-six
saved ordinary/fringe specified structural comparisons reproduce their earlier
output within 1e-10, including three intervention types. This comparison is
useful economic regression evidence despite the waived API/target compatibility
requirement. The final installed suite passes 901 assertions, zero failures,
errors or skips, and twelve expected existing warnings (913 reported outcomes).
The three initial comparison failures caused by new array attributes were
repaired only in the independent test inputs; all numerical tolerances remain
unchanged. The failed run and 424-assertion focused rerun are retained.

Twelve compiled one-firm fixtures pass after preserving length-one arrays in
common data vectors. Twenty-four two-draw native fixtures pass posterior
summary, reconstruction and combined ownership/cost simulations: all 96
market-year draws are valid, baseline reproduction is below 2e-7 and absolute
normalized residuals below 2e-10. These are specified Fixed_param draws, not
posterior recovery. A fresh read-only R review passes 63 controls with no
defects. An initial summary finding was withdrawn after its independent oracle
indexing was corrected; explicit draw/component loops match means and weighted
medians within 1e-12 for both demands. See
native_prior_flags.log, native_target_comparison.log,
native_compiled_nesting.log, compiled_stackelberg_validation.csv and the
executable research scripts.

The final source-package check completes with Status: OK, zero errors,
warnings or notes. Its full test output reports PASS 901, FAIL 0, SKIP 0,
WARN 12 (expected existing simulation warnings). Installation replay uses the
unchanged audited native code: eight canonical source fingerprints agree,
the final R/database update retains the identical native DLL checksum, and
vignette rebuild checks pass. Coordination also passes a fresh source build
and check at b21b1ec. Runner review is HOLD pending checkpoint provenance,
partial-draw reporting, paired welfare failures and producer-surplus repair.

## Empirically estimated

No new-model synthetic recovery or PNB posterior has been run yet. The protocol
and priors were fixed before fitting in recovery_spec.md.

## Still unresolved

Production implementation and numerical integration gates are complete; the
source-package checks pass; empirical gates remain in progress.
Global leader-equilibrium existence/uniqueness, real-data identification,
posterior-partition versus threshold mapping sensitivity and repeated-seed recovery remain research
questions. The Bayesian demand-quality grouping retains the existing
firm/group hierarchy; the multiproduct recovery DGP therefore uses equal
normalized quality within each firm. General heterogeneous product-level
qualities are supported in deterministic calibration. The inherited single-market
CES normalization also ties the outside share to price units when its market
intercept is fixed at zero; the mean-value convention is retained and this
restriction is reported. All CES programs now use the physical outside-revenue
conversion without a second epsilon adjustment. A separate Bayesian
quality identifier is not part of this initial extension. CES input markets
and no-outside calibration are outside the initial
supported constructor domain, with explicit rejection. No-outside Logit lacks
an inverse common price level; fixed-expenditure CES does identify that level
but requires a different valuation calibration and boundary optimality audit.
