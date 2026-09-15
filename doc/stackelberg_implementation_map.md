# Stackelberg implementation map

2026-09-13. Economics is specified in
[the audit](stackelberg_math_audit.md). This memo precedes production edits.
Mathematical sign-off requires passing the independent research oracle and
resolving substantive independent review findings. The four structural cases
keep base conduct separate from timing/competitive-fringe structure.

## Repository boundaries and working copies

antitrustBayes owns inference; coordination owns deterministic Stackelberg.
antitrust remains an unmodified dependency providing S4 demand, quantities,
elasticities, surplus and normalization. PriceLeadership and PriceLeadershipBLP
retain their coordinated-effects economics and existing dispatch.
antitrust also owns an existing linear/log-linear quantity Stackelberg class
and constructor. Do not redefine that S4 class or its registry entry;
coordination::stackelberg is a separate namespaced constructor for the new
StackelbergLogit/StackelbergCES classes. The reviewed reconnaissance report
is preserved alongside the audit artifacts.

The active workspace allows research-note writes, but its Git metadata is
mounted read-only, including during an authorized branch-creation attempt.
Use local clones in /tmp/stackelberg-coordination and /tmp/stackelberg-bayes,
with isolated worker branches/worktrees. Preserve initial commits 6935013
(coordination) and 023617d (Bayes) as comparison references. Integrate reviewed
commits in those clones and copy the Bayes changes into the active workspace.
Do not mutate the sibling coordination checkout outside the writable roots.
Local commits are deliverable; no remote push is authorized by this memo.

## Deterministic class and API contract

Use two public S4 demand subclasses, StackelbergLogit extends antitrust::Logit
and StackelbergCES extends antitrust::CES, with a common virtual conduct
class/method contract if dispatch benefits. Additional state:

* conduct: bertrand or cournot;
* leadersPre, leadersPost: character firm IDs;
* firmOwnerPre, firmOwnerPost: product-to-firm ID vectors;
* diagnostics: calibration and equilibrium information.

Ownership matrices are still stored for inherited demand/output methods.
The constructor accepts firm ID ownership vectors; reject partial-control
matrices instead of silently interpreting them as disjoint firms. Every owned
product receives its firm's role. A merger may use different leader sets;
changed ownership requires explicitly supplied leadersPost in coordination.
Unchanged ownership/cost-only shocks preserve the existing specified roles.
The Bayesian caller computes its documented leader-union default and passes
it explicitly. No size/cutoff estimation occurs in coordination.

Parallel constructor interface:

    stackelberg(prices, shares, margins, ownerPre, ownerPost=ownerPre,
      leadersPre, leadersPost=NULL, demand=c("logit","ces"),
      conduct=c("bertrand","cournot"), output=TRUE,
      insideSize=1, normIndex=NA, priceOutside=...,
      mcDelta=rep(0,n), subset=rep(TRUE,n), priceStart=prices,
      weights=rep(1,n), alpha=NULL, gamma=NULL,
      control.slopes=list(...), control.equ=list(...), labels=...)

shares are unconditional quantity shares for Logit and unconditional revenue
shares for CES, as in antitrust. They must sum below one in the initial API.
alpha, when supplied, is a positive absolute Logit coefficient; storage uses
the package sign convention. gamma, when supplied, is the structural CES
elasticity. Supplied demand parameters bypass margin calibration, but recover
mean values/costs using the same baseline shares and hard-role FOCs. Calibration
uses all finite positive observed margins with positive weights; partial NA
margin vectors are supported. All-missing calibration without a fixed demand
parameter is rejected. Preserve observed margin provenance when supplying
fixed demand and document any parent-S4 validity handling.

calcSlopes calibrates only the appropriate demand parameter. calcMargins uses
the audit formulas and current firm shares. calcMC returns persistent constant
pre costs or pre costs*(1+mcDelta); it must never reinvert post FOCs as costs.
calcPrices solves the verified leader and follower FOCs. Solving the joint
closed-form price equations is equivalent to solving the nested game at an
interior root, including Cournot represented in prices. It is permitted as an
efficient structural equilibrium algorithm: do not claim prices are Cournot
choice variables. Standalone follower best responses expose actual choices
z=p for Bertrand and z=q for Cournot. The generic Jacobian route supplies an
independent diagnostic and optional equilibrium check.

Expose stackelberg_response (analytic and generic routes),
stackelberg_residuals, stackelberg_followers and stackelberg_simulate.
The latter accepts ownership, costs, subset and leadersPost, updates firm-ID
state, preserves calibrated demand/cost primitives, solves, and returns the
same S4 model with pre/post outputs. Bayes delegates directly to that API.
Do not mislabel a Stackelberg draw as a registered ordinary antitrust model
just to enter its counterfactual registry. Reuse inherited welfare methods
where their demand conventions apply.

Calibration diagnostics include parameter, observed/implied margins, residuals,
weights, maximum follower/leader FOCs, independently started baseline price
reproduction error, identities, solver status and inadmissible branches.
Counterfactual diagnostics retain both ownership and role specifications.
Tariffs are represented as explicitly specified proportional marginal-cost
shocks, with combined ownership/cost changes supported.

## Bayesian structure and data contract

Append structure and cutoff_scale arguments to avoid shifting old positional
arguments. If structure is omitted, use_cutoff=TRUE maps to fringe and FALSE
maps to simultaneous. Explicit Stackelberg supports only Bertrand/Cournot;
reject auctions/MonCom/bargaining. Resolve contradictory explicit settings
with a clear validation error. Retain 0.2 for legacy fringe; Stackelberg's
prespecified default scale is 0.01. Record structure and scale in data,
parameter objects, model metadata and reconstruction. Do not estimate scale.

Firm-market-year data is constructed from group and market-year IDs:
N_firm_mky, firm_mky_shareIn, firm_mky_revShare, firm_mky_market_year,
firm_mky_group, obs_to_firm_mky. Sort unique firms within each cell by the
chosen conditional size metric; retain rank and cell start/count arrays.
Sum each sector statistic exactly once per unique firm. Product observations
retain individual price, demand share and margin likelihood rows. CES converts
outside quantity share to revenue share with observed outside price exactly
as existing code does.

The audit rejects weighted-sector plug-ins. Stan marginalizes an ordered
random-boundary mixture over positive hard partitions per market-year,
conditioned on a nonempty leader core. Direct ordered-k supports a fixed count
prior and common-market or period scope; see both signed addenda. Reuse cutoff
beta prior, alpha/gamma hierarchy, outside-share hierarchy and Student-t joint
share/margin likelihood. Stackelberg has no mu_b_fringe sampled parameter or
quality shift. The user subsequently waived backward compatibility. Native specifications 7 and 8
normalize every active prior natively in all four Stan programs, apply the
physical outside-revenue conversion to every CES structure, and remove
unverified role clamps. The default fringe transition scale remains .2 and
is prespecified and tunable. Old fits require refitting; no external evidence
normalization path remains. Structural meanings remain distinct.
No deterministic demand/equilibrium solver enters Stan or Bayes.

Factor likelihood evaluation into reusable Stan functions for both model and
generated quantities. Sum component conditional-margin densities across a
whole market-year before log_sum_exp. Generated predictions use component
posterior weights; expose those weights and component IDs to reconstruction.
Predictive log densities must respect this joint partition: market-year/market/
cluster aggregate quantities are authoritative. If row factors are exposed,
use a valid sequential conditional decomposition and document order dependence
of individual factors. Preserve bridge-sampling evaluation of the exact model.

## Posterior summaries and deterministic reconstruction

cutoff_summary returns market-year hard k(c) PMF, mean/median/interval, core
share, cutoff posterior and firm leader probabilities. It separately reports
posterior realized mixture k; neither is described as a discrete Stan parameter.
Stackelberg reconstruction uses supplied posterior demand, samples a realized
partition from its conditional component weights, recovers hard-partition costs
through coordination::stackelberg, and delegates the requested intervention.
Provide role_mapping="mixture" (default), "threshold" and fixed/adjacent-role
sensitivity options. Preserve RNG reproducibility; report identities and
mapping in each simulated draw. No fractional-role deterministic equilibrium
has been justified. Near-cutoff sensitivity uses neighboring hard partitions.

## Implementation ownership and diagnostic sequence

1. Lead: math oracle, note, smooth-role design and architecture.
2. Luna A: read-only class/path/test inventory.
3. Luna B: approved Logit API, methods and tests in /tmp/stackelberg-logit.
4. Lead review and executed Logit test gate before CES production coding.
5. Luna C: approved CES equations in a new worktree based on integrated Logit.
6. Lead economic review and executed CES gate; stabilize deterministic API.
7. Luna D: Bayes structure/data/Stan/summaries/reconstruction in its own
   worktree based on stable contracts. Lead retains mathematical review.
8. Fresh Luna E: independent suites, old-model comparisons, derivatives,
   pathological cutoffs, multiproduct aggregation and repeated-seed recovery.
9. Lead: resolve substantive defects, repeated-seed Monte Carlo and mapping
   sensitivity, then eight-model PNB smoke diagnostics with fixed priors.

The full test matrix includes each demand/conduct combination, multiple leader
counts, multiproduct roles, missing margins, baseline reproduction from
perturbed starts, merger/cost/combined shocks, analytic-vs-generic reactions,
independent reduced-profit gradients/optimization, boundary/nesting, ownership
and explicit role changes, old fringe/simultaneous outputs, Stan compilation
and exact likelihood/GQ agreement, proper-prior bridge sampling and valid
posterior simulation. Monte Carlo records multiple seeds and failed recovery;
the PNB gate does not substitute for synthetic identification checks.

No existing PriceLeadership file should change except package metadata needed
to expose additional classes. No mergerBayes application code changes precede
the final package and PNB gates.

## Nonempty-core support revision before fitting

The user subsequently restricted new Bayesian Stackelberg estimation to positive k. The [nonempty-core addendum](stackelberg/nonempty_core_addendum.md) derives the conditional public-boundary probabilities and governs estimation/simulation summaries. All deterministic hard-game FOCs above remain unchanged; k=0 remains a mathematical and deterministic nesting diagnostic. No prior or empirical result was tuned in response to a fit.
