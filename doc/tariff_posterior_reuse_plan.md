# Tariff counterfactuals using an existing posterior

Status: implementation contract, 2026-09-18. The initial cost-transformation
route and reusable posterior-model API have been implemented and tested in
isolated package installations. Integration and executed evidence are recorded
in `tariff_posterior_reuse_results.md`. Existing fits are reused only when the
equivalence conditions below are satisfied; heterogeneous post tariffs remain
a future extension.

## Mathematical contract

Initially restrict the extension to supported output-market Logit/CES games
with constant marginal production costs: simultaneous Bertrand, Cournot and
MonCom, and exact BF/CF/BS/CS. Auction, bargaining, input markets, quotas,
capacity constraints and smooth fringe require separate audits.
Check the recorded market side explicitly: current bayes_est() defaults Logit
to input markets, so a default Logit fit must not automatically enter this
output-tariff route. PNB deposit-rate fits are not output tariff models.

Use trade's product-vector convention. Consumer price is p_j, the tariff
fraction of that price is tau_j, and seller retention is r_j = 1 - tau_j > 0.
This is not the conventional producer-value tariff rate t_j: in that convention
r_j = 1 / (1 + t_j), so tau_j = t_j / (1 + t_j).

For baseline firm f, require a common retention r_f across its owned products:

\[
\pi_f=\sum_{j\in f}(r_f p_j-c_j)q_j
     =r_f\sum_{j\in f}(p_j-\kappa_j)q_j,
\qquad \kappa_j=c_j/r_f.
\]

Positive multiplication of each firm's profit leaves its optimal choices
unchanged. Followers' best responses and their derivatives are unchanged at
the corresponding effective costs. Each leader's reduced profit is likewise
multiplied by its own constant r_f. Therefore the same noncooperative game,
including simultaneous leaders, can be estimated in effective-cost units.
This argument is independent of the supported demand family.

The margin supplied to the existing likelihood must be the net-revenue margin:

\[
m_j^{net}=(r_f p_j-c_j)/(r_f p_j)=(p_j-\kappa_j)/p_j.
\]

For each saved posterior draw, recover production costs algebraically:

\[
c_j^{pre}=r_j^{pre}\kappa_j^{pre}.
\]

For a proportional production-cost shock d_j and a subsequent tariff, use:

\[
c_j^{post}=c_j^{pre}(1+d_j),\qquad
\kappa_j^{post}=\kappa_j^{pre}(1+d_j)
                 \frac{r_j^{pre}}{r_j^{post}}.
\]

The cost-only counterfactual route is valid only if retention is also common
within every POST-policy firm. It may differ between firms. A merger across
different pre-policy tariff rates is compatible with this route if the
post-policy rates are uniform within the resulting firm; product-level cost
recovery continues to use each product's own baseline retention.

## Assumptions required to avoid re-estimation

| Requirement | Consequence if absent |
|---|---|
| Baseline tariffs are known, finite, constant with respect to firm choices, and have positive retention. | Unknown or endogenous wedges need additional inference or economics. |
| Every baseline firm has one tariff rate across all its products in each observed market/period. | Existing untaxed margin predictions generally cease to be valid. |
| Recorded margins and their error model represent net-revenue margins. | Rescaling recovered costs cannot correct a likelihood fitted to incompatible margins. |
| Demand uses the recorded consumer prices, inclusive of the modeled tariff. | A price-definition mismatch requires correcting the estimation inputs/specification. |
| The existing fit's recovered costs can be interpreted as effective costs. No separate restriction was imposed on physical production costs. | A physical-cost prior, cost-shifter equation or intertemporal physical-cost restriction must be translated consistently. |
| Demand quality, demand parameters, outside-price normalization, market size/budget and role interpretation are unchanged by the tariff. | Entry, quality changes, income feedback or other omitted mechanisms require an explicit extension. |
| Tariff effects enter the supported game through the stated revenue wedge; constant production costs remain the cost primitive. | Other policy mechanisms cannot be represented by the same algebraic recovery. |

Tariffs may differ across baseline firms and periods, provided these
requirements hold in every cell and the fitted model imposes no incompatible
physical-cost linkage across those cells. The existing persistent quality
structure need not change.

If the recorded margin is instead M_j=(p_j-c_j)/p_j, its noiseless conversion
is m_j^{net}=(M_j-tau_j)/(1-tau_j). If it is profit divided by consumer
expenditure, a_j=(r_j p_j-c_j)/p_j, then m_j^{net}=a_j/r_j. These identities
do not establish equivalence of a noisy likelihood on a different reporting
scale. Initial reuse rejects such historical fits unless their actual
likelihood already makes the necessary transformation. There is no automatic
margin conversion after fitting.

No re-estimation also requires a usable posterior under the package's existing
diagnostic policy. Tariff promotion does not repair a failed sampler or
authorize previously withheld BMA results.

When baseline tariffs are zero, production and effective costs coincide and
the margin-basis distinctions above disappear for the stated profit model.
Introducing a new tariff therefore needs no baseline tariff correction.
For nonzero baseline tariffs, checking the actual dataset's accounting
definition is the first application gate, not a metadata formality.

## Package responsibilities and public API

* antitrust owns the common StructuralFit lifecycle and structural cost state.
* coordination owns the demand/conduct game, ownership, firm roles, follower
  responses, leader conditions and equilibrium diagnostics.
* trade owns tariff normalization, production-versus-effective cost bookkeeping,
  tariff policy state, government revenue and the welfare convention.
* antitrustBayes supplies reusable structural/quality/k draw plans and aggregates
  counterfactuals. It acquires no hard dependency on trade and contains no
  tariff equations.

Prefer extending trade::as_trade_fit() with methods for registered native
coordination models, returning TradeFit through the common StructuralFit
lifecycle. Audit the registry/dispatch adaptation before choosing internal
tariff subclasses. Do not build a second public fitted-model lifecycle merely
because the source coordination classes currently differ from AntitrustFit.

Implemented promotion contract in the tested trade extension:

```r
tariff_fit <- trade::as_trade_fit(
  structural_model,
  policy = "tariff",
  tariffPre = tau_pre,
  cost_basis = "effective",
  margin_basis = "net_revenue"
)
trade::simulate(tariff_fit, tariffPost = tau_post)
```

Explicit provenance is mandatory for nonzero baseline tariffs. Promotion must
record the tariff convention, input cost basis, margin basis, baseline
retention, recovered production costs, original effective costs, model
identity, demand normalization and pre/post role sets. A historical fit
without margin provenance requires the caller to declare that interpretation;
unknown provenance must never silently mean net revenue.
Current bayes_market() inputs do not store tariff or cost-basis provenance,
and bayes_est() documents the ordinary output Lerner margin. The new adapter
therefore cannot infer baseline tariff compatibility from existing metadata.
Missing provenance, incompatible margin semantics and unsupported structural
models should be distinct error categories.

Existing trade promotion preserves the source game and treats policy changes
relative to tariffPre. That is not by itself proof that copied cost slots are
physical production costs. Audit compatibility with that existing convention;
the new explicit cost contract must not silently reinterpret old promoted
objects. Keep physical costs as authoritative stored state, and distinguish
temporary effective costs used by the solver. Never compound tariff ratios
through previously adjusted costs when a scenario is repeated or resumed.

## Staged implementation

### 1. Establish input provenance and exact equivalence

Write the short tariff math/architecture specification in trade. Define the
supported consumer-price tariff convention, output-only initial scope, margin
and cost bases, welfare jurisdiction, and typed errors for unsupported cases.
Test independent payoff and derivative identities before production edits.
Check existing source-fit costs and parameter-loading paths for the stated
effective-cost interpretation. Do not alter Stan or recompute a posterior.

### 2. Implement the cost-only route first

Add trade promotion adapters and tariff/cost state for supported coordination
models. Recover production costs from saved effective costs without demand
calibration or an equilibrium solve. Verify the source baseline's FOCs under
the recovered tariff payoff.

For each scenario, derive post effective costs from the original production
cost state and current policy. Check tariff uniformity using the POST ownership
partition. Reuse public coordination counterfactual APIs with the caller's
ownership and leader/core choices. Fail explicitly if the new firm has
heterogeneous tariffs. All-zero and unchanged-tariff policies must reproduce
the corresponding existing structural game.

Compute tariff revenue as sum_j tau_j p_j q_j. Compute producer profit using
net revenue and physical costs, not effective-cost profit. Reuse the audited
consumer-surplus normalization. Total/national welfare must identify whose
producer surplus and government revenue are included; transfers must not be
double counted. Do not automatically rebate tariff revenue into demand.
For nonzero baseline tariffs, equivalence concerns equilibrium choices,
prices and shares, not equality of physical profit to the source object's
effective-cost profit report. An unchanged-policy scenario has zero welfare
change under consistent pre/post accounting, while baseline tariff revenue
can remain nonzero.

### 3. Connect reusable posterior plans and run paired benchmarks

Use existing posterior_policy_draws() plans, including market-specific
conditional quality and exact-k states. Promote each reconstructed structural
draw through trade's public API. Hold its realized demand quality, outside
share, roles and physical costs fixed across comparable scenarios. Apply only
the specified ownership, tariff, cost and explicit role changes. Do not draw
new quality after policy, choose MAP k, or refit for each scenario.

Compare a no-cost merger, tariff change, production-cost shock and combined
intervention on the same draw plans. Retain failures and draw indices. Record
construction/solve/assembly time, solver calls and available evaluations, peak
memory where practical, and price/share/welfare differences. Promotion must
add zero equilibrium solves; unchanged-policy checks use the existing
baseline rather than solving a dummy post state.

### 4. Add heterogeneous POST tariffs only after the initial route passes

If needed, expose narrowly scoped payoff-derivative extension points in
coordination. trade supplies the revenue wedge and its derivatives; the common
game solver handles follower Nash conditions, the implicit response Jacobian,
and simultaneous noncooperative leader conditions. Do not copy Stackelberg
equations into trade or let tariff subclasses inherit inappropriate untaxed
analytic response formulas.

A valid baseline posterior remains reusable when the NEW policy or merger
creates within-firm tariff heterogeneity. This stage changes counterfactual
payoffs, not baseline estimation. Validate against direct optimization and
fully solved follower finite differences before enabling it.

### 5. Defer estimation extensions

Heterogeneous tariffs within baseline firms, incompatible reported margins,
unknown tariffs, and physical-cost restrictions are separate estimation
extensions. Known tariffs would enter structural mean predictions and/or
consistent measurement transformations; no new Bayesian inference engine is
required. Initial implementation reports that a baseline refit or likelihood
extension is needed instead of presenting cost post-processing as a cure.

## Validation gates and separable commits

1. **Math/API specification and independent references:** firm-uniform payoff,
   FOC and response identities for single/multiproduct firms, one/several
   leaders, Logit/CES and differing tariffs across firms. Include a
   heterogeneous-within-firm counterexample. Validate tariff-rate conversion.
2. **trade promotion and policy state:** zero/unchanged tariffs, baseline
   production-cost recovery, zero promotion solve count, immutable original
   structural state, repeated/resumed scenarios, and unsupported provenance
   errors. Add public registry/lifecycle hooks only where needed.
3. **Counterfactual and welfare:** merger, tariff, cost and combined shocks for
   B/C/MC/BF/CF/BS/CS. Match cost-only solutions to independently coded tariff
   payoffs within existing solver tolerances. Reject within-post-firm
   heterogeneity in the first implementation. Check government revenue and
   physical producer profit independently.
4. **Posterior integration and benchmark:** sampled/collapsed quality, multiple
   markets/periods, missing margins, active rho and multiple k states. Preserve
   plan fingerprints, quality/k realizations and posterior indices. Check
   fixed-state Stan/bridge targets and original stanfit draws are unchanged.
   Run paired 200-draw policy fixtures, recording failures and timing rather
   than changing tolerances to obtain agreement.
5. **Optional generic payoff solver:** separate commit/gate for heterogeneous
   POST tariffs, implicit derivatives and leader optimization. Do not block
   the useful firm-uniform route on this larger extension.

Use Luna xhigh workers for bounded implementation and fresh independent
verification after the specification is approved; assign disjoint files or
isolated worktrees and integrate sequentially. Root retains mathematical and
architecture ownership. Run one heavy R/build job at a time, checkpointing
results on the durable workspace because of WSL instability. Commit source,
documentation and compact CSV summaries; retain raw fits/builds locally.

Completion of the initial route means existing VALID fits can evaluate
supported tariff/merger/cost scenarios under the explicit assumptions. It
does not establish that existing datasets use the necessary margin definition,
or that all future tariff models can avoid baseline re-estimation.

## Independent plan review

A read-only Luna xhigh review confirmed the payoff-scaling argument and the
productwise pre/post retention mapping, including mergers across differing
pre-policy tariffs followed by uniform post-policy tariffs. It identified
missing current tariff/cost provenance, the Logit input-market default, native
coordination-to-TradeFit registry/dispatch adaptation, and explicit government
revenue accounting as integration requirements. Those findings are reflected
above. This is mathematical/API review of the plan, not production tariff
validation or a completed implementation benchmark.
