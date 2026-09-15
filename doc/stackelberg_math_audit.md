# Noncooperative Stackelberg mathematical audit

Research lead: independent derivation, 2026-09-13. The audit concerns
simultaneous **noncooperative** first-stage leaders and a simultaneous Nash
second-stage follower sector. It does not concern coordinated PriceLeadership.
No production equations were edited before executing the research oracle.
**Phase 1 sign-off:** approved for production on 2026-09-13 after the lead
oracles and fresh independent review passed. The architecture memo is also
complete. Approval covers the stated interior model and ordered mixture,
with global uniqueness and empirical identification explicitly unproved.

## 1. Definitions, units, signs, and scope

Products are j, firms are f, and ownership is a disjoint partition of products.
L and F partition firms, not product rows. Each firm controls all its products
at its assigned stage. Other leaders' actions are held fixed in a leader's
optimization. Marginal costs c_j are constant; capacity constraints, entry,
partial ownership, and nonlinear costs are outside these derivations.

For output markets profit is sum_j (p_j-c_j) q_j. For input markets the
primitive is marginal value v_j and profit is sum_j (v_j-p_j) q_j. Reported
input margins are positive markdowns (v-p)/p. They need not be below one.
Let a=|alpha|>0. antitrust stores beta=-a for output Logit, beta=+a for input
Logit. s_j=q_j/M, S_f=sum_owned s_j and s0=1-sum_inside s_j are
**unconditional quantity shares**. Conditional inside shares must first be
multiplied by 1-s0. In particular, the cutoff can use conditional firm size
without using that conditional size in a supply FOC.

Logit: s_j=exp(delta_j+beta*(p_j-p0))/(1+sum_inside exp(...)). The outside
utility normalization is zero at fixed price p0. Mean utilities are therefore
delta_j=log(s_j/s0)-beta*(p_j-p0), exactly antitrust's Logit convention.
For derivations set M=1; quantities are multiplied by M in the API.

CES uses antitrust's fixed nominal expenditure E, elasticity gamma>1, fixed
outside price p0>0, and outside valuation normalized to one:

    D = p0^(1-gamma) + sum_j A_j p_j^(1-gamma)
    r_j = A_j p_j^(1-gamma)/D
    q_j = E r_j/p_j,   r0 = p0^(1-gamma)/D
    R_f = sum_owned r_j,   H = 1+(gamma-1) r0.

Here r and R are unconditional **revenue shares**, not quantity shares.
Conditional quantity shares t_j and outside quantity share s0 imply
r0=s0*p0/[s0*p0+(1-s0)*sum_j t_j*p_j]. Calibration uses
A_j=(r_j/r0)*(p_j/p0)^(gamma-1), which is the existing CES meanval formula.
E=insideSize/(1-r0), where CES insideSize is observed inside revenue.

The legacy Bayesian CES code bounds its outside quantity probability away
from zero/one and applies the same epsilon bounding again after converting it
to revenue share. That second operation differs from the physical identity
above by order 1e-6. New Stackelberg likelihoods retain the existing quantity
hierarchy but use the exact revenue conversion, with no second epsilon shift,
so estimation and deterministic reconstruction use identical physical shares.
The initial audit retained legacy simultaneous/fringe conversions. The user
subsequently waived compatibility; native specifications 7 and 8 use the physical
conversion above in every CES program. A revenue
share that is not representable in the open unit interval is an unsupported
numerical domain, not grounds for clipping a new Stackelberg FOC.

The first supported CES implementation is output-market CES. Existing Bayesian
CES rejects input markets; changing gamma's sign or merely negating a margin
would not establish an input-market inverse supply system. CES input support
must not be inferred from the output equations. Positive outside shares are
required for the initial constructor. In fixed-total-quantity Logit without
an outside good, relative-price normalization does not identify the common
inverse-price level. This objection does not apply mechanically to CES:
fixed nominal expenditure E identifies its inverse-price scale even without
outside expenditure. The initial API nevertheless excludes that CES boundary
because the outside-valuation calibration above is undefined there, and a
single positive-cost firm with all expenditure has no finite interior optimum.

## 2. Logit Bertrand: follower and leader derivations

For output demand J_jk = -a s_j (1{j=k}-s_k). A firm's ordinary price
gradient is s_k[1-a mu_k+a sum_owned mu_j s_j], where mu=p-c.
Subtracting the normalized gradients for any two owned products gives
mu_j=mu_k. Thus the follower dollar markup is

    mu_f = 1/[a(1-S_f)].                                      (LB-F)

This includes the single-product case by S_f=s_j. Holding costs fixed,
every follower product has the same price differential dp_f=dmu_f.
Write T=sum_inside s_j dp_j and d_f=1-S_f+S_f^2. Differentiating the
follower equilibrium gives

    dS_f = -a S_f(dp_f-T)
    dp_f = dS_f/[a(1-S_f)^2] = (S_f/d_f) T.

With T_L=sum_leader_products s_j dp_j, it follows that

    B_F = sum_f_in_F S_f^2/d_f
    phi_F = 1/(1-B_F)
    T = phi_F T_L
    dp_j_follower/dp_k_leader = (S_f/d_f) phi_F s_k.             (LB-R)

Consequently each leader's own reduced gradient is
s_k[1-a mu_k+a phi_F sum_owned mu_j s_j]. Products owned by the same
leader again have common dollar markups:

    mu_l = 1/[a(1-phi_F S_l)].                                 (LB-L)

No other leader's profit is present in this equation. It holds for any
number of simultaneous leaders. Firm composition affects product shares and
individual reaction columns s_k; the sector multiplier only needs firm shares.

Since d_f-S_f=(1-S_f)^2>0, B_F<sum_F S_f. With s0>0,
B_F+S_l<1 and all denominators are positive. If F is empty phi_F=1:
all-leader pricing is ordinary multiproduct Bertrand. If followers vanish,
B_F vanishes quadratically in their shares. With L empty use ordinary Nash
follower equations, not an artificial first stage.

## 3. Logit Cournot: independent inverse-demand derivation

Output inverse demand is p_j=p0+[delta_j-log(s_j/s0)]/a, so
K_jk=dp_j/ds_k=-(1{j=k}/s_j+1/s0)/a. A firm's output quantity
gradient is p_k-c_k-(1+S_f/s0)/a. Therefore

    mu_f = (1+S_f/s0)/a, common across follower products.       (LC-F)

Let t_f=S_f/s0. At a follower equilibrium
s_j=s0 exp(delta_j-a(c_j-p0)-1-t_f), and summing gives
t_f exp(t_f)=sum_owned exp(delta_j-a(c_j-p0)-1). The right side
does not depend on leader quantities. Thus every follower's t_f and every
individual ratio s_j/s0 are constant along its reaction equilibrium.
Let V=s0+S_F and hold other leaders fixed. Adding-up implies

    ds0/ds_k_leader = -s0/V
    ds_j_follower/ds_k_leader = -s_j/V.                        (LC-R)

On the follower manifold a leader's own inverse-demand derivative is
-(1{j=k}/s_j+1/V)/a. Its product quantity gradient therefore gives

    mu_l = [1+S_l/(s0+S_F)]/a.                                (LC-L)

These are dollar markups, common within every firm even if its products have
different prices, costs and qualities. Responses in quantities depend on
individual composition s_j, whereas leader and follower markup multipliers
only need aggregate firm shares. With F empty V=s0 and Cournot is nested;
as S_F approaches zero the same limit obtains. Finite markups require
s0>0 for followers and V>0 for leaders. s0=0 is not a valid interior
full-vector Logit Cournot inverse-demand game.

For input Logit inverse supply has K_jk=+(1{j=k}/s_j+1/s0)/a and the
quantity gradient is v_k-p_k-(1+S_f/s0)/a. Price gradients have the
opposite direct price term and upward-sloping demand. Repeating both
calculations gives exactly (LB-F/L) and (LC-F/L) for **positive markdowns**,
and the same reaction matrices. Recover v=p+mu, rather than c=p-mu.
This result is a sign derivation, not an assumption that inputs are outputs.

## 4. CES Bertrand

CES demand derivative is
J_jk=(q_j/p_k)[-gamma 1{j=k}+(gamma-1)r_k]. Divide a follower price
gradient by r_k/p_k. It is
1-gamma m_k+(gamma-1) sum_owned r_j m_j, where m_j=(p_j-c_j)/p_j.
Owned products therefore have common **Lerner margins**, not dollar markups:

    t_f = gamma-(gamma-1)R_f
    m_f = 1/t_f.                                             (CB-F)

All follower prices satisfy p_j=c_j/(1-m_f); with positive costs their
log-price differentials are identical, du_f=dlog p_j=dm_f/(1-m_f).
Write U=sum_inside r_j dlog p_j. Then
dR_f=-(gamma-1)R_f(du_f-U) and
du_f=(gamma-1)dR_f/[t_f(t_f-1)]. Solving gives

    a_f = (gamma-1)^2 R_f/[t_f(t_f-1)+(gamma-1)^2 R_f]
    B_F_CES = sum_F R_f a_f
    phi = 1/(1-B_F_CES)
    du_f = a_f phi sum_leader_products r_k dlog p_k
    dp_j_follower/dp_k_leader = p_j a_f phi r_k/p_k.            (CB-R)

An algebraically identical implementation avoids cancellation in t_f-1 when
gamma is near one. With h=gamma-1 it is
a_f=h*R_f/[1-R_f+h*(1-R_f+R_f^2)]. Ten thousand randomized comparisons
against the original expression, including gamma-1 down to exp(-10), agree
within 3e-11; the stable expression is preferred in production.

The reduced own leader gradient is
1-gamma m_k+(gamma-1)phi sum_owned r_j m_j. Therefore

    m_l = 1/[gamma-(gamma-1)phi R_l].                          (CB-L)

For gamma>1, r0>0, 0<R_f<1, we have t_f>1 and 0<a_f<1.
B_F_CES<sum_F R_f and phi R_l<1. Every predicted margin is in (0,1)
and finite. Composition appears in product reaction coefficients and
counterfactual demand, although firm shares suffice for these FOC multipliers.
F empty and negligible followers nest ordinary CES Bertrand.

## 5. CES Cournot

Inverting the CES demand Jacobian, rather than changing symbols in Logit,
gives

    K_jk = -p_j/(gamma q_j) 1{j=k}
           - [(gamma-1)/(gamma H)] p_j p_k/E.

Put b=(gamma-1)/(gamma H). A follower's quantity gradient is
p_k-c_k-p_k/gamma-b p_k R_f, hence

    m_f = 1/gamma+b R_f.                                    (CC-F)

The above equals antitrust's simultaneous CESCournot formula after replacing
its alpha=r0/(1-r0). Dollar markups generally differ by product.

Follower p_j=c_j/(1-m_f) implies fixed within-firm product ratios in both
prices and quantities. Write x_f=dlog q_j for its common log-quantity
differential, W=sum_inside r_j dlog q_j, and A_f=1-m_f. Inverse demand gives

    dlog p_j_follower = -x_f/gamma-b W
    dR_f = R_f[(gamma-1)x_f/gamma-b W]
    dr0 = -(gamma-1)r0 W/H
    dH = -(gamma-1)^2 r0 W/H.

Differentiate m_f=1/gamma+bR_f and impose dlog p_f=dm_f/A_f.
Algebraically collecting coefficients gives

    a_f = gamma b [A_f-bR_f+R_f (gamma-1)^2 r0/H^2]
          /[A_f+(gamma-1)bR_f]
    C_F = sum_F R_f a_f
    x_f = -a_f W
    W = sum_leader_products r_k dlog q_k/(1+C_F)
    dq_j_follower/dq_k_leader = -q_j a_f p_k/[E(1+C_F)].        (CC-R)

Thus a leader's effective inverse-demand derivative replaces b by
b/(1+C_F). The noncooperative leader FOC is

    m_l = 1/gamma + b R_l/(1+C_F).                            (CC-L)

This compact closed form is valid on an interior, regular follower
equilibrium branch with positive costs. It is not an unverified elasticity
approximation. Note that a_f need not be positive for every possible market:
do not impose universally negative follower responses. The independent review
prompted the following strengthened result, checked algebraically by the lead.
Put h=gamma-1, x=r0 and y=R_f. A_f=b(H-y)>0 and
A_f+h*b*y=b[H+(h-1)y]>0. The numerator N of a_f+1 after
multiplication by its positive denominator H^2[H+(h-1)y] is affine in y.
At y=0 it is H^2(H+h)>0, and at y=1-x it is
x[h^3+3h^2+2h+1+(h^3+h^2+2h)x+h^2*x^2]>0.
Thus a_f>-1, 1+C_F>1-SUM_F R_f>=r0+SUM_L R_f>0. Consequently
0<m_l<1/gamma+b<=1 on the interior output domain. The homogeneous
follower-response system has within-firm common log changes, then
x_f=-a_f W and (1+C_F)W=0, establishing local follower-Jacobian
nonsingularity. Retain explicit conditioning/admissibility checks numerically,
especially near zero outside revenue. Ordinary Cournot is nested
when F is empty. As followers shrink, C_F tends to zero.

An alternative implementation, used as a diagnostic oracle in all four games,
is R=-solve(dF_F/dz_F,dF_F/dz_L), evaluated at solved followers. For any
leader and product k its reduced derivative is
partial_pi_l/partial_z_k + sum_F partial_pi_l/partial_z_j R_jk.
Follower FOCs must be profit gradients with costs held fixed during
differentiation; differentiating recalibrated costs would give a false zero.

## 6. Existence, optimality, and identification gates

These are interior stationarity results. An invertible follower Jacobian
establishes a local differentiable reaction branch by the implicit-function
theorem. Negative definite own reduced-profit Hessians verify strict local
leader maxima; follower own Hessians also need checking. A small FOC residual
alone is not a global existence or uniqueness theorem for the leader game.
Reject nonfinite/negative actions, invalid outside shares, singular response
systems, inadmissible margins and failed solver residuals. Do not invent
cross-conduct price rankings. Numerical tests below prove neither global
uniqueness for arbitrary markets nor empirical identification of the cutoff.

Calibration moments can be noisy. Logit x_j=h_f/p_j and lambda=1/a give
lambda_hat=sum w_j x_j m_obs_j / sum w_j x_j^2, using all usable margins.
Default product WLS weights are one, matching antitrust's exposed weights.
Costs must be p-h/a (values p+h/a for inputs), not p*(1-m_obs).
CES fits log inverse margins in one dimension gamma>1 with the verified
role-specific formulas, then costs p*(1-m_pred). Demand normalization exactly
fits the observed shares. Fixed-demand reconstruction accepts a posterior
a or gamma directly and again inverts the hard-role FOCs for costs.

## 7. Numerical verification and corrections

Executed independent base-R/numDeriv/nleqslv research oracle:
[math_check.R](stackelberg/math_check.R),
full log (local evidence: `stackelberg/math_check.log`),
[derivative table](stackelberg/derivative_validation.csv), and
[direct optimization table](stackelberg/optimization_validation.csv).
The oracle constructs demand from observed shares, inverts the proposed FOCs
for constant costs, solves followers using raw profit gradients, and compares
rank-one analytic reactions against both solved-follower finite differences
and implicit differentiation of raw product profit FOCs. Leaders optimize
their own reduced profits with other leaders held fixed.

| Model | Maximum reaction FD error | Maximum implicit error | Maximum leader gradient | Maximum direct optimum action error |
|---|---:|---:|---:|---:|
| Logit Bertrand output | 1.16e-11 | 4.26e-12 | 1.08e-12 | 2.65e-10 |
| Logit Bertrand input | 7.72e-12 | 4.48e-12 | 8.66e-13 | 3.10e-10 |
| Logit Cournot output | 2.23e-11 | 5.25e-11 | 3.02e-11 | 8.94e-9 |
| Logit Cournot input | 5.77e-11 | 4.21e-11 | 3.24e-11 | 3.37e-9 |
| CES Bertrand output | 1.05e-11 | 6.14e-12 | 6.64e-13 | 3.03e-9 |
| CES Cournot output | 1.10e-8 | 4.19e-11 | 6.30e-9 | 1.20e-9 |

Each row covers seven scenarios: single-product one leader; two simultaneous
leaders; multiproduct leader; multiproduct follower; all leaders; followers
with shares 1e-4--3e-4; and no leaders. There are 42 hard-game scenarios.
The largest follower FOC residual is 4.45e-16. Initial raw-Foc solving produced
an unacceptable small-follower finite-difference discrepancy; normalizing
solver residuals by quantity/price fixed precision without changing equations
or loosening the audit gates. Direct optimization is bounded and local; it is
independent objective validation, not proof against every distant deviation.
The oracle also optimizes every follower against fixed rival actions and
checks negative definite own-profit Hessians for both stages. The maximum
follower action error is 9.55e-8; all tested own Hessians are strictly negative
(the least negative maximum eigenvalue is -5.87e-5).

**Correction status:** all four provisional Logit formulas are established
algebraically and pass numerical tests in the stated interior domains.
CES formulas above were derived independently from the package convention.
The provisional smooth weighting requires replacement as follows.
An expanded [boundary oracle](stackelberg/boundary_check.R) and
[validation table](stackelberg/boundary_validation.csv) also pass at gamma
1.1, 1.5, 2, 6 and 15 and at follower shares shrinking from 1e-2 to 1e-8.
The largest relative reaction error is below 8.5e-7. Log quantities and
Richardson extrapolation avoid loss of relative precision for tiny responses;
the original acceptance gates were retained.
The [fresh independent review](stackelberg/independent_review.md) found no
substantive error in the hard-role equations or mixture. Its separately coded
[primitive CES inversion](stackelberg/independent_ces_check.R) was executed
by the lead; output (local evidence: `stackelberg/independent_ces_check.log`) shows maximum
reaction discrepancy 1.59e-8 and leader gradient 6.71e-9 at gamma=1.5,2,6,
with six products and both multiproduct leader and follower firms.

An additional [dominant-follower Cournot check](stackelberg/ces_negative_response_validation.R)
uses revenue shares (.03,.82,.08,.04) and outside share .03 at gamma=1.1,
1.5,2,6. The dominant follower's quantity response is positive in the first
three cases: its a_f is negative. The [table](stackelberg/ces_negative_response_validation.csv)
shows maximum solved-follower FD discrepancy 2.2e-8, IFT discrepancy 1.1e-8,
leader gradient 1.5e-10 and direct local optimum relative action error 5.7e-8.
All tested reduced leader curvatures are negative. This confirms the sign
warning above without asserting an economic ranking or global uniqueness.

## 8. Smooth boundary: coherent differentiable relaxation

Replacing F membership by weights in B_F or S_F and interpolating leader and
follower margins does not describe the original hard-role game. A firm with
fractional weight then contributes to a follower response system that also
enters its own purported leader derivative. No timing, payoff or equilibrium
system has been specified that makes that self-response valid. Even a
one-firm output Bertrand market, where both hard timing assignments have
identical monopoly FOCs, obtains a spurious different leader component from
B=w*S^2/(1-S+S^2). Nonlinear expected markups also cannot be obtained by
inserting expected membership into a nonlinear denominator. Treating that
construction as established Stackelberg economics is unjustified.

Use a **random common size boundary** instead. Fix the size metric a_f to
observed conditional firm quantity shares for Logit and conditional firm
revenue shares for CES (the existing cutoff convention). For each market-year,
let the publicly realized boundary C=c_m+tau*epsilon, where epsilon has the
standard logistic distribution and tau>0 is prespecified. Firm f leads iff
a_f>C. This gives exactly w_f=P(follows|c)=inv_logit((c-a_f)/tau), with
correlated role assignments respecting firm size ordering. Every realized
component is an actual noncooperative hard-role game. Marginalize these games
in the likelihood; do not solve a fictitious fractional-role game.

Order firm sizes a_1>=...>=a_n. There are only n+1 possible partitions,
with k largest firms leading. Define v_i=inv_logit((a_i-c)/tau).
Their probabilities are pi_0=1-v_1, pi_k=v_k-v_(k+1) for 0<k<n,
and pi_n=v_n. Ties have zero-probability intermediate partitions; preserve
ties instead of imposing an arbitrary economically active tie break.
The margin likelihood for a market-year is

    log_sum_exp_k [log(pi_k)+sum_observed_rows log p(margin_j | share_j,k,theta)].

Share likelihoods and their correlation-based conditional margin adjustment
are reused unchanged; Stackelberg receives no mu_b_fringe quality shift.
The same partition equations must generate predictions and log likelihoods.
Use market-year joint likelihoods (or valid sequential conditional factors)
for predictive comparisons: an arbitrary sum of independent per-row mixture
likelihoods wrongly gives each product a different market partition.
All partitions sum firm statistics once, never once per owned product.
Logistic interval probabilities need stable log-CDF differences in Stan.

This model is smooth in c and theta, coherent at every tau>0, and tends to
the original hard threshold away from ties as tau approaches zero. It adds
one prespecified scale, no extra estimated weak parameter. This is a model of
uncertain realized timing, not fractional firms and not a sampled discrete k.

## 9. Transition scale and posterior-to-simulation mapping

The legacy tau=0.2 means the 10--90% follower-probability transition spans
4.394*tau=0.879 share units. At c=.03 a .21-share firm still has follower
probability .289, and a .001-share firm only .536. This cannot sharply
describe a 1--5 firm strategic core in a 39-bank market. Preserve 0.2 for
legacy fringe calls. Expose cutoff_scale for all structures, defaulting to
0.01 only for explicitly requested Stackelberg. Prespecify sensitivity at
0.0025, 0.005, 0.01 and 0.02; report failures rather than choosing the scale
that gives an attractive PNB result. Cutoff beta-prior hyperparameters remain
application-selectable; no package default prior on k=1--5 is imposed.
Executed [scale_check.R](stackelberg/scale_check.R) records
[probability geometry](stackelberg/scale_geometry.csv) and
[prior sensitivity](stackelberg/scale_prior_sensitivity.csv), without fitting
PNB. At c=.03, hard k is 6; the randomized k mean ranges from 6.55 at
tau=.0025 to 19.25 at tau=.2. An unbounded logistic boundary can realize
C<0 (all leaders): at tau=.01 the 95% realized-k interval reaches 39.
This tail is an explicit modeling implication requiring sensitivity reporting,
not evidence that a .01 scale is empirically identified or optimal.

Report the requested hard posterior summary k_hard(c)=sum_f 1{a_f>c},
its PMF, mean/median/interval and total core share market-year by market-year,
alongside cutoff and role probabilities. Label k_hard as a summary of the
continuous c posterior. Also report the separate posterior realized-partition
probabilities proportional to pi_k times component margin likelihood.
These are different quantities under a finite transition scale.
The size metric is a fixed baseline classification covariate when deriving
unilateral profit deviations. Roles must not be reassigned as a deviating
leader/follower changes its current share; that would define a different
endogenous-timing game and invalidate the stated FOCs.

The coherent simulation mapping samples a realized hard partition using those
posterior component probabilities, recovers costs with that partition's FOCs,
and delegates the deterministic draw to coordination. Each reconstructed
baseline then reproduces observed prices/shares. Compare this with simply
thresholding at c (a sensitivity mapping that discards boundary noise), and
with adjacent partitions for firms near c. There is no justified deterministic
fractional-role analogue to compare as an equilibrium. The mixture is the
coherent soft-role *stochastic* analogue; report both interventions under the
same mapping. The Bayesian bridge defaults to preserving any pre-merger leader
status within the merged owner, with leadersPost explicitly overridable.
coordination requires explicit post roles for changed ownership and does not
decide a merger's timing assignments on the caller's behalf.

## 10. Evidence labels and remaining gates

Analytically established: hard-role FOCs, local IFT reactions, common dollar
Logit markups/common CES Lerner margins, the ordered mixture's probability
model, and the scale calculation. Numerically verified: the 42 synthetic
stationary equilibria and local leader optimization tables above. Empirically
estimated: none at this stage. Still unresolved: global equilibrium uniqueness,
finite-sample timing/cutoff identification and model discrimination, Monte
Carlo recovery, counterfactual mapping sensitivity, and PNB estimation.
Production implementation, full test suites, repeated-seed recovery and PNB
diagnostics are subsequent gates and must not be described as completed here.

Following user steering, [ordered_k_addendum.md](stackelberg/ordered_k_addendum.md)
also specifies a direct prespecified P(k) option using the same exact hard-game
mixture. A panel-wide k_m is marginalized over all years jointly. This option
samples neither discrete k nor continuous cutoff; its actual posterior count
probabilities differ conceptually from the hard-k summary of cutoff draws.
Existing fringe equations remain unchanged. The pre-fit Bayesian single-market
normalization restriction and compatible recovery DGP are documented separately
in [recovery_spec.md](stackelberg/recovery_spec.md).

The fresh deterministic implementation review subsequently found a numerical
issue, not a correction to these FOCs. Logit cost primitives may be negative,
so the target c+mu(p) can be negative during a price-root step even when the
positive-price equilibrium exists. Solve 1-[c+sigma*mu(p)]/p=0 rather than
log(p/[c+sigma*mu(p)])=0. For Bertrand followers the algebraic normalized raw
gradient is 1+beta*(mu_j-sum_owned(mu_i*s_i)); solving this avoids apparent
zeros caused by tiny q_j. Keep raw gradients in diagnostics and check the
normalized equation. The repaired implementation passes 120 randomized
baseline cases, the original independent derivative/calibration/policy gates,
and its package check. One follower deviation has an independently verified
negative unconstrained best-response price; positive-price corner equilibria
are outside the current deterministic solver. See the complete
repair evidence (local evidence: `stackelberg/results.md#independent-review-numerical-repair`).

## Nonempty-core support revision before fitting

The user subsequently restricted new Bayesian Stackelberg estimation to positive k. The [nonempty-core addendum](stackelberg/nonempty_core_addendum.md) derives the conditional public-boundary probabilities and governs estimation/simulation summaries. All deterministic hard-game FOCs above remain unchanged; k=0 remains a mathematical and deterministic nesting diagnostic. No prior or empirical result was tuned in response to a fit.
