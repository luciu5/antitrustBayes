# Cournot Logit policy failures: an independent interior-equilibrium oracle

This is a diagnostic for the existing Cournot/fringe and Cournot/Stackelberg
FOCs. It changes neither estimation nor the production solver. It distinguishes
an interior solution with nonpositive prices from a missed positive root.

Let the signed price coefficient be beta, with beta=-a in output markets and
beta=+a in input markets, a>0. Write Logit utility relative to the outside good as

\[
u_j=d_j+\beta(p_j-p_0),\qquad s_j/s_0=\exp(u_j).
\]

The cost primitive c_j is marginal cost in output markets and marginal value in
input markets. The audited inverse FOC in both cases is

\[
p_j=c_j-\frac{h_f}{\beta},\qquad j\in f.
\]

Define A_j=exp(d_j+beta(c_j-p_0)), A_f=sum_{j in f} A_j and x_f=S_f/s_0.
Common dollar markups within each firm imply

\[
x_f=A_f\exp(-h_f),\qquad
s_0=\frac{1}{1+\sum_f x_f}.
\]

For a Cournot strategic-core firm, h_f=1+x_f. Hence

\[
x_f e^{x_f}=A_f/e,\qquad x_f=W(A_f/e).
\]

For a monopolistically competitive fringe firm, h_f=1 and x_f=A_f/e.
The positive real branch of W is unique: log(x)+x is strictly increasing from
minus infinity to plus infinity on x>0. Aggregation includes every owned
product exactly once.

For Cournot-Stackelberg, first calculate each follower's
x_f=W(A_f/e). With D=1+sum_{f in F}x_f, the audited leader FOC is

\[
h_\ell=1+\frac{x_\ell}{D},\qquad
x_\ell=D\,W\!\left(\frac{A_\ell}{eD}\right).
\]

Other leaders are noncooperative. Each leader's equation is unique conditional
on the follower set, and the normalization recovers the equilibrium shares.
This uses the previously derived leader multiplier, not joint-profit
maximization or fractional roles. With no followers D=1, it nests simultaneous
Cournot. These calculations determine the unique interior FOC solution with
real prices. They do not solve price-constrained boundary/KKT problems.

For observed-share anchoring, one can reconstruct primitives without a new
calibration or a demand-normalization choice. At the observed prices,

\[
A_j=\frac{s_j^{obs}}{s_0^{draw}}\exp(h_j^{pre}),\qquad
c_j=p_j^{obs}+\frac{h_j^{pre}}{\beta}.
\]

Use the saved structural draw, exact-k realization, and actual pre/post
ownership and role mapping. The same primitives are held fixed across the
merger. The reference solves log(x)+x=log(A)-1 through a scalar monotone root;
it does not call either production nonlinear equilibrium solver.

## Executed checks

Thirty synthetic checks cover output and input signs, a multiproduct leader,
a multiproduct follower, and two simultaneous leaders. Direct Logit shares
agree within 1.11e-16. Finite differences of follower profits and leader profits
with fully solved follower responses have maximum gradient magnitude 2.43e-11.

The observed-anchor PNB check includes all 200 CF and 200 CS draws, including
every saved failure. It uses the physical coefficient a_market/price_sd and
the saved draw permutation, outside shares and k assignments. Prices and
shares are independently substituted back into the audited FOCs.

| model | saved successes with positive oracle prices | saved failures with negative oracle prices | maximum FOC residual | maximum successful price-change difference |
|---|---:|---:|---:|---:|
| CF | 164/164 | 36/36 | 2.16e-12 | 1.33e-13 |
| CS | 195/195 | 5/5 | 4.36e-13 | 5.98e-13 |

All negative CF rates occur at the target bank; all negative CS rates occur
at the acquirer. The target is ranked third; at k=2 it joins the acquirer's
strategic owner. CF's target previously had the fringe multiplier h=1, whereas
CS's target was already a strategic follower. Ownership and role propagation
therefore change the implied markdown differently across the two games.

The production solvers search p=exp(z)>0. For these 41 draws, the unique interior
FOC solution lies outside that search domain. Increasing solver iterations or
loosening residual tolerances cannot create a strictly positive interior root.
This establishes the positivity obstruction for the observed-anchor draws; the
failed objects still do not identify the particular internal call that first
raised the error. Latent-anchor failures are not classified by this check.

A negative deposit rate can represent a charge to depositors. Allowing real
Logit counterfactual prices or imposing a zero-rate boundary is a separate
domain/economic decision. Simply clipping a negative rate to zero would not
solve the coupled equilibrium or its boundary conditions. No such change is
made here. The two PNB Bertrand HMC divergences remain a separate diagnostic.

The [400-row oracle ledger](market_quality_pnb_cournot_boundary.csv) retains
each draw's sign classification, k, minimum-rate bank and numerical errors.
A fresh Luna reviewer approved the mathematics, physical slope conversion,
saved ownership/role mapping, structured scalar calculation, and fixed-state
validation. Its 400-row summary reformats the root fixture and is not a second
independent 400-state calculation. An earlier raw-Stan-CSV probe used chain rows
as permuted draw indices; its mismatched results are preserved as untrusted.
The approved reference uses RStan's saved permutation. Earlier native target
and factorization reviews remain separate evidence.

The [synthetic reference](market_quality_cournot_boundary_reference.py) and
[saved-draw reference](market_quality_cournot_boundary_pnb_reference.R) are
standalone scripts. The latter requires the locally retained fits and policy
plans and verifies the stored merger scenario. It extracts only a_market and
s0_prob, reads the saved k realization, and retains all 400 draws. The full
ledger maxima include failed states; a joint R aggregate with missing successful
price errors would omit them and is not used for the all-draw residual bound.
