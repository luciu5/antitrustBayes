# Explicit real-rate domain for input Logit Cournot counterfactuals

This extension targets the diagnosed CF/CS simulation failures. It does not
change Stan targets, demand coefficients, costs/values, ownership, role mapping,
or the existing interior conduct equations. It expands the deterministic action
outcome domain when explicitly requested. The positive-rate path remains the
default. It is not a zero-floor or constrained-equilibrium model.

## Approved initial scope

`price_domain = "real"` is supported initially for input-market Logit,
Cournot conduct, and native exact core-fringe or Stackelberg simulation.
Unsupported output/CES/Bertrand combinations must fail explicitly. Observed
prices supplied to estimation/calibration remain positive; this option concerns
solved latent pre-policy and post-policy rates. The common deterministic
implementation belongs in coordination, not in Stan or antitrustBayes.

## Exact interior solution

Let beta > 0 be the physical input Logit coefficient, d_j the stored mean value,
p_0 the outside normalization, and v_j the preserved unit procurement value.
For any finite real rate,

\[
u_j=d_j+\beta(p_j-p_0),\qquad s_j/s_0=e^{u_j}.
\]

Write the verified markdown FOC as

\[
p_j=v_j-h_f/\beta,\quad
A_j=e^{d_j+\beta(v_j-p_0)},\quad A_f=\sum_{j\in f}A_j.
\]

Firm aggregation includes every active owned product exactly once. With
x_f=S_f/s_0, the equilibrium identity is x_f=A_f e^{-h_f}.

In CF, strategic-core h_f=1+x_f and MonCom fringe h_f=1. Therefore

\[
x_f=W(A_f/e)\quad(f\in\mathrm{core}),\qquad
x_f=A_f/e\quad(f\in\mathrm{fringe}).
\]

In CS, first solve every strategic follower x_f=W(A_f/e). Set
D=1+sum_{f in F}x_f. Simultaneous noncooperative leaders satisfy

\[
x_\ell=D W\{A_\ell/(eD)\},\qquad
h_\ell=1+x_\ell/D,\qquad h_f=1+x_f\ (f\in F).
\]

The common-dollar-markdown property preserves product-level demand and values.
Recover s_0=1/(1+sum_f x_f), product quantities from Logit demand and
p_j=v_j-h_{f(j)}/beta. The equation x+log(x)=L has exactly one positive
solution for finite L, since its derivative 1+1/x is positive. A log-domain
scalar root supplies W(exp(L)) without computing a large exponential first.
No K^M enumeration, joint-leader maximization, fractional role or new demand
normalization is introduced.

For CS, follower quantities respond affinely to leader total quantity: writing
M for market size, q_f=(x_f/D)(M-Q_L), with product composition fixed by its
within-firm A_j weights. Thus a leader's verified reduced-profit derivative
uses s_0+S_F, as in the existing model. Input follower and reduced leader profit
are strictly concave in their own positive quantity vectors when an outside
quantity remains: the inverse Logit own-quantity derivative contributes the
negative diagonal 1/(beta*q_j) term and additional negative firm-total terms.
For own portfolio quantity q and firm total Q, the Hessian is

\[
-\frac{1}{\beta}\left[\operatorname{diag}(1/q_j)
 +\left(\frac{2}{H}+\frac{Q}{H^2}\right)11^\prime\right],
\]

where H is the remaining outside quantity for a Cournot follower, and the
remaining quantity after all leaders for a reduced CS leader. Both are
positive; this matrix is negative definite. Market size simply rescales these
quantity coordinates. This preserves the interior-game interpretation; it does not establish an
alternative equilibrium under a rate floor.

## Reporting and verification requirements

Dollar markdown is computed directly as h/beta. Its price-normalized version
is h/(beta*p): undefined at zero, and signed at negative rates. A zero rate
must return NA for that ratio while retaining the finite dollar markdown.
Producer surplus must use dollar markdown times quantity; no 0*Inf arithmetic.
Revenue/expenditure ratios with nonpositive rates must not be presented as
probability shares. Quantity shares, level-rate changes and the existing Logit
utility-based welfare formulas remain meaningful. Percentage changes from a
zero latent pre-policy rate are undefined and must be NA with a flag.

The production real-domain solver must check the original FOCs at the returned
prices with the existing tolerances. Tests must independently verify product
shares, firm aggregation, costs/values preserved through pre/post, residuals,
quantity best responses/leader reduced profits, and welfare. Successful
positive-domain reference cases must agree numerically. Preserve every prior
failure; classify latent pre versus post separately. The PNB comparison must
reuse the exact saved posterior/quality/k plan and costs, not refit or resample.

The option has an economic interpretation only if negative rates can represent
net depositor charges in the application. It is an explicit extension, not a
claim that historical recorded nominal rates could be negative.
