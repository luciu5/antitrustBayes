# Candidate: exact firm-level Logit quantity equilibrium

Status: analytically derived under the assumptions below; independent numerical
and package-equivalence verification pending. No production solver change is
approved by this note. The lifecycle/extraction units remain first priorities.

Consider one market, MNL outside utility normalized to zero, absolute slope
\(a>0\), constant product marginal costs \(c_j\), full firm ownership, no
capacity constraints and an exogenous market size. Let \(s_0>0\),
\(S_f=\sum_{j\in f}s_j\), and normalized mean utilities \(\delta_j\).
Define

\[
A_f=\sum_{j\in f}\exp(\delta_j-a c_j),\qquad
\omega_{j|f}=\exp(\delta_j-a c_j)/A_f.
\]

Inverse output prices satisfy
\(p_j=(\delta_j-\log(s_j/s_0))/a\). Follower quantity FOCs give
\(p_j-c_j=(1+S_f/s_0)/a\), hence, with \(u_f=S_f/s_0\),

\[
u_f e^{u_f}=A_f/e,\qquad u_f=W_0(A_f/e).
\]

The \(u_f\)'s do not depend on the leader quantities. If
\(U_F=\sum_{f\in F}u_f\) and \(Q_L=\sum_{\ell\in L}S_\ell\), then

\[
s_0=(1-Q_L)/(1+U_F),\qquad S_f=u_f s_0.
\]

For simultaneous, noncooperative quantity leaders the verified FOC is
\(p_j-c_j=[1+S_\ell/(s_0+S_F)]/a\). Set
\(v_\ell=S_\ell/(1-Q_L)\). Since \(s_0+S_F=1-Q_L\),

\[
v_\ell e^{v_\ell}=\frac{A_\ell}{e(1+U_F)},\qquad
v_\ell=W_0\!\left(\frac{A_\ell}{e(1+U_F)}\right).
\]

Writing \(V_L=\sum_{\ell\in L}v_\ell\) gives the complete equilibrium:

\[
s_0=\frac{1}{(1+U_F)(1+V_L)},\quad
S_f=\frac{u_f}{(1+U_F)(1+V_L)},\quad
S_\ell=\frac{v_\ell}{1+V_L}.
\]

Within each firm \(s_j=\omega_{j|f}S_f\). Prices are recovered with common
dollar margins \((1+u_f)/a\) for followers and \((1+v_\ell)/a\) for leaders.
With no followers, this nests ordinary simultaneous Cournot. Ordinary
Cournot also follows by treating every firm as a follower with no leaders.

Each firm's quantity profit has negative-definite Hessian
\(-[\operatorname{diag}(1/s_j)+(2/R+S_f/R^2)11']/a\), where \(R=s_0\)
for followers holding other quantities fixed and \(R=1-Q_L\) for leaders
after solving followers. Thus these stationary allocations are best responses;
the positive Lambert equations and normalization determine a unique interior
allocation. Product price positivity is an additional admissibility condition
for package paths that require positive prices.

For input markets with constant product values \(v_j\) and utility
\(\delta_j+a w_j\), replace \(\delta_j-a c_j\) by
\(\delta_j+a v_j\); recover \(w_j=v_j-\mu_f\). Positive procurement
prices must be checked explicitly.

Implementation must work in log-A coordinates to avoid overflow, retain the
package's actual normalization and cost state, and verify every result against
the current demand and FOC methods. Fractional/cross ownership, an inside
normalization good, quantity-dependent costs, nests, capacity constraints and
CES are outside this derivation. Generic paths remain necessary.

A firm-level implementation could avoid both dense ownership algebra and a
product-sized nonlinear solve. Before acceptance require independent direct
optimization, package prices/shares/welfare agreement, multiproduct cases,
mergers/cost shocks, input cases, role changes and existing solver tolerances.
