# Approved CES implementation contract

Production CES work starts only after the Logit gate. This contract specializes
the signed-off audit; it does not replace that derivation.

## Demand and state

StackelbergCES extends antitrust::CES, with the same firm ownership, role,
conduct and diagnostics fields as StackelbergLogit. Reuse the common constructor
and public followers/response/residuals/simulation functions; dispatch by demand
class. Keep output-only CES, gamma>1, positive outside price and an explicit
outside revenue share. Constructor shares are unconditional revenue shares.
Do not redefine antitrust's legacy Stackelberg class or coordinated leadership.

Store antitrust CES slopes exactly: alpha=1/shareInside-1,
gamma=gamma, meanval_j=(r_j/r0)*(p_j/p0)^(gamma-1),
mktSize=insideSize/(1-r0). insideSize is inside expenditure and E=mktSize.
For forward demand, D=p0^(1-gamma)+sum A_j*p_j^(1-gamma),
r_j=A_j*p_j^(1-gamma)/D, q_j=E*r_j/p_j. Reuse inherited demand/output
methods for equilibrium shares, quantities, surplus and elasticities.

The generic inverse-demand oracle may solve one positive scalar w:

    p_j=w*(A_j/q_j)^(1/gamma)
    p0^(1-gamma)*w^gamma + w*sum(A_j^(1/gamma)*q_j^((gamma-1)/gamma)) = E.

The left side increases strictly from zero, so bracket and solve without an
arbitrary upper quantity bound. Any positive finite inside quantities admit a
positive outside expenditure. This is the same CES system, not a new quantity
normalization. Use log arithmetic when price/quantity ranges require it.

## Exact firm FOCs and analytic responses

Let h=gamma-1, R_f=sum_owned r_j, H=1+h*r0, b=h/(gamma*H).
All products within a firm have a common Lerner margin, generally different
dollar margins. Roles remain fixed during derivatives and deviations.

Bertrand:

    t_f=1+h*(1-R_f)
    a_f=h*R_f/[1-R_f+h*(1-R_f+R_f^2)]
    B=sum_followers R_f*a_f; phi=1/(1-B)
    m_follower=1/t_f
    m_leader=1/[gamma-h*phi*R_leader]
    dp_jF/dp_kL=p_j*a_firm(j)*phi*r_k/p_k.

The stable a_f above equals h^2*R_f/[t_f*(t_f-1)+h^2*R_f]; verify this
identity numerically before production use. For multiproduct followers, the
price response is proportional to each follower product price.

Cournot:

    a_f=h*(H-2*R_f+gamma*h*r0*R_f/H)/[H*(H+(gamma-2)*R_f)]
    C=sum_followers R_f*a_f
    m_follower=1/gamma+b*R_f
    m_leader=1/gamma+b*R_leader/(1+C)
    dq_jF/dq_kL=-q_j*a_firm(j)*p_k/[E*(1+C)]
    K_jk=dp_j/dq_k=-p_j/(gamma*q_j)*I(j=k)-b*p_j*p_k/E.

The stable Cournot coefficient is independently verified in the audit. It can
be negative. Never assume follower quantities always fall with leader output.
For interior positive outside share a_f>-1 and
1+C>1-sum_followers R_f=r0+sum_leaders R_f>0.
Use conditioning diagnostics near very small outside share.

## Calibration and solver gates

At each trial gamma compute hard-partition margins from observed revenue
shares; minimize sum usable w_j*(log(1/m_obs_j)-log(1/m_pred_j))^2.
Default weights one, partial missing margins supported; all missing requires
fixed gamma. Preserve observed margins in diagnostics. A fixed gamma bypasses
calibration. Search gamma>1 with transparent bounds/adaptive expansion or
explicit controls and report boundary optima. Never silently claim a bounded
optimum is globally interior. Costs are p_j*(1-m_pred_j), then proportional
mcDelta; do not invert noisy moments or recalibrate after shocks.

The joint price representation solves p*(1-m(p))-c=0 for either conduct and
the exact given roles. It is an efficient representation of the quantity game
for Cournot, not a change of strategic actions. Expose actual follower pricing
or quantity equilibria for fixed leader actions. Return quantities and prices
with leader actions preserved and actual follower FOC residuals.

Generic response uses raw profit gradients with fixed costs and demand:
F_j=q_j+sum_owned (p_i-c_i)*dq_i/dp_j for Bertrand;
F_j=p_j-c_j+sum_owned q_i*dp_i/dq_j for Cournot.
Implicit R=-solve(dF_F/dz_F,dF_F/dz_L) with z=p or q, respectively.
Every leader own total profit gradient uses its own owned products, holding
other leaders fixed. Implement an optional equilibrium route through this
generic response and fully solved followers as an independent check.

Tests must cover gamma on both sides of two, negative Cournot reaction
coefficients, multiproduct roles, all/no followers, very small followers,
outside normalization not one, margin calibration from several noisy moments,
fixed gamma/all-missing moments, ownership-only/cost-only/combined shocks,
explicit post roles, active products if supported, perturbed baseline starts,
analytic/generic reactions and equilibrium solutions, finite differences of
fully solved followers, and one-leader direct reduced-profit optimization.
Keep tight gates; failure is a defect or documented conditioning limitation.
