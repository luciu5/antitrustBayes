# Exact ordered-k Stan contract

Status: lead implementation contract, 2026-09-15. This specializes
`exact_ordered_k_design.md` to the existing Stan programs. The equations come
from `../stackelberg_math_audit.md`; this note does not introduce a continuous
role game.

## Partition and ranking

An explicit exact fringe or Stackelberg fit has one partition unit per market.
Every unit uses the user-supplied `k_support` and positional `k_prior`.
Fringe admits nonnegative k and interprets k=0 as exact monopolistic
competition. Stackelberg requires k>=1. A support value may not exceed the
number of distinct firms in any selected market.

Firm roles use one predetermined rank for the whole observed market panel.
The observable ranking metric is the equal-period mean across all observed
periods in the market of its product-summed conditional inside share (Logit)
or conditional inside revenue share (CES), treating an absent firm-period as
zero. These are the size coordinates available before HMC when
the outside share is latent. Using a posterior outside share in the ranking
would violate the fixed-ranking requirement. The rank is descending in this
metric with the character firm label as the deterministic tie break. Every
firm-period record inherits the firm-market rank; entry and exit do not rerank
the firms that are present in a period.

Observation-level prices, responses, masks, market/year IDs, covariates, and
outside-price fields retain caller row order. `obs_to_firm_mky` and
`obs_to_firm_market` map those rows to the auxiliary tables. Any likelihood
traversal permutation is stored separately and applied symmetrically to all
row-indexed quantities.

## Conditional inverse-margin predictions

All Logit shares below are unconditional at the evaluated continuous state:

    S_f = firm_shareIn_f * (1 - s0).

Let h_f satisfy the common firm dollar markup `mu_f=h_f/alpha`. The predicted
inverse Lerner for product i owned by f is `alpha * price_i / h_f`.

For exact Bertrand fringe BF(k):

    h_f = 1/(1-S_f)  if rank_f <= k,
          1          otherwise.

For exact Cournot fringe CF(k):

    h_f = 1+S_f/s0  if rank_f <= k,
          1         otherwise.

For Bertrand Stackelberg BS(k), followers have rank>k. Define

    B_F = sum_{rank_g>k} S_g^2/(1-S_g+S_g^2),
    phi = 1/(1-B_F).

Then

    h_f = 1/(1-S_f)       for a follower,
          1/(1-phi*S_f)   for a leader.

For Cournot Stackelberg CS(k), define

    S_F = sum_{rank_g>k} S_g.

Then

    h_f = 1+S_f/s0          for a follower,
          1+S_f/(s0+S_F)    for a leader.

For CES, first convert the sampled outside quantity share to the physical
outside revenue share using the package normalization:

    R0 = inv_logit(logit(s0) + log(priceOutside) - log(insidePriceIndex)),
    R_f = firm_revShare_f * (1-R0),
    gamma > 1, h = gamma-1, H = 1+h*R0, b = h/(gamma*H).

Stan predicts inverse Lerner `1/m_f`. For exact fringe:

    BF core:    m_f = 1/(gamma-h*R_f)
    CF core:    m_f = 1/gamma+b*R_f
    BF/CF fringe: m_f = 1/gamma.

For CES Bertrand Stackelberg, followers have

    t_f = gamma-h*R_f,
    a_f = h^2*R_f/[t_f*(t_f-1)+h^2*R_f],
    B_F = sum_{followers} R_f*a_f,
    phi = 1/(1-B_F),
    m_follower = 1/(gamma-h*R_f),
    m_leader   = 1/(gamma-h*phi*R_f).

For CES Cournot Stackelberg, follower `a_f` and `C_F` use the algebraically
verified stable expression in the math audit/current Stackelberg program:

    a_f = h*(H-2*R_f+gamma*h*R0*R_f/H)
          / [H*(H+(gamma-2)*R_f)],
    C_F = sum_{followers} R_f*a_f,
    m_follower = 1/gamma+b*R_f,
    m_leader   = 1/gamma+b*R_f/(1+C_F).

Membership tests use the inherited firm-market rank, not a record's position
in the firm-period table. This is necessary when a higher-ranked firm is absent
in a period.

## Whole-market mixture

The share equation and its likelihood are common across k and are evaluated
once per observation outside the mixture. For market m and component k:

    score[m,k] = log(k_prior[k])
                 + sum_{i in market m, margin observed}
                     log p(y_margin_i | y_share_i, theta, k).

The model contribution is `log_sum_exp(score[m,])`, once for each market.
There is no Cartesian enumeration across markets and no observation-level
mixture. Implementation may use a market-specific observation index so each
component visits only rows in its market; its unavoidable cost is
`sum_m O(K_m*N_margin,m)` and its state count is `sum_m K_m`.

The conditional margin density must retain the existing correlated Student-t
factorization. Exact ordered-k models have no `mu_b_fringe` quality shift, so
the share residual used by that conditional density is identical across k.
Zero prior probabilities contribute `negative_infinity()` and are never
renormalized to a different user prior.

## Generated quantities and accounting

For every market/component, return its k, supplied prior probability, complete
component log score, and posterior probability from `softmax(score)`. Return
per market posterior mean, median, and modal k; posterior core probability for
each firm-market record; and posterior strategic-core share under the fixed
ranking metric.

Posterior predicted margins average component predictions with the conditional
posterior component probabilities. Whole-market `log_lik_market` equals the
share likelihood plus the exact margin mixture. If row-level diagnostic
`log_lik` fields are retained, compute ordered telescoping increments using one
common row permutation and test that they sum exactly to the market value.
They are accounting allocations, not rowwise latent-k likelihoods.

The public `posterior_k()` helper aggregates generated component probabilities
over posterior draws without presenting k as a sampled Stan parameter.

## Required deterministic gates before sampling

1. R and compiled Stan conditional predictions agree for Logit/CES x BF/CF/BS/CS.
2. BF(0)=CF(0)=MC to machine tolerance.
3. BF(N)/CF(N) equal ordinary simultaneous multiproduct conduct.
4. BS(N)/CS(N) equal ordinary simultaneous conduct where established by the audit.
5. Whole-market scores and probabilities match an independent brute-force R oracle.
6. Panel ranking is invariant to observation order and period-specific rank reversals.
7. Multiproduct firms enter follower/core sums exactly once per firm-period.
8. Likelihood and generated-quantity accounting telescope to machine tolerance.

Only after these gates and the full package suite pass should multi-seed HMC
recovery resume.
