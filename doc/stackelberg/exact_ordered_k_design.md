# Exact ordered strategic-core marginalization

Status: implementation specification, 2026-09-14. This supersedes the smooth-cutoff preference in the earlier implementation map. The deterministic economics and the Stackelberg math audit remain unchanged.

## Models and support

For observed firms ranked by firm market share, define BF(k) and CF(k) as simultaneous multiproduct Bertrand or Cournot conduct for the top k firms and product-level monopolistic competition for every remaining firm. Their support is k=0,...,Kmax; k=0 is exactly the existing monopolistic-competition equation. Define BS(k) and CS(k) as the audited noncooperative Stackelberg games with the top k firms as leaders and all remaining firms as strategic followers. Their support is k=1,...,Kmax. No Stackelberg k=0 state is admitted.

`ordered_k` is the preferred partition mode for explicit `structure="fringe"` and `structure="stackelberg"`. The existing smooth fringe cutoff remains available as an explicit legacy/robustness mode. Legacy calls that select the old cutoff without an explicit structure retain their prior meaning. Simultaneous models have no partition.

The public API is `k_support` plus `k_prior`. `k_support` is an explicit integer vector: it may include zero for fringe and must be strictly positive for Stackelberg. `k_prior` has the same length, contains finite nonnegative values, and sums to one within numerical tolerance; omitting it gives a uniform prior over the supplied support. The package never hard-codes the PNB 0:5 or 1:5 supports. A supplied support must be feasible in every selected market and is rejected rather than silently truncated or renormalized. Users can pass the same support/prior to BF and CF, and the same positive support/prior to BS and CS. Simultaneous models ignore both arguments. Store both vectors verbatim in fit metadata and provenance.

## Firm ranking and panel scope

One latent k_m and one firm ranking apply to the complete market likelihood, including all selected years. Ranking is firm-level, never product-level. Because the outside share is latent during estimation, the predetermined observable metric is conditional inside quantity share (Logit) or conditional inside revenue share (CES); using posterior unconditional shares would make roles parameter-dependent during HMC. For a panel, construct one unique firm-market table, sum products within each firm-period, give an absent firm-period zero, and average across all observed periods with equal period weight. Rank descending by that market aggregate, with an explicit deterministic firm-label tie break. Store `obs_to_firm_market` and the resulting firm-market rank; every product and year of the firm inherits that rank. Entry/exit merely makes a ranked firm absent in a period and does not rerank the other firms for that period. Record the aggregation and tie rule in fit metadata. Conditional and unconditional shares induce the same ranking within a period; the conditional panel metric is the identified prespecified extension when latent outside shares vary by period.

This replaces the current ordered-k behavior that shares only the count while ranking firms separately in each market-year. Period-specific ordered k is not part of the preferred exact specification. A separately named robustness mode may retain it, but it must not be called market-level marginalization.

## Exact conduct predictions

Use unconditional firm shares S_f and outside share s0 for Logit. For a core firm, BF uses h_f=1/(1-S_f), and CF uses h_f=1+S_f/s0. Every fringe product in BF(0:K) or CF(0:K) uses h=1, the exact monopolistic-competition margin equation, regardless of portfolio. Predicted dollar markup is h/alpha and predicted inverse Lerner is alpha*p/h.

For CES, let R_f be unconditional firm revenue share, R0 the audited outside revenue share, gamma>1, h=gamma-1, H=1+h*R0, and b=h/(gamma*H). A BF core firm has Lerner margin 1/[gamma-h*R_f]; a CF core firm has margin 1/gamma+b*R_f. Every CES fringe product has margin 1/gamma. These are the existing simultaneous multiproduct and monopolistic equations under the package's CES convention.

BS/CS retain the signed-off leader/follower equations and firm aggregation exactly. No formula is interpolated across roles.

The exact BF/CF likelihood does not introduce the legacy `mu_b_fringe` demand-quality shift. That shift remains specific to the smooth-cutoff robustness model. Consequently the share likelihood is common across k in every exact ordered-k model; only the conditional margin likelihood changes with k. With correlated share/margin residuals, use the existing conditional margin density evaluated using the common share residual.

## Market likelihood and generated quantities

For market m and continuous state theta, sum every observed margin contribution in every period/product under one k before mixing:

    score[m,k] = log_pi[m,k] + sum_{i in market m, margin observed} log p(y_margin_i | y_share_i, theta, k)
    target += share_loglik[m] + log_sum_exp(score[m,])

Computing share likelihood outside the mixture is algebraically identical to mixing the complete joint market likelihood because it is common across k. Never apply `log_sum_exp` observation-by-observation. Multiple markets use one independent `log_sum_exp` per market, giving O(sum_m K_m), not a Cartesian enumeration.

Generated quantities return, for every market/component, the supplied prior probability, component log score, normalized posterior probability, k, expected k, modal/median k, core-share summaries, and firm leader/core probabilities. A public `posterior_k(fit)` helper returns one row per market and supplied k with prior and posterior probabilities; it may also attach mean/median/credible summaries without obscuring the probability table. Whole-market log-likelihood/LOO contributions use the same mixture. Product-row likelihood fields may be diagnostic allocations only and must telescope exactly to the whole-market contribution; they cannot imply row-specific k states.

## Simulation mapping

For deterministic posterior simulation, sample or select one k_m from each draw's posterior component probabilities and use the same firm-market ranking and k across all periods and scenarios for that draw, unless the counterfactual explicitly supplies a new role rule. BF/CF reconstruct simultaneous strategic-core plus exact MonCom-fringe objects/equations. BS/CS reconstruct the existing coordination Stackelberg object. Sensitivity reports include posterior-mixture simulation, modal k, posterior expected outcomes by k, and adjacent states for firms near the Kmax/prior boundary.

## Required gates

Tests must cover k=0 fringe=MC exactly; k=N fringe=ordinary simultaneous Bertrand/Cournot; k=N Stackelberg=ordinary simultaneous conduct where established; one market with several periods and changing per-period share order but invariant market ranking; several markets with separate O(MK) mixtures; multiproduct firms counted once; common prior reuse and Stackelberg conditioning off k=0; exact whole-market mixture versus brute-force R; generated quantities versus likelihood; synthetic recovery for BF/CF/BS/CS with more than one seed; model discrimination; ownership, cost, and combined counterfactual recovery. The existing smooth-fringe tests remain as robustness/legacy tests.
