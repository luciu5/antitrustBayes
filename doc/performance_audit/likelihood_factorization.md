# Likelihood factorization gates for vectorization and threading

This specification preserves the normalized baseline density. It grants no
permission to change conduct, covariance, priors, support or sampler settings.
Each candidate requires native fixed-parameter target/gradient agreement and
measured benefit before acceptance.

## Sampled quality, simultaneous and smooth fringe

Conditional on the sampled quality coordinates and common parameters, the
share terms are observation additive. Observed margins contribute either an
independent margin density or the current conditional margin-given-share
density. Missing margins contribute no margin term. The exact current
Normal/Student-t scale convention, including the t standard-deviation
conversion and rho-dependent conditional degrees of freedom, is retained.

A task may own a disjoint set of observations and return the sum of their
share and observed-margin densities. Parameter priors, truncation constants,
Jacobian terms and HMT penalties are evaluated once outside tasks. Vectorized
Normal/t evaluation is considered first; grouping cannot change observation
matching or treat a missing margin as an observed zero.

## Exact ordered-k, sampled quality

For market m, one k applies to all periods and products in its partition.
Its margin contribution is

\[
\log L_{m,margin}=\operatorname{logsumexp}_{k\in K_m}
\left[\log\pi_{mk}+\sum_{i\in m,\,margin\ observed}
\log p(y_i^m\mid y_i^s,\theta,k)\right].
\]

A task owns complete markets, including their complete mixtures. The common
sampled-quality share density is evaluated once, either as an independent
observation sum or within the complete market task when valid. Never place
log_sum_exp around individual product likelihoods. Predetermined firm ranks
and aggregation use the baseline maps unchanged.

## Collapsed Gaussian quality

Global persistent quality is integrated out, so the joint share density is
not an observation or market sum when a quality group crosses markets. Keep
the existing structured global share kernel outside any observation/market
reduction. With rho fixed to zero, margins do not condition on the globally
correlated share residuals, and complete-market margin mixtures remain
additive; only that separated margin term may be threaded. The share density
is evaluated once per posterior state.

With active rho, integrating quality changes the conditional margin covariance.
No market reduction is permitted unless an explicit valid block factorization
has been established for the selected supported configuration. Existing
single-market exact-k joint conditioning supplies no useful market concurrency.
Do not silently replace this density with independent conditional margins.

## Evidence and deployment

Use a real threaded Stan build and explicit runtime controls; linked TBB or
an unused R grainsize alone provides no threaded target execution. Compare
single-thread and multi-thread normalized targets, chain diagnostics, resource
use and wall time over frozen N/G/market counts and grainsizes. Floating-point
sum order may change draws, so require target/gradient tolerance and posterior
equivalence rather than identical HMC arrays for likelihood reductions.
Threading remains explicit for small fixtures where overhead may dominate.
