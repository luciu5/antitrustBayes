# Likelihood threading: implementation contract

This specifies the next measured unit after the Helmert gate. It does not claim
that any current sampling call uses within-chain threads.

## Preserve a serial reference

Keep priors, truncation constants, coordinate Jacobians, HMT terms, role support
and all economic predictors outside the likelihood reduction exactly as before.
First compare a serial vectorized likelihood with the existing normalized target.
Only then add threaded dispatch of that same likelihood function. Do not combine
these changes with a generated-quantities split or a covariance rewrite.

For sampled quality in simultaneous/smooth-fringe kernels, conditional share
and margin densities are additive by observation. Build the existing predicted
log inverse-margin mean for observed margins with the same stabilizing constants.
Vectorize the Normal/Student-t density calls where supported. Missing margins
contribute no margin term; share-only and margin-only likelihood modes retain
their existing meaning. Preserve all conditional rho means, scales and degrees
of freedom, including the package's variance-standardized t convention.

For the serial vectorization unit, construct the observation share-scale
vector from the existing market-indexed scales and call one `normal_lpdf`
or `student_t_lpdf` on the share vector. For margins, first construct vectors
of the current means and scales on `margin_idx` only, including the current
rho adjustments, and call one matching density. An empty margin index must
contribute exactly zero. Keep the collapsed joint-Gaussian branch unchanged.
The identity being used is the sum of normalized scalar log densities equals
the normalized vectorized density; the predictors and conditional
distributions themselves remain unchanged. This unit precedes and is
committed separately from reduction dispatch.

Use an observation index slice with explicit shared predictors and scales for
`reduce_sum`. A one-thread serial path is the numerical reference. Never send
the collapsed share likelihood through an observation reduction: it remains
one global structured evaluation. Initially restrict observation threading to
the sampled-quality path; broaden only if a separately profiled factorization
justifies it.

## Exact ordered-k kernels

The reduction slice contains market IDs. One task owns that market's full
likelihood across its periods and all admissible k states, ending in exactly
one `log_sum_exp(log_prior_k + market_margin_log_lik_k)`. Share likelihood stays
outside the finite-state mixture and is computed once. Do not distribute one
market's observations to independent k mixtures.

Collapsed global quality with rho off has a global share term plus separable
market margin mixtures, so only the latter may be reduced by market. Active
rho configurations must follow their currently supported joint covariance
factorization. The supported single-market collapsed/rho branch stays serial;
do not expand its statistical support to create threading opportunities.

## Runtime and compatibility

Verify compile-time Stan threading and the actual installed rstan runtime
control. R `threads_per_chain` and `grainsize` must reach a consumed Stan/runtime
input rather than disappear from sampling arguments. Preserve caller options
and thread state on success and failure. Default within-chain threads remain
one; explicit threads do not change priors, sampler settings or solver controls.
Distinguish observation grainsize from market grainsize and validate each in
its own units. Handle single-observation/market and zero observed margins.

Every Stan source-hash migration requires its own target/gradient/raw-schema
gate and explicit bridge allow-list entry. Preserve rotation-era private flags.

## Gate and benchmark

Compare old scalar, new vectorized serial, and new threaded normalized targets
and gradients at nonzero fixed coordinates for both demands, both sampled share
likelihoods, Gaussian/t margins, missing margins, rho, pooling/rotation options,
and all BF/CF/BS/CS exact-k models. Compare full GQ at identical raw posterior
rows. Numerical summation order may cause small floating-point differences;
constants or systematic density differences are defects.

Benchmark N=100, 1,000 and 10,000 with several prespecified grainsizes, holding
compiler flags, data, parameter points and HMC settings fixed. Use safe group
counts for the old dense Helmert reference and state their configuration.
Measure fixed target/gradient times before small paired HMC runs. For exact k,
include multiple market counts and periods; a one-market PNB fixture is the
expected overhead control. Retain failures and regressions. Integrate only
paths whose measured gains justify their complexity.
