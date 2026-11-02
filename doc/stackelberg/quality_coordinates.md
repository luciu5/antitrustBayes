# Equivalent firm-quality coordinates (native specification 8)

This is a computational change proposed after the initial posterior recovery
runs and longer dense-metric pilot failed prespecified sampling diagnostics.
It changes neither the economic model nor the prior on structural quantities.
Sampling improvement remains unverified until completed posterior fits pass diagnostics.

## Analytically established

Let there be G firm groups, sigma>0, and P=I-11'/G. Previously, raw coordinates
u had a normalized standard multivariate Normal prior, and demand quality was
b=sigma P u. The four likelihoods use u only through b. The discarded direction
t=1'u/sqrt(G) is therefore an independent standard Normal absent from the
likelihood, and its normalized density integrates to one.

Choose the G by (G-1) orthonormal Helmert matrix Q, with Q'Q=I, Q'1=0,
and QQ'=P. For column k, the first k entries are 1/sqrt(k(k+1)), entry k+1
is -k/sqrt(k(k+1)), and the remaining entries are zero. Define

    beta = sigma Q'u,     b = Q beta,
    beta | sigma ~ Normal_(G-1)(0, sigma^2 I).

Both representations give Cov(b|sigma)=sigma^2 P, including the original
component variance sigma^2(1-1/G). No component-variance rescaling is applied.
The prior on sigma remains the same normalized positive-truncated
Normal(1, prior_sigma_meanval_strat).

For the augmented inverse transformation u=Q beta/sigma+t1/sqrt(G),
the absolute Jacobian is sigma^(-(G-1)); retaining sigma makes the derivative
matrix block triangular. Consequently,

    log phi_G(u) - (G-1) log(sigma)
      = log phi_(G-1)(beta; 0, sigma^2 I) + log phi(t).

Integrating t preserves the joint prior on b and sigma, the likelihood,
posterior structural quantities, and normalized model evidence. Stan uses
`target += normal_lpdf(b_group_centered | 0, sigma_b_strat)` and
`b_quality_all = Q * b_group_centered`. The Normal density already contains
the sigma normalization: no separate Jacobian term and no second sigma
factor belong in the Stan target or quality expression.

For G=1, Q is 1 by 0, beta is zero-length, and b=0. Its zero-dimensional
Normal contributes zero log density. The scalar sigma retains its normalized
prior; the old independent u direction integrates to one.

## Numerically verified

Independent algebraic checks cover 48 combinations of G=1,3,5,8,
sigma=.02,.1,1,5 and t=-1.3,0,.8. Maximum quality identity error is
2.36e-16, density/Jacobian error 1.42e-14, and covariance error 3.33e-16.
A separate reviewer re-executed these checks and inspected all four
likelihoods for use of the discarded direction.

Against saved native-specification-7 evaluations, 288 mapped compiled cases
verify the native-specification-8 target identity above, product-quality
identity, and generated log-likelihood identity. Maximum discrepancies are
1.14e-13, 1.67e-16, and 3.55e-15 respectively. These include 12 single-group
cases and all four Stan graphs.

Independent reference checks additionally pass 180 prior/flag cases,
96 whole-unit likelihood cases, 20 finite-difference gradient checks,
24 two-draw reconstruction and combined ownership/cost panels,
12 singleton panels, and 24 all-leader nesting comparisons.
The maximum relative gradient error is 8.56e-8 against the unchanged 2e-7
gate. Detailed maxima and unchanged tolerances are in
[native8_gate_summary.csv](native8_gate_summary.csv).

These are deterministic mathematical, compiled-density, and simulation
plumbing checks. They establish no posterior recovery, policy coverage,
model discrimination, or PNB result. Full package checks and independent
implementation review are separate release gates. The final normal source-package
check passes Status: OK, with 971 passing assertions, 12 expected warnings,
zero failures and zero skips. Both vignettes rebuild successfully. The fresh
independent implementation review also passes. The subsequent unchanged-data
computational pilot is separate from these implementation checks.

## Fit provenance

Native specification 8 has different internal parameter coordinates and
dimension from specification 7. Serialized old fits require refitting for
current comparison and simulation APIs. New computational experiments use
new frozen source directories, strict specification-8 metadata and installed
code fingerprints; old failed runs and checkpoints remain local evidence.
No prior, DGP, sampling acceptance gate, or economic equation is changed
to obtain a passing result.
