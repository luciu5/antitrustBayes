# Joint market/outside-share rotation diagnosis

Status: the existing shifted rotation is an exact representation of the intended
posterior, but executed evidence confirms that it creates poor boundary geometry.
The proposed joint orthogonal rotation is also exact, removes the divergence and
treedepth failure on the controlled fixture, and passes one of two strict HMC
gates; the other misses only the R-hat threshold. On the exact Supreme Court
year cross-section, however, it produces 12 and 9 divergences in two fixed-seed
runs. The candidate materially improves the old failure but does not fix it.
This is conclusion B. The upper Supreme Court ladder is not authorized.

## Reproduced facts from the committed Supreme Court run

The source for this section is `mergerBayes` commit
`082c68d51b75bbf55cd1c449ecf1ed6f6771ade5`. The compact committed CSVs, rather
than recollection, establish the following pattern.

* The one-market case panel was clean for all four conduct models: zero
  divergences and zero treedepth hits, maximum R-hat 1.0034, and minimum BFMI
  0.904.
* The six-market year cross-section without outside-share pooling was much
  healthier. MC and Auction had zero divergences, Bertrand had zero or one
  depending on alpha pooling, and Cournot had one or fourteen.
* With `s0_pooling="auto"` (market-level pooling) and the shifted rotation, the
  s0-only fits had 39, 45, 41, and 64 divergences for Auction, MC, Bertrand, and
  Cournot. Bertrand and Cournot also had 879 and 1,285 treedepth hits.
* The both-pooled cross-section was also unusable: 5--76 divergences, as many as
  1,230 treedepth hits, and Cournot R-hat 1.196.
* The 22-cell full panel was unusable. Its four models had 5--599 divergences or
  151--1,159 treedepth hits; Auction R-hat was 1.314 and Bertrand bridge-run SD
  was about 123.

The raw failed shifted-rotation Supreme Court fits were not committed, so their
divergent draws cannot be reconstructed locally. The corresponding bridge
weights and model probabilities must not be interpreted.

## Existing shifted target

For market-level, noncentered pooling the legacy hierarchy is

\[
b_m=\mu_b+\sigma_bu_m,
\qquad
\eta_m=\mu_s+\sigma_sz_m,
\qquad
u_m,z_m\stackrel{iid}{\sim}N(0,1).
\]

Here \(\eta_m\) is the base outside-share coordinate. Production code retains
the existing bounded inverse-logit map when reconstructing \(s_{0m}\).

The historical shifted rotation samples \((u_m,w_m)\), where

\[
w_m=z_m+\frac{b_m}{\sigma_s},
\qquad
z_m=w_m-\frac{b_m}{\sigma_s},
\qquad
\eta_m=\mu_s+\sigma_sw_m-b_m.
\]

The model evaluates the original standard-normal density on the reconstructed
\(z_m\). This is target-equivalent to the legacy hierarchy. With scale
parameters retained, the inverse map is block triangular with unit determinant,
so no Jacobian term is missing. All normalized prior and likelihood constants
are retained.

The problem is geometry rather than target correctness. When \(\sigma_s\) is
small, HMC must maintain cancellation between \(w_m\) and
\(b_m/\sigma_s\), even when the reconstructed economic state is ordinary.

## Joint orthogonal coordinates

Define

\[
S=\sqrt{\sigma_b^2+\sigma_s^2}
\]

and sample independent standard-normal coordinates \(v_m,h_m\). Reconstruct

\[
u_m=\frac{\sigma_bv_m-\sigma_sh_m}{S},
\qquad
z_m=\frac{\sigma_sv_m+\sigma_bh_m}{S}.
\]

Then use the unchanged hierarchy

\[
b_m=\mu_b+\sigma_bu_m,
\qquad
\eta_m=\mu_s+\sigma_sz_m.
\]

Conditional on the two scales, this is an orthogonal map with determinant one.
The full scale-dependent transformation also has determinant one because the
scale coordinates are retained and the Jacobian is block triangular. Isotropic
standard-normal density is invariant under the rotation, so placing independent
standard-normal densities on \(v_m,h_m\) reproduces exactly the original joint
prior on \(u_m,z_m\). The coefficients \(\sigma_b/S\) and \(\sigma_s/S\) are
bounded, and no division by either individual hierarchy scale occurs. Stan's
positive scale transforms keep \(S>0\) in the parameter space; `hypot` is used
for stable evaluation when both scales are small.

The identified sum simplifies to

\[
b_m+\eta_m=\mu_b+\mu_s+Sv_m.
\]

Thus \(v_m\) follows the share-likelihood direction and \(h_m\) follows the
orthogonal direction that the priors and supply likelihood use to separate
\(b_m\) from \(s_{0m}\).

## Fixed-state and gradient evidence

`doc/s0_joint_rotation_reference.R` implements an independent dense Gaussian
oracle with unequal markets, market-centered persistent quality, active rho,
missing margins, and the production bounded outside-share link. It covers
ordinary scales, either scale near zero, and both scales small.

`doc/s0_joint_rotation_native_validation.R` evaluates the compiled Stan targets
and generated quantities for all four kernels: Logit and CES, ordinary supply
and exact ordered-k Stackelberg. It also compares Stan automatic gradients with
finite differences.

The compact results are:

| Evidence | Checks | Maximum error | Result |
|---|---:|---:|---|
| Independent dense target/state oracle | 68 | 8.57e-10 | all pass |
| Independent chain-rule gradients | 8 | 3.46e-8 | all pass |
| Native Stan target/generated quantities | 136 | 5.56e-17 | all pass |
| Native Stan autodiff vs finite differences | 8 | 2.11e-9 | all pass |

The native checks include the exact-k component probabilities and priors. No
conduct equation, quality covariance, rho term, missing-margin rule, bounded
link, CES normalization, or generated structural quantity changed.

## Controlled HMC evidence

The local fixture has 45 observations in six markets and uses simultaneous
Bertrand, joint Gaussian share-margin residuals, collapsed market quality,
active rho, no alpha pool, and the same priors and sampler settings in every
coordinate comparison. Each principal run used four chains, 800 iterations,
400 warmup iterations, `adapt_delta=.95`, and `max_treedepth=13`.

| Coordinate case | Seed | Divergences | Treedepth | Max R-hat | Min bulk ESS | Seconds | Strict gate |
|---|---:|---:|---:|---:|---:|---:|---|
| no s0 pool, shifted | 19091 | 0 | 0 | 1.0067 | 673 | 28 | pass |
| market pool, legacy NC | 19091 | 0 | 0 | 1.0110 | 534 | 44 | fail by 0.0010 R-hat |
| market pool, legacy NC | 29091 | 0 | 0 | 1.0065 | 567 | 51 | pass |
| market pool, shifted NC | 19091 | 17 | 0 | 1.0135 | 316 | 208 | fail |
| market pool, shifted NC | 29091 | 20 | 0 | 1.0468 | 124 | 103 | fail |
| market pool, shifted centered | 19091 | 99 | 0 | 1.0748 | 52 | 82 | fail |
| market pool, joint NC | 19091 | 0 | 0 | 1.0107 | 661 | 46 | fail by 0.0007 R-hat |
| market pool, joint NC | 29091 | 0 | 0 | 1.0061 | 682 | 47 | pass |

The short 200/100 run is retained separately as a failed-budget experiment; it
was not used to declare sampler usability.

In the failing shifted run at seed 29091, the 20 divergent draws had median
\(\sigma_s=0.024\), compared with 0.097 in nondivergent draws. Their median
\(\max_m|b_m/\sigma_s|\) was 23.2, compared with 6.2 in nondivergent draws.
The reconstructed \(z_m\) was not unusually large in divergent draws. This is
consistent with the sampler failing where a large cancellation is required; it
does not rely on correlation alone to establish the mechanism.

The economic correlation between \(b_m\) and reconstructed
\(\operatorname{logit}(s_{0m})\) remains about -0.95 under both legacy and joint
fits, as it must for an exact reparameterization. Under the joint fit, the
posterior correlation between the new raw sampling coordinates \(v_m,h_m\) is
0.016. Across 25 economic posterior quantities, the largest legacy-versus-joint
mean difference was 1.29 combined Monte Carlo standard errors and 0.063 pooled
posterior standard deviations. All prespecified equivalence checks passed.

## Implementation and compatibility

The implementation adds an explicit experimental value:

```r
market_s0_parameterization = "joint"
```

It is accepted only for multiple markets with noncentered market-level
outside-share pooling. Codes 1 (`legacy`) and 2 (`rotated`) retain their exact
historical meanings; code 3 selects the joint rotation. Using a new code avoids
silently changing saved-fit interpretation. Package defaults are not changed.

All four generated Stan bindings were rebuilt. Bridge restoration contains an
explicit forward allow-list from the previously validated kernels to the new
kernel hashes. The installed model-code hashes match the allow-list, and tests
verify that saved code-1 and code-2 private data retain their coordinate flags.
There is no backward migration from code 3. Statistical model specification
version 10 remains appropriate because this adds sampler coordinates without
changing the normalized target or economic model.

## Exact Supreme Court gate

The exact data were subsequently supplied in the local `mergerBayes/data`
directory. The fixture reconstruction matches the committed cross-section: six
markets with 39, 6, 31, 85, 8, and 7 observations, totaling 176. Its stable R
object SHA-256 prefix is `6d4065fe2be5`.

The exact joint-coordinate Bertrand runs used active rho, collapsed market
quality, market-level s0 pooling, no alpha pooling, four chains, 2,000
iterations, 1,000 warmup iterations, `adapt_delta=.99`, and
`max_treedepth=13`. Two chains ran concurrently to limit local resource use.

| Seed | Seconds | Divergences | Treedepth | BFMI | Max R-hat | Min bulk ESS | Min tail ESS |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 19091 | 1,954 | 12 | 0 | 0.724 | 1.0041 | 559 | 492 |
| 29091 | 2,121 | 9 | 0 | 0.652 | 1.0060 | 528 | 741 |

The joint rotation therefore removes the old treedepth failure and improves
mixing and divergence counts, but it fails the prespecified zero-divergence
gate reproducibly. Healthy R-hat, ESS, and BFMI do not override this result.

The exact joint sampling coordinates are nearly uncorrelated (-0.017 and
-0.013), and the economic market-effect/outside-share correlation is only about
-0.15 on this fixture. The remaining divergences are not a recurrence of the
original `b_market / sigma_s0_pool` cancellation.

A separately labeled half-budget diagnostic reproduced eight divergences. In
those draws, the largest standardized mean differences relative to
nondivergent draws were:

* rho: -1.07 standard deviations (divergent median -0.44 versus 0.05);
* residual share scale: -0.94 standard deviations (0.24 versus 0.43);
* quality scale: +0.35 standard deviations (1.47 versus 1.44).

The divergent median of
`sigma_quality / sigma_share` was 6.28 versus 3.45. The effective conditional
share scale `sigma_share * sqrt(1-rho^2)` had medians 0.19 versus 0.41. These
associations are based on eight rare events and do not establish causality, but
they point away from the s0 hierarchy and toward the collapsed quality/share
variance decomposition and its interaction with rho.

## Next exact geometry experiment

The next candidate should align the two share-variance components without
changing their priors. Let

\[
T=\sqrt{\tau^2+\sigma_s^2},
\qquad
\lambda=\frac{\sigma_s^2}{T^2},
\]

and sample `log(T)` and `logit(lambda)`, reconstructing

\[
\sigma_s=T\sqrt{\lambda},
\qquad
\tau=T\sqrt{1-\lambda}.
\]

This follows the total centered-share variance and the weakly identified
quality-versus-noise decomposition directly. It must remain an exact coordinate
change: evaluate the existing priors on reconstructed \(\sigma_s,\tau\), retain
every normalizer, and include the correct Jacobian. Independent derivation gives,
for \(\ell=\log\sigma_s\), \(x=\log T\), and
\(y=\operatorname{logit}(\lambda)\),

\[
\left|\frac{\partial(\ell,\tau)}{\partial(x,y)}\right|
=\frac{\tau}{2}.
\]

Thus an implementation with unconstrained `x,y` must evaluate the old density
at the reconstructed \((\ell,\tau)\) and add \(\log\tau-\log 2\) exactly once.
The existing normal prior on \(\ell\) and normalized lower-truncated Normal
prior on \(\tau\) remain unchanged. The rho interaction also requires an audit
because the joint Gaussian likelihood uses the conditional
slope \(\rho\sigma_m/\sigma_s\) and conditional variance
\(\sigma_m^2(1-\rho^2)\). No implementation should begin until the full
Jacobian and the current parameter measure are independently verified.

The benchmark runner remains resumable, writes only compact CSV evidence, and
does not save large fit objects. The historical invocation is retained for
reproducibility:

```text
S0_BENCH_DATA=/path/to/sc_year_cross_section.rds \
S0_BENCH_OUTPUT=/path/to/s0_joint_rotation_sc_bertrand.csv \
S0_BENCH_CASES=market_legacy_nc,market_rotated_nc,market_joint_nc \
S0_BENCH_SEEDS=19091,29091 \
S0_BENCH_CHAINS=4 S0_BENCH_CORES=4 \
S0_BENCH_ITER=2000 S0_BENCH_WARMUP=1000 \
S0_BENCH_ADAPT_DELTA=0.99 S0_BENCH_MAX_TREEDEPTH=13 \
S0_BENCH_SUPPLY=bertrand \
Rscript doc/s0_joint_rotation_hmc_benchmark.R
```

The RDS must contain the selected rows with columns `prices`, `quantities`,
`margins`, `market`, `year`, and `group`. A Cournot run changes only
`S0_BENCH_SUPPLY=cournot`. Alpha pooling can be activated separately with
`S0_BENCH_ALPHA_POOL=true`; it was false for this s0-only gate.

## Decision

**B. Diagnosis confirmed, but the candidate does not fully fix the exact-data
sampler failure.** Keep the exact joint implementation and its evidence as an
experimental coordinate, but do not make it the default and do not submit the
full upper ladder. Audit and test the total-scale/variance-fraction coordinate
next. The existing likelihood, priors, rho, collapsed quality model, and sampler
gate remain unchanged.

## Follow-up: collapsed quality/share scale and rho

The prescribed follow-up has now been executed. The historical conditional
collapsed kernel agrees with a dense joint-Gaussian oracle in 120/120 density
checks and 72/72 finite-difference gradient checks, including residual share
scales down to 0.001 and \(|\rho|\) up to 0.94. The kernel was not an
algebraic defect.

An exact total-scale/variance-fraction coordinate was then implemented:

\[
T=\sqrt{\tau^2+\sigma_s^2},\qquad
\lambda=\sigma_s^2/T^2.
\]

The implementation retains the old priors on \(\tau\) and
\(\log\sigma_s\) and includes the exact Jacobian. Reference and native Stan
target/gradient checks all pass. It reduced the frozen Supreme Court result
from 12 to six divergences but did not meet the zero-divergence gate.

Finally, the same joint density was evaluated margin-first, conditioning share
residuals on observed margin residuals. This replaces
\(\rho\sigma_m/\sigma_s\) with
\(\rho\sigma_s/\sigma_m\) and uses the existing structured centered-Gaussian
kernel. It agrees with the dense oracle, the old native target, and all 1,760
checked aggregate generated-quantity elements. On the exact fixture it produced
10 and six divergences in two independent sets of 4,000 post-warmup draws. It
therefore does not repair the remaining geometry.

The final conclusion remains **B**. The original market-\(s_0\) cancellation
was real and the joint rotation repairs it, but a second, weaker ridge remains
in the quality-versus-residual-share decomposition with active rho. Neither
the total/fraction coordinates nor an equivalent reversal of Gaussian
conditioning eliminates that ridge on the real fixture. No prior, likelihood,
diagnostic threshold, or sampler setting was changed to obtain these results.
