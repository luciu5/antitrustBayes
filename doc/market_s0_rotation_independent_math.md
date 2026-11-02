# Independent derivation: exact market/outside-share coordinates

This memo is an independent audit of the four Stan programs at repository
HEAD `12b26c0`. It concerns only a change of coordinates. The production
programs were not edited. The paired numerical reference is
`/tmp/market_s0_rotation_reference.R`.

## What the live programs actually compute

All four programs define

\[
 h(u)=\epsilon +(1-2\epsilon)\operatorname{logit}^{-1}(u),
 \qquad \epsilon=10^{-6},
\]

in `bounded_inv_logit`, then set `s0_prob = h(u)` and use

\[
 g(u)=\operatorname{logit}(h(u))
\]

in the share mean. Here `u` is the base outside-share location plus the
market-year covariate term. Thus the effective logit outside share is `g(u)`,
not `u`. This distinction is negligible around ordinary values but is not
zero in the tails and must be retained for exact target equality.

For Logit, the share mean is

\[
 g(u_{mt})+b_m+b_j-\alpha_m p_{jmt}+\cdots .
\]

For CES, the common outside-share part is also `g(u)+b_m`; the source adds the
CES price and outside-price terms and uses

\[
 R_{0,mt}=\operatorname{logit}^{-1}\{g(u_{mt})+
 \log p_{0,mt}-\log \bar p_{mt}\}.
\]

The CES files have no `is_cournot`, `log_kappa`, or Cournot alpha rotation.
Their positive scale is `exp(mu_log_a)` (with its existing market hierarchy).

The existing `gamma_market` (Logit) and `b_market_rot` (CES) values are
transformed quantities. They are `b_market + mean(g(u))` for markets with
more than one market-year and `b_market` for one-year markets. They are not
parameters, do not enter the share target, and are not used by the current
margin calculations. The proposed rotation must not introduce a second
meaning for these names.

## Exact affine rotation

Let `eta` denote the latent *base* outside-share location, before `X_s0`
terms, and let

\[
 c=\eta+b.
\]

The reconstruction is

\[
 \eta=c-b,
 \qquad u=c-b+X_{s0}\beta_{s0},
 \qquad s_0=h(u).
\]

The exact share contribution is therefore

\[
 g(u)+b
 =c+X_{s0}\beta_{s0}+\underbrace{[g(u)-u]}_{\delta(u)}.
\]

The correction `delta(u)` is part of the target. Dropping it and replacing the
share mean by `c + X_s0 beta_s0` changes the model, especially in the tails.
The implementation should expose the effective derived quantity

\[
 c^{\mathrm{eff}}_{mt}=g(u_{mt})+b_m
\]

separately from the raw affine coordinate `c`; they are equal only when the
bounded-link correction is zero.

`X_s0` is added after the base rotation. This preserves the source convention
that market-level pooling is on one base location per market, with the
market-year covariate adjustment applied afterward. Rotating `eta+X beta`
instead would induce a different hierarchy and would not preserve that
information structure.

## Jacobians and the noncentered hierarchy

For an unpooled or centered latent `eta`, the map

\[
 (\eta,b)\longmapsto(c,b)=(\eta+b,b)
\]

has triangular Jacobian with determinant one. The old prior is evaluated at
`eta = c-b`; all normal density constants remain present.

The old noncentered hierarchy is

\[
 \eta=\mu+\sigma z,\qquad z\sim N(0,1).
\]

If the sampler replaces `z` by the actual-level coordinate

\[
 c=\mu+\sigma z+b,
\]

then

\[
 z=\frac{c-b-\mu}{\sigma},\qquad
 \log\left|\frac{\partial z}{\partial c}\right|=-\log\sigma
\]

per latent coordinate. Consequently the new target contains

\[
 \log\phi\!\left(\frac{c-b-\mu}{\sigma}\right)-\log\sigma
 =\log N(c-b\mid\mu,\sigma).
\]

For `n` free latent coordinates this is `-n log(sigma)`. Omitting this term
does not preserve the density with respect to the new `c` measure and will
change bridge targets.

There is an exact unit-Jacobian option that preserves the noncentered raw
style. Define

\[
 w=z+\frac{b}{\sigma},
 \qquad z=w-\frac{b}{\sigma},
 \qquad c=\mu+\sigma w.
\]

The map `(z,b) -> (w,b)` is triangular with determinant one. The target uses
`log phi(w-b/sigma)` with no scale Jacobian, and reconstructs

\[
 \eta=\mu+\sigma w-b,
 \qquad c=\eta+b=\mu+\sigma w.
\]

This coordinate is a standardized shifted version of `c`, not the literal
level `c`; its derived `c` is shown above. It preserves both the original
noncentered standard-normal style and the exact model. The actual-level and
shifted-raw choices are mathematically equivalent after their respective
Jacobians are included.

Recommended sampler design: use direct actual-level `c` for the unpooled and
centered branches, and use shifted raw `w=z+b/sigma` for noncentered branches.
The latter is likelihood-aligned while retaining the native standard-normal
noncentered block and unit determinant. Keep the actual-level `c` construction
in the reference tests as a guard: if production samples it directly, it must
add `-n log(sigma)` exactly once.

The full map remains triangular when `b=mu_b+sigma_b b_raw`, because the
market intercept is retained as a separate coordinate. The existing Stan
constraint Jacobian for positive `sigma` is unchanged and cancels in paired
old/new comparisons; the additional `-n log(sigma)` is the Jacobian of the
new actual-level rotation itself.

## Each supported s0 structure

The number of rotated coordinates is the number of free s0 base effects, not
the number of observations.

| Existing source structure | Base latent | Rotated coordinate | Reconstruction |
|---|---:|---|---|
| Unpooled market-year | `eta_mky` | `c_mky=eta_mky+b_m` | `eta_mky=c_mky-b_m` |
| Hierarchical market-year, centered | `eta_mky ~ N(mu,sigma)` | `c_mky` | `eta_mky=c_mky-b_m` |
| Hierarchical market-year, noncentered | `z_mky` | `w_mky=z_mky+b_m/sigma` (unit-Jacobian choice) or actual `c_mky` | `eta_mky=mu+sigma*w_mky-b_m` or `c_mky-b_m` |
| Hierarchical market, centered | `eta_m ~ N(mu,sigma)` | `c_m=eta_m+b_m` | expand `eta_m=c_m-b_m` to each year |
| Hierarchical market, noncentered | `z_m` | `w_m=z_m+b_m/sigma` or actual `c_m` | expand reconstructed market base to years |
| Single-market market pool | direct `mu_s0_pool` with `b=0` | no-op (`c=eta`) | source direct prior unchanged |

For market-level pooling, a repeated `c_m` must not be represented as `M*T`
independent coordinates. Its repeated market-years are deterministic copies,
and the Jacobian counts `M` free effects. The market-year pool counts `M*T`
free effects. A one-period market still has a valid algebraic rotation but
has no repeated-period ridge to improve.

The single-market branches set `b_market=0`. Unpooled and market-year cases
therefore have no `eta+b` separation to repair. The single-market market-pool
branch uses a direct proper prior on `mu_s0_pool` and has no between-market
scale or standardized market effect; it is a no-op rotation.

## Cournot interaction

Only the Logit programs contain the special Cournot alpha/s0 construction.
With no alpha market pool, they compute

\[
 \bar\ell_s=\frac{1}{MT}\sum_{mt}\log s_{0,mt},\qquad
 \mu_{\log a}=\log\kappa-\bar\ell_s,\qquad
 a=\exp(\mu_{\log a}),
 \qquad \alpha=a/\mathrm{price\_sd}.
\]

With alpha pooling they compute, per market,

\[
 \log a_m=\log\kappa_m-
 \frac{1}{T_m}\sum_t\log s_{0,mt},
 \qquad \alpha_m=\exp(\log a_m)/\mathrm{price\_sd},
\]

and retain the existing `mu_log_a`/`sigma_log_a` prior on the derived
`log_a_m`. The rotation leaves `log_kappa`, `log_kappa_m`, and every other
Cournot coordinate untouched. It reconstructs the same `u`, `s0_prob`, and
therefore the same `a`, `alpha`, and FOC inputs. Ordered-k and Stackelberg
margin scores continue to receive the reconstructed `s0_prob` and alpha.

This is why a Bertrand-style replacement of the Cournot coordinates with a
new alpha/s0 combination would be wrong: the existing kappa rotation already
encodes a separate identification relationship. Only the base s0/intercept
coordinate is changed.

CES supports a Cournot conduct formula in its margin function, but its current
Stan parameterization does not have the Logit `log_kappa` construction. The
same raw affine s0 rotation is valid for CES because its share mean contains
the same `g(u)+b` term; its `R0`, gamma, and CES FOC/margin formulas are
evaluated after reconstructing `s0_prob` exactly.

## Density and gradient validation

The reference script uses normalized `dnorm`/Student-t terms, normalized
half-normal hyperpriors, the source logistic-normal s0 hyperprior values
(`s0_prior_mean=-1.446786403764666`, `s0_prior_sd=0.456234014550216`), and
the normalized low-rank Gaussian quality marginal. It evaluates a synthetic
panel with products, multiple periods, nonzero `X_s0`, missing margins, fixed
`k` states, Logit and CES margin functions, Logit Cournot alpha construction,
both share paths, and both sampled and collapsed Gaussian quality.

The corrected run produced 1,000 paired target checks and 28 finite-difference
gradient chain-rule checks, all passed. The maximum absolute target error was
`1.136868e-13`; the maximum gradient error was `7.835175e-08`. The actual-level
noncentered checks compare

\[
 \log p_{\rm legacy}+\log|\partial\theta_{\rm legacy}/\partial c|
 \quad\text{to}\quad
 \log p_{\rm actual\ c},
\]

while centered and shifted-raw checks have unit Jacobian. The expected-failure
control is retained: replacing `g(u)+b` by raw `c+X beta` at extreme `u` had
absolute log-target error `409.4559`, and is recorded as
`expected_failed_control_raw_c_without_bounded_link_correction`, with absolute
error `405.1276`.

The first exploratory run is preserved under the unique attempt-001 names
`/tmp/market_s0_rotation_reference_attempt001_failures.csv` and
`/tmp/market_s0_rotation_reference_attempt001_failures.log` (the earlier
`initial_failures` copies are retained as well). It exposed and was
corrected for two bookkeeping errors in the reference itself: applying the
actual-level noncentered Jacobian twice, and using `b` when recovering the old
noncentered raw coordinate. Those artifacts are retained as failed controls;
they are not evidence against the derivation.

The final machine-readable results are at
`/tmp/market_s0_rotation_reference_results.csv`.
