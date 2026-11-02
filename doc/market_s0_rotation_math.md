# Market/outside-share coordinate rotation: mathematical audit

Status: derivation and native Stan equivalence gate passed independent review;
experimental production implementation integrated locally. Legacy remains the
default. Frozen paired HMC performance evidence is pending.

Independent reference evidence: 1,000 paired density checks passed with maximum
absolute difference `1.14e-13`; 28 finite-difference chain-rule checks passed
with maximum error `7.84e-8`. A separate root reference gave inverse-map error
`5.55e-17`, determinant error `8.23e-11`, and gradient-chain relative error
`1.21e-9`. The subsequent independent native audit passed 468 cases, including the full
176-case multivariate matrix and 256 single-market no-op cases. Maximum target
difference was `2.84e-14`, and maximum AD/finite-difference relative error was
`6.86e-9`. See `doc/market_s0_rotation_native_validation.md`. This establishes
compiled density/coordinate equivalence; it does not establish HMC improvement.

The initial reference run is retained as a failed attempt. Root review found
double-counted actual-level scale Jacobians, incorrect reconstruction of an old
raw effect, and reference share-mean/sign errors. They were corrected before
approval. The deliberately invalid bounded-link shortcut remains an expected
failure control. No production diagnostic was weakened.
Baseline: `12b26c0ee715319396aa8c185a756847dd80963b` (native statistical
specification 9). This experiment changes coordinates only.

## Actual baseline coordinates and outside-share link

Write \(v_u\) for the outside-share **base location**, indexed by market-year
for unpooled/market-year pooling and by market for market pooling. The existing
market effect is noncentered:

\[
b_m=\mu_b+\sigma_b z_{b,m},\qquad z_{b,m}\sim N(0,1).
\]

For one market, the package instead fixes \(b_m=0\) and omits that hierarchy.
The outside base location is either a directly sampled level \(v_u\), a
centered hierarchical level \(v_u\sim N(\mu_0,\sigma_0^2)\), or
\(v_u=\mu_0+\sigma_0 z_{0,u}\), with \(z_{0,u}\sim N(0,1)\).
The package adds optional observed covariates afterwards:

\[
t_{mt}=v_{u(mt)}+X_{0,mt}\beta_0.
\]

The actual production link is **not** plain inverse logit. It is

\[
B(t)=\epsilon+(1-2\epsilon)\operatorname{logit}^{-1}(t),
\qquad \epsilon=10^{-6},\qquad s_{0,mt}=B(t_{mt}).
\]

Thus the effective outside-share logit in the share equation is
\(h(t)=\operatorname{logit}(B(t))\), rather than \(t\). Removing this
protection would change the baseline model. All four Stan programs must retain
it exactly.

## Corrected rotation preserving existing noncentered safeguards

For directly sampled outside levels (including centered pooling), sample

\[
q_u=v_u+b_{m(u)},\qquad v_u=q_u-b_{m(u)}.
\]

For a noncentered outside hierarchy, retain the standardized scale by sampling

\[
w_u=z_{0,u}+\frac{b_{m(u)}}{\sigma_0},\qquad
z_{0,u}=w_u-\frac{b_{m(u)}}{\sigma_0}.
\]

Its implied likelihood-aligned base level is

\[
q_u=\mu_0+\sigma_0w_u=v_u+b_{m(u)},\qquad
v_u=\mu_0+\sigma_0w_u-b_{m(u)}.
\]

This is a noncentered representation of the rotated level. It preserves the
existing standardized hierarchy instead of silently switching a requested
noncentered hierarchy to centered level coordinates. Whether this improves
geometry is an empirical question, especially near small \(\sigma_0\).

Unpooled/year pooling has one free rotated coordinate per market-year. Market
pooling has one per market; adding free period-specific coordinates there would
change the model and is prohibited. Covariate coefficients remain unchanged.
For the single-market normalization \(b=0\), the rotation is a documented
no-op: retain the original sampling coordinates and dimensions.

## Jacobian and normalized target equality

Keep the original market-effect coordinates and all hyperparameters. The direct
map \((v,b)\leftrightarrow(q,b)\) is triangular with diagonal one. The
standardized map \((z_0,\mu_0,\sigma_0,z_b,\mu_b,\sigma_b)\leftrightarrow
(w,\mu_0,\sigma_0,z_b,\mu_b,\sigma_b)\) is also triangular with diagonal
one: dependence on \(b/\sigma_0\) occurs only off the diagonal. Therefore

\[
\left|\det\frac{\partial\theta_{legacy}}
 {\partial\theta_{rotated}}\right|=1.
\]

All existing bounded/positive parameter transformations are unchanged, so the
same equality holds in the corresponding Stan unconstrained coordinates.
Evaluate the original proper priors on the reconstructed values:

\[
\begin{aligned}
\text{unpooled: }&\quad \log N(q-b\mid\mu_{fixed},\sigma_{fixed}^2),\\
\text{centered pool: }&\quad \log N(q-b\mid\mu_0,\sigma_0^2),\\
\text{noncentered pool: }&\quad \log N(w-b/\sigma_0\mid0,1).
\end{aligned}
\]

No new independent prior on \(q\) or \(w\) is permissible. Keep every prior
normalizer, including the existing positive-truncated scale priors, and all
likelihood constants. Native log targets at mapped states must agree directly,
without an additional scale Jacobian, under this standardized-shift design.

A distinct design replacing the old raw \(z_0\) with actual level \(q\)
would have determinant \(\sigma_0^{-d}\) in the inverse map and need
\(-d\log\sigma_0\). That design is not the initial implementation. Treat
omitting that term as a failed control, not an equivalent implementation.

## Exact share equation and bounded-link correction

The effective sum requested in the motivating model remains available as a
derived economic quantity:

\[
c^{effective}_{mt}=h(t_{mt})+b_m.
\]

For the sampled rotated base level the exact share intercept is

\[
h(t_{mt})+b_m
=q_{u(mt)}+X_{0,mt}\beta_0+
\underbrace{[h(t_{mt})-t_{mt}]}_{\text{existing bounded-link correction}}.
\]

The correction must never be dropped, including near the outside-share bounds.
The transparent implementation can simply retain the original expression
\(h(t_{mt})+b_m\), after reconstructing \(t\) from the rotated coordinates.
This changes the sampler's inputs while keeping the existing likelihood
calculation. In the interior the sampled base sum aligns closely with the
effective likelihood sum; exact effective-logit coordinates would require a
nonlinear map with state-dependent bounded support and a separate Jacobian.
That alternative is not assumed equivalent to the unbounded affine proposal.

## Logit conduct and Cournot composition

Bertrand keeps its existing price-coefficient coordinates. Cournot currently
uses, for the common coefficient,

\[
\mu_{\log a}=\log\kappa-\overline{\log s_0},
\]

or, with market coefficient pooling,

\[
\log a_m=\log\kappa_m-\overline{\log s_0}_m.
\]

The original prior is on the implied \(\mu_{\log a}\)/\(\log a_m\).
Keep \(\log\kappa\), its market vector, and their prior evaluations unchanged.
At mapped outside states the same bounded \(s_0\) and hence the same
\(a_m\) and \(\alpha_m=a_m/price\_sd\) are reconstructed. The composition
remains triangular with unit determinant because the kappa coordinates are
unchanged. The coordinate change does not undo that Cournot rotation.

All supply equations continue to use reconstructed \(s_0\); no conduct FOC,
firm share, role statistic, or ordered-k prior is altered.

## CES treatment

The current CES programs have the same additive \(h(t)+b_m\) share intercept,
so the same coordinate change is justified. They store a positive standardized
slope \(a_m\) and structural elasticity
\(\gamma_m=1+a_m/log\_price\_sd\); there is no Logit Cournot kappa coordinate
in those programs. Keep their coefficient hierarchy unchanged. Keep the existing
quantity-to-revenue outside-share transformation:

\[
R_{0,mt}=\operatorname{logit}^{-1}\!\left(
 h(t_{mt})+\log p_{0,mt}-\log P^{inside}_{mt}\right).
\]

The CES price normalization, outside-price term in mean demand, and all margin
equations remain unchanged.

## Quality likelihood, rho, and ordered k

Both sampled-quality likelihood families and the collapsed Gaussian block use
the same reconstructed mean and unchanged quality priors/covariance. Rho,
missing-margin matching, and supported-combination restrictions are unchanged.
The share density still occurs once; complete market margin states are still
marginalized with the original normalized ordered-k probabilities.

## Gradient identity and diagnostic gate

Let \(J\) be the Jacobian of the inverse mapping from new to old unconstrained
coordinates. Since the log determinant is zero,

\[
\log p_{rot}(u)=\log p_{legacy}(f(u)),\qquad
\nabla_u\log p_{rot}=J^T\nabla\log p_{legacy}.
\]

Independent R references, finite differences, and production Stan automatic
gradients must verify this identity before HMC comparisons. Stress the bounded
link tails, covariates, small pooled scales, all hierarchy modes, and Cournot
composition; retain failed checks. Algebraic target equality does not establish
sampler improvement or empirical policy equivalence.

## Second-stage panel decomposition: not implemented

An orthogonal mean/deviation basis could separate market-level and time-level
share information. A literal average plus orthonormal deviations has a scale
Jacobian; using \(\sqrt{T_m}\bar\eta_m\) as the mean coordinate instead gives
an orthonormal transformation. Market pooling permits no new free deviations.
This requires a separate explicit derivation and experiment after the simple
rotation; it is not part of the initial production change.
# Geometry caveat for shifted noncentered coordinates

The shifted standardized coordinate w=z+b/sigma_0 preserves the original
noncentered hierarchy and its normalized density under a unit-Jacobian map.
It is nevertheless not independent of the hierarchy scale under the induced
prior. If sigma_0 approaches zero while b remains nonzero, w becomes large.
Thus exact equivalence does not guarantee improved geometry for pooled modes.
The primary unpooled experiment cannot establish performance for automatically
selected, noncentered market pooling. Those configurations require their own
paired sampler evidence before any general default replacement.

An eventual scale-dependent orthogonal rotation of the two Gaussian market
effects could avoid this particular division by a small scale, but is a separate
coordinate experiment. The panel mean/deviation decomposition requested as a
second stage is likewise deferred. Neither is implemented to rescue a failed
prespecified primary comparison.
