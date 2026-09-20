# Collapsed-quality scale rotation

## Scope

This experiment changes sampler coordinates only. It does not change the
Gaussian share-margin likelihood, the centered quality prior, the conduct
equations, or any prior family. The implementation is available through

```r
quality_scale_parameterization = "total_fraction"
```

and is initially restricted to Gaussian shares, analytically collapsed
quality, and a common share residual scale. The historical coordinates remain
the default under `"legacy"`.

## Historical coordinates and measure

Write

\[
\ell=\log\sigma_s,
\qquad
\tau=\sigma_b>0.
\]

The normalized prior density is defined with respect to

\[
d\ell\,d\tau.
\]

The package retains its existing Normal prior on \(\ell\) and its normalized
lower-truncated Normal prior on \(\tau\). No prior is placed directly on the
new coordinates.

For a centered quality block with projector

\[
P=I-\frac{11'}G,
\]

the one-observation-per-group share covariance is

\[
\tau^2P+\sigma_s^2I.
\]

On the \(G-1\) centered contrast directions, its eigenvalue is
\(\tau^2+\sigma_s^2\). This motivates the **contrast total scale**

\[
T=\sqrt{\tau^2+\sigma_s^2}
\]

and share-noise fraction

\[
\lambda=\frac{\sigma_s^2}{T^2}.
\]

Because quality is centered, \(T\) is not generally the marginal standard
deviation of an individual product. With repeated periods it also does not
describe the entire covariance spectrum. It is a sampler coordinate aligned
with the weakly identified variance decomposition, not a new economic
parameter.

## Exact coordinate map

The unconstrained coordinates are

\[
x=\log T,
\qquad
y=\operatorname{logit}(\lambda).
\]

Reconstruction uses stable log-scale expressions:

\[
\ell=x+\frac12\log\{\operatorname{logit}^{-1}(y)\},
\]

\[
\log\tau=x+\frac12\log\{\operatorname{logit}^{-1}(-y)\}.
\]

The inverse map is

\[
x=\frac12\log\{e^{2\ell}+\tau^2\},
\qquad
y=2(\ell-\log\tau).
\]

The Jacobian under the historical measure is

\[
\left|\frac{\partial(\ell,\tau)}{\partial(x,y)}\right|
=\frac{\tau}{2}.
\]

Consequently, an unconstrained implementation in \((x,y)\) evaluates every
historical normalized prior at reconstructed \((\ell,\tau)\) and adds

\[
\log\tau-\log2
\]

exactly once.

The production Stan programs preserve their historical raw parameter schema
for saved-fit bridge compatibility. Under `total_fraction`, the constrained raw
field historically named `sigma_b_strat` stores \(T\), while the unconstrained
raw field historically named `log_sigma_share` stores \(y\). Since Stan adds
\(\log T\) for the positive raw coordinate, the model block adds

\[
\log\tau-\log2-\log T.
\]

The resulting unconstrained target therefore contains precisely
\(\log\tau-\log2\). Economically named transformed quantities
`sigma_b_strat_economic`, `log_sigma_share_economic`, and `sigma_share` expose
the reconstructed historical scales. Posterior-quality extraction uses the
economic alias.

## Active rho and the boundary

For share residuals \(r_s\), observed-margin residuals \(r_m\), matching
matrix \(E\), quality incidence \(Z\), and centered projector \(P\), the
collapsed joint covariance is

\[
\Sigma=
\begin{pmatrix}
\tau^2ZPZ'+\sigma_s^2I & \rho\sigma_s\sigma_mE\\
\rho\sigma_s\sigma_mE' & \sigma_m^2I
\end{pmatrix}.
\]

The historical structured evaluator conditioned margins on shares and contained
the algebraic ratio \(\rho\sigma_m/\sigma_s\). An independent dense covariance
oracle tested sparse and repeated-panel fixtures, missing margins,
\(|\rho|\le .94\), \(\tau\in\{.2,1.5\}\), and
\(\sigma_s\in[.001,.5]\). All 120 density comparisons and 72 gradient
comparisons passed. The maximum log-density error was \(3.91\times10^{-9}\).
That evidence established that the historical evaluator was correct; it did
not establish favorable sampler geometry.

The exact Supreme Court scale-coordinate run still produced 6 divergences in
4,000 post-warmup draws, compared with 12 under the historical scale
coordinates. This triggered the prespecified second-stage experiment using the
exactly equivalent margin-first factorization:

\[
r_m\sim N(0,\sigma_m^2I),
\]

\[
r_s\mid r_m\sim N\left(
\rho\frac{\sigma_s}{\sigma_m}Er_m,
\tau^2ZPZ'+\operatorname{diag}(v_i)
\right),
\]

where \(v_i=\sigma_s^2(1-\rho^2)\) for a row with an observed margin and
\(v_i=\sigma_s^2\) otherwise. This avoids division by \(\sigma_s\) and is
evaluated in \(O(N+G)\) with the existing centered-quality sufficient-statistic
kernel. A separate dense oracle again passed 120/120 density comparisons and
72/72 finite-difference gradient comparisons, including missing margins,
\(|\rho|=.94\), and \(\sigma_s=.001\). The production implementation changed
only after those checks passed.

## Validation rule

Mapped legacy and total/fraction native Stan targets differ by \(-\log2\) in
their respective unconstrained coordinates, exactly as required by the
coordinate-map determinant. Economic transformed quantities must be identical.
Bridge sampling remains valid because each target retains all normalizing and
Jacobian constants.
