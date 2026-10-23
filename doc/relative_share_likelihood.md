# Normalized relative-share Gaussian calculation

This contract concerns the relative-share likelihood, called Option B. It does
not assert posterior recovery or sampler performance. The absolute-share
likelihood, Option A, retains its separate target and observation measure.

Let `cell` identify each observed market-year equilibrium. For cell c with J_c
observations, choose an orthonormal J_c by (J_c - 1) contrast basis. Assemble
these into Q, and let E contain the cell indicators. Then Q'E = 0 and Q'Q = I.
The observed demand coordinates are Q'log(r), with dimension N - C. Their
density is relative to ordinary Euclidean measure on those coordinates.

Z maps observations to the existing product-quality indices. The existing
centering projector P removes one mean per quality block. Each cell lies
within one quality block, so Q'ZPZ'Q = Q'ZZ'Q. The proper quality prior remains
N(0, tau^2 P) on its centered subspace; no unidentified intercept parameter or
improper intercept prior is added.

For Gaussian shares and margins, first condition the share errors on the
observed margin residuals. Only rows with observed margins receive a correction:

```
u_i = log(r_i) - v_i - rho * sigma_s_i / sigma_m_i * margin_residual_i
D_i = sigma_s_i^2 * (1 - rho^2)
```

Other rows use `u_i = log(r_i) - v_i` and `D_i = sigma_s_i^2`. The full joint
log likelihood adds the normalized Gaussian density of every observed margin.
No missing margin is imputed. For share-only conditioning, use the unconditional
share variances and omit the margin correction. The mean v includes the current
price-sign, units, outside-price and deterministic conduct terms, excluding any
cell-constant level. Both input- and output-market signs are retained.

Set A = D + tau^2 ZZ' and S = E'A^-1E. The contrast determinant and quadratic are

```
log|Q'AQ| = log|A| + log|S| - sum_c log(J_c)
c_hat = S^-1 E'A^-1 u
quadratic = (u - E c_hat)' A^-1 (u - E c_hat)
```

Thus the normalized contrast log density is minus one half of
`(N-C)*log(2*pi) + log|Q'AQ| + quadratic`. Profiling c_hat is only a numerical
identity. Omitting its determinant correction changes model evidence.

## Structured evaluation and numerical stability

Write w_i = 1/D_i, t_g = sum_{i:g} w_i and B_gc = sum_{i:g,c} w_i. Then

```
h_g = tau^2 / (1 + tau^2 t_g)
a_g = t_g / (1 + tau^2 t_g)
log|A| = sum_i log(D_i) + sum_g log1p(tau^2 t_g)
```

S is a positive diagonal ridge plus a graph Laplacian. Its ridge at cell c is
`sum_g a_g B_gc/t_g`; its edge between c and d is `sum_g h_g B_gc B_gd`.
Groups with no active share observations contribute no ridge or edge. The R
implementation stores observation/group sufficient statistics and a cell-sized
factor; it forms neither observation nor quality covariance matrices.

For ordinary conditioning S uses Cholesky. When its small ridge would be lost
against much larger Laplacian entries, QR of positive square-root ridge/edge
rows retains that direction. This changes no variances or target and adds no
jitter. Log determinants use the factor's diagonal. Reciprocal-sum and log-scale
calculations avoid unnecessary overflow in the group shrinkage terms.

The Stan evaluator additionally rotates the cell factor into a data-determined
orthonormal basis: one constant column for each connected component of cells
sharing quality groups, followed by within-component Helmert contrasts. Edge
rows are exactly zero in constant columns. The constant right-hand sides are
computed analytically as sums of group means times `t_g/(1+tau^2*t_g)`, avoiding
cancellation of large weighted residuals. The orthogonal map has determinant
magnitude one, so this improves arithmetic without changing the density,
variances or Jacobians. Isolated cells and tau zero retain their exact limits.

Before the solve, u is centered arithmetically within cells. The quadratic is
evaluated as positive weighted within-quality-group deviations plus
`a_g * weighted_group_mean^2`; it is not a difference of large quadratic forms.
Cells with one observation have no share contrast and contribute exactly zero
share log density. Their original rows, margins, quality indices and prior
centering membership remain meaningful and are retained by the caller.

## Conditional quality recovery

With H = diag(h), the auxiliary uncentered quality posterior is Gaussian with

```
mean = H (Z' W u - B c_hat)
covariance = H + H B S^-1 B' H
```

The actual centered quality has mean P times that mean and covariance P times
that covariance times P. Conditional variances use positive diagonal and
cell-factor terms, without forming a quality covariance matrix. An unused or
unidentified quality direction retains its proper prior uncertainty.

An exact conditional draw can instead use the Matheron identity. Draw centered
prior quality b0 and independent error e0 with covariance D, then return

```
b0 + tau^2 P Z' R_A (u - Z b0 - e0)
R_A = A^-1 - A^-1 E S^-1 E' A^-1
```

R_A is applied through the same structured cell solve. This preserves posterior
quality correlations and centering. Tau zero returns exactly zero quality.
Derived cell intercepts must use this same quality realization, rather than
quality zero or a plug-in conditional mean, because log-sum-exp is nonlinear.

## Exact-k and comparison contracts

One ordered k is shared across the periods of its market. Its normalized joint
likelihood is the common marginal contrast density times a prior-weighted
mixture of conditional margin densities. Equivalently, each component uses the
margin-first Gaussian factorization above. The original k prior enters once.
Conditional-k fits use singleton prior1 and restore the declared original prior
only in `bayes_combine_conditional_k()`. The existing restriction on global
collapsed quality with active rho and multi-market exact-k remains.

Evidence comparisons require one observation basis throughout. Relative and
absolute shares use different reference measures; raw evidence levels must not
be combined. Within one basis, sampler backend is computational provenance and
can vary by model or k component. Predictive scoring uses complete market-year
cells and whole shared-k units; observation-level relative-share LOO is invalid.
Existing collapsed-quality scoring restrictions remain.

## Internal R interface

- `.relative_share_gaussian(residual, variance, quality_index, quality_block,
  cell, tau, return_moments=FALSE)` returns `log_density` and `dimension`.
  With moments requested it also returns `conditional_mean` and
  `conditional_variance`, one value per original quality index.
- `.relative_share_quality_draw(..., n_draws=1L, seed=NULL)` returns a draws by
  quality matrix. A supplied seed preserves the caller's global RNG state.
- `.relative_share_precision(x, variance, quality_index, quality_block, cell,
  tau)` applies R_A to a vector or an N-row matrix. It supports predictive
  conditional means and covariance calculations without duplicating the
  factorization mathematics.
- `.relative_share_basis(cell)` constructs the orthonormal basis for small
  deterministic checks; production likelihood evaluation does not need it.

`quality_index` has one integer index per observation. `quality_block` has one
label per quality index, including indices with no active share contrast.
`cell` has one label per observation. `variance` is one positive variance or N
positive variances; tau is one nonnegative economic quality SD across blocks.
The margin-first correction and margin log density are the caller's responsibility.
