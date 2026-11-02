# Market-specific persistent quality: mathematical specification

This is a substantive prior change for multiple markets. The historical global
option remains a distinct model. No demand, supply, ownership, or welfare
equation changes with quality scope.

## Quality coordinates and normalization

Let observations retain their existing economic group and market identifiers.
Under market scope the quality unit is `(market, group)`, persistent across
periods. For market m, an orthonormal Helmert matrix C_m has G_m rows and G_m-1
columns, C_m'C_m=I and C_m C_m'=I-11'/G_m. Define

\[
b_m=C_m\beta_m,\qquad \beta_m\mid\tau\sim N(0,\tau^2 I).
\]

The beta blocks are independent conditional on the common tau. The prior is a
normalized density on orthonormal contrast coordinates, rather than a density
on an unconstrained G_m-dimensional vector. There is no extra coordinate
Jacobian. A one-group block has no coordinates and b_m=0. With one market this
is exactly the historical global Helmert prior, including its normalization.
The existing proper prior on tau is retained.

## Joint Gaussian integration

Write r_s=y_s-mu_s and r_m(k)=y_m-mu_m(k), retaining the package's actual
log-share and log-inverse-margin outcomes and standardization. A margin row
is matched to one share row i using margin_idx; missing margins introduce no
innovation row. For common residual scales define

\[
\lambda=\rho\sigma_m/\sigma_s,
\qquad v_m=\sigma_m^2(1-\rho^2).
\]

The following augmented observations have independent Gaussian noise
conditional on quality:

\[
z_i=r_{s,i},\quad x_i=1,\quad v_i=\sigma_s^2;
\qquad
z_{N+q}=r_{m,q}(k)-\lambda r_{s,i(q)},\quad
x_{N+q}=-\lambda,\quad v_{N+q}=v_m.
\]

Here z_a=x_a b_{g(a),m(a)}+e_a. The triangular transformation of observed
residuals has determinant one. Quality affects shares alone before this
transformation. This representation therefore preserves the original
observation-level share/margin covariance and every density constant.

For one market, form group sufficient statistics

\[
a_g=\sum_{a:g(a)=g}x_a^2/v_a,\quad
c_g=\sum_{a:g(a)=g}x_a z_a/v_a,\quad
W=\sum_a(z_a-x_a c_{g(a)}/a_{g(a)})^2/v_a.
\]

Terms with a_g=0 use zero fitted mean. Use stable incremental weighted
residual sums in production. Put w_g=(1+tau^2 a_g)^{-1}. Exact integration
of the centered contrast prior gives

\[
q=W+\sum_{g:a_g>0}(c_g^2/a_g)w_g
 +\tau^2\frac{(\sum_g c_g w_g)^2}{\sum_g w_g},
\]

\[
\log|\Sigma_m|=\sum_a\log v_a+
\sum_g\log(1+\tau^2a_g)+\log\left(\frac{\sum_gw_g}{G_m}\right),
\]

\[
\log L_m(k)=-\tfrac12\{(N_m+Q_m)\log(2\pi)+\log|\Sigma_m|+q\}.
\]

These determinant and precision identities retain the centering correction.
For tau=0 they reduce to the original Gaussian observation density. For
G_m=1 the quality contribution vanishes for every tau.

## Exact ordered-k factorization

Conditional on all common parameters, market-quality priors and residual
blocks are independent across markets. Consequently

\[
p(D\mid\theta)=\prod_m\sum_k\pi_{mk}L_m(k).
\]

One k is shared by the whole market across periods. Common shrinkage tau and
other common parameters do not invalidate this conditional factorization.
After averaging over common posterior parameters, market roles can still be
dependent. Joint posterior restrictions therefore average the product of
market probabilities within each structural draw; multiplying separately
averaged market probabilities is generally incorrect.
Global quality does invalidate it when quality connects markets and rho is
active; that legacy collapsed exact-k combination remains explicitly
unsupported. No fractional role or observation-level k mixture is introduced.

Share-only sufficient statistics can be accumulated once per market. Each
component adds only its matched margin innovations. Subtracting that market's
share-only log density from each component and adding it once to the target
is algebraically identical and preserves existing conditional component-score
semantics. Marginalization sums log_sum_exp contributions across markets.
Production stores only incidence maps and sufficient statistics, never a dense
global observation covariance or the Cartesian product of k states.

## Posterior quality and policy realizations

For a realized k_m at a given common-parameter draw, define

\[
d_g=\tau^2/(1+\tau^2a_g),\qquad u_g=d_gc_g.
\]

The conditional centered posterior has

\[
E[b_g]=u_g-d_g\frac{\sum_hu_h}{\sum_hd_h},\qquad
\operatorname{Cov}(b)=\operatorname{diag}(d)-dd'/\sum_hd_h.
\]

Draw independent a_g~N(u_g,d_g), then project
b=a-d sum(a)/sum(d). This produces the exact conditional Gaussian without
a dense covariance factorization. Singleton and tau=0 blocks are identically
zero. For exact k, first draw each market's k from its conditional component
probabilities, then draw quality conditional on that same state. Averaging
over k yields the proper mixture; a MAP k is not a replacement.

A reusable policy draw plan stores the structural draw, market-specific k
realizations, and conditional quality realization together. It holds that
quality fixed through pre- and post-policy equilibrium. This affects no
existing cost recovery, normalization, or welfare formula.

## Scope and identification

Market scope allows persistent local quality, without a common national firm
effect. A future g_j+c_jm extension would be a separate richer-information
model and reintroduce cross-market dependence. Sparse single-period markets
need not separately identify tau and share-error variance. Collapsing preserves
that uncertainty; it does not create variance-component identification.

For a balanced panel with G groups observed T times each and known share
locations, the share-only covariance has eigenvalues

\[
 \sigma_s^2+T\tau^2\quad(G-1\text{ directions}),\qquad
 \sigma_s^2\quad(GT-G+1\text{ directions}).
\]

Repeated periods therefore add noise-only directions that can distinguish
transitory noise from persistent quality. In a one-period cross section there
is just one such direction; its mean also depends on the outside-share and
intercept specification. These fixed-parameter covariance identities do not
establish nuisance-variance identification in the complete posterior.

Conditional on structural parameters and k, the information for an individual
quality loading before imposing the centered constraint is

\[
 a_g=\frac{n_{g,\mathrm{missing}}}{\sigma_s^2}
      +\frac{n_{g,\mathrm{matched}}}{\sigma_s^2(1-\rho^2)}.
\]

The dimensionless quantity tau^2*a_g compares data curvature with prior
curvature. With equal observation counts and rho=0 it becomes
T*tau^2/sigma_s^2; all matched margins replace sigma_s^2 by
sigma_s^2*(1-rho^2). Unbalanced centering couples quality contrasts, and
unknown structural parameters add posterior dependence. This is an
information diagnostic, not a universal sampled-versus-collapsed runtime
threshold. More markets with the same sparse panel design do not increase
each local group's observation count.

Implementation and defaults remain gated on independently executed density,
gradient, sampler, bridge, policy, and PNB validation. This derivation alone
does not establish empirical performance.

## Prior-predictive interpretation

The retained baseline scale prior is positive-truncated Normal(1, 1).
Prespecified simulations (seed 1729, 20,000 draws) give an approximate central
95% interval of 0.01 to 110 for a pair's relative share odds
`exp(b_g-b_h)`, across representative markets with 2, 6, and 39 groups. This is
a broad prior under the new local-quality interpretation; no hyperprior has
been changed. Because a centered pair difference has conditional variance
`2*tau^2`, its odds expectation is `E[exp(tau^2)]`, which is infinite under this
scale prior. The distribution is proper, but that moment does not exist. Report
quantiles rather than a misleading finite Monte Carlo odds mean. The initial
diagnostic's empirical means and the corrected output are both retained.
