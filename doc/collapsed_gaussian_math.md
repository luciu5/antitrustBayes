# Collapsed Gaussian share likelihood: approved derivation

This note is the signed mathematical specification for the collapsed-Gaussian experiment. It
is based on the current specification-8 Stan graphs (logit_supply.stan,
ces_supply.stan, and their Stackelberg counterparts), whose quality block is
the centered orthonormal Helmert representation. It does not change a
production model.

## 1. Current coordinates and observation map

Let \(G\) be the number of global product groups and \(N\) the number of
observations. The Stan transformed-data matrix \(C\), called Q in the Stan
source, is \(G\times(G-1)\). Its \(r\)-th column is

\[
C_{jr}=\begin{cases}
1/\sqrt{r(r+1)},&j\le r,\\
-r/\sqrt{r(r+1)},&j=r+1,\\
0,&j>r+1.
\end{cases}
\]

For \(G>1\),

\[
C'C=I_{G-1},\qquad C'\mathbf 1_G=0,\qquad CC'=P_G,
\quad P_G=I_G-\mathbf1_G\mathbf1_G'/G.
\]

For \(G=1\), \(C\) has one row and zero columns. Products are mapped to groups
by the fixed one-hot matrix \(Z\in\{0,1\}^{N\times G}\),

\[
Z_{ig}=1\{\mathrm{group}(i)=g\},\qquad H=ZC.
\]

Thus the product-level quality in the Stan code, Q multiplied by the sampled
centered group vector and then indexed by group, is exactly \(H\beta\), where
\(\beta\) is that sampled vector. This distinction matters for multiproduct
firms: several rows of \(Z\), and therefore \(H\), are identical for products
of one firm. The map is not an \(N\times N\) identity matrix in that case.

For the collapsed experiment write the share equation as

\[
y=\mu(\theta)+H\beta+e_s,\qquad
e_s\sim N_N(0,\sigma_s^2I_N),\qquad
\beta\mid\tau\sim N_{G-1}(0,\tau^2I_{G-1}).
\]

Here \(y_i=\log(s_{i,\mathrm{inside}})\). The vector \(\mu(\theta)\) contains
the current package's non-quality terms: the logit of outside share, market
effects, the signed standardized-price coefficient, and any other term active
in the chosen specification. Nothing in the demand or conduct equations is
redefined. In the current Stan graphs \(\tau\) is `sigma_b_strat` and
\(\sigma_s\) is `sigma_share_market[m]`; the first experiment uses one common
\(\sigma_s\), as requested.

The market-intercept, outside-share, price-sign, and price-standardization
normalizations remain exactly those in Stan. In particular, for a single
market b_market is fixed to zero; for multiple markets it remains a separate
effect. A quality draw is centered:

\[
\mathbf1_G'C\beta=0,\qquad
\operatorname{Var}(ZC\beta\mid\tau)=\tau^2ZP_GZ'.
\]

The same construction applies to Logit and CES. Only the formula for
\(\mu(\theta)\) differs between their demand equations.

## 2. Expanded and collapsed share densities

Conditional on \(\beta\), the fully normalized expanded share density is

\[
p(y\mid\beta,\theta,\sigma_s)
=(2\pi\sigma_s^2)^{-N/2}
\exp\left[-\frac{(y-\mu-H\beta)'(y-\mu-H\beta)}
{2\sigma_s^2}\right].
\]

The normalized quality-coordinate density in the current parameterization is

\[
p(\beta\mid\tau)
=(2\pi\tau^2)^{-(G-1)/2}
\exp\left[-\frac{\beta'\beta}{2\tau^2}\right].
\]

There is no additional quality Jacobian in this expression. The current
parameter is already the centered contrast coordinate with its Normal scale;
the old raw-coordinate derivation in doc/stackelberg/quality_coordinates.md
establishes the equivalence, including the discarded mean direction.

Integrating \(\beta\) gives

\[
y\mid\theta,\tau,\sigma_s\sim N_N(\mu,\Sigma_s),\qquad
\boxed{\Sigma_s=\sigma_s^2I_N+\tau^2HH'
=\sigma_s^2I_N+\tau^2ZP_GZ'.}
\]

Therefore the collapsed log density is

\[
\log p(y\mid\theta,\tau,\sigma_s)
=-\frac N2\log(2\pi)-\frac12\log|\Sigma_s|
-\frac12(y-\mu)'\Sigma_s^{-1}(y-\mu).
\]

This is an exact identity, not a new approximation. Let \(r=y-\mu\),

\[
V=(\tau^{-2}I+\sigma_s^{-2}H'H)^{-1},\qquad
m=VH'r/\sigma_s^2.
\]

Completing the square gives

\[
\beta\mid y,\theta,\tau,\sigma_s\sim N_{G-1}(m,V).
\]

For \(G=1\), the contrast dimension is zero, \(H\) has zero columns, the
quality contribution vanishes, and the collapsed density is
\(N_N(\mu,\sigma_s^2I)\). This is the correct zero-dimensional limit.

## 3. Numerically stable evaluation

An \(N\times N\) Cholesky is retained only as a small-\(N\) reference. It is
not the production path. The one-hot firm incidence and the centering
projection yield an exact linear-time likelihood.

The more general form needed for the rho calculation has pseudo-observations
\(z_l\), diagonal residual variances \(v_l>0\), group indices \(g(l)\), and
one nonzero group loading \(x_l\) per row. Let \(X\) be the sparse matrix with
\(X_{l,g(l)}=x_l\), let \(D=\operatorname{diag}(v_l)\), and let
\(u=X\mathbf 1_G=(x_1,\ldots,x_L)'\). Then

\[
K=D+\tau^2XP_GX'
 =\underbrace{D+\tau^2XX'}_{V}
 -\frac{\tau^2}{G}uu'.
\]

The matrix \(V\) is block diagonal by group and each block is diagonal plus
rank one. Define

\[
a_g=\sum_{l:g(l)=g}\frac{x_l^2}{v_l},\qquad
c_g=\sum_{l:g(l)=g}\frac{x_lz_l}{v_l},\qquad
e=\sum_l\frac{z_l^2}{v_l},
\]

\[
w_g=\frac{1}{1+\tau^2a_g},\qquad W=\sum_gw_g.
\]

The block determinant lemma and a final Sherman--Morrison downdate give

\[
\boxed{
\log|K|=
\sum_l\log v_l+
\sum_g\log(1+\tau^2a_g)+
\log\left(\frac{W}{G}\right).}
\]

The last factor is evaluated as \(W/G\), rather than as
\(1-(\tau^2/G)\sum_ga_g/(1+\tau^2a_g)\), to avoid cancellation. The exact
quadratic form is

\[
q_0=e-\tau^2\sum_gc_g^2w_g,\qquad
b=\sum_gc_gw_g,
\]

The subtraction in that algebraic expression is **not** the production
evaluation: it loses precision as share noise approaches zero. Define the
weighted within-group regression residual sum

\[
T=\sum_l\frac{(z_l-x_l\bar b_{g(l)})^2}{v_l},\qquad
\bar b_g=c_g/a_g\quad(a_g>0).
\]

Groups with \(a_g=0\) have \(\bar b_g=0\). Since
\(e=T+\sum_ga_g\bar b_g^2\), an exactly equivalent positive-term formula is

\[
\boxed{q_0=T+\sum_ga_g\bar b_g^2w_g.}
\]

Production accumulates \(T\) with a weighted regression update. For a new
row in a group whose previous statistics are \(A>0,C\), write
\(h=x^2/v\) and update

\[
T\leftarrow T+\frac{A}{A+h}\frac{(z-xC/A)^2}{v},\quad
a_g\leftarrow A+h,\quad c_g\leftarrow C+xz/v.
\]

If \(A=0\), the first nonzero-loading row has zero within residual; a
zero-loading row instead contributes \(z^2/v\). This update permits zero and
negative loadings without dividing by \(x\). The rho pseudo-rows use this
same rule and extend the cached share statistics. The independent R reference
uses a separate two-pass residual calculation. No negative quadratic is
clamped and no statistical approximation is introduced.

\[
\boxed{z'K^{-1}z=q_0+\frac{\tau^2b^2}{W}.}
\]

Consequently

\[
\boxed{
\log p(z)=-\frac12\left[
L\log(2\pi)+\log|K|+z'K^{-1}z
\right].}
\]

This exact density requires \(O(L+G)\) work and \(O(G)\) group-statistic memory. The current Stan wrapper also stores the linear-size pseudo-observation vectors and automatic-differentiation tape, so total production storage is \(O(L+G)\), rather than quadratic. It
constructs neither \(K\) nor a dense contrast matrix. The \(G=1\) limit is
stable: centering removes the sole quality coordinate and the expression
reduces exactly to the diagonal Gaussian density.

For the share-only likelihood, set \(L=N\), \(x_i=1\),
\(v_i=\sigma_s^2\), and \(z_i=y_i-\mu_i\). The dense Cholesky formula in
Section 2 and the structured formula above are algebraically identical. The
production gate requires numerical agreement with the dense oracle.

## 4. Optional jointly Gaussian rho block

Let \(q\) margins be observed, with margin residual

\[
y_m=\nu(\theta,k)+e_m,\qquad e_m\sim N_q(0,\sigma_m^2I_q).
\]

The vector \(\nu\) is the current conduct prediction for the relevant exact
ordered-\(k\) state. Let \(R\in\mathbb R^{N\times q}\) encode which share
observation is paired with each margin. With one margin for rows
\(i_1,\ldots,i_q\), \(R_{i_a,a}=1\) and all other entries are zero. This
handles missing margins without inventing observations; a general matching
matrix is also allowed.

Before quality is integrated, assume

\[
\operatorname{Cov}(e_s,e_m)=\rho\sigma_s\sigma_mR,\qquad |\rho|<1.
\]

Quality is independent of both disturbances. Integrating quality changes only
the share-share block:

\[
\boxed{\Sigma=
\begin{pmatrix}
\Sigma_s & B\\
B' & \sigma_m^2I_q
\end{pmatrix},\qquad B=\rho\sigma_s\sigma_mR.}
\]

This is the exact collapsed joint Gaussian density. The conditional factor
useful for computation is

\[
y_m\mid y_s,\theta\sim
N_q\left(\nu+B'\Sigma_s^{-1}(y_s-\mu),
\sigma_m^2I_q-B'\Sigma_s^{-1}B\right).
\]

The production implementation does not form that Schur complement. For each
observed margin define

\[
r_i^s=y_i^s-\mu_i,\qquad
r_i^m=y_i^m-\nu_i,\qquad
\lambda=\rho\frac{\sigma_m}{\sigma_s},
\]

and use the unit-determinant triangular transformation

\[
t_i=r_i^m-\lambda r_i^s.
\]

Conditional on group quality, a share contributes the pseudo-row

\[
z_i=r_i^s,\qquad x_i=1,\qquad v_i=\sigma_s^2,
\]

and an observed transformed margin contributes

\[
z_{N+i}=t_i,\qquad x_{N+i}=-\lambda,\qquad
v_{N+i}=\sigma_m^2(1-\rho^2),
\]

with the same group index as its matched share. These pseudo-row residuals are
independent conditional on quality, so the augmented vector has exactly the
structured form in Section 3. Its normalized joint density is evaluated in
\(O(N+q+G)\), including missing margins, without constructing the joint or
conditional covariance. At \(\rho=0\), margin loadings are zero and the
formula reduces to the collapsed share density times independent Gaussian
margin densities.

If \(\tau=0\), this reduces to the familiar rowwise conditional mean
\(\nu+\rho(\sigma_m/\sigma_s)R'(y_s-\mu)\). With nonzero quality variance,
\(\Sigma_s^{-1}\) generally makes a margin's conditional mean depend on
residuals from other products of the same firm and, through the quality
covariance, on other observations. Reusing the old rowwise conditional
expression after collapsing would therefore be a different model.

For a one-to-one matching matrix, \(|\rho|<1\), \(\sigma_s>0\), and
\(\sigma_m>0\) imply positive definiteness. Since
\(\Sigma_s\succeq\sigma_s^2I\),
\(R'\Sigma_s^{-1}R\preceq\sigma_s^{-2}I_q\), and the Schur complement is

\[
\sigma_m^2[I_q-\rho^2\sigma_s^2R'\Sigma_s^{-1}R]\succ0.
\]

For arbitrary \(R\), the exact condition is that this Schur complement be
positive definite. A conservative sufficient condition is

\[
|\rho|<
\{\sigma_s^2\lambda_{\max}(R'\Sigma_s^{-1}R)\}^{-1/2}.
\]

If \(q=0\), the margin block is empty and the density is the share density.

The current quality centering is global across every group in the fitted data.
After collapse, \(\Sigma_s\) therefore need not be block diagonal by market.
With \(\rho\ne0\), the Schur complement can correlate conditional margin
errors across markets. An exact-k model with a separate \(k_m\) for each market
would then require a joint sum over all market-state combinations. Replacing
that joint sum with separate market-level sums would change the probability
model. The initial implementation therefore supports collapsed nonzero
\(\rho\) for exact-k models only when there is one market. The default
experiment has \(\rho=0\), for which the margin likelihood remains
market-factorized and the existing \(O(MK)\) marginalization is exact.

A Gaussian-share/Student-t-margin hybrid is coherent when \(\rho=0\): its
joint density is the product of the collapsed Gaussian share density and the
existing independent fixed-df Student-t margin densities. Nonzero \(\rho\) for
that hybrid is outside this specification and must be rejected.

## 5. Quality posterior with and without rho

Let

\[
x=\begin{pmatrix}y_s-\mu\\y_m-\nu\end{pmatrix},\qquad
X=\begin{pmatrix}H\\0_{q\times(G-1)}\end{pmatrix},\qquad
K=\begin{pmatrix}\sigma_s^2I_N&B\\B'&\sigma_m^2I_q\end{pmatrix}.
\]

Conditional on \(\beta\), the expanded joint observation model is
\(x\mid\beta\sim N(X\beta,K)\). Hence the exact conditional posterior is

\[
\boxed{
V_\beta=(\tau^{-2}I+X'K^{-1}X)^{-1},\qquad
m_\beta=V_\beta X'K^{-1}x,\qquad
\beta\mid y_s,y_m,\theta\sim N(m_\beta,V_\beta).}
\]

For \(\rho=0\), this reduces to the share-only formulas in Section 2; margins
then provide no information about \(\beta\) conditional on \(\theta\). Recover
group and product quality moments by

\[
E[b_{\mathrm{group}}\mid y]=Cm_\beta,\qquad
\operatorname{Var}(b_{\mathrm{group}}\mid y)=CV_\beta C',
\]
\[
E[b_{\mathrm{product}}\mid y]=ZCm_\beta,\qquad
\operatorname{Var}(b_{\mathrm{product}}\mid y)=ZCV_\beta C'Z'.
\]

These formulas allow posterior quality draws after a collapsed fit. Sample
\(\beta^*\sim N(m_\beta,V_\beta)\), then return \(C\beta^*\) for group
qualities and \(ZC\beta^*\) for product qualities. Do not add a second scale
or Jacobian.

The production reconstruction can avoid a dense \(G\times G\) contrast
calculation. Using the structured pseudo-row statistics from Section 3, first
ignore the centering constraint and define

\[
d_g=(\tau^{-2}+a_g)^{-1},\qquad
\tilde m_g=d_gc_g,\qquad D=\operatorname{diag}(d_g).
\]

The centered posterior is the independent Gaussian
\(N(\tilde m,D)\) conditioned on \(\mathbf1'b=0\):

\[
\boxed{
m_b=\tilde m-D\mathbf1
\frac{\mathbf1'\tilde m}{\mathbf1'D\mathbf1},\qquad
V_b=D-D\mathbf1(\mathbf1'D\mathbf1)^{-1}\mathbf1'D.}
\]

An exact draw requires only independent scalar Normal draws. Draw
\(u_g\sim N(\tilde m_g,d_g)\), then project

\[
b^*=u-D\mathbf1\frac{\mathbf1'u}{\mathbf1'D\mathbf1}.
\]

This is \(O(N+G)\), has exact zero sum, and is algebraically equivalent to the
Helmert-coordinate posterior. With rho, use the augmented pseudo-row
statistics. For exact k, the joint-conditioned quality posterior is a finite
mixture over k-specific margin residuals; at rho zero it is common across k.

## 6. Exact-k compatibility and identification limits

The share likelihood has no \(k\) argument. For exact ordered-\(k\) conduct,
one evaluates it once and combines it with the complete market-level margin
mixture:

\[
\log L_M=\log L_{\mathrm{share}}+
\operatorname{logsumexp}_k\{\log\pi_{Mk}+\log L_{\mathrm{margin},Mk}\}.
\]

The collapsed density is therefore identical across BF/CF/BS/CS states with
the same continuous parameters. It must never be recomputed once per \(k\),
and products of per-product mixtures would be incorrect.

For the one-market rho extension, compute the share sufficient statistics and
share log density once. The transformed margin rows add to \(a_g,c_g,T\), and
their residual terms depend on the k-specific conduct mean. For component
\(k\), evaluate the structured joint density from the augmented sufficient
statistics and set

\[
\log L_{m\mid s,k}=\log L_{s,m,k}-\log L_s.
\]

The market mixture then uses these conditional scores. This preserves the
same factorization while avoiding a repeated pass over the share rows for each
k.

With one observation per firm in one period, the data generally do not identify
\(\tau^2\) separately from \(\sigma_s^2\). The covariance is identified through
the composite \(\tau^2ZP_GZ'+\sigma_s^2I\), while the decomposition is largely
prior-sensitive. More periods and/or repeated markets can identify additional
within-group covariance patterns and pool common variance components; they do
not require a different economic model. The experiment must therefore report
\(\tau/\sigma_s\) prior sensitivity and compare policy posteriors, rather than
claiming that either variance is precisely learned in a sparse cross section.

This covariance formula is also the nested information-ladder architecture.
With a persistent group observed at several dates, its rows of \(Z\) repeat and
the off-diagonal covariance \(\tau^2 ZP_GZ'\) supplies information about
persistent quality that is absent from a one-period cross section.  Additional
markets supply replication for common scale parameters and can use the
package's existing alpha and outside-share pooling switches.  No demand or
conduct equation changes across these cases.  The initial implementation
retains the existing global Helmert normalization exactly; market-block
centering would define a distinct statistical specification and is not an
automatic consequence of adding markets.

When \(\rho\ne0\) and exact \(k\) is marginalized, the conditional quality
posterior is a finite mixture. Conditional on component \(k\), use the
Gaussian formula in Section 5 with margin mean \(\nu_k\), then weight or sample
components using their posterior probabilities. At \(\rho=0\), the quality
posterior is common across every \(k\).

## 7. Status of the reference checks

The accompanying base-R reference script checks the following deterministically:
analytic expanded/collapsed equality, one-dimensional direct integration, a
two-dimensional low-dimensional integration, a multiproduct \(Z\ne I\) map,
centered Helmert identities, conditional posterior moments, and the full
\(\rho\) covariance with missing margin rows. Those are mathematical/numerical
checks only. They do not establish HMC improvement, posterior recovery, or
policy similarity.

An independent lead calculation separately matched the normalized
one-contrast multiproduct integral to \(3.33\times10^{-16}\) and the rho-block
joint/conditional factorization to \(4.44\times10^{-16}\). These checks approve
the algebra for production implementation; they are not sampler evidence.
