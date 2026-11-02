# Sampled-quality contrast matrix: exact structured alternative

This is an audit derivation, not a production change or a new prior.
The existing sampled-quality kernels declare a dense orthonormal contrast
matrix Q with G rows and G-1 columns. Its kth column is

\[
Q_{jk}=\begin{cases}
[k(k+1)]^{-1/2},&j\leq k,\\
-k[k(k+1)]^{-1/2},&j=k+1,\\
0,&j>k+1.
\end{cases}
\]

The actual sampled coordinates remain beta = `b_group_centered`, with the
existing Gaussian prior and scale. For v_k = beta_k/sqrt(k(k+1)),

\[
b_1=\sum_{k=1}^{G-1}v_k,\qquad
b_j=\sum_{k=j}^{G-1}v_k-(j-1)v_{j-1},\quad 2\leq j\leq G.
\]

A descending recurrence computes these suffix sums using O(G) work and
storage. Initialize a=0; for j=G,...,2 compute b_j=a-(j-1)v_(j-1), then
update a=a+v_(j-1). Finally b_1=a. For G=1 the quality effect is zero and
there are no sampled contrast coordinates.

This computes the same Q beta in exact arithmetic. It changes neither sampled
coordinates nor the density measure, so it introduces no prior transformation
or Jacobian. Numerical checks must cover outputs, centering, derivatives,
normalized native targets, and reconstructed demand. Different summation
orders can introduce small floating-point differences.

At G=10,000, dense double storage for Q alone is 799,920,000 bytes. Native
setup and gradient costs must be measured rather than inferred solely from
this size. Sampling benchmarks at N=10,000 must state G as well: a panel with
100 persistent firms and 10,000 observations has a very different contrast
cost from 10,000 single-product firms. The collapsed-quality kernels allocate
zero contrast columns; this derivation does not change their existing O(N+G)
likelihood or quality distribution.

No optimization is accepted on this derivation alone. Independently compare
the recurrence with dense Q on resource-bounded reference fixtures, then
benchmark the native before/after paths with identical parameters and settings.
