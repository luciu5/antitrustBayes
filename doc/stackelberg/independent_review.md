# Independent Stackelberg mathematics review

Date: 2026-09-13  
Scope: `doc/stackelberg_math_audit.md`, `doc/stackelberg/math_check.R`, the
recorded oracle outputs, and the ordered random-boundary proposal.  No
production files were edited.  The independent executable check is
`/tmp/independent_ces_check.R`.

## Finding

I found no substantive algebra error in the supplied hard-role equations.
The Logit Bertrand/Cournot signs, common dollar markups, multiproduct firm
aggregation, CES Bertrand common Lerner margins, and CES Cournot rank-one
inverse Jacobian and response all follow from the stated primitives.  The
CES Cournot check below independently inverts the primitive CES demand and
uses numerical derivatives of that inverse, rather than importing the
oracle's Jacobian formula.

The CES Cournot response coefficient `a_f` can be negative for a very large
follower firm and a small outside revenue share, as the audit correctly
warns.  This does not by itself imply a singular response.  In fact, on the
stated interior domain one can strengthen the denominator gate.  Write
`h = gamma - 1`, `x = r0`, `y = R_f`, `H = 1 + h*x`, and
`b = h/(gamma*H)`.  Then

```
A_f = 1 - m_f = b*(H-y) > 0
D_f = A_f + h*b*y = b*(H + (h-1)*y) > 0
a_f = h*(H - 2*y + gamma*h*x*y/H) / (H*{H + (h-1)*y}).
```

For `0 < x < 1` and `0 <= y <= 1-x`, `a_f > -1`.  One direct proof is that
the numerator of `a_f + 1`, after multiplication by the positive `H`, is
affine in `y`; at `y=0` it is `H*(H+h)`, and at `y=1-x` it is

```
x * [h^3 + 3*h^2 + 2*h + 1
     + (h^3 + h^2 + 2*h)*x + h^2*x^2] > 0.
```

Consequently `C_F = sum_f R_f*a_f > -sum_f R_f > -1`, so
`1 + C_F > 0` whenever the outside revenue share is strictly positive.  The
homogeneous follower-Foc nullspace argument also gives nonsingularity: a
null vector must have a common log-quantity change within each follower
firm, then `x_f = -a_f W` and `(1+C_F)W=0`, forcing every change to zero.
This is a stronger interior result than the audit's conservative
“require ... nonsingular and `1+C_F` nonzero” gate.  Keeping an explicit
finite/regularity check is still sensible numerically, especially near
`r0=0`, where `1+C_F` can be very small and reactions become ill-conditioned.
It remains correct that `a_f` itself need not be positive and that leader
margins/costs must be checked for admissibility.

## Follow-up: no-outside CES boundary

The amended audit and `doc/stackelberg/ces_contract.md` correctly distinguish
the Logit and CES boundaries. With CES outside coefficient zero, the
primitive inverse remains

```
p_j = w*(A_j/q_j)^(1/gamma),
w * sum_j A_j^(1/gamma) q_j^((gamma-1)/gamma) = E,
```

so fixed nominal expenditure identifies `w=E/T` even without outside
expenditure. The initial constructor can nevertheless reject this boundary:
the usual `A_j=(r_j/r0)*(p_j/p0)^(gamma-1)` calibration is undefined at
`r0=0`, and a single positive-cost firm with all expenditure has profit
`E-c*q`, whose finite interior optimum does not exist. This is a valid API
identification/equilibrium gate, rather than a claim that CES inverse demand
itself loses its price scale. The contract's explicit positive outside
revenue-share requirement and its inverse-demand scalar equation are
consistent with this distinction.

## Independent CES Cournot computation

The script uses six products, ownership `A,A,B,B,C,D`, multiproduct leader
`A`, multiproduct follower `B`, and single-product leader/follower `C,D`.
It constructs CES primitives from target prices and revenue shares, recovers
inverse prices by solving the scalar CES price-index equation, calibrates
constant costs from the proposed leader/follower FOCs, solves follower raw
profit gradients with leader quantities fixed, and compares the closed-form
follower response with central finite differences.  The inverse Jacobian is
also compared directly with `numDeriv::jacobian(inv_p, q)`, where `inv_p`
comes from the CES primitive.

| gamma | follower FOC | response error | leader reduced gradient | inverse-Jacobian error | inverse round trip |
|---:|---:|---:|---:|---:|---:|
| 1.5 | 0.0e+00 | 1.58e-08 | 6.70e-09 | 1.43e-08 | 2.78e-17 |
| 2.0 | 1.78e-15 | 2.48e-11 | 1.10e-10 | 1.09e-08 | 1.53e-16 |
| 6.0 | 5.55e-16 | 1.30e-11 | 5.98e-11 | 2.27e-08 | 5.55e-17 |

The finite-difference errors are consistent with the selected step and
root tolerance.  The check passed its gates.  It is a local stationarity
and reaction check, not a global equilibrium or uniqueness result.

The supplied R oracle also passed when rerun to `/tmp/stackelberg-oracle-check`.
Its FOC, implicit-reaction, leader-gradient, and bounded local-optimization
tables are useful evidence.  The expanded oracle described by the research
lead adds follower optimization/Hessian checks and stress cases; those are
appropriate additions.

## Ordered random-boundary likelihood

The proposed probability model is coherent.  With descending firm sizes
`a_1 >= ... >= a_n`, a common `C=c+tau*epsilon`, and `epsilon` standard
logistic, `v_i = P(a_i>C) = plogis((a_i-c)/tau)` is descending.  The only
realized role partitions are nested top-`k` leader sets, with

```
pi_0 = 1-v_1,
pi_k = v_k-v_(k+1),  0<k<n,
pi_n = v_n.
```

These are nonnegative and telescope to one.  Exact ties produce zero mass
for the skipped intermediate partitions and should be retained as ties,
without arbitrary firm tie-breaking.  The mixture must not clip `C` to
`[0,1]`; its logistic tails are exactly what produce all-follower and
all-leader mass.

The margin component must be formed for the whole market-year before
`log_sum_exp`: one common boundary gives one common partition for all
products.  Per-product mixture likelihoods would assign different role
partitions within one market and are incoherent.  Share terms can be reused
if they are genuinely role-independent; with no observed margins, the
component scores coincide and the cutoff is prior-identified.  If a
conditional margin model uses the share residual (`rho`), the same
conditional factor must be used in every component and in generated
quantities.  Posterior component weights differ from prior `pi_k` after
conditioning on margins.

The implementation needs a numerically stable log-CDF difference for
intermediate `pi_k`, including exact ties and very large/small logits.  It
also needs to gather product rows through `obs_order` when using
`mky_start/mky_n`; those arrays describe market-year-sorted rows while the
observations may retain their original interleaving.  Required edge tests
are ties, one firm, zero and all leaders, all margins missing, row
permutations, extreme `c/tau`, and equality of generated mixture scores to
the model likelihood.

## Validation limits and sign-off recommendation

The oracle's costs are intentionally generated from the candidate FOCs, so
it cannot establish the FOC formulas by itself.  Its independent raw-profit
gradients, primitive inverse-demand checks, implicit derivatives, and local
optimization/Hessian checks substantially reduce that circularity.  The
remaining claims must stay bounded:

* interior IFT reactions do not prove global follower or leader existence;
* local negative Hessians do not prove global optimality or uniqueness;
* near-zero outside share and near-singular response denominators need
  numerical rejection/diagnostics;
* cutoff/timing identification, posterior recovery, and PNB behavior remain
  empirical gates;
* a single-cross-section cutoff cannot identify a continuous `c` within a
  gap between observed firm sizes without prior/model information.

Subject to those gates and the ordered-likelihood implementation checks
above, I recommend mathematical sign-off for the approved hard-role
equations.  I do not recommend describing the result as a global existence,
uniqueness, or empirical-identification theorem.
