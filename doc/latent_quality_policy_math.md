# Latent-quality counterfactuals: normalization and cost specification

Status: implementation integrated locally; independent mathematical and native
reference checks passed. Frozen posterior policy benchmarks remain a separate
gate. No simulation default has changed.
This is a substantive simulation extension, separate from the exact sampler
rotation. The existing observed-share anchor remains directly comparable.

Root re-execution verified actual refactor antitrust Logit/CES share and Logit
cost predictions. An independent native test passed 24 checks: Logit/CES x
Bertrand/Cournot x Stackelberg/fringe x merger/cost/combined interventions.
Maximum reference share error was `2.22e-16`; solved pre/post FOC residuals were
below `4.78e-13`; costs, demand values, and market-size invariants agreed
exactly with the stated interventions. No coordination source change was
needed. These are tiny reference fixtures, not posterior policy evidence.

Failed reference attempts are retained locally. A low-alpha Cournot-fringe
fixture correctly failed positive-cost interiority; a separate analytically
interior fixture was then used for the support check. A CES test initially
compared default quantity shares against revenue shares and was corrected to
request revenue explicitly. Neither failure led to relaxed solver tolerances.

## Quality posterior and structural shares

For sampled quality, use the quality vector from the same Stan draw. For
collapsed Gaussian quality, draw it conditionally after fitting. Reuse one
global group-quality vector across all markets, periods, and policy scenarios
for that common-parameter draw. Do not independently redraw quality by market
or policy. With supported active Gaussian rho, recovery must condition on
matched margins; unavailable full conditional recovery must error rather than
silently use the share-only posterior.

Let ell_j denote the fitted complete log conditional-share location, including
that quality draw. It generally does not add up after exponentiation. The
initial latent simulation uses an explicit outside-share-preserving projection:

\[
t_j^L = \frac{\exp(\ell_j)}{\sum_h\exp(\ell_h)},\qquad
s_j^L=(1-s_0)t_j^L.
\]

Calculate through log-sum-exp. This propagates relative quality uncertainty and
retains the posterior outside-share draw at the observed reference prices.
Common utility level is reanchored; the construction does not claim that the
current transformed-outcome likelihood is already a simplex likelihood.
Market-common intercept terms cancel from t. Their posterior role in outside
share and other structural parameters remains unchanged.

The corresponding quantity odds relative to the outside good are

\[
z_j=\frac{s_j^L}{s_0},\qquad
\log z_j=\log t_j^L-\operatorname{logit}(s_0).
\]

Equivalently, for arbitrary inside scores x_j, the level shift is
`-logit(s0) - log_sum_exp(x)`. Using positive `logit(s0)` swaps the intended
outside and inside shares. Independent package predictions must test this
sign; an internal algebra test alone is insufficient.

## Logit normalization

At the antitrust storage boundary use its signed slope a_p, negative for
output demand and positive for the supported input convention. If its utility
at prices is `meanval_j + a_p * (p_j - p0)`, recover

\[
meanval_j^L=\log z_j-a_p(p_j-p_0).
\]

Verify the actual package normalization with `calcShares()`; do not introduce
a competing outside-price convention. The Bayesian structural alpha remains
the positive magnitude. For output Bertrand, dollar markup is
`mu_j = 1 / (alpha * (1-S_f))`, Lerner margin is `mu_j / p_j`, and cost is
`p_j - mu_j`. These units must remain distinct.

For input procurement profit is
\[
\pi_f=\sum_{j\in f}(v_j-p_j)q_j(p),\qquad
\frac{\partial\pi_f}{\partial p_i}
=-q_i+\sum_{j\in f}(v_j-p_j)\frac{\partial q_j}{\partial p_i}.
\]
Thus the positive gap is value minus procurement price. For input Logit,
the antitrust convention uses a positive price slope and
reports the positive dollar markdown
`mu_j = v_j - p_j` from `calcMargins(level = TRUE)`. The stored structural
value or procurement-price primitive is therefore `v_j = p_j + mu_j`, and the
latent pre-policy residual is `mc_j - p_j - mu_j`. The same primitive is held
fixed when demand is replaced, and the native input solver is used for both
pre- and post-policy equilibria. This support is limited to input Logit;
input CES remains unsupported pending a separate normalization audit.

## CES normalization

For the implemented output CES convention, quantity odds satisfy
`q_j/q0 = A_j * (p_j/p0)^(-gamma)`. Therefore

\[
\log A_j^L=\log z_j+\gamma\log(p_j/p_0).
\]

Construct latent revenue shares using actual prices and latent quantity shares:

\[
R_j^L=\frac{p_js_j^L}{p_0s_0+\sum_hp_hs_h^L},\qquad
R_0^L=\frac{p_0s_0}{p_0s_0+\sum_hp_hs_h^L}.
\]

The observed inside-price index generally cannot be reused for this conversion.
Input CES requires an explicit sign/normalization audit before support.

In the antitrust input interface, a `costs` counterfactual is a proportional
shock to the stored value/procurement-price primitive. It is a value shock,
not a tariff or a claim about tariff pass-through, so policy results using
that field should be described as value shocks.

## Fixed costs and baseline equilibrium

For the principal comparison, first construct the original observed-share
structural object and retain its cost primitive. Replace its demand state with
the latent-quality state, solve the pre-policy equilibrium at those same
costs, and solve post-policy with the same quality and the stated ownership or
cost intervention. Pre-policy prices can differ from observed prices. Report
that discrepancy and equilibrium residuals explicitly.

Do not accidentally recalculate costs after replacing shares: re-inverting
FOCs at observed prices would change costs together with quality and define a
different sensitivity experiment. Nonconstant cost closures must be preserved
as closures, not only evaluated at observed quantities. Supported deterministic
methods and cost primitives are audited before choosing the implementation.

For exact-k conduct, draw k with the same common-parameter draw and retain the
predetermined observed-share firm ranking. Latent shares affect demand and
equilibrium; they do not redefine the estimation partition. Quality and k
plans must be reusable across interventions without refitting.

## Validation and interpretation

Require independent antitrust/coordination share predictions at the reference
prices, valid quantity/revenue adding up, solved pre/post FOC residuals,
unchanged captured costs absent a shock, draw/quality/k identity across paired
scenarios, multiproduct and panel label mapping, and explicit missing-margin
and supported-rho handling. Compare policy distributions against the anchor
with common draw indices and scenario seeds.

This is a specified posterior pushforward using an outside-preserving
normalization and the existing inferred-cost convention. It is not an exact
reparameterization of the old counterfactual map, nor a claim that all absolute
quality-level uncertainty has been retained. A fully generative compositional
share/latent-cost likelihood is a separate research extension.
