# Market generation, posterior prediction, and inference contract

This package targets sparse structural IO data, including one market with one
observed margin. CmdStanR HMC is the exact-inference default. Pathfinder is an
explicit approximate-draw option; its optimization is not a structural MAP
estimate. `role_mapping = "map"` remains the conditional most-probable ordered-k
role mapping for a selected draw, not a continuous-parameter MAP.

## Calibrated output-Logit market

For one market/year with constant normalized prices \(p_j=1\), `safe_sd()`
returns 1 and standardized prices are zero. Let \(r_j>0\) be conditional
inside shares, \(\sum_j r_j=1\), and center quality so \(\sum_j b_j=0\). The
current Stan share equation is \(\log r_j=\operatorname{logit}(s_0)+b_j\) at
this baseline. Consequently,

\[
\operatorname{logit}(s_0)=J^{-1}\sum_j\log r_j,\qquad
b_j=\log r_j-J^{-1}\sum_i\log r_i,\qquad
s_j=(1-s_0)r_j.
\]

This is analytically derived from the live Stan parameterization and
numerically checked for the screened design. Constant prices intentionally
leave the demand slope identified through supply/margins rather than
cross-sectional price movement. A prespecified anchor margin then determines
the model-specific slope through a public deterministic-model margin oracle.
For auction conduct the oracle must use the fitted model's conditional-margin
convention and value index. Fresh equilibrium solves from perturbed starts,
positive costs, and fixed-role checks validate each generated cell. The
coordination routes additionally check FOC residuals; direct FOC residuals for
the MC/B/C/A public routes remain unknown (`NA`). Failed seeds are counted,
never silently replaced. Structural truth retains all margins and costs;
estimation data hide unobserved margins with `margin_observed`.

The focused test covers all eight conduct models on one \(J=5\) geometry. The
[deterministic sweep](calibrated_market_sweep.R) attempted ten \(J=7\) seeds
for each model (80 cells), with Dirichlet shape 2.5 and an explicit anchor
interval \((0.08,0.16)\). The [recorded results](calibrated_market_sweep.csv)
show zero structural failures in that bounded design. This does not establish
a general admissible interval. The public `antitrust` API exposes no direct
FOC-residual function for MC/B/C/A; those residuals remain `NA`, even though
their perturbed equilibrium price/share solves passed.

## Posterior predictive target

Stan fits `log_shareIn` and `log_marginNormInv`. Generated quantities such as
`pred_shareIn` and `pred_marginNormInv` are structural conditional means, not
replicated observations. `posterior_predict()` draws fresh errors from the
selected likelihood: Gaussian or variance-normalized Student-t, with market
scales and an active joint correlation when specified. With correlation, a
replicated margin error is conditioned on the *new replicated* share error.
Collapsed Gaussian quality and intercepts are drawn conditionally on observed
data; exact ordered-k roles are drawn before component-specific margin means
are computed. An entire response equation omitted by the fitted likelihood is
`NA` in the replication. Within an active margin equation, rows whose margins
were hidden from estimation receive posterior predictive imputations; the
observed-margin mask accompanies the output.

Independent log-share errors need not yield a simplex. The literal
`log_share_rep` is the fitted-scale predictive variable; any normalized share
vector is explicitly a derived summary. Policy counterfactuals never receive
observational residual noise.

## CmdStanR and evidence

The four existing Stan kernels are compiled from installed package sources
into a writable cache keyed by source hash and CmdStan version. A fitted
object records kernel/source/inference provenance and materializes its draws
so it survives ordinary R serialization. For bridge evaluation of a saved fit,
the package reconstructs the original Stan JSON data from the retained Stan
data and validates its recorded hash. LOO and bridge evidence require an HMC
fit and a matching current Stan source; failed evidence is reported as
unavailable, not imputed. Posterior model probabilities remain `NA` if any
candidate's bridge log marginal likelihood is unavailable. Pathfinder records
the pre-resampling Pareto-k
diagnostic when present. The tested 4- and 16-path runs had \(k\) about 1.8
and 1.7, respectively, so they are exploratory and do not support reliable
posterior or model-evidence claims. The short HMC pilot had divergences; its
speed comparison is engineering evidence, not sampler certification.
