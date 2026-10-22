# Default-model history and promotion discipline

## Single-market Logit utility correction (2026-10-21)

Version `0.2.0.9007`, specification 12, removes the incorrect zero inside-market
intercept in the two Logit kernels. The single-market intercept uses the existing
proper population-mean `Normal(0,1)` prior directly; no hierarchy SD is introduced.
Outside-good utility, centered product quality, price standardization, supply
equations and outside-share pooling are unchanged. There is no legacy switch.
Existing fits must be rerun before current comparison or simulation. CES kernels
are unchanged. The calibrated-market API independently accepts `outside_share`
(default `.20`) instead of deriving it from conditional-share geometry.
See the [normalization contract](single_market_utility_normalization.md) and
[executed correction checks](validation/single_market_intercept_results.md).

## Hierarchy SD prior defaults (2026-10-12)

Version 0.2.0.9003 made the active alpha and outside-share hierarchy SD
defaults Gamma(shape=2, rate=5). In version 0.2.0.9004, the user also approved
Gamma defaults for market-intercept SD (shape 2, rate 5) and the optional
market residual-scale heterogeneity SDs (shape 2, rate 7). The latter remains
disabled unless requested. These means are close to the historical
Half-Normal means; the Gamma density vanishes at zero. All four families
remain user-selectable, with proper normalization and active prior metadata.

Persistent-quality SD retains its positive-truncated Normal centered at 1;
share and margin error SDs retain Lognormal priors. Their modes already lie
away from zero. Fixed scales on signed-coefficient priors are hyperparameters,
not sampled SDs, so replacing them with Gamma would change the coefficient
support rather than an SD prior.

The [paired known-truth comparison](validation/hierarchical_prior_recovery_results.md)
tested alpha/outside-share Gamma versus Half-Normal while explicitly holding
the intercept Half-Normal and residual hierarchy disabled. It is not evidence
of HMC improvement from these additional default changes. Their verification
uses deterministic API, provenance, Stan compilation and normalized-density
checks. Existing fits and banking libraries keep their frozen priors; adopting
new defaults requires a newly installed package and separately labeled fits.

Verification on 2026-10-12 passed 1,500 assertions in nine focused files,
installation into an isolated library, and stanc plus C++ compilation for
all four kernels. Ninety-six full-kernel log-density differences for
intercept, residual and combined prior changes matched R's normalized
Gamma/Half-Normal densities, with maximum absolute error 5.33e-14. These
checks used default and custom hyperparameters and no HMC. Independent
review found no substantive issues. Raw verification artifacts are local at
`/tmp/sd-gamma-defaults-zaZJRh`.

## Likelihood and quality defaults

The historical specification used joint Gaussian share/margin disturbances and
explicit quality coordinates. Student-t residuals were subsequently adopted
amid divergence concerns. That history motivates investigation; it does not
prove that the Normal likelihood caused those divergences.

The controlled collapsed-Gaussian experiment dated 2026-09-16 held the Gaussian
statistical model fixed between sampled and analytically collapsed quality.
All seven collapsed Gaussian primary fits passed the frozen diagnostic gate;
none of the seven sampled Gaussian or seven sampled Student-t fits passed.
Poor BFMI and mixing, rather than many divergences, were the main reproduced
failure: each sampled group had only one divergence in total. Thus these
results support a sampled-quality geometry problem, not a universal claim
about the cause of historical divergences. The experiment and retained failures
are documented in [collapsed Gaussian results](collapsed_gaussian_results.md).

Analytical collapse preserves the historical Gaussian prior/likelihood exactly.
Market-specific quality is a subsequent substantive prior change: separate
persistent firm-market effects are centered within relevant markets and share
a common quality scale. It removes the historical global centering dependence
across markets. In one market the prior is exactly unchanged. It also permits
exact market-level ordered-k marginalization with active Gaussian rho after
quality integration.

The intended routine specification is Gaussian shares and margins, analytically
collapsed market-specific quality, and active rho. Student-t sampled-quality
and Gaussian-share/Student-t-margin hybrid specifications remain available.
Scope alone does not ensure usable HMC or robust policy inference. Mathematical
validation and executed sampler comparisons support the new architecture;
every application still requires its own sampler and counterfactual checks.

Status: on 2026-09-17 the user approved promoting Gaussian shares/margins,
collapsed market quality and active rho to the routine defaults after reviewing
the rare PNB sampling warning and the separate negative-rate counterfactual
obstruction. Function signatures, help pages and tests are updated together.
The default margin family follows the share family when `margin_likelihood=NULL`,
so both resolve to Gaussian without requiring a degrees-of-freedom argument.
Historical serialized fits without scope metadata still denote global quality;
they are never reinterpreted as market quality.

This decision supersedes the original all-gates promotion hold; it does not
alter the frozen usability criteria or claim that all gates passed. The
remaining conduct/scaling HMC work is incomplete, and provisional PNB evidence
and BMA withholding retain their original labels. No prior, pooling,
market/outside-share coordinate, sampler, solver or economic default changes
with this promotion. See [validation contract](market_quality_validation.md)
and [current evidence](market_quality_hmc_results.md).

The historical model remains explicit:

```r
bayes_est(...,
  share_likelihood = "student_t",
  margin_likelihood = "student_t",
  quality_parameterization = "sampled",
  quality_scope = "global",
  use_rho = TRUE,
  df = 5
)
```

The first two-seed multi-market benchmark with active rho does not reproduce
the earlier collapsed-Gaussian usability success: both collapsed and sampled
Gaussian market-quality fits have divergences, and none of the six initial
Gaussian/Student-t fits passes the complete frozen gate. Collapsing improves
BFMI and effective counts but lowers minimum bulk ESS per wall second in this
comparison. These failures remain evidence; the earlier seven-fit result is
not transferred to this different configuration. See
[current sampler evidence](market_quality_hmc_results.md).

Later, separately labeled market/outside-share coordinate sensitivities
preserved the target and materially improved multi-market geometry. Both
collapsed fits at the doubled budget passed the frozen HMC gate, with zero
divergences. This improves the computational evidence without replacing the
original failed runs or completing the remaining conduct benchmarks.

Seven of eight archived-input PNB Gaussian fits passed the full sampler gate;
Bertrand had two divergences in 24,000 post-warmup draws. Complete canonical and
rich bridge calculations are retained as provisional evidence. Policy draws
also retain Cournot-fringe and Cournot-Stackelberg solver failures. A direct
interior-equilibrium oracle diagnoses every observed-anchor failure as an
implied negative deposit rate outside the positive-rate solver domain; this
counterfactual restriction is separate from quality collapse and HMC. These
are small or localized diagnostic problems, not a finding that Gaussian
residuals or market quality are invalid. Prespecified benchmark evidence remains
incomplete; complete-universe BMA is withheld. Default promotion is a separate,
explicitly approved package choice, not retrospective reclassification.

The promotion verification passed 1,909 package assertions and 32 exact
before/after Stan-data comparisons. The compiled kernels and normalized targets
of explicitly selected models did not change. Details and the unchanged failure
labels are retained in the [validation record](market_quality_validation.md).
