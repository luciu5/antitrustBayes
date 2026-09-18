# Default-model history and promotion discipline

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
