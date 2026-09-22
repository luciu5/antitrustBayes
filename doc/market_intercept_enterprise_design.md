# Enterprise market-intercept validation design

The benchmark asks one narrow question: does exact Gaussian integration of a
multi-market intercept make the intended Bayesian model reliably sampleable,
without altering its posterior target? The analytical density, gradient, and
all-four-kernel checks are already complete. This design adds broad HMC
evidence without making the Supreme Court application the certification DGP.

## Evidence hierarchy

1. **Analytical and fixed-target evidence.** The existing dense oracle and
   native Stan checks establish normalized target and gradient equivalence.
2. **Synthetic HMC robustness.** Logit and CES, all eight conduct cases,
   cross-sections and panels, exact ordered-k mixtures, active rho, persistent
   market quality, missing margins, uneven markets, and quality-dominant
   compositions are tested under frozen budgets. The normalized compositional
   perturbation is not an exact iid-Gaussian DGP for the package's log-share
   pseudo-likelihood, so these runs certify sampler behavior under coherent
   share data rather than exact Gaussian recovery.
3. **Sampled-versus-collapsed comparison.** B, C, BS, and CS are paired on the
   same datasets, DGP seeds, priors, sampler seeds, and HMC budgets. The sampled
   model uses the validated joint outside-share/intercept coordinates; the
   collapsed model removes that coordinate and uses the legacy outside-share
   hierarchy. These are representations of the same economic/statistical
   model.
4. **Supreme Court regression.** The six-market, 176-observation base-year
   cross-section and the base-year-minus-one through base-year-plus-two panel
   test the actual sparse, weakly identified application. They are not used to
   tune priors or sampler controls.

## Resource policy

No production HMC runs on the unstable WSL machine. Local validation is
limited to deterministic fixture construction, task isolation, manifest
checks, aggregation checks, and package tests. The Slurm array runs one fit per
process, uses four CPUs for four chains, pins BLAS/OpenMP to one thread, and
writes only compact CSV evidence. Raw `stanfit` objects are not retained.

Run the 32 collapsed synthetic gate tasks first. If a conduct/demand/panel cell
fails, preserve it and diagnose before running the corresponding sampled
comparison or expanding seeds. Run the second seed only after the first seed
passes or when a rare-event diagnostic requires replication. Model evidence
and bridge comparisons remain invalid until all constituent fits pass the
predeclared HMC gate.

## Interpretation

Synthetic success establishes broad numerical robustness under the supported
model architecture. Supreme Court failure after synthetic success would point
to application-specific weak identification or boundary geometry; it would
not falsify target equivalence. Conversely, one clean Supreme Court fit cannot
establish enterprise robustness. Both forms of evidence are required before
changing `market_intercept_parameterization` from opt-in to a default. The
rich synthetic FOC schedules do not replace deterministic equilibrium and
counterfactual tests in `coordination`.
