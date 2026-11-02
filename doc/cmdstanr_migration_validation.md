# CmdStanR migration and calibrated-market validation, 2026-09-25

The starting local `antitrustBayes` commit was
`438812ad729b2908dd38a324aeb429875dd44cf8`. No Stan kernel source was
changed. The older [pilot](cmdstanr_pilot/results.md) is retained as historical
engineering evidence; its timing scopes differ and its divergent short chains
do not certify a speed or inference advantage.

## Executed checks

* CmdStan 2.39.0 compiled all four current kernels from package sources into a
  writable temporary cache. A source/CmdStan-version keyed cache is used by
  the package, so these builds do not populate `inst/stan`.
* The full non-opt-in `devtools::test(reporter = "summary")` suite completed
  with no failed expectations. It reported twelve existing opt-in Stan skips
  and warnings from deterministic `antitrust` simulation calls.
* A source tarball built from a clean package-source staging directory passed
  `R CMD check --no-manual --no-vignettes` installation, code, documentation,
  examples, and installed-package tests. The check reported two vignette
  warnings because vignette building was explicitly disabled; it did not run
  the package vignettes. The source staging excluded the repository's
  historical 6.6 GB `doc/` archive, which `.Rbuildignore` already excludes
  from the distributed package.
* [The deterministic sweep](calibrated_market_sweep.csv) attempted every cell
  in a 10-seed by 8-conduct design at `J = 7`, Dirichlet shape `2.5`, and
  anchor-margin interval `(0.08, 0.16)`. It recorded 80 successes and zero
  structural failures without redrawing cells. This is a bounded screen, not
  proof of general admissibility. The public `antitrust` MC/B/C/A routes do
  not expose direct FOC residuals; these checks are `NA` in the artifact.
* Calibrated one-market `B` and `BF` HMC smokes used two chains, 50 warmup and
  50 retained iterations per chain, `J = 5`, a fixed `0.12` anchor margin,
  and one observed margin. Both fits had `price_sd = 1` and zero reported
  divergences. These chains are too short for convergence claims.
* Saved `B` and `BF` HMC fits were reloaded in a fresh session. Ten posterior
  predictive draws per fit yielded finite fitted log shares and finite margins
  for all five products, with precisely one margin marked observed. The `BF`
  output retained an exact-k component state. Prediction leaves the serialized
  CmdStanR fit unchanged despite its usual lazy draw cache. A saved `BF` fit
  also produced
  two valid policy models and two valid merger simulations after the legacy
  one-year plan-key compatibility fix.
* A longer two-chain, 500-retained-draw HMC engineering fit was serialized and
  reloaded. Bridge sampling using a writable CmdStan installation and restored
  original Stan JSON data returned finite log marginal likelihood `54.48978`
  with `validated` bridge status and a one-model posterior probability of 1.
  That fit had five divergences and adverse PSIS/WAIC warnings; the numerical
  bridge result validates the API route, not the model evidence.
* Pathfinder four- and sixteen-path exploratory runs completed quickly, but
  pre-resampling Pareto-k diagnostics were approximately `1.8` and `1.7`.
  These draws are not certified for posterior or policy inference. HMC remains
  the default; Pathfinder requires an explicit request and policy opt-in.

The R-side posterior predictive tests independently check Gaussian and
Student-t residual formulas, active correlation, seeded draws, missing-margin
imputation, collapsed quality states, exact-k component means, and empirical
exact-k mixture frequencies against known conditional probabilities. A
one-market/two-year regression checks global collapsed quality with active
rho and a shared exact-k state. Replicated
`log_shareIn` is the literal likelihood-scale variable; normalized shares are
a derived view.

## Limits and next gate

The current checks establish executable backend and API routes, not production
sampler quality for all eight models. Longer converged HMC runs, diagnostic
thresholds, and calibrated-market feasibility screens over broader geometries
remain necessary before using generated cells or model probabilities in a
paper-facing Monte Carlo. No downstream `mergerBayes` files were modified.
