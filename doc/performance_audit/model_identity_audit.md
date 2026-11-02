# Model identity audit

The parameterized simulation path already maps the five Bayesian conduct codes
to canonical `antitrust::model_spec()` entries in `.bayes_model_spec()`.
`as_antitrust_fit()` checks any supplied specification against that canonical
entry. No second registry is required for ordinary registered models.

Integer codes remain necessary in Stan data, so removing every integer mapping
would add indirection without reducing inference work. The display-name table
in `bayes_compare()` is a drift risk, not a measured performance bottleneck.
It can consume the same package-owned conduct mapping in a separate correctness
unit once the extraction changes to `R/bma.R` are reviewed.

Two deliberate exceptions must remain explicit: smooth-cutoff conduct and
Logit multiproduct bargaining use Bayesian conduct classes because their
equations differ from the canonical registered class. Changing their class
merely to reduce identity duplication would change economics. Exact fringe
and Stackelberg realizations retain their deterministic structural dispatch.

Trade extends `antitrust::StructuralFit`; antitrustBayes must not depend on
trade to initialize or simulate its own structural draws. Public baseline
lifecycle behavior should be consumed through antitrust.
