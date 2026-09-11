# antitrustBayes

Bayesian estimation of differentiated-products demand and supply models using Stan, with posterior counterfactual analysis via the [antitrust](https://github.com/luciu5/antitrust) R package.

## Overview

`antitrustBayes` provides three core entry points:

- **`bayes_est()`** — Estimate a hierarchical Bayesian Logit or CES demand + supply model using Stan MCMC. Supports Bertrand, monopolistic competition, Cournot, second-score auction, and (for applications where it is institutionally defensible) Nash bargaining.

- **`bayes_specify()`** — Create the same structural-fit interface from supplied parameter draws and observed market geometry, without running Stan.

- **`simulate()`** — Apply an `antitrust::counterfactual()` over the posterior distribution. Mergers, firm-specific cost shocks, exit, entry, quality changes, and sequential scenarios use `antitrust`'s common machinery. `bayes_sim()` remains a convenience wrapper for legacy ownership arguments.

Additional utilities:

- **`bayes_compare()`** — Compare models using WAIC/LOO-IC with stacking weights and optional bridge sampling. Logit/CES comparison is allowed when both fits use the same quantity-space outcomes and market side.
- **`bayes_avg()`** — Bayesian model averaging across supply model specifications.
- **`bayes_prior_manifest()`** and **`bank_model_priors()`** — Reproducible proper parameter priors and candidate-set probabilities. The banking baseline assigns 0.30 each to Bertrand, Cournot, and monopolistic competition and 0.10 to the auction model; equal-weight and core-three sensitivities are available.
- **`conduct_distance()`**, **`conduct_information_matrix()`**, and **`conduct_validation()`** — Theory-guided conduct distinguishability and targeted Monte Carlo validation.
- **`merger_decision()`** and **`decision_regret()`** — Translate posterior merger effects into threshold decisions and decision regret.

## Installation

The `refactor` branch currently depends on the matching development stack.
Install the known-compatible `antitrust` revision before installing
`antitrustBayes`; do not rely on either repository's default branch:

```r
# Install from GitHub (requires remotes)
remotes::install_github(
  "luciu5/antitrust@103c2491607d7f164f5e4e00fc21d6ea79331138"
)
remotes::install_github("luciu5/antitrustBayes@refactor")
```

The dependency-stack CI uses the same pinned `antitrust` SHA by default and
records both resolved repository SHAs. Its manual dispatch accepts another
`antitrust` branch or SHA for compatibility testing.

## Quick Start

```r
library(antitrustBayes)

# 1. Estimate the Bayesian model
fit <- bayes_est(
  prices  = dat$prices,
  quantities = dat$quantities,
  margins = dat$margins,
  market  = dat$market,
  year    = dat$year,
  group   = dat$group,
  supply_model = "bertrand",
  prior_manifest = "baseline",
  model_prior = bank_model_priors()[["bertrand"]]
)

# 2. Describe and simulate the merger over posterior draws
scenario <- antitrust::counterfactual(ownership = post_merger_ownership)
result <- simulate(fit, scenario, n_draws = 500)

# 3. Summarize results with credible intervals
summary(result)
```

For a banking model comparison, pass the candidate-set probabilities explicitly
so that the institutionally motivated auction downweighting is visible and
reproducible:

```r
fits <- lapply(names(bank_model_priors()), function(sm) {
  bayes_est(..., supply_model = sm,
            model_prior = bank_model_priors()[[sm]])
})
comp <- bayes_compare(
  Bertrand = fits[[1]], Cournot = fits[[2]], Auction = fits[[3]],
  MonCom = fits[[4]], prior_weights = bank_model_priors()
)
```

The primary share outcome is the log conditional inside quantity share,
`log_shareIn`. By default, `use_rho = TRUE` couples it to the log inverse-margin
outcome with a bivariate Student-t likelihood factored as
`p(share) * p(margin | share)`. The correlation is estimated only for a joint
fit with at least one observed margin; otherwise it is structurally fixed at
zero and is not sampled. Shared proper priors are required for bridge-sampling
model probabilities. Use `likelihood = "share"` and `likelihood = "margin"`
fits with `bayes_evidence_decomposition()` when you want to report how much
evidence comes from each block.

## Dependencies

- [rstan](https://mc-stan.org/rstan/) (>= 2.26.0)
- [antitrust](https://github.com/luciu5/antitrust) (>= 0.99.38)
- [loo](https://mc-stan.org/loo/) (>= 2.5.0)
