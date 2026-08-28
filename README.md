# antitrustBayes

Bayesian estimation of differentiated products demand and supply models using Stan, with merger simulation via the [antitrust](https://github.com/luciu5/antitrust) R package.

## Overview

`antitrustBayes` provides two core functions:

- **`bayes_est()`** — Estimate a hierarchical Bayesian Logit or CES demand + supply model using Stan MCMC. Supports four supply models (Bertrand, Monopolistic, Cournot, Second-Score Auction).

- **`bayes_sim()`** — Simulate merger effects over posterior draws. Each draw produces a full merger simulation via `antitrust::sim()`, yielding posterior distributions of price changes, share changes, and all other `antitrust` package metrics.

Additional utilities:

- **`bayes_compare()`** — Compare models using WAIC/LOO-IC with stacking weights.
- **`bayes_avg()`** — Bayesian model averaging across supply model specifications.

## Installation

```r
# Install from GitHub (requires devtools)
devtools::install_github("luciu5/antitrustBayes")
```

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
  supply_model = "bertrand"
)

# 2. Simulate merger effects
result <- bayes_sim(
  fit       = fit,
  ownerPost = post_merger_ownership,
  n_draws   = 500
)

# 3. Summarize results with credible intervals
summary(result)
```

## Dependencies

- [rstan](https://mc-stan.org/rstan/) (>= 2.26.0)
- [antitrust](https://github.com/luciu5/antitrust) (>= 0.99.27)
- [loo](https://mc-stan.org/loo/) (>= 2.5.0)
