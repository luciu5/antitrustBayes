# antitrustBayes

Bayesian estimation of differentiated-products demand and supply models using CmdStanR, with posterior counterfactual analysis via the [antitrust](https://github.com/luciu5/antitrust) R package.

## Overview

`antitrustBayes` provides three core entry points:

- **`bayes_est()`** — Estimate a hierarchical Bayesian Logit or CES demand + supply model using CmdStanR HMC. An explicit Pathfinder option provides exploratory approximate draws. Supports Bertrand, monopolistic competition, Cournot, second-score auction, and (for applications where it is institutionally defensible) Nash bargaining.

- **`bayes_specify()`** — Create the same structural-fit interface from supplied parameter draws and observed market geometry, without running Stan.

- **`simulate()`** — Apply an `antitrust::counterfactual()` over the posterior distribution. Mergers, firm-specific cost shocks, exit, entry, quality changes, and sequential scenarios use `antitrust`'s common machinery. `bayes_sim()` remains a convenience wrapper for legacy ownership arguments.

Additional utilities:

- **`bayes_calibrated_market()`** — Generate a one-market output-Logit design from Dirichlet inside shares, normalized prices, and one observed margin; retain full structural truth for oracle checks.
- **`posterior_predict()`** — Draw replicated likelihood-scale shares and margins from a fitted posterior, with the original observation mask retained.

- **`bayes_compare()`** — Compare models using bridge sampling, with WAIC/LOO-IC and stacking where the likelihood supports them. Collapsed-quality fits use bridge evidence; ordinary pointwise LOO/WAIC is unavailable. Logit/CES comparison is allowed when both fits use the same quantity-space outcomes and market side.
- **`bayes_avg()`** — Bayesian model averaging across supply model specifications.
- **`bayes_conditional_k_arguments()`**, **`bayes_est_conditional_k()`**, and
  **`bayes_combine_conditional_k()`** — Separate fixed-k computation for one
  market's Bertrand/Cournot fringe or Stackelberg hypothesis. Components keep
  the original continuous priors; their evidence and predictions are mixed
  with the original k prior and resulting posterior k weights.
- **`bayes_prior_manifest()`** and **`bank_model_priors()`** — Reproducible proper parameter priors and candidate-set probabilities. The banking baseline assigns 0.30 each to Bertrand, Cournot, and monopolistic competition and 0.10 to the auction model; equal-weight and core-three sensitivities are available.
- **`conduct_distance()`**, **`conduct_information_matrix()`**, and **`conduct_validation()`** — Theory-guided conduct distinguishability and targeted Monte Carlo validation.
- **`merger_decision()`** and **`decision_regret()`** — Translate posterior merger effects into threshold decisions and decision regret.

## Installation

The package depends on CmdStanR, an installed CmdStan toolchain, and the matching development stack.
Install the known-compatible `antitrust` revision before installing
`antitrustBayes`:

```r
# Install from GitHub (requires remotes)
remotes::install_github(
  "luciu5/antitrust@103c2491607d7f164f5e4e00fc21d6ea79331138"
)
remotes::install_github(
  "luciu5/antitrustBayes"
)
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
`log_shareIn`. New fits default to Gaussian share and margin residuals,
analytically collapsed persistent quality centered within each market, and
`use_rho = TRUE`. Quality persists across periods for each firm-market group.
The share-margin residual correlation is estimated only for a joint fit with
at least one observed margin; otherwise it is structurally fixed at zero and
is not sampled. The historical specification remains available with
`share_likelihood = "student_t"`, `margin_likelihood = "student_t"`,
`quality_parameterization = "sampled"`, and `quality_scope = "global"`.
Shared proper priors are required for bridge-sampling
model probabilities. Use `likelihood = "share"` and `likelihood = "margin"`
fits with `bayes_evidence_decomposition()` when you want to report how much
evidence comes from each block.

## Market-intercept SD prior

For multiple markets, select the prior through `bayes_est(..., priors = ...)`:

```r
priors <- list(
  prior_sigma_b_market_family = "gamma",
  prior_sigma_b_market_shape = 2,
  prior_sigma_b_market_rate = 5
)
# fit <- bayes_est(..., priors = priors)
```

Gamma uses shape and **rate** and is the default from version 0.2.0.9004.
Select the historical alternative with `prior_sigma_b_market_family = "half_normal"` and
`prior_sigma_b_market_scale = 0.5`. Both choices use normalized densities for
bridge evidence and are recorded in resolved fit priors. The setting applies
to Logit/CES, every conduct family, and sampled/collapsed intercepts. It is
inactive for a single market. Keep the scientific prior consistent across
candidate models; changing it requires new evidence under that prior.

## Alpha and outside-share hierarchy SD priors

From version 0.2.0.9003, active log-price-effect and outside-share hierarchy SDs
default to Gamma(shape=2, rate=5). Both means are 0.4, close to the historical
HalfNormal(0, 0.5) mean 0.399, with less probability near negligible heterogeneity.
The log-price effect means log-alpha for Logit and log of the standardized
price effect for CES. Outside-share SD is on the logit scale.

```r
priors <- list(
  prior_sigma_alpha_family = "gamma",
  prior_sigma_alpha_shape = 2,
  prior_sigma_alpha_rate = 10, # Mean 0.2 for a smaller heterogeneity scale.
  prior_sigma_s0_pool_family = "gamma",
  prior_sigma_s0_pool_shape = 2,
  prior_sigma_s0_pool_rate = 5
)
# fit <- bayes_est(..., priors = priors)
```

To reproduce the historical defaults, explicitly select both `"half_normal"`
families and set `prior_sigma_alpha = 0.5` and `prior_sigma_s0_pool = 0.5`.
These existing unsuffixed fields control only the Half-Normal scales; changing
them does not change a Gamma rate. Family names also accept codes 1 and 2.
Alpha SD is absent without multi-market pooling. Outside-share SD is active
for multi-market market pooling or market-year pooling, including a one-market
panel. Both choices preserve normalized evidence and active prior provenance.
Use a frozen historical library for existing fits; rebuild and refit under a
new run label when adopting the new priors. Banking libraries are not rebuilt
by these package source changes.

## Optional residual-heterogeneity SD priors

When `use_market_resid_scale = TRUE`, the SDs of market deviations in log
share-error and log margin-error scales default to Gamma(shape=2, rate=7).
The mean 2/7 = .286 is close to the historical HalfNormal(0, .35) mean .279;
the mode is .143. This changes the heterogeneity prior without enabling the
optional hierarchy, which remains disabled by default.

```r
priors <- list(
  prior_sigma_resid_market_family = "gamma",
  prior_sigma_resid_market_shape = 2,
  prior_sigma_resid_market_rate = 7
)
# Supply with use_market_resid_scale = TRUE in a supported likelihood configuration.
```

To select the historical prior, set `prior_sigma_resid_market_family =
"half_normal"` and `prior_sigma_resid_market = .35`. Both share and margin
heterogeneity SDs use this specification. Densities are normalized, and
inactive families/hyperparameters do not change the model-evidence identity.

All sampled SDs that previously had zero-centered Half-Normal defaults now
default to Gamma. Persistent-quality SD keeps its positive-truncated Normal
centered at 1, expressing quality dispersion on its existing scale. Share and
margin error SDs keep Lognormal priors. These priors already have modes away
from zero. Fixed Normal prior scales for signed coefficients are not sampled
SD parameters and remain unchanged. Select all four `*_family` fields
explicitly as `"half_normal"` to reproduce historical hierarchy priors.

## Separate fixed-k computation

For one market with several periods, a marginalized ordered-k fit can mix
poorly between continuous posterior regions. An alternate calculation fits
each k separately, then combines its normalized evidence:

```r
calls <- bayes_conditional_k_arguments(
  prices = dat$prices, quantities = dat$quantities, margins = dat$margins,
  market = dat$market, year = dat$year, group = dat$group,
  demand_model = "logit", market_side = "output",
  structure = "fringe", supply_model = "cournot",
  k_support = 1:6, k_prior = rep(1/6, 6),
  stan_options = list(chains = 4, iter = 4000, warmup = 2000, cores = 4)
)
# No HMC above. Execute/schedule each ordinary bayes_est call explicitly:
# fits <- lapply(calls, function(args) do.call(bayes_est, args))
# Review every component's diagnostics and obtain its normalized bridge logml.
# mixture <- bayes_combine_conditional_k(
#   logml = log_z, k_support = attr(calls, "k_support"),
#   k_prior = attr(calls, "k_prior")
# )
```

`bayes_est_conditional_k()` executes the same component calls sequentially and
retains failed components as error conditions. The argument builder permits
independent scheduling and explicit component seeds/output paths. The mixer
requires the complete original support; it never renormalizes over missing
evidence. Use its posterior k weights to mix conditional policy draws, keeping
one draw plan per component for paired counterfactuals. The components represent
one conduct hypothesis; its top-level model prior is unchanged.

This route applies to all four rich families. It rejects multiple markets:
separate market k values and shared parameters require joint combinations.

## Dependencies

- [CmdStanR](https://mc-stan.org/cmdstanr/) (>= 0.9.0) and a working [CmdStan](https://mc-stan.org/users/interfaces/cmdstan) installation
- [antitrust](https://github.com/luciu5/antitrust) (>= 0.99.38)
- [loo](https://mc-stan.org/loo/) (>= 2.5.0)
