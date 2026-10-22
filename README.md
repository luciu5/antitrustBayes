# antitrustBayes

Bayesian estimation of differentiated-products demand and supply models using CmdStanR, with posterior counterfactual analysis via the [antitrust](https://github.com/luciu5/antitrust) R package.

CmdStanR remains the default. An optional Nutpie backend evaluates the same
installed Stan programs through BridgeStan; it is an alternative NUTS sampler,
with diagonal or low-rank adaptation. It does not change the economic model.

For development, current validation commands, and research-document status,
start with the [documentation index](doc/README.md).

## Overview

`antitrustBayes` provides three core entry points:

- **`bayes_est()`** — Estimate a hierarchical Bayesian Logit or CES demand + supply model using CmdStanR HMC (default) or optional Nutpie NUTS. CmdStanR's explicit Pathfinder option provides exploratory approximate draws. Supports Bertrand, monopolistic competition, Cournot, second-score auction, and (for applications where it is institutionally defensible) Nash bargaining.

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

Install CmdStanR and a CmdStan toolchain following the
[CmdStanR installation instructions](https://mc-stan.org/cmdstanr/).
Install compatible revisions of `antitrust` and `coordination` before
installing `antitrustBayes`. The declared package minima are
`antitrust >= 0.99.38`, `coordination >= 0.0.0.9003`, and `cmdstanr >= 0.9.0`;
these are dependency requirements, not a tested combination of arbitrary
revisions. Check the selected R library before development:

```r
stopifnot(
  packageVersion("antitrust") >= "0.99.38",
  packageVersion("coordination") >= "0.0.0.9003",
  packageVersion("cmdstanr") >= "0.9.0"
)
```

With compatible dependencies already installed, install the package from GitHub
without upgrading the selected dependency versions:

```r
remotes::install_github("luciu5/antitrustBayes", upgrade = "never")
```

The [dependency-stack CI](.github/workflows/qa.yml) defaults to `antitrust`
SHA `103c2491607d7f164f5e4e00fc21d6ea79331138` (version 0.99.39) and
records resolved `antitrust` and `antitrustBayes` repository SHAs. It does
not pin `coordination` or the complete dependency stack. In the sources
inspected on October 13, 2026, `coordination` 0.0.0.9006 requires
`antitrust >= 0.99.41`, so it cannot be combined with that default CI pin.
Some integration tests also require `coordination >= 0.0.0.9006` and skip
with older versions. CI manual dispatch accepts another `antitrust` branch
or SHA for explicit compatibility testing. See the documentation index for
local environment checks and validation tiers.

## Optional Nutpie backend

The [executed validation](doc/validation/nutpie_backend_results.md) covers 34
sampler fits in 12 fixtures, all four Stan kernels and pooled C/CF/CS cases.
Scored posterior/policy comparisons agree within the declared Monte Carlo
tolerances. Nutpie diagonal is a validated same-target alternative in these
fixtures and is worth trying for difficult hierarchical C/CF/CS fits. CmdStanR
remains the default and mature fallback. Diagnostic concerns and CES policy-tail
uncertainty remain explicit; these results do not establish general superiority.
Low-rank remains an expert/experimental option and did not consistently improve
on diagonal adaptation in the executed studies. The separate
[qualification pilot](doc/validation/nutpie_certification_pilot_results.md)
also found residual pooled mixing concerns and CES divergences on both samplers.

A CmdStan-only installation does not need `nutpieR`, BridgeStan's R client,
Python, or Rust. Request the optional sampler explicitly:

```r
fit <- bayes_est(
  prices = dat$prices, quantities = dat$quantities, margins = dat$margins,
  market = dat$market, year = dat$year, group = dat$group,
  supply_model = "cournot", backend = "nutpie",
  stan_options = list(chains = 4, iter = 2000, warmup = 1000, cores = 4,
    seed = 604, control = list(adapt_delta = .95, max_treedepth = 13)),
  nutpie_options = list(adaptation = "diag")
)
plan <- posterior_policy_draws(fit, n_draws = 200, seed = 701)
# Reuse plan as draw_plan in both merger and tariff simulate() calls.
```

`backend = "cmdstanr"` gives the original route. `inference = "pathfinder"`
remains CmdStan-only; Nutpie initially supports `inference = "hmc"`.
`stan_options` supplies chains, total iterations, warmup, cores, seed, refresh,
initial values and warmup retention; `control$adapt_delta` and
`control$max_treedepth` translate to Nutpie target acceptance and maximum depth.
Nutpie requires positive warmup, no thinning, and one thread within a chain.
Unsupported controls fail explicitly. Its adaptation choice belongs in
`nutpie_options`: `"diag"` (default) or `"low_rank"`. Low-rank
`mass_matrix_gamma` and `mass_matrix_eigval_cutoff` are native Nutpie controls;
CmdStan `dense_e` is a different adaptation strategy and is never translated
implicitly. Native diagnostics and the resolved controls are retained.

Expert debugging can retain divergence states, gradients, unconstrained
coordinates and warmup mass-matrix traces. These options can substantially
increase fit size and are off by default:

```r
debug_options <- list(adaptation = "diag", store_divergences = TRUE,
  store_mass_matrix = TRUE, store_unconstrained = TRUE, store_gradient = TRUE)
# Pass nutpie_options = debug_options and
# stan_options = list(..., save_warmup = TRUE) to bayes_est().
# Inspect fit@stanfit$raw_diagnostics and fit@stanfit$warmup_diagnostics.
```

`mindepth`, `extra_doublings` and `max_energy_error` are additional native
expert controls. The package does not tune them automatically. Python flow
adaptation and `draw_diag` are not exposed by this interface.

The adapter requires `nutpieR >= 1.9.1`. Prefer a compatible prebuilt package
from the [official R-universe repository](https://andytimm.r-universe.dev):

```r
install.packages("nutpieR",
  repos = c("https://andytimm.r-universe.dev", "https://cloud.r-project.org"))
```

Check the installed version. If that repository does not yet provide the
required version for your platform, the audited source revision is:

```r
remotes::install_github("andytimm/nutpieR",
  ref = "8e88a975d5e0b8ed3c76de84d22adb33de601a1e", upgrade = "never")
```

Source installation needs Rust >= 1.85; supported prebuilt binaries avoid that
requirement. Runtime Stan compilation needs a C++ toolchain. No Python or PyMC
is used. Set `options(antitrustBayes.nutpie_cache = "/writable/cache")` for
compiled artifacts; installed package files are never compiled in place.
Bridge evidence also needs a BridgeStan R client matching the compiled target
ABI. The executed stack used nutpieR 1.9.1 and BridgeStan 2.8.0; its client can
be installed separately with
`remotes::install_github("roualdes/bridgestan", ref = "v2.8.0", subdir = "R", upgrade = "never")`.
Other toolchains must satisfy the runtime ABI, coordinate and density checks.
Saved fits retain ordinary draw arrays, diagnostics, source/data SHA-256,
model specification, resolved priors, versions and controls. Posterior and
policy operations on a reloaded fit use those arrays without live native
pointers. Evidence evaluation additionally requires the optional `bridgestan`
R client. A verified compatible saved library is reused; if it is missing or
unusable, the exact matching installed Stan kernel and saved data are rebuilt
in the external cache through the compatible optional toolchain. Changed
source/data fail closed. Reconstruction records its own artifact provenance
without rewriting the saved fit; the new binary need not match the old byte hash.
CmdStanR's CSV precision remains unchanged (eight significant digits in the
validated toolchain); CmdStan users can request `stan_options = list(sig_figs = 18)`
for strict algebraic checks. Nutpie retains native double arrays.

The [evidence experiment](doc/validation/nutpie_evidence_results.md) passed
native target, bridge-uncertainty and model-probability checks. Its historical
zero-divergence qualification criterion remains recorded. Public `bayes_compare()`
now supports Nutpie and mixed CmdStan/Nutpie evidence under the same scientific
contracts. Sampler choice is computational provenance, retained model by model;
it must not be selected to obtain a preferred policy result or model probability.

The [public closeout validation](doc/validation/nutpie_closeout_summary.md)
completed all 42 tasks, including fixed-k CF/CS evidence, a full conditional-k
family with mixed backends, and exact-source reconstruction for all four kernels.
Conditional merger/tariff means agree within the declared precision bounds.
Some SD, tail and threshold comparisons remain precision-inconclusive at the
fixed draw budget; their uncertainty and realization attrition remain explicit.
This does not establish universal sampler quality or precise policy-tail equivalence.

Finite evidence, bridge repetition completeness and sampler diagnostics are
reported separately. Every attempted repetition is retained; finite estimates
survive other failed repetitions, with SD unavailable when fewer than two
repetitions succeed. Every positive-prior candidate still requires finite
evidence before model probabilities exist. `bayes_weights(comp, "bridge")`
requires complete bridge probabilities and never substitutes stacking or
renormalizes surviving candidates; request `"stacking"` explicitly for that
predictive method. A finite estimate alone does not establish accuracy.

Predictive scoring uses the common `log_lik` interface where the existing
likelihood permits LOO/WAIC; collapsed-quality restrictions remain. Conditional-k
call generation forwards `backend` and `nutpie_options`. Components may use
different samplers, retain every original k state and restore the original k
prior exactly once. No `backend_by_k` API is needed: edit separately scheduled
component argument lists when different samplers are appropriate.

See the [evidence contract](doc/validation/nutpie_evidence_design.md) and
[paired validation runner](doc/validation/nutpie_backend_validation.R).
The runner uses fixed package fixtures, records unavailable tasks and raw
warnings, separates compile time, and requires explicit opt-in for sampling.
The executed timing/ESS tables describe these fixtures and matched draw budgets;
they do not establish a general efficiency or robustness advantage.

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

## Logit utility normalization

Every Logit market has a free inside-market utility intercept, including one
market observed in one or several periods. Outside-good utility remains zero,
product qualities stay centered, and prices retain their standardization. A
single-market intercept has a proper `Normal(0, 1)` prior with no hierarchy SD;
multiple markets retain the existing hierarchical specification.

Version `0.2.0.9007` corrects the former single-market zero-intercept restriction,
which incorrectly tied the outside share to the geometric mean of conditional
inside shares. Earlier fits require refitting under model specification 12.
The intercept and outside-share log odds enter conditional demand as a sum;
supply information and priors determine their separate identification. See the
[normalization and identification contract](doc/single_market_utility_normalization.md).

`bayes_calibrated_market(..., outside_share = 0.20)` now accepts the outside
share as an independent baseline primitive (default 20%), while preserving the
chosen conditional-share geometry. Subsequent equilibrium shocks may change
the outside share.

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

## HMC controls

Supply sampling controls through `bayes_est(..., stan_options = list(control =
list(...)))`. Partial control lists preserve the package's `adapt_delta = .95`
and `max_treedepth = 13L` defaults unless those fields are explicitly overridden.
For example, a pooled Cournot fit can use a dense metric:

```r
fit <- bayes_est(
  prices = dat$prices, quantities = dat$quantities, margins = dat$margins,
  market = dat$market, year = dat$year, group = dat$group,
  demand_model = "logit", market_side = "output", supply_model = "cournot",
  use_alpha_market_pool = TRUE,
  stan_options = list(control = list(metric = "dense_e"))
)
```

A dense metric can help when posterior parameters are strongly correlated.
Compare divergences, R-hat, ESS, BFMI and elapsed time for the application;
it can also take longer. Both `metric = "diag_e"` and the existing
`use_centered_alpha` argument remain available. These controls affect the
sampler while preserving the model and its priors.
Setting a control field to `NULL` delegates that field to CmdStan's default;
`control = NULL` delegates all controls. An empty `control = list()` retains
the package defaults.

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
- [coordination](https://github.com/luciu5/coordination) (>= 0.0.0.9003; the selected revision may require a newer `antitrust`)
- [loo](https://mc-stan.org/loo/) (>= 2.5.0)
