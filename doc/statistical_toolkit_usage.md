# Experimental policy simulation and model universe

Estimated-fit policy simulations default to the fitted posterior state through
`demand_anchoring = "latent_quality"`. Observed-share calibration is the explicit
alternative `demand_anchoring = "observed"`. Specification 14 and its inference
contracts are unchanged; changing the simulation default does not require a refit.
`BayesSpecifiedFit` has supplied primitives rather than an estimated quality
posterior and retains its supplied-state/calibration default.

## Reuse posterior latent states across policies

```r
plan <- posterior_policy_draws(fit, n_draws = 100, seed = 20260916)

merger <- antitrust::counterfactual(ownership = owner_post)
cost_shock <- antitrust::counterfactual(costs = cost_change)
# cost_change is a named product vector of proportional cost changes;
# for example, rep(0.10, n_products) with the product labels as names.

observed <- bayes_sim(
  fit, ownerPost = merger, draw_plan = plan,
  demand_anchoring = "observed"
)
latent <- bayes_sim(
  fit, ownerPost = merger, draw_plan = plan
)
latent_cost <- bayes_sim(
  fit, ownerPost = cost_shock, draw_plan = plan,
  demand_anchoring = "latent_quality"
)
```

The plan selects common posterior draws, one persistent quality vector per quality block and
draw, and the conditional market-level ordered-k state where applicable.
Reuse the plan to hold these states fixed across policies and periods.
Input Logit is supported using the native positive dollar markdown and held-fixed
value primitive. Input CES and unsupported conduct/likelihood combinations fail explicitly.
Second-score Logit auctions use their cost-based value index, with the correct
input/output sign, rather than a posted-price demand inversion. Supported
bargaining routes retain their existing firm-share equations and cost/value
recovery. Neither route silently falls back to observed-share anchoring.
Relative Gaussian fits support active-rho quality reconstruction. For input Logit, an
antitrust `costs` counterfactual changes the value primitive; it is not a tariff
pass-through calculation.

The latent baseline projects relative fitted quality to shares while retaining
the outside-share draw at reference prices. It preserves the costs recovered
under the observed-share anchor, then solves its own baseline equilibrium.
Baseline prices can therefore differ from recorded prices. Relative-share fits derive one intercept per market-year from the same quality
state. Their observation measure is the orthonormal log-ratio coordinate space;
absolute-share fits retain the raw log-share working likelihood.

## Compare disjoint economic hypotheses through bayes_compare

```r
fits <- list(MC = fit_mc, B = fit_b, C = fit_c,
             BF = fit_bf, CF = fit_cf, BS = fit_bs, CS = fit_cs)

universe <- bayes_model_universe(
  list(MC = "MC", B = "B", C = "C",
       BF = list(fit = "BF", k = "interior"),
       CF = list(fit = "CF", k = "interior"),
       BS = list(fit = "BS", k = "interior"),
       CS = list(fit = "CS", k = "interior")),
  prior_weights = rep(1 / 7, 7)
)
comparison <- bayes_compare(fits, universe = universe, bridge = TRUE)
bayes_dimensions(comparison, dimension = "conduct")
bayes_dimensions(comparison, dimension = "strategic_core")
```

The equal weights in this example are a user choice. The API requires explicit
scientific weights and does not redistribute them after inspecting results.
Interior means `0 < k < n_market` in every market. Full exact-family fits can
include endpoints; restricted evidence is obtained using the normalized
encompassing-prior identity. Shared endpoints are represented by MC, B, and C.
The source fits must satisfy the common data/statistical specification contract.

The universe currently provides evidence and dimension probabilities. Policy
simulation for a restricted family additionally requires continuous posterior
draws to be reweighted by their restriction probabilities. An unrestricted
draw plan with k alone restricted is insufficient. This additional conditional
policy-BMA facility is not implemented by the evidence interface.

## Hybrid likelihood

```r
fit_hybrid <- bayes_est(
  prices = p, quantities = q, margins = m,
  market = market, group = owner,
  share_likelihood = "gaussian",
  quality_parameterization = "collapsed",
  margin_likelihood = "student_t",
  df = 5,
  use_rho = FALSE
)
```

All seven frozen small hybrid benchmark fits passed the prespecified HMC gate.
This supports availability of the option; it does not yet establish a routine
default or a robustness advantage in a margin-outlier application.
