# Explicit rate-domain simulation

The Bayesian posterior is fitted once and reused. `price_domain` is a
deterministic reconstruction/simulation option and is excluded from fitted
model, Stan data, and posterior draw-plan identities. A saved plan therefore
realizes the same structural draw, quality vector, and market-level k under
both domains. Neither domain redraws quality between policy states.

## Responsibilities

| Layer | Responsibility |
|---|---|
| coordination | Domain validation; native input Logit Cournot CF/CS equilibrium; original-FOC residual checks; dollar markdowns; typed positive-domain obstruction; revenue diagnostic |
| antitrust | Existing quantity, producer-surplus, and utility-based welfare methods; economic counterfactual definitions |
| antitrustBayes | Explicit API forwarding; preserved posterior realization; latent pre-equilibrium; stage-specific failure records; result aggregation |

No Stan or economic counterfactual specification changes are needed.

## API

`simulate()`, `bayes_sim()`, `bayes_core_fringe_model()`, and
`bayes_stackelberg_model()` accept `price_domain = c("positive", "real")`.
The first value remains the default. Real mode requires coordination >=
0.0.0.9004 and initially supports input-market Logit Cournot on native exact
fringe or Stackelberg routes. Other combinations are rejected before posterior
extraction. Ordinary simultaneous, smooth fringe, CES, output markets, and
Bertrand retain their existing domains.

The deterministic model stores the selected domain in its existing solver
control and diagnostics; no new S4 slot is required. An old serialized model
without this field uses positive mode. Observed calibration inputs remain
positive. Real mode changes the feasible equilibrium domain explicitly; it
does not model a zero floor or establish a constrained equilibrium.

## Failure and result contracts

Every requested draw retains its position, including result-assembly failures.
The existing failure message remains available. Native simulation metadata
also records the stage, condition category, posterior index, domain, and
minimum rate/product when supplied by the deterministic condition.

`latent_pre_equilibrium` and `post_equilibrium` distinguish a changed latent
baseline from the policy solve. Other stages are `baseline_construction`,
`scenario_construction`, and `result_assembly`.
An unsuccessful parallel worker is retained with stage `worker_process`.

`price_level_change` remains meaningful at finite real rates. Relative changes
from a zero baseline are NA and flagged by
`undefined_relative_price_change`. Ratios from negative baselines retain
their algebraic value and are flagged by `signed_relative_price_change`.
`policy_distribution_conditioning` states whether all requested draws or only
successful draws contributed. Failed draws are retained and never replaced
by zero or a different posterior state.

Native quantity shares remain probabilities. Revenue/expenditure ratios at
nonpositive active rates are NA with an explicit diagnostic. Dollar markdowns
are evaluated directly, so producer surplus remains finite at zero; normalized
markdowns are undefined at zero. The inherited welfare formulas are preserved.
For real Cournot Stackelberg, residuals normalized by a zero implied rate are
also NA, with `undefined_normalized_zero` identifying the affected products.
Aggregate normalized maxima are NA when their portfolio contains an undefined
term. The original raw residual and strict price FOC remain the acceptance
criteria; neither gate is changed by this reporting rule.

## Source map

- `R/price_domain.R`: compatibility validation and failure-condition metadata.
- `R/stackelberg_bridge.R`: public native reconstruction and domain forwarding.
- `R/latent_quality_simulation.R`: domain-aware solved baseline validation.
- `R/simulate.R`: per-draw stage capture and complete result-position retention.
- coordination native model files and a shared scalar-W helper: deterministic
  implementation, with a separate positive-domain solver path.

Raw fits, simulation objects, and reference matrices stay in ignored local run
directories. Only compact source notes, reference code, and curated validation
tables are candidates for Git.
