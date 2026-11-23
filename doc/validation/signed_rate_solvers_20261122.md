# Signed rates and counterfactual solver comparisons, 2026-11-22

These comparisons use the local development sources of antitrust, coordination,
and antitrustBayes. Eligible input Logit and BLP simulations accept finite
signed rates by default. A positive-rate restriction can be selected explicitly.
CES uses positive prices. These choices apply to counterfactual simulation.

## Matched deterministic timings

Times are milliseconds per call, using one thread, 20 calls per batch, and
medians of three paired batches with alternating order. Each pair used
identical fitted primitives, ownership, and policy changes.

| Calculation | Reference | New | Maximum absolute price difference |
| --- | ---: | ---: | ---: |
| Output Logit Stackelberg Cournot, 6 products | 18.15 | 1.70 | 1.8e-15 |
| Output Logit Stackelberg Cournot, 24 products | 70.00 | 1.90 | 3.6e-15 |
| Output Logit bargaining | 1.75 | 1.10 | 4.5e-9 |
| Output CES bargaining | 13.30 | 1.95 | 1.2e-9 |
| Input BLP Bertrand | 17.85 | 3.45 | 4.2e-9 |
| Output BLP bargaining | 4.40 | 2.60 | 4.4e-9 |
| Input BLP bargaining | 4.60 | 3.75 | 6.0e-9 |

The Stackelberg rows compare the previous product-price root with the checked
firm-level solution on an outside-good Logit fixture with two leaders,
multiproduct firms and a 5% cost reduction. They time the equilibrium solve.
The remaining rows time a complete single-market `antitrust::simulate()` call
on a fitted model. For BLP Bertrand the shock produces negative equilibrium
rates. The BLP bargaining cases use unequal probability weights, changing
bargaining power and retention, and a product exit. No HMC or fit construction
is in these timings.

The Bayesian Logit bargaining model has its own firm-level margin equation.
Three matched supplied draws on each market side completed the same merger
under the ordinary solver and `solver = "zeta"`. All six were valid and the
maximum price difference was 5.4e-9. The upstream flat-Logit bargaining timings
measure a different margin equation.

## Verification

| Suite | Passed expectations | Warnings | Skipped tests | Failures or errors |
| --- | ---: | ---: | ---: | ---: |
| antitrust, extended | 1,738 | 148 | 45 | 0 |
| antitrustBayes, ordinary | 8,509 | 6 | 17 | 0 |
| coordination, ordinary | 743 | 15 | 0 | 0 |
| vertical, lifecycle dispatch and linked exit | 118 | 6 | 0 | 0 |

The regressions exercise input and output equilibria, explicit positive-rate
restrictions, product exits, alternative Stackelberg solves, and resuming
policies on native and downstream models. The Bayesian comparison used supplied
draws. These runs measure deterministic counterfactual performance only.

The local scripts, console timings and package test logs are retained in
`doc/statistical_toolkit_runs/signed-games-20261122/`.
