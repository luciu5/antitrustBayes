# Stackelberg extension reconnaissance

Read-only inventory for the requested coordination/antitrustBayes extension.
No production files were changed.

## Repository state and local instructions

- `coordination`: `/mnt/d/Projects/coordination`, branch `main`, clean at `6935013` (`Install exact antitrust checkout in CI`).
- `antitrustBayes`: `/mnt/d/Projects/antitrustBayes`, branch `main`, `023617d` (`Update installation instructions after branch merge`). The only current untracked files are documentation under `doc/` from the lead work.
- `antitrust`: `/mnt/d/Projects/antitrust`, branch `refactor`, clean at `3ab5f37` (`Fix nested Logit shares for inactive nests`).
- `/mnt/d/Projects/antitrust/AGENTS.md` requires preserving legacy S4 economics as the behavioral source of truth, keeping calibrate/update/respecify/counterfactual/simulate semantics distinct, preserving uncommitted changes, and reporting tests honestly. No repository-specific `AGENTS.md` was found in the coordination or antitrustBayes trees.

## coordination: current structural path

`coordination/NAMESPACE` imports antitrust's `calcMC`, `calcMargins`, `calcPrices`, `calcQuantities`, `calcShares`, `calcSlopes`, ownership conversion, and simulation-related methods, and exports its own `PriceLeadership`, `PriceLeadershipBLP`, `ple`, `ple.blp`, `calcMargins`, `calcPrices`, `calcSlopes`, and `ownerToMatrix` interfaces.

The existing classes are in `coordination/R/PriceLeadershipClasses.R:14` and `:135`. They extend antitrust demand classes and store `supermarkupPre/Post`, `timingParam`, `coalitionPre/Post`, `bindingFirm`, and `slackValues`. They implement the coordinated-effects model with supermarkups and incentive-compatibility constraints. Their constructors and calibration paths are in `coordination/R/PriceLeadershipFunctions.R` (constructor around `:114`, `calcSlopes` around `:513`, margins around `:1007`, prices around `:1107`) and `PriceLeadershipBLPFunctions.R`. These must remain untouched in economic meaning.

Current coordination tests and fixtures:

- `tests/testthat/helper-fixtures.R`: `fixture_ple` (three-product coalition plus one fringe product), `fixture_bertrand`, and deterministic BLP fixture.
- `tests/testthat/test-price-leadership-invariants.R`: fringe Bertrand best response, coalition supermarkup relation, slack/timing and finite-output checks.
- `tests/testthat/test-price-leadership-parity.R`: known numerical parity and export checks.
- `tests/testthat/test-price-leadership-blp.R`: BLP behavior and dispatch.
- `tests/testthat/test-s4-dispatch.R`: loading coordination must not redefine ordinary antitrust `Logit`, `Bertrand`, or `Cournot` methods.
- `tests/testthat/test-antitrust-boundary.R`, `test-grim-trigger.R`, and the boundary/parity tests protect the package boundary.

Implementation constraint: add separate Stackelberg classes/methods (the lead's map proposes `StackelbergLogit` and `StackelbergCES`) and dispatch only on those classes. Do not add methods that alter ordinary antitrust classes or reuse `PriceLeadership`.

## antitrust: legacy demand and quantity Stackelberg

Relevant definitions are in:

- `R/BertrandRUMClasses.R:128`: `Logit` class and slots (`prices`, `margins`, `priceStart`, `normIndex`, `shareInside`, `weights`, `priceOutside`, `mktElast`, `insideSize`, `mktSize`); `:564` contains `CESCournot`.
- `R/LogitCournotClasses.R:3`: `LogitCournot` extends `Logit`.
- `R/LogitFunctions.R:241`: `logit()` constructor; `:314`: `logit.cournot()`.
- `R/CESFunctions.R:185`: `ces()`; `:418`: `ces.cournot()`.
- `R/CournotFunctions.R:334`: legacy generic `stackelberg()` constructor.
- `R/CournotClasses.R:177-230`: legacy `Stackelberg` class, extending `Cournot`, with plant/product-level `isLeaderPre/Post` matrices and quantity-dependent `dmcfunPre/Post` slots.
- `R/OutputMethods.R:60`: simultaneous Cournot quantities; `:157-250`: legacy Stackelberg quantity calculation using leader indicators and quantity-dependent costs; `R/PricesMethods.R:49`: Cournot inverse demand and `:82`: Logit price solver.
- `R/MarginsMethods.R:34`: Bertrand price FOCs and output/input signs; `:92-136`: `LogitCournot`; `:139-182`: `CESCournot`.
- `R/CostMethods.R:174`: Bertrand MC recovery; `:293`: Cournot MC; `:423-452`: legacy Stackelberg `calcdMC`.
- `R/ParamsMethods.R:45`: Cournot calibration; `:230`: legacy Stackelberg calibration.

The legacy `Stackelberg` is linear/log-linear quantity conduct with product/plant role matrices and quantity-dependent cost functions. It is not the requested firm-level Logit/CES price/quantity game. It is exported and used by old tests, so do not redefine it or change its registry semantics. There is also a likely legacy validity typo around `R/CournotClasses.R:211` (`dmfunPost` versus `dmcfunPost`); flag for separate maintenance rather than changing it during this extension.

The modern architecture is in `R/ModelRegistry.R` and `R/ModelArchitecture.R`. The registry currently has `linear::stackelberg` and `loglin::stackelberg`, but no Logit/CES Stackelberg. `calibrate()` only accepts the legacy linear/log-linear Stackelberg path and rejects post-scenario fields during calibration. Simulation has legacy Stackelberg handling around `ModelArchitecture.R:735-796`; old registry/parity tests (`tests/testthat/test-model-registry.R`, `test-legacy-model-parity.R`) must remain unchanged.

## antitrustBayes: input, cutoff, Stan and reconstruction paths

Data preparation:

- `R/helpers.R:193-407`, `prepare_stan_data()`: validates prices/quantities/margins, encodes market/year/group, calculates conditional product shares, and computes `firm_shareIn` by summing product shares within firm-market-year before repeating that aggregate on each product row (`:268-278`).
- `R/helpers.R:419-649`, `prepare_ces_stan_data()`: output CES only; computes quantity shares, revenue shares, inside price index, and repeated `firm_revShare` per product row (`:501-523`); outside price/revenue share are included in the Stan data.
- `tests/testthat/test-helpers.R` checks data fields, dimensions, margin transforms, input validation, and model-spec version.

The repeated firm aggregate is valid for a product's own FOC response, but summing it over follower rows double-counts multiproduct firms. A Stackelberg data path needs one firm-market-year record and an observation-to-record map (e.g. `N_firm_mky`, `firm_mky_shareIn`/`firm_mky_revShare`, `firm_mky_market_year`, `firm_mky_group`, `obs_to_firm_mky`) and tests with multiproduct firms.

Estimation/API:

- `R/estimate.R:218-245` defines `bayes_est()` with current `demand_model`, `market_side`, `supply_model`, and `use_cutoff`; supply codes are mapped around `:265-272` to Bertrand, MonCom, Cournot, auction, bargaining. `prepare_stan_data()`/`prepare_ces_stan_data()` are called around `:393-451`. `use_cutoff` is explicitly ignored for MonCom (`:288-292`).
- `R/params.R:1-130`: `bayes_params()` stores alpha/gamma, outside share, cutoff share, cutoff scale (default `0.2`), and supply model. `bayes_specify()` is around `:211-269`.
- `R/params.R:275-351`: `build_conduct_object()` reconstructs a custom `BayesConductLogit` or `BayesConductCES`; it converts conditional to unconditional shares and recovers Logit/CES mean values using the existing outside normalization. `R/params.R:353-545` routes registered no-cutoff models through antitrust, while cutoff fits remain custom because the blend is not a registered antitrust model.
- `R/AllClasses.R:1-110`: `BayesFit`/`BayesLogit`/`BayesSpecifiedFit` and custom `BayesConductLogit`/`BayesConductCES` classes. The custom classes currently carry supply model, cutoff and bargaining slots, but no Stackelberg leader partition.

Cutoff/conduct equations:

- `R/conduct.R:1-158`: `conduct_weight()` is a logistic weight with orientation 0 for large/leader-like shares and 1 for small/fringe-like shares; `conduct_state()` aggregates product shares to firms; custom Logit/CES `calcMargins` methods blend the existing simultaneous/MonCom/Cournot/auction/bargaining schedules row-wise.
- `inst/stan/logit_supply.stan`: fixed logistic scale `.2` around `:29-35`; `logit_inv_markup()` around `:38-62` converts repeated conditional firm share to unconditional share and blends old supply equations. The likelihood and generated quantities around `:364-469` duplicate the old equations and include `mu_b_fringe` as a demand-quality shift. This shift is fringe-specific and should not be copied automatically into Stackelberg.
- `inst/stan/ces_supply.stan`: analogous fixed `.2` scale (`:28-35`), repeated conditional firm revenue shares, old schedule blend, and duplicated likelihood/generated-quantity equations (`:332-440`).
- `R/simulate.R:222-259`: `.estimated_draw_params()` retrieves posterior cutoff but hard-codes `cutoff_scale = 0.2` at `:255`; `:348-620` builds a custom Bayes object and delegates current counterfactuals to `antitrust::simulate()`. It has no cutoff-to-hard-leader mapping or coordination reconstruction.

Compatibility requirements: append/introduce a new explicit structure argument while retaining old `use_cutoff` behavior (`TRUE` means existing fringe path when structure is omitted; `FALSE` means existing simultaneous path). Preserve old `.2` fringe behavior for reproducible workflows; expose a separate prespecified Stackelberg transition scale rather than estimating another weakly identified parameter. Validate Stackelberg support only for Bertrand/Cournot; current CES estimation is output-only (`prepare_ces_stan_data()` rejects input CES around `:437-441`). Keep generated quantities algebraically identical to the likelihood.

Posterior/simulation tests already present:

- `tests/testthat/test-params.R`: custom object construction, antitrust class routing, cutoff-neutral baseline and cost-shock/FOC recovery.
- `tests/testthat/test-conduct.R`: cutoff behavior, including recomputation from post-merger firm shares.
- `tests/testthat/test-model-spec.R`: Stan source assertions (`firm_shareIn`, `firm_revShare`, `R0_prob`, Student-t likelihood, log-lik fields).
- `tests/testthat/test-simulate-dispatch.R`: specified/estimated simulation and optional Stan tests (`ANTITRUSTBAYES_RUN_STAN_TESTS=true`).
- `tests/testthat/test-bma.R`: bridge sampling/model-weight fallback.
- `tests/testthat/test-validation.R`: recovery/validation utilities.

## PNB data and scripts

- Data: `data/pnb.rda`, loader/metadata in `R/data.R`. This is a single market/year Philadelphia National Bank case with 39 banks and input deposits; attributes include acquirer `298816`, target `889018`, and `insideSize`.
- `data-raw/make_pnb.R:1-34`: builds the cleaned event from an external Supreme Court data object, drops invalid/low margins, and orders firms by quantity.
- `data-raw/run_pnb.R`, `run_pnb_all.R`, `run_pnb_logit_all.R`, `run_pnb_ces_all.R`: standard Bertrand/Cournot/auction/MonCom estimation and merger simulation.
- `data-raw/run_cutoff.R:5-73`: existing `use_cutoff=TRUE` estimation, bridge summaries, and merger simulation.
- `vignettes/pnb-merger.Rmd:36+` and simulation around `:137-140`: documented PNB-Girard workflow.

Read-only inspection of `pnb.rda` shows 39 observations, 39 unique groups, one product per group, quantity sum about `4845.345`, `insideSize` `4845345000`, and observed margins approximately `0.083` to `1.840` (six at least one). Thus PNB smoke fits do not exercise multiproduct aggregation; synthetic multiproduct tests are essential. No new PNB fit was run in this reconnaissance.

## Compatibility and test matrix proposal

Coordination tests should cover each `demand × conduct` pair (`Logit`/`CES` × `Bertrand`/`Cournot`) with 3–5 single-product firms; one and multiple simultaneous leaders; several followers; no followers/all leaders; followers approaching zero share; multiproduct leader and follower; pre-equilibrium reproduction; FOC residuals; finite/positive margins; ownership-only merger, cost/tariff shock, and combined shocks; distinct pre/post role partitions; and firm-level role expansion to all products. Compare follower-response and leader derivatives against finite differences of fully solved equilibria, direct optimization/brute-force one-leader benchmarks, and generic-vs-closed Logit Jacobians. Add nesting tests against ordinary simultaneous Bertrand/Cournot where applicable and ensure leaders maximize individually rather than joint profit. CES tests must exercise outside revenue share, outside price normalization, gamma restrictions, and implicit Jacobian failures.

Bayes tests should retain old simultaneous outputs and old `use_cutoff=TRUE` fringe output, reject unsupported Stackelberg combinations, verify old positional calls, verify unique firm-market-year aggregation/no double counting, check role-weight direction and scale metadata/sensitivity, ensure likelihood/generated quantities use the same equations, and test bridge compatibility. Add synthetic posterior summaries market by market: role probabilities, hard `k` from each cutoff draw, mean/median/interval, strategic-core share, and cutoff recovery. Test posterior-draw reconstruction into coordination under hard threshold plus sensitivity near boundary firms, ownership/cost/combined shocks, valid diagnostics, and multiple-seed recovery/discrimination/merger/pass-through simulations. Run legacy antitrust Stackelberg linear/log-linear parity and full existing suites to verify no registry or dispatch changes.

PNB is a later smoke diagnostic only after package tests: compare simultaneous/fringe/Stackelberg Bertrand and Cournot alongside auction/MonCom, reporting cutoff, implied-k, core share, fit/marginal likelihood, sampler diagnostics, margin fit, merger effect, and cutoff/k sensitivity. No application prior or tuning should be changed opportunistically.

## Findings classification

- **Observed in code:** class/method locations, API mappings, repeated firm aggregates, fixed `.2` cutoff, custom reconstruction path, legacy Stackelberg scope, existing tests, and PNB data shape.
- **Implementation constraints inferred from compatibility:** preserve legacy classes/registry/dispatch and use unique firm-market-year records for follower-sector sums.
- **Not evaluated here:** Stackelberg economic formulas, smooth-role coherence, numerical derivative checks, new deterministic solver, Stan fit/recovery, and PNB estimation. Those belong to the lead's signed math/architecture phases.
