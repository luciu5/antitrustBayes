# Collapsed-Gaussian architecture

This memo records the reviewed implementation map and approved architecture for
the collapsed Gaussian share-quality experiment. The detailed reconnaissance
below refers to native specification 8, the baseline inspected before changes.

## Approved public contract

`bayes_est()` gains `share_likelihood = c("student_t", "gaussian")` and
`quality_parameterization = c("sampled", "collapsed")`. It also gains an
explicit optional `margin_likelihood`; when omitted it matches the share
family. The production default remains Student-t shares and margins with
sampled quality until the benchmark is complete.

| mode | share | quality | margin |
|---|---|---|---|
| A | Gaussian | sampled | Gaussian |
| B | Gaussian | collapsed | Gaussian |
| C | Student-t | sampled | Student-t |
| D, optional | Gaussian | collapsed | fixed-df Student-t |

Collapsed quality is supported only with Gaussian shares. Mode D requires
`use_rho = FALSE`. The reference implementation initially rejects collapsed
quality with market-specific residual scales. Gaussian collapsed nonzero rho
uses the full Schur-complement likelihood; exact-k plus nonzero rho is initially
restricted to one market because global quality centering otherwise couples
the markets' conditional margin likelihoods. Every unsupported combination
receives an explicit R validation error.

The resolved families and parameterization are stored in Stan data, fit
metadata, diagnostics, and prior-normalization provenance. This remains one
statistical architecture for one-period, panel, and multiple-market data;
additional data use the same economic equations and existing pooling flags.

The public outside-share policy also has an information-aware selector,
`s0_pooling = c("auto", "none", "market", "market_year")`.  `auto` resolves
to one directly parameterized outside share in a one-market panel and to the
existing hierarchical market-level outside shares when several markets are
present.  `none` retains independent market-year parameters with their common
fixed prior, while `market_year` selects the existing global market-year
hierarchy.  Explicit legacy pooling flags map to the same modes, so a benchmark
can freeze the old unpooled specification without relying on the data-dependent
default.  A one-market market-level mode does not estimate a between-market
variance from one unit: it gives the common market logit outside share the
proper direct prior and lets `X_s0` carry any prespecified period variation.
This single-market special case removes the former unidentified hyper-scale
and standardized market effect. Its direct outside-share prior equals the
unpooled single-cell prior, rather than the wider marginal prior induced by
the former one-unit hierarchy. This statistical change is confined to the
single-market market-pooling mode; the frozen A/B/C experiment explicitly uses
unpooled controls and is unaffected.

## Nested information ladder

The implementation is organized as one likelihood with increasingly
informative designs, rather than a PNB-specific model and separate panel
models.  The observation-to-quality design matrix is always built from the
supplied persistent `group` identifier.  Repeating a group in later periods or
markets repeats the corresponding row of that design matrix and therefore
adds the covariance implied by one persistent quality effect.  The same Stan
program evaluates every rung.

| rung | data added | statistical information added | optional pooling |
|---|---|---|---|
| 0 | one market, one period | composite cross-sectional share dispersion | one directly parameterized outside share |
| 1 | later periods in the same market | repeated observations of persistent group quality, separating persistent quality from transitory share error | one direct market outside share plus optional `X_s0` variation |
| 2 | additional markets | replication for common variance scales and market intercepts | hierarchical market price coefficients and market outside shares |
| 3 | markets by periods | both within-group persistence and across-market replication | hierarchical market outside shares plus `X_s0`; market-year latent effects remain explicit opt-in |

At rung 0, `use_market_resid_scale`, `use_alpha_market_pool`,
`use_s0_market_pool`, and `use_s0_market_year_pool` are explicitly false in the
controlled A/B/C benchmark.  That is a frozen experimental setting and does
not rely on the public `auto` preset.  Higher rungs use the same likelihood and
may activate its pooling dimensions; they do not change the demand mean,
conduct equations, or treatment of ordered k.  Exact k remains a market-panel
state, so every period in a market uses the same predetermined firm ranking and
k state.

The first production implementation deliberately reproduces the package's
current *global* Helmert quality coordinates.  Thus `group` labels are global:
the same label in two markets means the same persistent entity, and unrelated
market-specific firms should have distinct labels.  Global centering also
induces the small cross-group covariance implied by the existing sum-to-zero
normalization.  Replacing it with separately centered market blocks would be a
different statistical model; it is not silently selected at a higher ladder
rung.  If later evidence favors market-block quality populations, the same
design-matrix interface can accept a block contrast basis without changing
the economic model or creating a separate Stan program.

## Target and generated-quantity contract

In collapsed mode the normalized multivariate share density is evaluated once
per posterior state. At rho zero, exact BF/CF/BS/CS retain one complete
market-level `log_sum_exp` per market. At nonzero rho the Gaussian margin target
uses the exact conditional multivariate normal in the math note.

Collapsed shares have no honest observation-wise likelihood factor after
quality integration. Stan exposes the full joint share log density and retains
the exact-k market margin/component outputs. LOO/WAIC is rejected for collapsed
fits until a valid held-out-group or held-out-market calculation exists. Bridge
sampling remains supported because the target retains every determinant,
Gaussian constant, and normalized prior.

The bridge target is a proper density for the package's recorded transformed
outcomes.  In Gaussian-rho mode it is evaluated as the normalized share
marginal times the normalized conditional margin density, which is the same as
the normalized joint Gaussian.  Exact-k weights are categorical probabilities
that sum to one within every market.  Fixed-df Student-t terms retain their
normalizing constants.  Positive-scale and bounded-correlation priors retain
their truncation or transformation normalizers.  The existing HMT soft penalty
is the explicit exception: it is not a normalized probability model, and
bridge comparison continues to reject `use_hmt = TRUE`.  Comparisons also
require the same transformed observations; this contract does not claim that
evidence values from different outcome transformations are interchangeable.

`posterior_quality()` reconstructs conditional contrast and centered group
quality draws after fitting. HMC never samples these coordinates in collapsed
mode.
For active Gaussian rho, automatic recovery conditions on both shares and
matched margins using the same sparse augmented-row representation. Exact-k
rho recovery currently requires explicit share-only conditioning: component
margin predictions cannot be replaced by their posterior average when
conditioning quality. This restriction affects recovery, not the supported
one-market exact-k rho fitting target.

## Structured production likelihood

Production Stan code must not construct or factor an \(N\times N\) collapsed
covariance. The dense Cholesky calculation is an independent small-\(N\)
oracle only. The production density uses the one-nonzero-per-row group design,
diagonal residual variance, within-group diagonal-plus-rank-one blocks, and the
single global rank-one correction induced by the centered quality constraint.
It accumulates per-group sufficient statistics and evaluates the exact log
determinant and quadratic form in \(O(N+G)\) time. Group statistics use \(O(G)\) memory; the current Stan wrapper and differentiation tape also retain linear-size observation vectors, so total production storage is \(O(N+G)\).

The rho path uses the same function. A unit-determinant triangular
transformation turns every matched margin into a pseudo-observation with a
known diagonal residual variance and a scalar loading on the same group
quality. Thus the joint Gaussian calculation costs \(O(N+N_{margin}+G)\) and
does not form a dense conditional margin covariance. Exact-k rho scores use
the k-specific margin mean in these sufficient statistics; the initial
one-market restriction avoids an incorrect product of market mixtures under
the package's global centering constraint.

## Current likelihood and quality coordinates

The four active Stan programs all use the native specification-8 centered
Helmert coordinates. In `logit_supply.stan`, transformed data construct an
`N_group x (N_group - 1)` Helmert basis `Q` (lines 136--156), parameters contain
`b_group_centered` and `sigma_b_strat` (lines 159--194), and the transformed
quality is `b_quality_all = Q * b_group_centered` (lines 196--216). The exact
same contract appears in `ces_supply.stan` (lines 134--154 and 156--208),
`logit_stackelberg.stan` (lines 182--216 and 218--272), and
`ces_stackelberg.stan` (lines 190--223 and 225--274). The prior is
`normal_lpdf(b_group_centered | 0, sigma_b_strat)`; the positive prior on
`sigma_b_strat` is explicit in the model blocks. The supporting derivation and
identity checks are in [quality_coordinates.md](stackelberg/quality_coordinates.md),
especially lines 10--44.

The Stan object named `Q` is the *quality contrast basis*, not an observation-
to-firm incidence matrix. A collapsed implementation therefore needs an
observation design matrix `A` (one-hot `group[i]`) and the current basis `H`
(the Stan `Q`):

$$
b_{obs}=A H b_{group\_centered},\qquad
\operatorname{Cov}(b_{obs}\mid\sigma_b)=\sigma_b^2 AHH'A'.
$$

For `N_group = 1`, `H` has zero columns and quality is identically zero. The
current package's exact covariance is thus based on the contrast projection
`HH' = I - 11'/G`, with the separate market intercepts retained. It must not be
replaced by an iid `A A'` covariance.

R prepares shares by normalizing quantities within each market-year
(`R/helpers.R:323--329` for Logit and `:664--677` for CES), aggregates firm
shares for supply FOCs (`:331--341`, `:679--686`), and computes
`log_shareIn` and `clr_shareIn` in `compositional_share_data()` (`:58--75`).
The live likelihood uses `log_shareIn`; `clr_shareIn` is passed in the Stan
data but is not used by the current share target. Prices are standardized in
`R/helpers.R:343--346` (Logit) and `:688--691` (CES).

The current expanded Student-t share blocks are:

* `logit_supply.stan:372--387` and `:422--452` use independent product-level
  Student-t share terms conditional on the explicitly sampled quality vector.
* `ces_supply.stan:342--358` and `:393--424` do the analogous CES calculation.
* `logit_stackelberg.stan:423--431` and `:518--532` calculate the share target
  once before the exact-k margin mixture.
* `ces_stackelberg.stan:402--412` and `:509--525` do the same for CES.

The collapsed Gaussian block must therefore be a multivariate normal block
over all share observations (or an exactly equivalent low-rank evaluation). It
cannot be implemented as a product of independent `normal_lpdf` calls while
claiming to integrate out the shared quality effects. It also cannot use
`reduce_sum` over observations in its reference form unless the covariance
structure is first handled analytically.

## R estimation and data path

`R/estimate.R` is the public fitting entry point. The current arguments are
documented at lines 8--134 and defined at lines 235--268. It currently has
`likelihood = c("joint", "share", "margin")`, `use_rho`, `df`, and the exact-k
arguments, but no explicit `share_likelihood` or `quality_parameterization`.
The Logit and CES preparations are dispatched at lines 471--557; Stan sampling
is called at lines 566--570. Fit metadata is assembled at lines 572--605 and
must record the new likelihood and quality modes so that comparisons cannot
silently combine incompatible targets.

`R/helpers.R:196--511` (`prepare_stan_data()`) and
`:523--860` (`prepare_ces_stan_data()`) construct the full Stan data list,
factor encodings, market-year maps, outside-price normalization, observed
margin mask, and quality-group identifiers. Both paths append exact-k data
through `.stack_append_data()` when needed. A collapsed implementation will
need either an additional covariance/design block in these preparation
functions or a shared helper used by both, with strict validation for options
that cannot be collapsed.

The user-facing fit remains the S4 `BayesLogit` class (`R/AllClasses.R:48--58`)
inside virtual `BayesFit` (`:27--46`). There is no separate Gaussian fit class
needed merely for this experiment; mode and version belong in metadata and
Stan data. `R/methods.R:154--205` assumes a `BayesLogit` with a regular Stan
fit and currently prints `mu_log_a`, `sigma_share`, and `sigma_margin`.

## Rho and missing margins

The current rho path is a factored bivariate Student-t density. In
`logit_stackelberg.stan:57--70` and `ces_stackelberg.stan:62--75`,
`stack_margin_density()` conditions the margin residual on the matching share
residual. The same calculation is in the smooth programs at
`logit_supply.stan:389--406` and `ces_supply.stan:360--377`.

`rho_active` is set only when `use_rho == 1`, `likelihood_mode == 1`, and at
least one `margin_obs` is true (`logit_stackelberg.stan:188--215`, with the
same structure in the other three programs). `rho_raw` is consequently a
zero-length vector when inactive. The margin mask is prepared in
`R/helpers.R:347--358` and `:607--616`; masked rows remain in the share
denominator and share likelihood, while their margin values are replaced by a
finite dummy and excluded from the margin target.

For a Gaussian collapsed share block with `rho = 0`, the implementation can
retain the current independent margin block while replacing only the share
block. For `rho != 0`, the current per-row factorization is no longer a direct
collapsed equivalent: integrating shared quality creates cross-observation
share covariance, and the share-margin matching matrix plus missing-margin
selection must be represented in a joint Gaussian covariance or an equivalent
conditional Schur-complement calculation. A product of current
`stack_margin_density()` calls would omit those induced dependencies. A
Gaussian-share plus Student-t-margin option is likewise straightforward only
with `rho = 0`; exact nonzero-rho integration would require a different joint
calculation. This is a compatibility hazard, not a reason to alter the current
Student-t path.

## Exact-k organization and generated quantities

The exact ordered-k data contract is built in `R/stackelberg_bayes.R:245--342`
and `:344--422`. Firms are first aggregated to one row per firm-market-year;
the market ranking is computed from equal-period means at lines 260--294, and
the observation-to-firm maps are returned at lines 297--340. The partition
arrays and positional k priors are constructed at lines 344--421. This is the
existing multiproduct/double-counting protection and should be reused.

`logit_stackelberg.stan:433--485` and `ces_stackelberg.stan:414--476` evaluate
the complete market-level margin score for every k and apply one
`log_sum_exp(score)`. The share target is outside that loop at
`logit_stackelberg.stan:423--431` and `ces_stackelberg.stan:402--412`. Thus the
collapsed share density should be evaluated exactly once per posterior state,
then added to the market-level k mixture. The exact-k implementation already
has the intended `O(MK)` structure and should not be changed to k-specific
share evaluations.

Generated quantities in both exact programs expose
`posterior_component_prob`, `prior_component_prob`, `component_log_score`,
firm role probabilities, core-share summaries, and market/period/cluster log
likelihoods (`logit_stackelberg.stan:488--514` and
`ces_stackelberg.stan:479--505`). The k mixture is recomputed in GQ at
`logit_stackelberg.stan:534--663` and `ces_stackelberg.stan:527--664`, with
prefix accounting for observation-level `log_lik_margin`. A collapsed joint
share density is not decomposable into the current per-observation
`log_lik_share_marginal` values. GQ will need a clearly documented market-level
collapsed log-likelihood (and an explicit policy for observation-level LOO),
while preserving `posterior_component_prob` and the exact market mixture.

`R/posterior_k.R` extracts and summarizes the generated component probabilities
(see lines 269--384). `R/restriction_comparison.R:227--349` compares exact-k
restrictions using the encompassing-prior identity and requires the same
continuous-parameter prior across restrictions. These paths should continue to
work because k is orthogonal to the common share likelihood, but their tests
should verify that the collapsed share contribution is not duplicated by k.

## LOO, bridge, and simulation hazards

`R/bma.R:1--61` enforces identical observations, transformed outcomes,
margin masks, and shared priors across model comparisons. `bayes_compare()` at
`:64--130` requires generated `log_lik`, supports observation/period/market/
cluster units, and calls `bridgesampling::bridge_sampler()` at `:252--265`.
The documentation at `:120--130` explicitly requires normalized Stan targets
and rejects the non-normalized HMT bridge case. The collapsed Gaussian
implementation must retain multivariate normal determinant and constant terms;
otherwise bridge evidence and Gaussian-vs-Student comparisons are invalid.
Because the k state is market-level, current comparison logic already disallows
inappropriate observation/period LOO for market-shared exact-k mixtures. A
collapsed share likelihood makes observation-level LOO even less literal; the
default validation should prefer market or natural whole-panel units unless a
proper leave-one-out marginal density is implemented.

`R/stackelberg_bridge.R` reconstructs deterministic coordination draws. It
extracts posterior demand parameters and exact-k generated probabilities at
lines 606--730, then selects a role realization and assembles a coordination
object in `.stack_bridge_model_draw()` around lines 329--411. Legacy and exact
fringe routes are handled separately around lines 446--570. `R/simulate.R`
delegates fitted draws to these bridges (see lines 544--734) and stores
simulation diagnostics around lines 783--1118. The current bridge does not
extract a posterior firm-quality vector because quality is part of the Stan
draw; a future `posterior_quality()` helper must supply conditional Gaussian
quality draws for collapsed fits before any reconstruction that needs them.
`R/params.R` provides the Stan-free `bayes_params()`/`BayesSpecifiedFit` path
and should remain available for supplied point or draw matrices.

## Generated source and build files

`R/stanmodels.R:1--27` loads the four compiled modules and runs
`rstan::stanc_builder()` against `stan/<model>.stan` or `inst/stan/<model>.stan`
at package load. The generated C++ headers
`src/stanExports_logit_supply.h`, `src/stanExports_ces_supply.h`,
`src/stanExports_logit_stackelberg.h`, and
`src/stanExports_ces_stackelberg.h`, plus their small `.cc` wrappers, are
checked into the repository. `src/Makevars`, `src/Makevars.win`, `configure`,
and `configure.win` are rstantools-generated build plumbing. Any Stan source
change must be followed by the repository's normal rstantools regeneration and
compile checks; editing only `inst/stan` leaves installed/generated code out
of sync. Existing source tests inspect both Stan source and model dimensions,
so generated headers are a release artifact even though `stanmodels.R` rebuilds
model code during development.

## Recoverable Gaussian history

The original implementation is recoverable from reachable commit
`e0d33b0` (“Add antitrustBayes R package: Bayesian merger simulation with
Stan”). Its `inst/stan/logit_supply.stan` used direct product `shareIn` and
`normal_lpdf(shareIn | pred_s, sigma_share)` (`e0d33b0:inst/stan/logit_supply.stan:31--74`),
explicit iid `b_group_raw` with `b_quality_all = sigma_b_strat * b_group_raw`
(`:125--179`), and a normal margin conditional on the share residual. It did
not use the current centered Helmert normalization or log-share outcome.

The closest historical expanded Gaussian benchmark is the parent of commit
`78a2609` (“Harden logit_supply.stan geometry…”). Its Stan source used
`normal_lpdf(log_shareIn | pred_log_s, sigma_share)` and explicit
`b_group_raw`, with log-scale residual priors and current-style price/s0
structure. Commit `a851493` then moved the share block to centered/compositional
normal treatment while making margins robust Student-t; commit `962a130`
restored the current bivariate Student-t rho path. No separate old collapsed
Gaussian implementation was found in reachable history. `git fsck` found
unreachable coordination commits only; no additional antitrustBayes Gaussian
implementation was identified.

The historical source is useful as a benchmark fixture, but it is not an
automatic drop-in target: direct shares versus current log conditional shares,
uncentered iid qualities versus current Helmert contrasts, and old versus
current outside-share normalization are distinct statistical specifications.
The A-vs-B experiment should freeze whichever expanded Gaussian target is
chosen and compare it to its analytically integrated form before comparing
either with the current Student-t model.

## Existing tests and proposed additions

Relevant current tests are:

* `test-quality-coordinates.R` checks Helmert orthonormality, zero-sum quality,
  density/Jacobian identities, and that all four Stan graphs use the centered
  coordinate contract.
* `test-rho-studentt.R` checks the Student-t bivariate factorization, fixed and
  estimated df, rho activation/masking, generated log likelihoods, and Logit /
  CES parity.
* `test-exact-k-stan-source.R` checks one market-level exact-k mixture.
  `test-stackelberg-r-contract.R` checks unique firm-market-year aggregation,
  ranking, k support/prior validation, specified and estimated reconstruction,
  and exact-fringe simulation.
* `test-model-spec.R`, `test-bma.R`, and `test-s0-pooling.R` cover source
  contracts, comparison/bridge assumptions, and pooling dimensions.

Before any HMC benchmark, add independent small-matrix tests for:

1. expanded Gaussian density integrated over quality versus collapsed
   multivariate-normal density, including `G=1`, centered contrasts, repeated
   products owned by one firm, multiple market intercepts, and missing margins;
2. covariance and gradient checks for `sigma_b_strat`, `sigma_share`, s0,
   alpha/gamma, and rho where supported;
3. the structured \(O(N+G)\) sufficient-statistic evaluation versus a dense
   Cholesky reference, including rho pseudo-observations, log determinant, and
   all constants; benchmark the structured route at \(N=10^2,10^3,10^4\) for
   single-product and repeated-group panel layouts;
4. exact-k invariance of the share likelihood and exactly one whole-market
   k-mixture contribution for BF, CF, BS, and CS;
5. generated quantities and bridge-density normalization, with market-level
   LOO semantics made explicit;
6. conditional Gaussian quality posterior means/covariances versus direct
   conditioning formulas and expanded-model draws on tiny fixtures;
7. A/B/C API validation and metadata (`share_likelihood`,
   `quality_parameterization`, rho/missing-margin combinations), plus unchanged
   Student-t behavior;
8. the nested information ladder: one-market/one-period base, then additional
   periods/markets with the same economic equation and only the intended
   covariance/pooling additions.

The controlled benchmark should start with simultaneous Bertrand and the
prespecified no-rho/no-extra-hierarchy settings, then run a small Cournot,
fringe, and Stackelberg subset. Retain failed runs and report divergences,
treedepth, BFMI, R-hat, ESS/sec, posterior variance components, predictive
fit, bridge diagnostics, and policy simulations; no conclusion follows from
source inspection or one sampler run.

## Main compatibility hazards

1. Collapsing changes a product-factorized share target into a correlated
   multivariate target; `log_lik` and LOO units must not pretend otherwise.
2. The current `Q` name is a contrast basis, so an explicit incidence matrix is
   required to construct the covariance and avoid multiproduct errors.
3. Nonzero rho with collapsed quality is a joint Gaussian covariance problem,
   not a mechanical reuse of the current rowwise Student-t conditional.
4. Exact-k share likelihood must be common across k; only the market-level
   margin mixture varies with k.
5. Bridge sampling requires every Gaussian normalizing constant and a target
   matching the stored metadata; model comparison must reject incompatible
   share/quality modes.
6. Collapsed fits do not contain sampled quality coordinates, so simulation,
   posterior-quality summaries, and any fitted-demand reconstruction need a
   conditional recovery facility.
7. The old Gaussian commits use materially different outcomes and priors. They
   are recoverable reference code, not evidence that current A/B equivalence
   has already been established.
