# Prespecified exact ordered-k recovery and PNB diagnostic protocol

Status: pre-fit protocol, 2026-09-15. This protocol supersedes the cutoff-based
recovery design previously stored at this path. No result from an ordered-k
posterior fit was used to choose the design, priors, support, or diagnostics.
The legacy smooth-cutoff fringe remains a separate robustness model and is not
the data-generating process in this experiment.

## Estimand and finite-state prior

Each synthetic market has five firms and one fixed strategic-core count across
all eight periods. Firms are ranked once, before estimation, by the observable
rule used by the package: equal-period mean firm conditional quantity share for
Logit and equal-period mean firm conditional revenue share for CES, with an
ascending firm-label tie break. Every product and period of a firm inherits
that rank. The true core size is k=2.

The candidate supports and priors are fixed before fitting:

* BF and CF use k_support=0:5 and k_prior=rep(1/6,6). State zero is exact
  monopolistic competition.
* BS and CS use k_support=1:5 and k_prior=rep(1/5,5). Stackelberg has no empty
  leader state.
* BF and CF therefore share one prior, and BS and CS share the corresponding
  positive-core prior. These are application choices, not package defaults.

One k applies to the complete market likelihood. Stan evaluates the sum of all
period/product margin log densities conditional on k before `log_sum_exp`.
There is no period-level or observation-level mixture and no sampled discrete
parameter.

## Synthetic configurations

Cross Logit/CES demand with Bertrand/Cournot base conduct and seeds 1103, 2207,
and 3301. Use the single-product layout for seed 1103, a multiproduct leader for
2207, and a multiproduct follower for 3301. This gives twelve primary DGPs and
covers each required ownership layout without treating several layouts from
one random seed as independent replications.

Truth is |alpha|=3 for Logit, gamma=3.5 for CES, and outside price 1.3. The
Bayes-compatible outside quantity share is derived algebraically from the
single-market normalization used by the fitted model. Firm demand quality is
fixed across periods. Conditional shares vary over time, including changes in
the second/third firms' period-specific ordering, while the predetermined
whole-panel ordering and top-two core remain fixed. Product costs are obtained
by inverting the independently audited hard-role FOCs at the generated prices.

Margin measurements use Student-t noise with five degrees of freedom and
standard deviation .04 in log inverse margin. The primary experiment retains
exact demand shares and labels that choice explicitly. Because exact shares can
drive the share-residual scale toward zero, the sampler gate records this
geometry rather than hiding it. A prespecified normalized compositional-noise
stress experiment may follow, but it is reported separately and cannot replace
the primary result.

The generator must pass structural checks before fitting: normalized demand,
constant within-firm demand quality, firm-level roles, no product double
counting, positive costs, valid margins, fixed market ranking, and correct
Logit/CES BF/CF/BS/CS FOC predictions. `synthetic_dgp.R` is the independent
oracle and may not call the production Stan likelihood or deterministic solver.

## Candidate fits and priors

For each DGP, fit three candidates with the same demand family and base
conduct: exact Stackelberg, exact fringe, and simultaneous conduct. The
matching structure is the recovery candidate; the other two test model
discrimination. A separate exact-fringe DGP experiment reverses Stackelberg
and fringe roles so BF/CF recovery is not judged only as a misspecified
alternative to BS/CS.

Use the package baseline proper prior manifest for every candidate, with no
post-result changes. Use `use_rho=FALSE`, fixed df=5,
`use_market_resid_scale=FALSE`, `use_alpha_market_pool=FALSE`,
`use_s0_market_pool=FALSE`, and `use_s0_market_year_pool=FALSE`. The latter
settings fit one demand coefficient and one outside share for the synthetic
single-market panel. CES uses the same observed inside quantities as Logit;
the package converts them to revenue shares using price and outside price.

The primary HMC run uses two chains, 1000 iterations per chain, 500 warmup,
adapt_delta=.95, and max_treedepth=13. A short fixed-parameter or small-HMC run
is only a compilation/sampler smoke test and is never reported as recovery.
If the primary run fails diagnostics, any longer or differently parameterized
rerun is labeled as sensitivity evidence and does not replace the failure.

## Required recovery output

Retain one compact row per DGP/candidate with the seed and layout. Report:

* demand-parameter median, 95 percent interval, error, and truth coverage;
* outside-share median, interval, error, and coverage;
* P(k=r|D), posterior mean/median k, true-k probability, modal k, and posterior
  strategic-core share for exact structures;
* predicted-margin error by role and overall;
* divergences, maximum R-hat, minimum effective sample size, BFMI, and
  maximum-treedepth hits;
* elapsed time and any sampling or bridge-sampling failure;
* bridge marginal likelihood when its numerical diagnostics are usable, plus
  joint predictive fit without substituting one for the other.

Recovery is assessed across all three seeds. One visually favorable
replication, a diffuse posterior that happens to cover truth, or a sampler run
with unresolved diagnostics is not called successful. All failed outcomes are
kept in the table.

## Counterfactual recovery

After parameter and k recovery pass, evaluate a no-cost-change merger, a 10
percent marginal-cost shock, and their combination with the same deterministic
truth and posterior scenario definitions. Stackelberg uses the coordination
solver. Exact fringe uses a hard top-k strategic-core/MonCom-fringe solver;
smooth cutoff weights are not a substitute.

The main posterior simulation maps each continuous draw to one market-wide k
drawn from P(k|D,theta), cached across periods and pre/post scenarios. This is
an unbiased Monte Carlo representation of the finite mixture. Also report
fixed modal k, every supported fixed-k sensitivity, and, once implemented,
the exactly probability-weighted conditional outcome mean across k. The last
quantity is a Rao-Blackwell summary and must retain component failure rates;
it does not define a fractional-role game.

Compare true and posterior changes in product prices, shares, and welfare, plus
solver success and leader/follower or strategic/fringe FOC residuals. Do not
impose unproved price rankings as test expectations.

## PNB smoke diagnostic

PNB fitting begins only after package suites, compiled likelihood checks,
synthetic recovery, and exact-fringe simulation pass. Estimate the eight
prespecified models: simultaneous Bertrand, BF, BS, simultaneous Cournot, CF,
CS, auction, and MonCom. Use k_support=0:5 with a uniform prior for BF/CF and
k_support=1:5 with a uniform prior for BS/CS. These are explicit application
arguments. No Court-informed prior enters the smoke run.

Use the same proper prior manifest and sampling settings across candidates.
Report posterior k, core share, predicted margins, bridge and sampling
diagnostics, merger effects, fixed-k sensitivity, and simulation failure rates.
Do not change priors after inspecting results and do not modify the paper or
mergerBayes in this phase.

## Artifact and evidence policy

Compiled objects, stanfit objects, chain CSV files, and verbose logs remain in
temporary local storage and are deleted after compact summaries are extracted.
Only small source files and CSV summaries belong in Git. Each final result is
labeled analytically established, numerically verified, empirically estimated,
or unresolved.
