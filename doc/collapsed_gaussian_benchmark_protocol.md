# Prespecified collapsed-Gaussian benchmark protocol

Status: frozen before any collapsed-Gaussian HMC fit. Source and unit tests may
change to correct implementation defects; data, priors, sampler settings, and
acceptance gates do not change in response to posterior results. Failed and
interrupted runs remain in the ledger.

Before these HMC cells run, the production collapsed density must use the
exact sparse group-sufficient-statistic algorithm in
`collapsed_gaussian_math.md`; a full covariance Cholesky is permitted only in
the reference test. Numerical equality is required for plain and rho cases.
The optimized route is separately timed at \(N=100,1000,10000\), for both
one-product-per-firm and repeated-group multiproduct/panel layouts. These
algorithmic timings do not replace the HMC ESS-per-second benchmark below.

## Models

The three primary modes use identical economic equations and structural
priors:

* A: Gaussian shares, sampled centered qualities, Gaussian margins;
* B: Gaussian shares, collapsed centered qualities, Gaussian margins;
* C: Student-t shares, sampled centered qualities, fixed-df Student-t margins.

All primary fits set `use_rho=FALSE`, `use_market_resid_scale=FALSE`,
`use_alpha_market_pool=FALSE`, `use_s0_market_pool=FALSE`,
`use_s0_market_year_pool=FALSE`, `use_hmt=FALSE`, and `df=5`. The package
baseline prior manifest is used without modification. A and B must be the same
normalized Gaussian probability model; they differ only in whether the quality
contrasts are sampled or integrated.

## Data and run plan

The sparse base is one market, one period, 24 single-product firms. Demand
qualities are centered in the package's Helmert coordinates. Prices, true
conditional shares, outside share, and conduct-consistent margins are generated
from a fixed independent oracle. The exact-share scenario uses those shares
without perturbation. The noisy-share scenario applies prespecified Gaussian
log perturbations with standard deviation 0.04 and renormalizes the complete
inside composition; it is never constructed by adding unnormalized share
errors. Margin log inverse errors are Gaussian with prespecified standard
deviation 0.04. A, B, and C fit exactly the same realized data; C changes only
the fitted residual family to the fixed-df Student-t robustness model.

The primary seed is 271828. A second noisy-share sparse replication uses seed
314159. The primary cells are:

1. Logit simultaneous Bertrand, exact and noisy shares, A/B/C, seed 271828;
2. Logit simultaneous Bertrand, noisy shares, A/B/C, seed 314159;
3. Logit simultaneous Cournot, noisy shares, A/B/C, seed 271828;
4. Logit exact BF with `k_support=0:5` and exact BS with `k_support=1:5`,
   uniform priors, noisy shares, A/B/C, seed 271828;
5. CES simultaneous Bertrand, noisy shares, A/B/C, seed 271828.

This is 21 primary fits, far below the earlier 144-fit recovery ladder. A
four-period version of the same oracle is a nested-information-ladder check for
mode B only; it is not substituted for a failed sparse fit.  A compact
multi-market, multi-period mode-B fixture additionally verifies that the same
likelihood and API accept higher ladder rungs, preserve repeated group rows in
the quality design, and use the requested alpha/outside-share pooling flags.
That structural check is not counted as another A/B/C benchmark cell.

Each primary fit uses four chains, 1,000 iterations per chain, 500 warmup,
`adapt_delta=.95`, `max_treedepth=13`, and the same fixed initialization policy.
The package default diagonal metric is retained. A short compile/sampling smoke
may precede these fits but is never reported as benchmark evidence.

## Fixed diagnostics and comparisons

A usable fit has zero divergences, zero maximum-treedepth hits, minimum BFMI at
least 0.30, maximum rank-normalized R-hat at most 1.01, and minimum bulk and
tail ESS at least 100 for the monitored structural and variance quantities.
The monitored set also includes all sampled quality contrasts in A/C and the
sampled log-scale/outside-share coordinates. Thus an apparently well-mixed
structural coefficient does not conceal poorly mixed explicit quality effects.
Failures remain failures; a later alternative sampler configuration, if any,
is a separately labeled sensitivity.

For every fit record elapsed time, gradient evaluations where available,
diagnostics, ESS and ESS/second for the demand coefficient, outside share,
quality scale, share scale, and policy summaries, posterior variance-component
correlations, predictive share/margin summaries, and any bridge failure. Bridge
sampling is attempted for the seed-271828 noisy simultaneous-Bertrand A/B/C
fits only, with identical repetitions and limits.
The fixed bridge settings are three normal-method repetitions, at most 1,000
iterations per repetition, one core, and effective-sample-size weighting.
Warnings, nonfinite estimates, or iteration-limit hits are retained as bridge
failures; they are not accepted because an estimate happens to be finite.

A and B posterior summaries are compared using Monte Carlo uncertainty; source
or target equality alone is not treated as empirical equivalence. B and C are
compared as different likelihood families. No numerical threshold is used to
declare policy equivalence mechanically; the CSV reports absolute and relative
differences in posterior medians and interval endpoints.

## Counterfactuals

Every usable fit is reused, without refitting, for 100 posterior draws under:

1. merger of the two largest observed firms with no cost change;
2. a 10 percent marginal-cost increase for every inside product;
3. the merger and cost increase together.

The output records product price and share changes, available consumer,
producer, and total-surplus summaries, aggregate price summaries, and solver
failure rates. Exact-k simulations use the fit's posterior component
probabilities and retain the same market-wide role draw within each pre/post
scenario.

## Prespecified sensitivities

After the primary ledger is complete, mode B for the seed-271828 noisy sparse
Bertrand case is rerun with the quality-scale prior width multiplied by 0.5 and
2. These are labeled variance-decomposition sensitivities, not replacements.
One mode-B run with `use_rho=TRUE` and Gaussian margins checks computational
viability of the optional correlation. Neither sensitivity changes the package
default or rescues a failed primary fit.

Large fit objects, chain files, and verbose logs stay under `/tmp`. Git receives
only source, compact CSV summaries, and the final evidence note.

## Implementation failures and continuation

The initial runner misplaced `metric` outside rstan's `control` list; all
21 setup failures are retained. The next attempt saved A and B posteriors,
then exposed a scalar-array bug in quality recovery. That parser was repaired
and tested against the actual saved fit. `--resume` checks exact equality of
the frozen plan and reuses saved posteriors; it does not rerun their HMC.
Per-fit timers are now persisted before postprocessing. The missing enclosing
B timer remains unavailable. Uniform native warmup-plus-sampling timings are
reported separately for all fits and used for consistent ESS comparisons;
the first A wall timer includes a cold package load and is not used to infer
a speed ratio. These are implementation/measurement repairs, not changes to
data, priors, seeds, sampler settings, or diagnostic gates.
