# Prespecified synthetic recovery and PNB diagnostic protocol

This protocol is written before new-model posterior fitting. Results will be
reported as recovery evidence or failures, not used to retune package priors.

Synthetic configurations: Logit/CES crossed with Bertrand/Cournot; five firms,
eight periods; include a six-product version with one multiproduct firm in
each role. Seeds 1103, 2207 and 3301. Demand truth |alpha|=3 or gamma=3.5;
fixed firm demand quality across periods. The general economic DGP has outside
quantity share .30. Before fitting, a Bayes-compatible mode derives outside
truth from the existing single-market share normalization (see the addendum
below), keeping observed prices/shares and all other truths unchanged. Construct
varying conditional quantity shares that move the second/third firm's size
across a market cutoff c=.20, infer price paths from the fixed normalized
demand system, and invert the independently audited hard-role FOCs for costs.
This gives a coherent constant-demand, varying-cost DGP without requiring
a production equilibrium solver to generate the baseline oracle.

Timing is generated with the audited common logistic boundary at prespecified
tau=.005. Roles use conditional firm quantity size for Logit and conditional
firm revenue size for CES. Draw a publicly realized hard partition per period;
keep truth identities, effective k and costs. Margin measurements have
Student-t noise with five degrees of freedom and standard deviation .04 in
log inverse margin. Invalid margins cause the DGP to fail; they are never
clipped or selectively dropped.
Shares can first be exact, with a separately labeled measurement-noise stress
scenario. Do not incorrectly call independent noisy log shares a normalized
multinomial DGP.

The exact-share experiment sets use_rho=FALSE in every candidate fit, retaining
the independent marginal Student-t margin equation that generated its noise.
Zero correlation in a bivariate Student-t is not independence, so the default
conditional df+1 equation would otherwise change this DGP's error law. PNB
retains its existing default joint correlation setting.

Synthetic application priors are fixed across all replications and candidate
structures: cutoff beta(2,8), other prior manifest entries baseline, with any
noise-scale overrides listed before the first fit. This prior is chosen for
the simulated cutoff domain, not supplied as a package default. Core-count
truth can vary by period; do not attempt to point-identify a hard cutoff inside
one fixed share gap using a single cross section. Panel crossings provide
information that a single realized role partition does not.

The synthetic single-market panel uses one true price coefficient and one true
outside quantity share across periods. It fits one global price coefficient
and the existing independent-period outside-share hierarchy:
use_alpha_market_pool=FALSE,
use_s0_market_pool=FALSE, use_s0_market_year_pool=FALSE, and
use_market_resid_scale=FALSE in every candidate. PNB keeps the package's
existing default pooling settings. The coordination truth constructor takes
observed inside quantity for Logit and observed inside expenditure for CES.
The public Bayesian simulation API retains its existing inside-quantity
insideSize convention for both demands and converts CES quantity to expenditure
using observed conditional quantity shares and prices. Supply observed inside
quantity to Bayes in both families; welfare comparisons then use identical
physical quantity/expenditure units after the conversion.

Fit each synthetic dataset under matching Stackelberg, simultaneous and fringe
structures, using the same demand family and base conduct. Initial sampling:
two chains, 1000 iterations, 500 warmup, adapt_delta=.95, max_treedepth=13;
record all diagnostics. If diagnostics fail, any longer rerun is explicitly
labeled and does not replace the failed run in the record. Model discrimination
uses bridge marginal likelihood when its diagnostics are usable and otherwise
reports the failure alongside joint predictive fit. Do not infer successful
discrimination from the visually closest margin curve.

For every seed/configuration report demand-parameter posterior median/interval,
truth coverage, cutoff median/interval, hard and realized k PMFs, role
probabilities, core-share posterior, margin fit, model odds, sampling diagnostics
and simulation success/error rates. For both a no-cost-change merger of two
specified firms and a +10% cost shock, compare true intervention effects with
posterior intervals for prices, shares and welfare. Include combined shocks.
Use identical post-role rules in truth and posterior simulations and record
them. Compare mixture, posterior-cutoff threshold and adjacent-role mappings.

Recovery is not declared successful from a single seed. Report coverage and
errors across seeds with the individual rows retained. Wide/confounded cutoff
or demand posteriors, low probability for true k, wrong model odds or failed
counterfactual coverage are unresolved research outcomes, not test tolerances
to relax. This small experiment is a diagnostic, not a Monte Carlo precision
study of frequentist coverage rates.

PNB smoke fitting starts only after package suites and derivative/aggregation
tests pass. Eight models: Logit Bertrand simultaneous/fringe/Stackelberg,
Logit Cournot simultaneous/fringe/Stackelberg, auction and MonCom. Use the
existing input-market data, baseline proper prior manifest and identical
sampling settings across candidates. Keep the legacy fringe scale .2; new
Stackelberg tau=.01, plus prespecified .0025/.005/.02 sensitivity runs. No
Court-informed prior is added to this smoke diagnostic. Report cutoff, hard
and realized k, core share, margins, bridge diagnostics, sampling diagnostics,
merger effects and near-cutoff mapping sensitivity. Record computational or
identification failures explicitly; do not rewrite a paper or mergerBayes.

## Pre-fit single-market normalization repair

This repair is algebraic and precedes every new fit. Existing Stan code fixes
b_market=0 in a single market and centers b_quality over unique firm IDs.
The general structural DGP with arbitrary outside .30 therefore violates the
statistical share normalization. It remains a labeled misspecification stress
case. Production legacy equations and all priors remain unchanged.

For fixed base conditional quantity shares qbar and prices pbar, Logit uses
deltaNorm=log(qbar)+alpha*(pbar-p0). The compatible outside truth is

    logit(s0)=mean_unique_firms(deltaNorm)+alpha*(p0-mean_all_prices).

CES uses logQNorm=log(qbar)+gamma*log(pbar/p0) and

    logit(s0)=mean_unique_firms(logQNorm)-gamma*mean_all_log_prices.

Both means over prices include all observed product-period rows, matching
existing data standardization; quality means count each firm once. Observed
price paths do not depend on the chosen constant outside share. Rebuild the
physical normalized qualities, shares, costs and margins with the derived s0.
The resulting true firm-quality coordinates have mean zero and exactly match
the existing role-free Stan share predictions. Run synthetic_dgp.R with its
second argument 'bayes' to reproduce the 72 checks for cutoff and fixed-k modes.
Derived Logit s0 is .15345070/.13756237/.13123608 for single/multiproduct-leader/
multiproduct-follower layouts; CES is .06747971/.05721621/.05553923.
No posterior result motivated this adjustment. Single-market CES share
normalization with nonunit outside prices has a remaining dependence on price
units through this restriction; it is reported rather than silently repaired
in legacy models.

Direct-k implementation and its additional bounded recovery experiment are
prespecified in ordered_k_addendum.md. They marginalize exact games under a
direct count prior; cutoff-induced exact mixtures remain the baseline recovery
experiment. No cutoff recovery claim applies to direct-k fits.

The later [nonempty-core addendum](nonempty_core_addendum.md) governs positive-k support, conditional cutoff probabilities, and cutoff-center summaries. It supersedes any earlier zero-leader Bayesian-support examples.
