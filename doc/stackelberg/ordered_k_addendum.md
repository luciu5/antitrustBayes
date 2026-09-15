# Exact ordered-k prior: specification after user steering

This addendum precedes Bayesian implementation and all new-model fits. The
hard-game economics and deterministic interface are unchanged. Both Bayesian
prior choices below marginalize exact hard games; neither blends firm FOCs.

## One shared engine, two prior choices

Keep structure='stackelberg' and append partition_prior='cutoff', k_prior=NULL,
k_scope='market' to the public Bayesian APIs. The cutoff choice retains the
already approved latent-boundary architecture, but evaluates the exact joint
mixture over k described in the audit. The new partition_prior='ordered_k'
choice uses prespecified probabilities directly and samples no cutoff.
Non-Stackelberg models reject ordered_k initially; legacy fringe is unchanged.

For ordered_k, k_prior is a named numeric vector whose names are integer counts
(e.g. c(`1`=.2,`2`=.2,`3`=.2,`4`=.2,`5`=.2)), or a named list of such vectors
by market. Weights must be finite, nonnegative and sum to one within numerical
tolerance. Omitted counts have zero mass. Reject unsupported positive-mass
counts rather than silently truncating or renormalizing a supplied prior.
NULL means uniform over all supported counts, including all leaders; zero is excluded by the nonempty-core revision.
A Court-informed 1:5 prior is an explicit application choice, never a default.

The default k_scope='market' means one k_m for all observed years of market m.
Support is 1:min_t N_mt so every state is valid in every period. Explicit
k_scope='period' permits independent k_mt and supports 1:N_mt by period. In
both cases the number is fixed within its partition unit; leaders are that
period's top k firms by the existing demand-specific conditional size metric.
Public firm-ID ordering resolves exact size ties deterministically. Unlike a
scalar cutoff, a rank-count prior can split a size tie; label a resulting
implied scalar-boundary interval as nonexistent, not as an estimated cutoff.

## Marginal likelihood and generated quantities

Let u denote a partition unit (market or period), and J_u all its observations.
For every supported k use the audited hard-partition margins, that period's
ownership/ranking, and the same demand/outside-share draw. The margin target is

    logsumexp_k [log pi_uk + sum_{j in J_u} ell_jk].

Share densities add using the unchanged role-free share equation. Existing
rho-conditioned Student-t margin densities and margin masks are reused.
The cutoff prior choice has period units and pi_uk(c,tau) from the common
logistic boundary; the ordered_k choice uses fixed supplied pi_uk. Shared
functions compute hard-game moments and component scores in model and GQ.

Use generic flattened partition-unit component arrays: starts/counts, k values,
market-year-to-unit mapping, fixed prior log probabilities and posterior GQ
probabilities. Firm-market-year aggregation stays unique and unchanged.
For ordered_k market units, sum every period's component density before the
logsumexp. Multiplying separately marginalized period densities is a different
model and is forbidden under market scope.

GQ posterior probabilities are softmax(component scores). Their draw average
is P(k_u=r|D); no discrete parameter was sampled. Predictions E[m] and E[1/m]
remain separately labeled. Firm leader probabilities use the period's firm
rank with the shared unit posterior component probabilities. Core shares vary
by period even when k_m is fixed. Cutoff summaries for ordered_k return no
cutoff posterior; optional implied share-gap intervals are geometric summaries
only. Cutoff-mode summaries retain both hard k(c) and realized-partition k.

Sequential row/period likelihood factors may account for the joint target, but
are not valid leave-one-row/period-out densities after integrating shared k_m.
Under ordered_k market scope reject observation and period LOO units; require
market or cluster aggregation, with cluster constant throughout each market.
Under period scope retain the audit's observation-LOO rejection. Full-target
bridge sampling remains valid for either proper prior.

## Reconstruction and counterfactuals

Default role_mapping='mixture' samples a supported k using posterior component
probabilities (prior probabilities for specified fits without likelihood data).
Cache draw/unit selections so a joint panel simulation uses the same sampled
k_m in all selected periods. Role_mapping='map' selects the highest-probability
component and is a deterministic sensitivity. Fixed-k and adjacent-k overrides
are explicit sensitivity options. role_mapping='threshold' requires an actual
cutoff and is available only for cutoff-prior fits; do not invent a cutoff
posterior for ordered_k fits.

Selected roles instantiate coordination's exact hard game, with draw-specific
demand, FOC-inverted baseline costs, and the documented Bayesian merged-owner
constituent-role union or explicit leadersPost override. No deterministic
solver is added to Stan or antitrustBayes. Bayesian CES market size remains
inside quantity and is converted to coordination expenditure exactly once.

## Verification and research runs

Add independent two-period tests where per-period preferred counts differ:
market-scope GQ/target must match a single whole-market R logsumexp and differ
from the product of period mixtures. Test fixed/sparse/market-specific priors,
zero/all leaders, ties, absent margins, interleaved rows and shared simulation
state. Legacy targets and the current cutoff-induced recovery protocol stay
unchanged apart from the separately documented pre-fit normalization repair.

Add a small direct-k recovery experiment with constant known k=2 across eight
periods, Logit/CES x Bertrand/Cournot, three seeds and a multiproduct follower
layout. Matching direct-k, simultaneous and fringe candidates use the same
already prespecified demand/noise priors and HMC settings. Do not claim cutoff
recovery for a model that has no sampled cutoff. PNB may compare the two direct-k
Stackelberg prior variants with an alternative eight-model comparison replacing
the cutoff-prior Stackelberg candidates, rather than double-counting economic
models in one model-prior universe. Default direct-k prior is unrestricted
uniform, not 1:5; any Court-informed run requires a separately declared
application prior before fitting.

Analytically established: the marginalized likelihood and posterior component
probability identities. Numerically/empirically verified: pending new gates.

The later [nonempty-core addendum](nonempty_core_addendum.md) governs positive-k support, conditional cutoff probabilities, and cutoff-center summaries. It supersedes any earlier zero-leader Bayesian-support examples.
