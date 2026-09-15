# Exact Bayesian ordered-boundary contract

Approved mathematical specification supplementary to the audit and map.
This is implementation guidance, not production code. Implement only after
the deterministic coordination API passes its gates.

## Structure resolution

Append structure=NULL and cutoff_scale=NULL to public APIs, preserving old
positional arguments. Omitting structure maps old use_cutoff to fringe or
simultaneous. An explicitly supplied FALSE cutoff flag conflicts with explicit
fringe/Stackelberg; an explicitly supplied TRUE conflicts with simultaneous.
TRUE plus explicit Stackelberg is consistent. Internally cutoff is active for
fringe and Stackelberg, but the quality shift exists only for fringe.
Integer structure code: simultaneous=0, fringe=1, Stackelberg=2.
Keep a readable string in provenance. Stackelberg supply codes are 1 and 3.
Resolve NULL scale to .2 for legacy fringe and .01 for Stackelberg; validate a
positive finite prespecified scalar. The parameter is never sampled.

The user waived backward compatibility during implementation. The final
contract is specification 7: normalize priors natively in every Stan program
and require previous fits to be refitted. Economic meanings of simultaneous,
fringe, and coordinated PriceLeadership remain distinct. There is no external
bridge-normalization mode. See evidence_normalization_contract.md and
nonempty_core_addendum.md for the latest governing amendments.

For Stackelberg, .bayes_model_spec returns a separate clearly named Bayesian
structure specification recording demand, base conduct, structure and
coordination solver route; it must not advertise the registered ordinary
antitrust model ID as the complete Stackelberg identity. Legacy specifications
and as_antitrust_fit routes retain their existing registry identities.

## Unique firm cells and sorting

Retain original product-row order. Make unique (market_year,group) rows,
aggregate shareIn and revShare once each, and sort within market-year by
decreasing cutoff metric, with a stable arbitrary secondary ID order only for
storage. Equal-size firms have zero mass for partitions splitting the tie.
Store N_firm_mky, firm_mky_shareIn, firm_mky_revShare,
firm_mky_market_year, firm_mky_group, obs_to_firm_mky,
firm_mky_rank, firm_mky_start, firm_mky_n. Starts/counts refer to sorted unique
firm rows. Existing obs_order/mky_start/mky_n refer to sorted product rows;
always gather original observations through obs_order[start+j-1]. Test fully
interleaved market-years, product permutations and multiproduct tie blocks.

Logit cutoff size is conditional inside quantity share; CES size is conditional
inside revenue share. Logit supply S_f=shareIn_f*(1-s0). CES supply
R_f=revShare_f*(1-R0), with R0 converted from outside quantity s0 using
the package's observed-price/outside-price expenditure index. Never sum repeated
product-row firm totals for a follower statistic.

For every CES structure use the exact physical conversion
R0=inv_logit(logit(s0_prob)+log(outside_price)-log(inside_price_index)),
without applying bounded_inv_logit a second time. The quantity s0_prob retains
its existing bounded hierarchy. This corrects the historical second epsilon transform and removes a physical
outside-revenue mismatch. Reject a nonrepresentable open-interval revenue
share instead of clipping the equations. Logit Cournot's existing kappa coordinate still
derives one global alpha when pooling is disabled; preserve that hierarchy.

## Stable partition probabilities

For descending sizes a_1,...,a_n and cutoff c, let x_i=(a_i-c)/tau.
Use log pi_0=log1m_inv_logit(x_1), log pi_n=log_inv_logit(x_n).
For an intermediate partition k, the stable identity is

    log pi_k = log_inv_logit(x_k) + log1m_inv_logit(x_(k+1))
               + log1m_exp((a_(k+1)-a_k)/tau).

The gap is computed from sizes directly, avoiding cancellation of c in two
large logits. For exact ties return negative_infinity() explicitly. This
follows from sigmoid(x)-sigmoid(y)=exp(x)*(1-exp(y-x))/[(1+exp(x))(1+exp(y))].
Do not subtract saturated CDFs, drop tiny nonzero components or clip the common
boundary to [0,1]. All-leader/all-follower tails are modeling implications.
Use log_sum_exp normalization for component posterior probabilities.

## Component FOCs

For each k=0,...,n, firms of rank <=k lead and others follow. Compute audit
LB-F/L or LC-F/L for Logit; CB-F/L or CC-F/L for CES. Every equation uses the
component's hard partition and unconditional shares. In CES Cournot a stable
equivalent to the audit coefficient avoids a small b cancellation:

    h=gamma-1; H=1+h*R0; b=h/(gamma*H)
    a_f=h*(H-2*R_f+gamma*h*R0*R_f/H)
        /[H*(H+(gamma-2)*R_f)]
    m_F=1/gamma+b*R_f
    m_L=1/gamma+b*R_l/(1+sum_F R_f*a_f).

This follows by A_f=b*(H-R_f), not a new approximation. The independent
proof in the audit bounds 1+C_F below by R0+sum_L R_l. Verify equivalence
numerically across gamma and outside shares before using this stable form.
Avoid descending Stan loop ranges when k=n: guard follower loops by k<n.

For Logit predicted log inverse normalized margin is log(p_j*alpha/h_f+1e-9),
exactly the existing numeric offset. For CES it is log(1/m_f+1e-9).
All coefficients must be computed through shared functions called by both
the model and generated quantities.

## Joint likelihood, predictions and generated quantities

Stackelberg demand share predictions have no role quality shift. Reuse the
existing Student-t share density and per-row margin density conditional on
share/rho (when active). Let ell_jk denote that conditional margin density
for observed margin row j in partition k; zero for masked rows/share-only.

For a whole market-year, the component score is
t_k=log pi_k+sum_rows ell_jk and the margin target increment is
log_sum_exp(t_k). Share likelihoods add independently using existing equations.
Only this joint component construction has one public timing partition shared
across firms/products. Reusing independent per-row mixtures is disallowed.

GQ exposes posterior component probabilities softmax(t), their starts/counts
in flattened arrays (sum_period(n_firm+1)), prior role weights, and conditional
posterior firm leader probabilities. Predicted inverse margins can be posterior
mixture means. Also expose predicted margins as means of m, since
1/E[1/m] is generally not E[m]. Distinguish fitted conditional component
prediction from an unconditional prediction in an unobserved market.

For legacy row log_lik compatibility, a valid sequential conditional margin
factor is L_j-L_(j-1), where L_j=log_sum_exp_k[log pi_k+
sum_r_through_j ell_rk], using one documented deterministic observation order.
This telescopes exactly to the market-year margin target. Masked rows have
increment zero. Row factors depend on ordering individually; authoritative
market-year/market/cluster scores do not. Require a constant loo_cluster within
each Stackelberg market-year so cluster densities do not split a shared latent
partition. Add share densities to row factors and preserve existing aggregate
GQ fields, with model-target/GQ equality and permutation tests.
Sequential row factors are an accounting decomposition, not valid
leave-one-product-out densities after marginalizing a common partition.
Reject loo_unit="observation" for comparisons involving Stackelberg and
direct callers to period/market/cluster units. Bridge sampling remains valid
because it evaluates the full marginalized log posterior.

When no margin is observed, all component densities coincide: c has only its
prior information. Zero/all leaders both reproduce ordinary simultaneous FOCs
in every demand/conduct family; that endpoint observational equivalence is
another source of timing non-identification, not a computational defect.

## Posterior summaries and simulation

Public bayes_cutoff_summary returns per market-year:
cutoff draws/intervals, hard k(c) PMF/mean/median/interval, conditional and
unconditional strategic-core share, firm identities and role probabilities.
Logit core shares are quantity shares; CES core shares are revenue shares.
Conditional shares refer to the inside sector; unconditional CES core shares
use 1-R0, while Logit uses 1-s0. The hard
role probabilities P(a_f>c|data) differ from logistic prior role weights and
posterior realized-partition probabilities; label all three. Realized-k
summaries use GQ component probabilities, and hard-k summaries use continuous
cutoff draws directly. Support specified-parameter fixtures without Stan for
unit tests. Proper priors and bridge sampling evaluate the same marginalized
target; no discrete parameter is introduced.

Estimated draw reconstruction obtains component posterior weights from GQ;
sample one partition per draw/unit using role_mapping="mixture". Threshold
mapping uses strict a_f>c and is a labeled sensitivity. Allow a fixed number of
leaders or an adjacent-partition offset for near-boundary sensitivity, with
firm IDs/partition/mapping in diagnostics. A specified fit without component
likelihood information uses prior pi conditional on supplied c; say so.

Convert the selected partition to firm IDs, unconditional baseline shares,
fixed posterior alpha/gamma and supplied outside normalization, then call
coordination::stackelberg. Fixed demand bypasses noisy margin calibration;
costs invert that partition's true FOCs and baseline is reproduced. Each
reconstruction anchors utilities/values to observed conditional product shares;
it carries posterior coefficient, outside-share and partition uncertainty, but
does not transfer the full noisy group-quality posterior into counterfactuals.
This is a stated limitation of the initial calibration bridge, and must appear
in diagnostics and empirical reports. Threshold/adjacent-k reconstructions
also re-invert costs under their selected game; they are alternative structural
reconstructions, not role changes holding the mixture draw's costs fixed. Do not
wrap it with an ordinary registered antitrust conduct specification. In the
Bayes simulation loop branch to coordination::stackelberg_simulate and collect
the existing inherited output/welfare fields. Simultaneous and fringe retain their economic meanings.

Preserve the public Bayesian market-size convention: bayes_market/bayes_sim
insideSize is total inside quantity for both demand families. For CES convert
it to coordination's inside expenditure as
insideSize*sum(observed_price*conditional_quantity_share). Coordination's
StackelbergCES constructor then sets E=inside_expenditure/(1-R0). Passing an
already computed expenditure to the Bayesian API would multiply by average
price twice. Document this distinction explicitly; both API layers keep these explicitly different units. Physical welfare comparisons supply
observed inside quantity to Bayes and observed inside revenue to coordination.

One-step antitrust::Counterfactual objects can be translated from ownership,
costs and active-product fields to the coordination API. Reject unsupported
environment fields or multi-step scenarios explicitly if not implemented.
Ownership changes are firm-ID vectors; leadersPost is explicitly overridable.
Default merged-owner role is union of constituents' pre leader status, and
cost shocks use fixed primitive*(1+change). Support ownership-only, cost-only
and simultaneous ownership+cost scenarios. Deterministic partition selection
can be performed before fork workers to preserve reproducible RNG if needed.
Keep simulation_route and actual role mapping in BayesSim provenance.

## Stable deterministic bridge and public reconstruction

The reviewed coordination API is version 0.0.0.9003, locally installed in
/tmp/stackelberg-lib. Require that minimum version for the new bridge.
coordination::stackelberg defaults outside price to zero for Logit and one
for CES; the Bayesian bridge supplies its fitted positive outside price.
control.equ=list(implicitCheck=FALSE) skips optional nested implicit-equilibrium
checks during draw reconstruction. It never skips baseline reproduction or
analytic FOC checks. Residual diagnostics contain both actual raw profit
gradients and dimensionless gradients (price gradients divided by quantity,
quantity gradients divided by price).

Add bayes_stackelberg_model(fit, draw=1L, market=NULL, year=NULL,
role_mapping='mixture', fixed_k=NULL, role_offset=0L, insideSize=NULL,
priceOutside=NULL, labels=NULL, ...), returning the native coordination
StackelbergLogit/StackelbergCES object for a single selected market-year.
Ambiguous multi-cell selection is rejected. It uses the same reconstruction
helper as simulation, records firm roles/count/mapping, and samples no new
structural demand parameter. Multi-year simulation uses the shared partition
cache specified in the ordered-k addendum. Fixed_k and adjacent-k overrides
are explicit sensitivities, not latent sampled parameters.

Legacy as_antitrust_fit remains the registered/custom-fringe adapter. Do not
wrap the native Stackelberg object with an ordinary registry specification;
direct callers use bayes_stackelberg_model and coordination's simulation API.
Document that separation; API and target compatibility with earlier fits is
not required in this development project.

The later [nonempty-core addendum](nonempty_core_addendum.md) governs positive-k support, conditional cutoff probabilities, and cutoff-center summaries. It supersedes any earlier zero-leader Bayesian-support examples.
