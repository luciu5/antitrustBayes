# Model universe: specification before implementation

Status (reviewed 2026-10-13): original design contract with a public implementation
in [model_universe.R](../R/model_universe.R), [usage examples](statistical_toolkit_usage.md),
and [focused tests](../tests/testthat/test-model-universe.R). The contract below
remains the original design; implementation does not establish broad model-selection
recovery or restricted-family policy BMA. See the [current index](README.md).

## Economic states and Bayesian hypotheses

For a fixed market with n firms, exact fringe endpoints satisfy
`BF(0) = CF(0) = MC`, and all-core endpoints satisfy `BF(n) = BS(n) = B`
and `CF(n) = CS(n) = C`. These are statements about the economic FOCs.
An evidence identity also requires matching likelihoods, outcomes, priors,
outside normalization, and nuisance integration. A smooth fringe fit with a
quality shift is not silently identified with an exact fringe endpoint.

The initial supported universe consists of separate atomic simultaneous
conduct hypotheses and strictly interior exact ordered-k families. Demand
families and likelihood specifications remain explicit dimensions. Shared
economic endpoints are represented once within a matching statistical model.
No arbitrary automatic reassignment of model prior mass is permitted.

## Restricted-family evidence

For an encompassing family F with prior probability P(A | F) for a retained
set of ordered-k states A, the normalized restricted hypothesis has evidence

\[
Z_A = Z_F\,\frac{P(A\mid D,F)}{P(A\mid F)}.
\]

Both masses must be evaluated with the fit's recorded k prior and complete
market-state posterior probabilities. For a multi-market restriction,

\[
P(A\mid D,F) = E_{\theta\mid D,F}
 \left[\prod_m P(A_m\mid D_m,\theta,F)\right].
\]

Multiplying draw-averaged market probabilities is incorrect because common
structural parameters induce posterior dependence. The existing
`bayes_restriction_compare()` implements the encompassing-prior identity;
reuse it. Removing endpoint rows without the evidence correction is wrong.

If the retained interior set has zero prior mass, reject the requested
hypothesis explicitly. The uncertainty of rare-set restriction estimates and
of bridge evidence must remain visible.

If evidence fails for a positive-prior hypothesis, do not silently normalize
posterior weights over only the successful rows and present them as a complete
universe posterior. Return an explicit incomplete-evidence status and withhold
universe posterior probabilities. Numerically underflowed restriction masses
are not proof that a state has exactly zero posterior probability.

## Intended interface

Introduce an explicit universe descriptor used by `bayes_compare()`, rather
than a competing comparison interface. The descriptor names the atomic
hypotheses, source fits, retained k restrictions, and scientific prior weights.
It may request `k = "interior"` for exact families. Atomic simultaneous fits
retain their evidence. Exact-family rows receive corrected restricted evidence
and recorded prior/posterior restriction masses. Ordinary comparison without
a descriptor remains available, with overlap diagnostics rather than an
unannounced change in priors.

Return an atomic hypothesis table and support aggregation by economic
dimensions: price conduct, quantity conduct, other conduct, simultaneous,
fringe, Stackelberg, and the presence of a strict strategic core. Auction is
kept in its own category; it is not mechanically relabeled price competition.
Do not infer that these categories cover unmodeled economic states.

Use market-level firm IDs and predetermined ranking metadata. A panel's
market-level k restriction refers to its existing partition unit and union of
ranked firms. Reject unsupported period-level reinterpretations. The same
comparison data/prior contract used by `bayes_compare()` continues to apply.

## Required controls

Test exact evidence correction against independent finite mixtures; common
parameter multi-market integration; zero-mass restrictions; endpoint detection
with multiproduct ownership; explicit model priors; no duplicate atomic
conduct states within one matching statistical model; likelihood/demand
differences not silently consolidated; generated dimensional probabilities
summing correctly; and unchanged ordinary comparison behavior when no
universe is specified. No HMC or refitting is required for algebraic fixtures.

The descriptor does not yet implement market-specific conduct or a learned
population distribution of k. Those are separate statistical extensions.

## Restricted posteriors for policy simulation

Restriction also changes the posterior of common structural parameters. To
simulate a retained hypothesis A from encompassing-family draws, reweight
draw r by `P(A | D, theta_r, F)`, calculated as the within-draw product across
markets. Then draw each market k conditionally inside its allowed set. Drawing
an interior k from the unrestricted common-parameter posterior is insufficient.
Record importance ESS and the retained sets; report unsupported or poorly
resolved restrictions rather than silently averaging full-family simulations.

The comparison result should remain a data-frame table compatible with the
existing `bayes_compare()`/`bayes_weights()` contract. Store the descriptor,
source comparison, evidence completeness, and restriction diagnostics as
attributes. Exclude coordinate choices and requested/resolved provenance from
statistical duplicate fingerprints; retain effective pooling and prior choices.
