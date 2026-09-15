# Nonempty Stackelberg core: pre-fit support revision

User steering on 2026-09-14 restricts new Stackelberg estimation to k>0.
This changes no published simultaneous/fringe model. Deterministic economics
and the mathematical nesting oracle retain k=0 as an explicit diagnostic.
New Bayesian components have k=1,...,N; market-shared counts have support
1:min_t N_mt. Supplied priors with positive mass outside this support fail.
The direct default is uniform over the supported positive counts, never 1:5.

Requiring a nonempty first-stage leader sector is an economically coherent
restriction. It removes the empty-leader game, but k=N remains ordinary
simultaneous conduct and observationally nests that comparator. Its posterior
mass must be reported rather than suppressed for model discrimination.

For the cutoff prior, let C=c+tau*logit(U) be the publicly realized common
boundary and x_f=(size_f-c)/tau. Condition this lottery on C<size_1, with
size_1 the largest size. The exact probabilities of the positive hard games are

    log pi_k^+ = log pi_k - log_sigmoid(x_1),  k=1,...,N.

Here pi_k are the audited unconditioned share-gap probabilities. Since
sum_{k=1}^N pi_k=P(C<size_1)=sigmoid(x_1), the conditional probabilities
sum to one. Dropping k=0 without this continuous-parameter-dependent
normalization changes the cutoff likelihood and is forbidden. The
normalization stays inside each whole-equilibrium mixture; it is not a
constant evidence correction. The fixed cutoff Beta prior remains unchanged.

A firm's conditional leader probability before margin data is

    P(f in L | c, C<size_1) = sigmoid(x_f)/sigmoid(x_1).

Every realized component remains a hard noncooperative game, and at least one
largest firm is a leader. Under cutoff-size ties, the lottery cannot split the
tie: if all firms tie, only k=N has mass. Direct ordered-k priors can still
split ties by the declared public firm-ID ordering.

For numerical implementation, use stable log-CDF expressions and subtract
log_sigmoid(x_1); do not divide saturated floating-point CDFs. A one-firm
market has pi_1^+=1 and no cutoff information in the role likelihood.

The nominal center cutoff count #{size_f>c} can be zero, even though the
conditioned game cannot. Preserve it only as a separately labeled uncapped
cutoff diagnostic. The nonempty hard-center summary and threshold simulation
mapping are max(1,#{size_f>c}). Actual marginalized posterior P(k) and role/core
summaries use the conditional hard-game probabilities, not this center mapping.
When c exceeds the largest size, the truncated boundary lottery has its mode
at the largest-size boundary; the threshold mapping selects the largest
firm/tie block (see tie handling below), not a nonexistent zero-leader game.

If the largest firms tie and the uncapped count is zero, a cutoff-compatible
threshold must select the entire largest tied block. Thus the precise
cutoff-center mapping is the smallest supported positive count not splitting
a size tie that is at least the uncapped count. With unique sizes this is
max(1,uncapped_count). Direct fixed-k/adjacent-rank sensitivities may split ties
and are explicitly rank interventions rather than realizations of a scalar
boundary.

The same conditioning is used by the independent synthetic cutoff DGP:
replace U by U*sigmoid(x_1) before constructing C. This revision precedes all
new posterior fits and keeps seeds, other truths, priors and sampling settings.
Empirical recovery and new compiled-Stan validation remain pending.
