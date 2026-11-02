# Non-overlapping Bayesian conduct hypotheses

Economic nesting alone does not establish equality of Bayesian hypotheses.
Shared endpoints can be represented once when demand, likelihood, support,
normalization, nuisance integration, and parameter priors also agree. The
proposed universe uses canonical MC, Bertrand, Cournot, and auction hypotheses,
plus strictly interior fringe/Stackelberg families: 0 < k < n in each market.
Auction retains its own conduct category. Scientific prior weights are explicit;
endpoint weights are not redistributed automatically after examining evidence.

Observed-likelihood equality alone is also insufficient for consolidating
causal models. An executed one-firm, six-product baseline check gave identical
normalized Bertrand/Cournot targets (Logit error `1.42e-14`, CES error zero).
Yet an independent deterministic divestiture splitting a three-product
monopolist into two owners produced different Bertrand/Cournot policy prices
(maximum dollar differences `0.4116` for Logit and `0.3170` for CES), starting
from matching baseline costs and prices. B and C therefore retain separate
conduct hypotheses even in an observationally indistinguishable monopoly.
Endpoint aliases additionally require a consistent post-policy role rule.
The preserved scripts and results are in the ignored
`statistical_toolkit_runs/20260916/model-universe-monopoly-audit` archive.

Let F be an existing exact ordered-k fit, A an admissible restriction, and
Z_F its normalized marginal likelihood. Its restricted prior is the original
joint prior conditioned on A. Then

\[
Z_A=Z_F\frac{P(A\mid D,F)}{P(A\mid F)}.
\]

Both restriction masses must correspond to the entire data-generating market
partition. For market-wise admissible sets A_m and a shared continuous draw
theta_r, let

\[
a_r=\prod_m\sum_{k\in A_m}P(k_m=k\mid D,\theta_r,F).
\]

Posterior restriction mass is the posterior average of a_r. Taking the product
of separately averaged market probabilities is generally incorrect: common
continuous parameters couple markets after fitting. The independent finite
example in `model_universe_reference.R` gives full evidence 6.25, restricted
evidence 8.5, prior event mass 0.25, and posterior event mass 0.34. The incorrect
product-of-means calculation gives 0.25. Production tests must exercise the
actual encompassing-prior helper, in addition to this independent identity.

For explicit hypothesis weights omega_h, universe probabilities are

\[
P(H_h\mid D)=
\frac{\omega_h Z_h}{\sum_r\omega_r Z_r}.
\]

Use log-space evidence and normalization. A failed positive-prior evidence
calculation makes the comparison incomplete; normalizing only successful rows
changes the intended universe. Conversely, a rounded-to-zero *linear* mass
does not invalidate a finite, correctly computed log evidence ratio. Retain an
underflow diagnostic and use the finite log-space calculation when available.

Scientific dimension probabilities are sums over disjoint hypotheses. Price
competition contains B/BF/BS; quantity competition contains C/CF/CS; MC and
auction remain separate. In a multi-market comparison, a globally restricted
interior family is an all-market event, not a market-varying conduct model.

Restricted policy inference requires a further conditioning step:

\[
p(\theta\mid D,A,F)\propto p(\theta\mid D,F)\,a(\theta).
\]

Thus posterior simulation must first reweight/resample continuous draws by a_r,
then draw each k_m conditional on its admissible set and reconstruct quality
conditional on that same continuous draw. Selecting or renormalizing k alone
while retaining equally weighted unrestricted continuous draws does not give
the restricted posterior. In the finite example, the two continuous states
receive probabilities 16/17 and 1/17 after restriction. Until that draw-plan
extension is implemented and validated, a model-universe evidence table must
not be advertised as a complete restricted-policy BMA simulator.

Analytically established: the encompassing-prior identity and conditional
continuous-state weights. Numerically verified: the independent finite
reference. Native bridge estimates, posterior recovery, and policy integration
require separately executed tests; these are not implied by this derivation.

Integration status: 71 production assertions and the existing BMA tests pass.
The first native public-interface run caught and resolved a named-list routing
defect. A subsequent three-repetition bridge attempt on a saved fit failed
because its native Stan model instance was unavailable after reload, including
under the unchanged baseline DLL. The comparison correctly remained incomplete.
The failed logs are preserved; restoration of a live, target-matched context
without refitting is being independently validated. This runtime issue does
not establish an improper posterior density.

A further fixed-state audit exposed active-hyperprior identity relying on
metadata and a redundant manifest-label identity. Active Stan hyperprior
values now govern the fingerprint; conflicting active provenance is rejected,
while inactive priors remain excluded. Named manifest labels are not a distinct
prior once their actual active hyperparameters are known. The initial failed
audit is retained. These are comparison guards, without changes to fitting
priors or conduct equations.
