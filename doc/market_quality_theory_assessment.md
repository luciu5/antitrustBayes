# Current theory and computational assessment

Assessment date: 2026-09-17. This note assesses the implemented specification;
it does not change the likelihood, priors, economic equations or defaults.

## What is established

The noncooperative Logit/CES Stackelberg FOCs have separate derivations and
independent follower-response/direct-optimization checks in
[the economic audit](stackelberg_math_audit.md). They describe interior
stationary equilibria on regular follower branches. The checks do not prove
global existence/uniqueness for arbitrary data or empirical identification of k.

For Gaussian disturbances, integrating the per-market orthonormal quality
coordinates is exact. Sampled and collapsed market-quality specifications have
the same marginal model for the common parameters. Their parameter-vector
dimensions differ, so their joint targets at arbitrary sampled qualities are
not supposed to agree before integration. Proper prior normalizers, Gaussian
constants and the relevant coordinate Jacobians are retained.

Conditional on common parameters, independent market-quality blocks yield

\[
 p(D\mid\theta)=\prod_m\sum_k\pi_{mk}L_m(k;\theta).
\]

This includes matched share/margin correlation and missing margins. The
production implementation evaluates whole-market mixtures without enumerating
the Cartesian product of roles. Independence here is conditional on common
parameters; integrating their posterior can induce dependence across markets.
The independent small Cartesian-product oracle and native gradients pass.

Market-specific quality is a substantive multi-market prior change. It assigns
separate persistent quality to the same firm in different markets, with a
common shrinkage scale. It is not the old global model under new coordinates.
Single-market scope changes are exact equivalences: the PNB market/global
Gaussian regression produced identical checked posterior draws, including the
same two Bertrand divergences.

## The share observation model needs a precise interpretation

The current share likelihood is a Gaussian or Student-t density for the vector
of recorded log conditional shares. It is not a probability density on the
share simplex. Its fitted exponentiated locations need not sum to one.
Normalized Gaussian factors and a fully specified bridge target do not by
themselves make it a compositional observation model. Marginal evidence and
model comparisons must retain the same declared outcome measure.

An inherited single-market restriction makes this distinction economically
important. For one period, one observation per quality group, no smooth-fringe
quality shift, and the existing zero market-intercept normalization, the Logit
predictor is

\[
 \mu_j=\eta+b_j+\beta(p_j-\bar p),\qquad
 \eta=\operatorname{logit}(s_0),\quad \sum_jb_j=0.
\]

Consequently

\[
 \bar\mu=\eta,\qquad
 \overline{\log r}\mid\theta\sim
 N\left(\eta,\frac{\sigma_s^2}{G}\right)
\]

under the marginal Gaussian share block. The centered quality covariance has
zero variance in the common direction, so tau cannot absorb a mismatch in
this average. Rho changes joint conditioning with margins, not this marginal
share identity. Here eta is the effective logit of the package's bounded
outside probability, rather than its unbounded raw parameter.

For observed conditional shares, Jensen's inequality gives

\[
 \overline{\log r}\leq-\log G.
\]

Thus the *share-marginal, small-noise limiting point* would imply
\(s_0\leq1/(G+1)\). The archived PNB inputs have G=39 and
\(\overline{\log r}=-5.2773166184\), yielding the diagnostic point
\(\operatorname{logit}^{-1}(\overline{\log r})=0.0050801741\).
This is not a posterior estimate, a universal bound on posterior s0, or proof
of the cause of the Cournot anomaly. Margins, finite noise and priors also
matter. The argument does not automatically extend to unbalanced panels,
repeated products per quality group, or free multi-market intercepts.

The implication is that centering plus the intercept restriction supplies
substantive information about outside demand. A free inside utility level
relative to the outside good is economically meaningful; fixing it cannot
always be described as a harmless normalization. Price standardization is
algebraically equivalent to raw prices only with the corresponding intercept
and prior transformation. These choices require an explicit identification
interpretation before treating small outside-share estimates as data-driven
economic findings.

The independent review raised a possible outside-price interface mismatch.
The exact Logit identity resolves the algebra without changing the outside
share:

\[
 \delta_j^{physical}=B+b_j-\beta(\bar p-p_0),\qquad
 \eta+\delta_j^{physical}+\beta(p_j-p_0)
 =\eta+B+b_j+\beta(p_j-\bar p).
\]

Thus a centered-price quality/intercept prior corresponds to a particular
physical utility-level prior. It does not require relabeling or converting
the physical s0 alone. Converting only s0 before simulation would also change
the outside share used in supply FOCs, and cannot be treated as a harmless
interface fix. The limitation is the implicit level prior/restriction, rather
than an established missing-term coding defect.

For CES the exact existing predictor is

\[
 \mu_j=\eta+B+b_j+\gamma\overline{\log p}
                 +\gamma\log(p_0/p_j),\qquad
 \log A_j=B+b_j+\gamma\overline{\log p}.
\]

The mean-log-price term must not be dropped when simplifying the code.
Multiplying both inside and outside prices by a common factor c leaves all
physical CES price ratios unchanged but adds gamma*log(c) to this predictor
if B is unchanged. Keeping the same physical s0 and gamma requires
B_new=B_old-gamma*log(c), with the corresponding prior transformation.
The single-market fixed B=0 path does not permit that remapping. The current
CES share model therefore has a price-unit convention that merits a separate
robustness audit. This is a verified predictor identity, not a completed native
target or posterior rescaling experiment. No target has been changed.

[The independent price-coordinate script](market_quality_price_coordinate_reference.R)
passes 33 algebra checks, including both Logit signs, three CES elasticities
and common price factors .01, 2 and 100; maximum error is 8.0e-15. These
checks establish the equations above and resolve the review's initial CES
simplification. They do not validate a new likelihood or identify the cause
of existing sampler failures.

An ILR/logistic-normal share robustness model would remove common shifts from
relative-share information. That changes the statistical model and could
weaken outside-share identification; independent outside-share information,
margin information or defensible level priors would then have to supply it.
It is not a computational repair to make silently in the present task.

## Supply measurement and counterfactual mapping

Supply predictions use observed firm shares, adjusted by inferred outside
share, in the audited FOCs. Shares also receive observation noise in the
demand block. This is a conditional structural-moment specification; it does
not integrate latent true shares through every supply equation as a complete
measurement-error equilibrium model would. Prespecified observed rankings
likewise condition k roles on the recorded size ordering.

The latent-quality policy path draws quality consistently with the structural
draw and realized k, normalizes relative inside scores while preserving the
outside-share draw at reference prices, holds the observed-share-recovered
cost primitive fixed, and solves pre/post equilibria. This is a documented
posterior mapping. It propagates relative quality uncertainty, but is not a
joint posterior over all latent true shares, costs and equilibrium prices.
The [policy math note](latent_quality_policy_math.md) states this distinction.
Its historical global-quality language is superseded by the scoped quality
maps for new fits; the normalization and fixed-cost interpretation remain.

Rho represents correlated share and margin disturbances. It is not by itself
a general identification strategy for endogenous prices. Sparse data may
inform total residual variation better than its allocation among persistent
quality, transitory share error and their correlation with margins.

## What the current sampler evidence establishes

All seven original small collapsed-Gaussian primary fixtures passed the gate.
The earlier two-market rotation experiment used global quality and rho off;
both rotated fits passed with zero divergences and roughly 12-fold wall-time
improvements. The new initial three-market/two-period comparison instead uses
legacy coordinates and active rho; none of its six fits passes the complete
frozen gate. Those failures do not test whether rotation works in the new
configuration.

The PNB summary now has seven of eight proposed-model fits passing the full
gate. Bertrand has two divergences, BFMI .863, maximum R-hat 1.00058 and
minimum bulk ESS 7,474. The other models have zero divergences, minimum BFMI
.912--.946 and maximum R-hat below 1.001. Single-market scope equivalence
means this improvement over historical Student-t runs cannot be attributed
to the new market-specific prior alone.

The two new matched Gaussian comparisons show the following costs:

| DGP seed | collapsed / sampled wall seconds | collapsed / sampled total leapfrog steps | native ms per leapfrog, collapsed / sampled |
|---|---:|---:|---:|
| 271828 | 125.25 / 67.15 | 1,087,392 / 1,097,684 | .114 / .059 |
| 314159 | 137.25 / 81.67 | 1,082,660 / 1,217,297 | .125 / .066 |

Leapfrog totals include warmup and sampling. Native time per leapfrog is an
amortized diagnostic, not an isolated gradient benchmark: it also includes
Stan adaptation and generated quantities. Collapsing reduces steps by only
about 1--11% here while the amortized step cost is about 1.9 times larger.
An O(N+G) algorithm can still have a higher constant cost than conditional
Gaussian evaluation. Removing quality coordinates has not removed the
remaining outside-share/intercept ridge in these legacy-coordinate fits.

## Remedies and the information ladder

1. Test the existing rotation on the exact new saved fixtures, changing only
   coordinates and preserving priors, economic initial states and HMC budget.
   Retain the original primary failures. This is the first computational step.
   The implemented rotation acts on the existing unbounded outside-share
   coordinate and retains its bounded-link correction. The effective share
   sum is also exposed for diagnostics; dropping that correction would change
   the target, particularly near the outside-share bounds.
2. Diagnose any remaining divergent states before choosing another change.
   Longer chains can improve Monte Carlo precision but do not establish that
   biased exploration from divergences has disappeared. Any sampler-setting
   sensitivity must be separately recorded rather than replacing a failed run.
3. Compare sampled and collapsed market quality under the same successful
   market/outside-share coordinates. Choose on structural/policy ESS per second,
   memory and reliability; fewer coordinates alone is not a speed guarantee.
4. Audit sensitivity of outside-share and policy inference to the inherited
   share-level restriction. A compositional alternative or latent-cost model
   is a distinct research specification requiring explicit approval/design.

More markets add local quality effects as well as observations. With fixed
period counts, the observation count per firm-market effect need not increase;
the quality dimension grows as sum_m(G_m-1). More periods for existing local
groups do inform each persistent effect. Consequently collapsing may remain
valuable in many sparse markets, while sampled quality can become more
competitive in well-measured panels. Common hyperparameters can also become
better informed across markets. No unconditional rule that collapsing is
always faster, or unnecessary with many markets, follows from these facts.

At fixed structural parameters, a balanced panel has G-1 share covariance
directions with variance sigma_s^2+T*tau^2 and GT-G+1 directions with variance
sigma_s^2. Repeated periods create transitory-noise directions; merely adding
markets with the same period count does not do so for each local quality.
The quality-information diagnostic tau^2*a_g accounts for matched margins
through a_g=n_missing/sigma_s^2+n_matched/[sigma_s^2*(1-rho^2)]. Neither
identity proves full-posterior variance identification or a runtime crossover.
[The small dense information reference](market_quality_information_reference.R)
passes 12 covariance/curvature checks, including missing margins and both rho
signs, with maximum absolute error 1.07e-14.

The current recommendation is to retain both exact Gaussian parameterizations
and hold the package-default change. The quality integration and role
factorization are mathematically well supported; the broader interpretation
of share levels, outside-share identification and the policy mapping still
requires the stated robustness work.
