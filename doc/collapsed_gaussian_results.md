# Collapsed-Gaussian experiment results

Date: 2026-09-16. Controlled, local experiment; the Student-t/sampled-quality
likelihood default is unchanged. Failed runs and unusable posteriors are
retained. Full posteriors and logs (about 191 MiB) are also copied to the
Git-ignored, package-build-excluded Windows-drive folder
`doc/collapsed_gaussian_runs/20260916`, with a SHA-256 manifest. This durable
local archive avoids relying on WSL `/tmp` and is not added to GitHub. The separate requested information-aware `s0_pooling="auto"` API is
implemented; the controlled A/B/C comparison freezes explicit unpooled flags.

## Analytically established

The package's global orthonormal Helmert coordinates satisfy
\(CC'=I-11'/G\). Integrating their proper Gaussian prior gives exactly

\[
y\mid\theta,\tau,\sigma_s
\sim N\!\left(\mu(\theta),\;\sigma_s^2I+\tau^2ZCC'Z'\right).
\]

A (Gaussian, sampled quality) and B (Gaussian, collapsed quality) are the same
statistical model. Their true structural, quality, policy, and marginal-evidence
distributions therefore agree when reconstruction is the same; this algebra
alone does not establish empirical Monte Carlo agreement. C uses the existing
fixed-df Student-t likelihood and is a different statistical model. The share
likelihood remains the package's density for recorded log conditional shares
on its unrestricted transformed-outcome measure. The compositional DGP does
not silently replace it with a newly normalized simplex likelihood.

The collapsed Gaussian density uses firm sufficient statistics, a centered
determinant correction, and a positive within-group quadratic calculation.
That share/joint block takes
\(O(N+G)\) time and total storage, including the linear-size Stan wrapper and
automatic-differentiation tape. It constructs neither an \(N\times N\)
covariance nor a dense quality-contrast matrix in collapsed mode. Ordered-k
margin evaluation retains its additional finite-state cost.

Gaussian rho is implemented through a unit-Jacobian triangular transformation
of matched share/margin residuals. Quality affects the share covariance only;
missing margins retain their actual observation matching. With rho zero,
shares contribute once and each market's complete margin likelihood is
marginalized over ordered k. Global quality centering can couple conditional
margin errors across markets when rho is active, so collapsed multi-market
exact-k rho is explicitly rejected. The conditional quality posterior is
Gaussian and can be recovered after fitting, with full joint conditioning for
supported Gaussian-rho fits.

All accepted bridge targets retain likelihood constants, proper prior
normalizers, and transformation Jacobians. The existing HMT soft penalty is
rejected for bridge sampling because its model-dependent normalization is not
established. A normalized target is necessary for evidence estimation; it does
not make an inaccurate HMC or bridge estimate reliable. Evidence here is for
the package's common transformed-outcome measure.

Derivations: [math note](collapsed_gaussian_math.md). API, compatibility,
centering, pooling, and simulation limits:
[architecture memo](collapsed_gaussian_architecture.md).

## Numerically verified

| validation | result |
|---|---:|
| Independent expanded/collapsed log density | maximum error \(7.77\times10^{-15}\) |
| Independent direct integration | maximum error \(3.77\times10^{-15}\) |
| Conditional quality mean/covariance | maximum error \(2.86\times10^{-15}\) |
| Structured vs dense, 192 ordinary fixtures | maximum error \(2.42\times10^{-13}\) |
| Structured vs original joint covariance, 36 rho fixtures | maximum error \(2.16\times10^{-12}\) |
| Production normalized target comparison | maximum error \(2.84\times10^{-13}\) |
| Production AD vs finite differences | maximum relative error \(3.59\times10^{-10}\) |
| BF/CF/BS/CS whole-market mixtures, Logit and CES | maximum error \(3.55\times10^{-15}\) |
| Production validation | 113/113 checks passed |
| Optional Gaussian collapsed shares / fixed-df t margins | 12/12 checks passed; maximum error \(4.66\times10^{-15}\) |
| Information-aware s0 dimensions/values, Logit and CES | 204/204 checks passed |

Fixtures include uneven multiproduct firms, persistent panel group effects,
missing margins, rho zero and nonzero, the one-firm centered boundary, and the
zero-quality-variance limit. Optional mixed-margin validation uses the same
collapsed share target and checks normalized target/generated-quantity
differences for Logit/CES standard and singleton exact Stackelberg models.
Its first harness omitted the explicit missing-margin mask and was rejected
before fitting; those failures are archived separately. No optional-D HMC
benchmark was added. Independent economic-equation comparisons over
1,000 random cases per cell have maximum error \(2.22\times10^{-16}\) for
Logit and \(8.88\times10^{-16}\) for CES. Economic equations and `coordination`
were not changed.

The first subtractive quadratic implementation failed small-share-noise
controls. The replacement uses nonnegative within-group regression residuals.
At noise scales down to \(10^{-8}\), its independent centered-subspace spectral
checks have maximum error below \(2\times10^{-14}\). Original failed controls
remain in the validation CSVs. Some arbitrary-residual dense stress comparisons
were themselves ill-conditioned: densities of order \(10^{12}\) made a fixed
absolute \(10^{-8}\) gate inappropriate for that oracle. The ordinary dense
agreement gate was retained; tiny-noise checks use an exact spectral oracle.
This is an oracle correction, not evidence that every original stress failure
came from production arithmetic.

Production target and gradient timing covers all 12 requested size/layout
cells, with no failures:

| N | target seconds/call | gradient seconds/call |
|---:|---:|---:|
| 100 | about 0.00009–0.00011 | about 0.00016–0.00018 |
| 1,000 | about 0.00072–0.00094 | about 0.00138–0.00175 |
| 10,000 | about 0.0078–0.0090 | about 0.0155–0.0165 |

These are likelihood/gradient timings, not large-N HMC forecasts. Both
one-product-per-firm and two-product/four-observation panel layouts are covered,
with rho zero/nonzero. Dense N=10,000 was deliberately skipped with its memory
reason recorded; production at that size uses the structured evaluator only.
The full installed non-opt-in package suite passed, as did the earlier targeted
Student-t/rho, outside-share-pooling, and simulation opt-in tests. Independent
red-team review found no unresolved normalization, centering, matching,
whole-market mixture, or s0-default defect.

Artifacts: [production validation](collapsed_gaussian_production_validation.csv),
[production scaling](collapsed_gaussian_production_scaling.csv),
[structured reference](collapsed_gaussian_structured_reference.R),
[small-noise controls](collapsed_gaussian_small_sigma_validation.csv),
[ladder validation](collapsed_gaussian_ladder_validation.csv),
[mixed-margin validation](collapsed_gaussian_mixed_margin_validation.csv),
[test output](collapsed_gaussian_tests.txt),
[independent review](collapsed_gaussian_redteam.md),
[final independent evidence check](collapsed_gaussian_final_validation.txt).

## Sampler evidence

The frozen primary plan contains 21 fits: seven designs × A/B/C, four chains,
1,000 iterations including 500 warmup, adapt_delta 0.95, maximum treedepth 13,
and diagonal metric. Designs include exact-share stress, compositional noisy
shares with two prespecified seeds, Cournot, exact BF, exact BS, and CES. The
continuous priors, data, and HMC budget are identical between A and B. No
post-result rescue tuning was performed.

The fit gate requires zero divergences/treedepth hits, BFMI at least 0.3, maximum
rank-normalized R-hat at most 1.01, and minimum bulk/tail ESS at least 100.
Sampled quality coordinates are included in the gate. Discarding their bad
mixing would hide the phenomenon this experiment tests.

| mode | usable / 7 | total divergences | BFMI range | maximum R-hat range | median native Stan seconds |
|---|---:|---:|---:|---:|---:|
| A: Gaussian sampled | 0 | 1 | 0.053–0.123 | 1.018–1.133 | 4.899 |
| B: Gaussian collapsed | 7 | 0 | 0.719–0.933 | 1.002–1.005 | 3.832 |
| C: Student-t sampled | 0 | 1 | 0.057–0.134 | 1.031–1.311 | 6.801 |

All fits had zero treedepth hits. B's minimum bulk ESS is 605–972. B's estimated
demand-parameter ESS per native Stan second is 1.15–2.23 times A (median 1.32)
and 1.25–3.33 times C (median 1.53). These ratios are provisional: A/C's failed
mixing makes their ESS estimates and posterior comparisons unreliable. Native
Stan time is the sum of recorded warmup and sampling times; wall time is also
retained, except the first saved B fit's enclosing timer was lost during a
postprocessing failure. Cold-start wall time is not used to infer speed.

Collapsing materially improved energy exploration and mixing in these small
fixtures, including the exact-share stress case. The historical pattern of
many Gaussian divergences was not reproduced: A only had one divergence in
seven fits. Therefore this experiment supports a quality-coordinate geometry
problem, but does not establish the cause of every historical divergence or
that a Gaussian likelihood is universally preferable. A well-behaved B fit of
the same Gaussian model shows that Gaussian shares themselves need not prevent
usable HMC on these controlled data.

Collapsed posterior median alpha is 5.977–6.028 (truth 6); CES median gamma is
3.506 (truth 3.5). This is small-fixture recovery evidence, not a general Monte
Carlo recovery claim. Cross-sectional tau is about 0.058–0.064 while sigma_share
is about 0.161–0.176: the persistent/transitory decomposition is not separately
established from one-period data.

Three separately prespecified B sensitivities all passed, with zero divergences:
quality-prior width 0.5× and 2×, and active Gaussian rho. The prior-width
experiments retain the same untruncated Normal location and HMC budget; the
positive-truncated prior's arithmetic mean can change with its width. Their alpha medians are
6.024 and 6.029, while tau/sigma_share remain weakly decomposed. This limited
width check suggests structural stability on this fixture; it does not identify
the nuisance variances. The rho fit has BFMI 0.772 and maximum R-hat 1.0067,
showing computational viability for this small single-market Gaussian case,
not for all information-ladder or exact-k combinations.

Bridge runs on the noisy simultaneous-Bertrand A/B/C fits completed three
normal-method repetitions without algorithmic failure. Median log evidence:
A 24.9684, B 25.1997, C 25.8807. A/B evidence is exactly equal in the underlying
model; their 0.2313 numerical discrepancy is unresolved Monte Carlo error/bias
with A failing HMC checks. A/C bridge completion is not credible model-selection
evidence. No default or model preference is based on these estimates.

The two B-only information-ladder checks use the same likelihood and freeze
explicit unpooled controls:

| additional rung | N / groups | divergences | BFMI | maximum R-hat | status |
|---|---:|---:|---:|---:|---|
| one market, four periods | 96 / 24 | 0 | 0.852 | 1.0045 | usable |
| two markets, two periods | 96 / 48 | 10 | 0.890 | 1.0055 | failed |

The usable four-period fit has posterior median tau 0.157 and sigma_share
0.0478, closer to the persistent-quality SD 0.15 and compositional log-score
noise SD 0.04 than the one-period decomposition. This is one panel-fixture
result, not a general identification or recovery claim.

The second run took about 505 wall seconds and 503 native Stan seconds;
its tail end overlapped local policy computations. It is not included in the
primary A/B/C speed ratios. No tuning or refit was used to rescue its
divergences. A post hoc saved-draw diagnostic found 6 divergences in chain 1,
4 in chain 2, and none in the other chains. Across all 2,000 retained draws,
market intercepts correlate about -0.997 with their matched outside-share
logits. The market-intercept raw coordinates are already noncentered; there
is no clear small-tau or small-share-noise clustering. These are descriptive
findings, not causal identification of the divergences.

For this unpooled specification, a possible follow-up is the same-model
triangular coordinate change

\[
\xi_{mt}=\operatorname{logit}(s_{0,mt})+b_{market,m},
\qquad
\operatorname{logit}(s_{0,mt})=\xi_{mt}-b_{market,m}.
\]

Holding the market-effect coordinates fixed, this change has determinant one.
The original outside-share priors must be evaluated on the implied logits;
no independent new prior on xi should be substituted. It could reduce the
observed share-equation coupling without changing the model, but has not been
implemented or tested. Hierarchical pooling variants require their own explicit
coordinate/Jacobian audit. See the
[saved-draw diagnostic](collapsed_gaussian_divergence_diagnostic.txt). This is evidence that collapsing does not solve all multi-market
geometry. Separate 204-check Fixed_param validation establishes the new public
s0 auto selection and dimensions; it does not establish that auto pooling
would cure this failed HMC run. See
[ladder HMC results](collapsed_gaussian_ladder_benchmark.csv).

Artifacts: [frozen protocol](collapsed_gaussian_benchmark_protocol.md),
[21-fit benchmark](collapsed_gaussian_benchmark.csv),
[ESS/time comparisons](collapsed_gaussian_efficiency.csv),
[prior/rho sensitivity](collapsed_gaussian_sensitivity.csv),
[bridge diagnostics](collapsed_gaussian_bridges.csv).

## Policy-counterfactual evidence

Only the seven primary B fits passed the frozen fit gate. Each was reused for
100 paired posterior draws under a top-two-firm merger, a +10% marginal-cost
shock, and their combination. All 2,100 simulations completed with zero solver
failures, and all product/share/welfare summaries have 100 finite draws. The
counterfactuals required no Bayesian refitting. A/C are explicitly recorded as
unusable and cannot support an empirical A/B or Gaussian/Student-t
merger-distribution comparison at this budget.

For illustration, posterior median changes in the observed-share-weighted
price index are:

| design | merger | +10% costs | combined |
|---|---:|---:|---:|
| noisy Logit Bertrand, seed 271828 | 0.0533% | 9.7153% | 9.7200% |
| noisy Logit Cournot, seed 271828 | 0.2186% | 6.1773% | 6.2765% |
| noisy Logit BF, seed 271828 | 0.0513% | 9.7818% | 9.7866% |
| noisy Logit BS, seed 271828 | 0.0554% | 9.7135% | 9.7182% |
| noisy CES Bertrand, seed 271828 | 0.1889% | 9.9992% | 10.2068% |

These are conditional synthetic-fixture estimates, not general conduct
rankings or PNB findings. The complete CSV retains product prices, shares,
consumer/producer/total welfare, elapsed time, failure rates, and 90% intervals:
[policy comparison](collapsed_gaussian_policy_comparison.csv). Policy ESS per
second is explicitly unmeasured: these random subsamples do not provide a
reliable ordered-chain ESS calculation.

The existing simulation path anchors deterministic reconstruction on observed
shares and common structural posterior draws. It does not directly draw latent
firm quality. Quality uncertainty affects common-parameter inference, and
`posterior_quality()` is available, but full latent-quality policy propagation
is not implemented by this existing reconstruction. This material limitation
is preserved rather than silently changing the simulation/economic model.
Raw collapsed `pred_shareIn` generated quantities are common-mean predictions;
the benchmark uses conditional quality draws for comparable fitted-share
summaries.

## Unresolved issues

* Empirical A/B policy equivalence and Gaussian/Student-t policy differences
  remain inconclusive because A/C did not pass the fit gate. No gate was relaxed.
* Original historical Gaussian divergences were not reproduced in number in the
  primary fits; the observed gain there is primarily BFMI/mixing. The additional
  two-market unpooled stress fit had 10 divergences and remains unusable.
* Sparse one-period quality/share variance decomposition remains weakly learned.
* This benchmark has only two noisy-share seeds and deliberately excludes an
  outlier experiment, the full many-market recovery ladder, and PNB tuning.
* Active-rho collapsed multi-market exact-k fitting and full conditional-quality
  recovery for active-rho finite-state mixtures are unsupported and rejected.
  Explicit share-only quality conditioning remains available and labeled.
* Collapsed market-specific residual-scale hierarchies remain unsupported in the
  initial contract. More markets/periods do not require a separate model, but
  identification is not guaranteed by activating hierarchy flags.
* Raw collapsed observation-level LOO/WAIC is rejected; no valid held-out-block
  method has been implemented.
* Policy ESS/second has not been measured from full ordered chains; a paired
  100-draw random subsample can assess solver behavior and distributions, not
  reliably estimate chain ESS.
* Full latent-quality deterministic policy reconstruction remains a separate
  integration issue, described above.

## Recommendation

Keep Student-t/sampled quality as the likelihood default for now. Expose
collapsed Gaussian as a promising explicit sparse-market mode, with optional
Gaussian rho and Student-t margins as documented supported alternatives.
The controlled results support a material computation/geometry improvement;
they do not yet support changing the routine default without usable expanded/
Student-t policy comparisons and outlier/predictive evidence. The requested s0
information-aware selector is a separate API change and is implemented with
explicit overrides and proper priors.
