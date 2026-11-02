# Option B relative-share prototype: executed results, 2026-10-22

**Research prototype only. Option A remains the package specification and default.** All work is isolated in `.scratch/option-b-relative-shares-20261022`; no maintained R/Stan source, package priors/defaults, mergerBayes code/data/results or Slurm job was changed by this experiment.

Option B behaved as the identification argument predicts: with rho fixed to zero, its MC outside-share posterior equals the outside-share prior, while the matched free-intercept target shifts it downward. The PNB model-specific merger distributions change. The one noiseless known-truth Cournot example does not demonstrate better recovery under B: both targets recover approximately 20%. Normalized B evidence was computed correctly for the tested simultaneous models, but this is not a production-ready alternative: sampled-quality mixing is weak, the noiseless fixture has divergences, and native correlated/collapsed, rich-family and panel routes remain unimplemented.

The package baseline is `17c0e8aac6416c98d6db98e68c8545b559445c30`; experiments ran with R version 4.6.1 (2026-06-24). [Target provenance](target_provenance.csv) records source, data, native-library and common A/B data hashes. [Runtime versions](runtime.csv) record the actual toolchain. Compact CSVs and the [target design](DESIGN.md) are stored beside this report; executable scripts, seals and raw fits remain in the local research directory named above.

## Executed design

- Eight fits / four matched A/B pairs: PNB MC/B/C and one seven-firm noiseless Cournot equilibrium with true s0=.20, physical alpha=2.4.
- PNB reuses the sealed 39-bank, 32-observed-margin M1T1 dataset from the corrected-specification runs. Seven missing margins remain missing; PNB input-market sign and physical prices are preserved.
- Both A and B use Gaussian errors, rho=0, sampled market-scoped symmetric zero-sum quality, legacy separate quality/noise coordinates, and identical priors on shared parameters. A additionally retains its N(0,1) free intercept. The production rho-active/collapsed fits are not the one-change control.
- CmdStan 2.39: four chains, 1000 warmup + 1000 retained per chain, diagonal metric, target acceptance .999, depth 13, seed 451107, one thread/chain, at most four sampling cores. No tuning or reruns in response to diagnostics.
- One cumulative four-hour deadline covers preparation, compilation, sampling, evidence and policy. All execution finished well within it. Initial model executable checks stopped before HMC until threads_per_chain=1 was supplied; failed attempts remain recorded.
- Each successful fit receives three normal bridge repetitions, maxiter 5000, seed 19091. All 24 estimates are finite. All 200 observed-share merger states per fit are valid.

## Target and numerical checks

`DESIGN.md` and `make_kernel.py` specify the target. Q has J−1 orthonormal columns spanning the zero-sum subspace. Quality is b=Q theta with the existing normalized N(0,tau²I) coordinate prior; its subspace volume Jacobian is one. The demand likelihood is a normalized N(Q′log(r) | Q′v, sigma_s²I), a density with respect to J−1 dimensional orthonormal coordinates. It has no independent share-level equation or sampled intercept. The outside-good utility remains zero, and a=−logit(s0)−logsumexp(v) calibrates the latent utility system. PNB uses v=b+alpha_std*x; output markets use the negative price sign.

- 54 independent deterministic checks passed: basis/projector/rank at J=3/7/39; s0=.10/.20/.40 invariance/reconstruction; normalized Gaussian identities; masked joint covariance/Schur factorization and quality integration; B/C/CF equilibria, physical costs and FOCs.
- Exact-source BridgeStan comparison passed at 200 posterior states across all 8 targets. Maximum saved sampler log-density error 9.202e−13; maximum constrained/unconstrained round-trip error 2.776e−17.
- At 20 common A/B parameter points, their normalized target difference equals exactly the omitted orthonormal share-level normal density plus the N(0,1) intercept prior, maximum error 3.638e−12.
- MC B log-target changes with the outside-share logit by exactly its normal-prior change, error 5.462e−14. This proves prior factorization; it is stronger than an MCMC mean comparison.
- All 16,000 B retained states pass derived outside-share/zero-sum quality reconstruction checks. Exact source/data/config seals and cached B compilation-source bytes match in `post_execution_integrity_proof.rds`; this audit is explicitly post-execution, not described as a precompile check.
- The first generic CmdStanR B bridge attempts failed to restore native model methods and are preserved. Successful B bridge evaluation uses its actual exact-source BridgeStan library, unconstrained coordinates, propto=FALSE and jacobian=TRUE. No production fit identity or certification gate was spoofed.

Joint/collapsed identities are algebra only, not a native implementation. With observed-margin selection matrix S the transformed cross-covariance is rho*sigma_s*sigma_m*Q′S. Margin errors conditional on contrasts have dense covariance sigma_m²[I−rho² S′QQ′S]; rowwise substitution into the old conditional normal is invalid. Integrating isotropic quality adds tau²I to the share-contrast covariance. These normalized identities passed deterministic masked-data checks.

## PNB outside shares

| Conduct | A mean [95%] | B mean [95%] |
| --- | --- | --- |
| MC | 11.58% [5.17, 21.39] | 20.00% [9.12, 36.07] |
| B | 11.58% [5.34, 21.37] | 19.64% [8.78, 35.71] |
| C | 11.29% [5.69, 20.17] | 18.94% [8.62, 34.44] |

The common prior has mean 20.00006%, median 19.0497% and 95% interval [8.7786%,36.5265%]. MC B mean 20.0028% differs from its exact prior mean by 0.0194 posterior MCSE. Its SD is 7.062%, compared with 4.197% under A. PNB supply contains limited outside-share information in these rho0 fits; B should not be described as estimating s0 precisely.

A samples mean intercepts −3.123/−3.130/−3.112 for MC/B/C. B derives mean calibration intercepts −4.481/−4.480/−4.545. The meanings differ: B calibrates its normalized latent composition, whereas A also fits the absolute log-share level. Quality remains symmetric and centered in both.

## Model probabilities: restricted exploratory universe

Only MC/B/C were fitted. Their original PNB prior conditioned on these candidates is 60%/20%/20%. These are not the seven-model paper posterior or full BMA. Raw A/B marginal-likelihood levels have different observation reference measures; no A-versus-B Bayes factor is computed.

| Conduct | Restricted prior | A posterior | B posterior |
| --- | --- | --- | --- |
| MC | 60% | 43.98% | 37.95% |
| B | 20% | 18.14% | 19.49% |
| C | 20% | 37.87% | 42.56% |

Within-B log evidence (median of 3) and repetition SD: MC −137.7905(.0824), B −137.3584(.0596), C −136.5773(.0400). Every repetition and initial failure status is retained. Repetition SD does not capture possible posterior sampling bias from the mixing warnings below.

## Common merger definition and bank aggregation

PNB–Girard ownership change is identical for A/B. The public parameter/specification policy path uses posterior physical alpha and s0, observed conditional bank shares, and input-market real-rate support. It does not use zero placeholder intercept fields from the experimental native record. The observed-share anchoring recalibrates baseline utilities; it is distinct from validating latent-quality policy anchoring. MC gives zero by its unchanged conduct rule.

Table units are percentage changes in deposit rates. Each draw averages equally over 39 banks; SD/interval describe that posterior distribution of the mean bank effect. This differs from pooled bank-and-draw SD in the maintained paper helper. Both statistics are retained in policy.csv. Every draw is finite, so both aggregations use all 200 states and 39 banks; simulation metadata retains the actual draw permutation as well as the selected original-state set.

| Conduct | A mean / SD [95%] | B mean / SD [95%] |
| --- | --- | --- |
| MC | 0.000 / 0.000 [0.000, 0.000] | 0.000 / 0.000 [0.000, 0.000] |
| B | -2.064 / 0.390 [-2.919, -1.430] | -1.798 / 0.395 [-2.622, -1.125] |
| C | -3.935 / 0.960 [-6.112, -2.541] | -2.977 / 0.817 [-4.751, -1.677] |

No full BMA policy, tariff distribution or latent-quality policy certification is inferred from this prototype.

## Computational diagnostics

| Case | Spec | Div | Depth | Max Rhat | Min bulk/tail ESS | Min BFMI | Fit sec | Bridge sec |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PNB MC | A | 0 | 0 | 1.0403 | 74/83 | 0.095 | 9.70 | 131.20 |
| PNB MC | B | 0 | 0 | 1.0427 | 118/247 | 0.100 | 7.20 | 0.48 |
| PNB B | A | 0 | 0 | 1.0097 | 112/205 | 0.099 | 8.02 | 2.98 |
| PNB B | B | 0 | 0 | 1.0599 | 92/346 | 0.117 | 7.98 | 0.42 |
| PNB C | A | 0 | 0 | 1.0432 | 65/47 | 0.086 | 8.30 | 2.96 |
| PNB C | B | 0 | 0 | 1.0473 | 91/166 | 0.076 | 8.77 | 0.42 |
| Ctruth20 C | A | 19 | 0 | 1.0315 | 187/233 | 0.263 | 11.83 | 2.42 |
| Ctruth20 C | B | 29 | 0 | 1.0388 | 154/124 | 0.186 | 14.10 | 0.33 |

Fit seconds include call overhead. First A bridge includes 131s model-method setup; B BridgeStan compilation 87.8s is recorded separately from its successful .33–.48s bridge calls. Failed generic B bridge calls had near-zero duration and are not reported as successful bridge times.

PNB structural alpha/s0 Rhat is at most 1.0044 with minimum bulk ESS 1143, but nuisance quality/noise scales mix much less well: overall Rhat up to 1.0599, minimum bulk ESS 65 and BFMI down to .076. All six PNB fits have zero divergences/depth hits. This is not a broad sampler-quality pass.

**The noiseless Cournot controls have 19 (A) and 29 (B) divergences despite .999 acceptance.** A outside posterior mean .200665,95% [.187099,.217772]; B .200734,[.187354,.215834]. Physical-alpha means 2.39806/2.39699 contain truth 2.4 in their central intervals. Merger mean 2.80614%/2.80572%. These one-fixture checks preserve the economic truth but do not establish improved inference or repeated-sample recovery. No retuning, failed-state removal or threshold-based evidence suppression was performed.

## Remaining work before considering a production specification

1. Implement the transformed correlated likelihood and exact collapsed-quality target/reconstruction in native code, retaining all determinants, normalization constants and missing-margin masks; deterministic algebra alone is insufficient.
2. Extend and validate ordered-k core/fringe and Stackelberg economics, full k support, native and conditional-k evidence, all seven PNB models and complete-universe BMA.
3. Define panel/repeated-product contrast blocks and pooling/derived-intercept contracts; one-period calibration does not establish a time-invariant intercept under equilibrium shocks.
4. Validate latent-quality policy reconstruction using the derived intercept once, in physical input/output price units, in addition to the observed-share anchoring shown here.
5. Resolve the sampled-quality geometry through principled exact parameterization/integration and perform bounded noisy known-truth repeated-sample comparisons before calling B superior or changing defaults.

## Artifacts

`manifest.rds`, `DESIGN.md`, `make_kernel.py`, `option_b.stan`, `prepare.R`, `run_v2.R`, `bridge_native.R`, `policy.R`, `final_audit.R` retain the executable experiment. `deterministic_tests.csv`, `deterministic_economics.csv`, `native_target_tests.csv`, `additional_target_tests.csv`, `post_execution_integrity_proof.rds`, `fit_results_final.csv`, `posterior.csv`, `mc_prior_posterior.csv`, `bridge_repetitions.csv`, `restricted_probabilities.csv`, `policy.csv`, `attempt_statuses.csv` retain compact results. The successful audit is `final_audit_v3.log` and sealed `final_audit.rds`; earlier failed `final_audit.log`/`final_audit_v2.log` remain as attempt evidence. All raw fits, failed attempts, warnings, chain CSVs, native libraries and sealed input copies remain in the ignored research directory. No Git ref, commit, push, dependency install or production artifact was changed by this prototype agent.
