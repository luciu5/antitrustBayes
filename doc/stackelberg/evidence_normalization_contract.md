# Native Stan prior normalization: audit and final contract

This audit precedes every new-model posterior fit. It preserves the prior
distributions and corrects their target normalizers. The user explicitly waived
backward compatibility for this active-development project. The final
implementation uses normalized priors directly in all Stan programs, with
model specification version 7; earlier fits must be refitted. The earlier correlation-prior note establishes
only that one particular common constant cancels; it is insufficient to
justify uncorrected cross-structure bridge odds.

Stan distribution statements (`~`) omit parameter-independent terms. Explicit
`target += ..._lpdf(...)` calls retain them. See the official
[Stan proportionality-constant documentation](https://mc-stan.org/docs/2_38/stan-users-guide/proportionality-constants.html).
The bridge-sampling authors likewise require retained density constants for
Stan evidence calculations ([Gronau et al.](https://arxiv.org/abs/1710.08162)).
The installed rstan `stan_fit.hpp:1120` calls `log_prob_propto`; the installed
bridgesampling `.stan_log_posterior` calls rstan::log_prob. Thus the historical version-6
sampling target was a prior kernel with full explicit share/margin densities.
Specification 7 retains the prior constants natively.

For an unconstrained Normal prior, the implementation is exactly
`target += normal_lpdf(theta | mu, sigma)`. For a parameter constrained positive
and assigned a Normal truncated to that support, it is

```stan
target += normal_lpdf(theta | mu, sigma) - normal_lccdf(0 | mu, sigma);
```

Vector priors subtract the support term once per active scalar. A parameter
bound alone supplies the constraint transform, not this probability normalizer.
The [Stan statements reference](https://mc-stan.org/docs/reference-manual/statements.html#log-probability-increment-vs.-distribution-statement)
explains retained constants and lower-bound truncation. The constants below
describe the historical kernel audit; the final API applies none externally.

Let C_M be the additive difference between the fully normalized prior density
and that compiled prior kernel. Then

    log p(D | M) = log Z_compiled(M) + C_M.

The same correction applies to every bridge repetition and leaves sampling,
bridge relative error estimates and posterior parameter draws unchanged.
Correct all candidate models before computing model odds. Matching active
correlation priors is insufficient: cutoff and fringe-quality dimensions vary.
For one market and baseline beta(3,100), sd_fringe=1, the fringe-minus-
simultaneous correction is -lbeta(3,100)-log(sqrt(2*pi)), approximately 12.233.
Stackelberg with a cutoff has the cutoff correction but no fringe-quality
correction. Direct ordered-k Stackelberg has neither sampled prior factor.
Its fixed discrete probabilities already sum to one; there is no k Jacobian.

## Analytically derived constants

Each row is per active scalar parameter/prior factor. Counts follow actual
parameter dimensions and flags, not the number of observation rows.

| Prior factor | Additive correction |
|---|---|
| Normal with fixed sd t | -log(t)-log(sqrt(2*pi)) |
| Conditional Normal with sampled sd | -log(sqrt(2*pi)); the sd term is already retained |
| Standard Normal | -log(sqrt(2*pi)) |
| Positive Normal(0,t) | log(2)-log(t)-log(sqrt(2*pi)) |
| Positive Normal(1,t) | -log(t)-log(sqrt(2*pi))-log Phi(1/t) |
| Beta(a,b) cutoff or bargaining unit | -lbeta(a,b) |
| Gamma(a,rate=b) degrees-of-freedom excess | a*log(b)-lgamma(a) |
| Manual bounded-correlation kernel | -log I_eta, as derived in legacy_prior_normalization.md |

All coordinate transformations declared through Stan parameter constraints
already have their Jacobians in log_prob(adjust_transform=TRUE). Do not add
those Jacobians again to C_M. Logit Cournot log-kappa translations conditional
on outside shares have determinant one. Centered/noncentered child priors must
be counted in their actual sampled coordinates. Positive prior support needs
its truncation correction even when the omitted density mass is small.

The independent lead reference
[stan_prior_normalization_reference.R](stan_prior_normalization_reference.R)
evaluates the full normalized constrained prior, reproduces the constraint
Jacobian independently, and compares it with the compiled target and generated
likelihood. Across 24 legacy graphs at three unconstrained points each, the
recovered offset agrees with the formula to 2.0e-13. See the complete
[table](legacy_prior_full_normalization_validation.csv) and
[log](stan_prior_normalization_reference.log). A further 180 compiled checks
cover centered/noncentered alpha and outside-share pools, outside regressions,
market residual scales, estimated degrees of freedom, bargaining, and
share-only/margin-only modes; see the
[flag table](legacy_prior_flag_normalization_validation.csv) and
[runner](stan_prior_flag_validation.R). All agree within 1e-9 (observed errors
are approximately 1e-13). The corrected new Stackelberg graphs remain a
separate required gate.
This is numerical density validation, not a fitted-model result.

## Final implementation contract

Use `target += ..._lpdf(...)` for every prior in all four supply/Stackelberg
programs. Standard-Normal, fixed/sampled-scale Normal, Beta and Gamma density
constants are retained directly. For each positive Normal prior add the
support normalization `-normal_lccdf(0 | mu, fixed_sd)` per active scalar.
The manual rho kernel subtracts

    log I_eta = lbeta(0.5,eta) + beta_lcdf(0.95^2 | 0.5,eta) - log(0.95).

No additional Jacobian is added for declared Stan parameter constraints;
Stan already supplies it. There is no external bridge offset, legacy evidence
mode or `bridge_normalization` tuning argument. Ordinary `bayes_compare`
compares full normalized-prior targets. Record specification version 7 and
normalization provenance. Prior distributions and their target gradients are
unchanged, while their log densities shift by the audited constants. The
separate CES physical outside-revenue conversion corrects a second epsilon
transform and can change the likelihood and its gradients; it is tested and
reported separately.

The optional heuristic HMT penalty has no audited likelihood/prior
normalization interpretation. Bridge comparisons reject active HMT penalties;
`bridge=FALSE` remains available for predictive summaries. The prespecified
synthetic and PNB cohorts use `use_hmt=FALSE`. Both outside-share pooling flags
together declare unused coordinates without a prior and are rejected.

Evidence is for the existing transformed share and inverse-margin working
likelihood. Dataset-scale Jacobians cancel only across identical observations
and transforms. Do not compare different margin masks or likelihood modes as
if they modeled the same dataset without an explicit decomposition.

New Stackelberg likelihoods marginalize positive hard k over each complete
partition unit. All cutoff conditioning and discrete prior probabilities
remain inside the likelihood. They are not removable normalization constants.

## Executed and required numerical gates

The historical version-6 normalized-prior reference recovered C_M from 72
compiled target/likelihood/Jacobian comparisons and a further 180 flag checks.
The archived `prior_kernel_offset_reference.R` validates that derivation;
it is research reference code, not a production correction path.

For every version-7 graph verify

    normalized R prior + generated likelihood + constraint Jacobian
      - compiled target = 0.

For simultaneous/fringe graphs also compare against the version-6 snapshots:
new-minus-old target must equal the derived constant plus any separately
verified likelihood correction. Logit gradients must agree,
and deterministic economic outputs must be unchanged at default scale.
The formerly required identical legacy target gate is superseded by this
constant-shift/gradient gate. Test custom prespecified fringe transition scales.
For Stackelberg additionally execute the independent whole-unit likelihood,
generated quantities, posterior probabilities, and finite-difference gradient
gates. Native compiled version-7 gates pass: 180 flag cases, 72 historical-target
comparisons, 96 whole-unit Stackelberg cases, 20 finite-difference gradient
checks, and 24 proper-target nesting comparisons. See results.md and the
executable tables. These are numerical density checks, not posterior fits.
