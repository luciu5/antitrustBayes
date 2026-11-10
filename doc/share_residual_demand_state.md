# Demand state and measurement error in specification 15

New Gaussian relative-share Logit/CES fits use `share_residual_model = "demand_shock"`.
`"auto"` selects it for relative Gaussian observations; absolute-share fits
retain `"combined"`. The latter explicitly selects the original statistical
model when fitting new data as well.

In orthonormal within-market-year log-share contrasts, the observation equation
contains structural demand, persistent centered quality, a period demand shock
and independent measurement error. Write period-demand SD as omega, measurement
SD as m, margin-error SD as sigma_m, and demand-margin correlation as rho_d.
The observed residual covariance is

\[
\begin{pmatrix}
\omega^2 + m^2 & \rho_d\omega\sigma_m \\
\rho_d\omega\sigma_m & \sigma_m^2
\end{pmatrix}.
\]

The previous share-residual prior now governs omega. The quality prior and
its total/fraction coordinate Jacobian are preserved. The correlation prior
governs rho_d. Independent measurement precision has a proper Normal prior
on its primitive log SD, default `Normal(log(.10), .5)`; no extra exponential
Jacobian belongs on that primitive log coordinate. Set the median and log SD
through `prior_sigma_measurement` and `prior_log_sigma_measurement_sd` in
`priors`. These defaults remain fixed across prior manifests.

`sigma_share` and `sigma_share_market` report total observed residual SDs.
`sigma_demand_shock` and its market counterpart report economic period-demand
SDs; `sigma_measurement` is the independent measurement SD. `rho` is the base
effective observed correlation; `rho_market` contains the market-specific
effective correlations when residual scales vary. `rho_demand_shock` is the
correlation with economic period demand. In the original combined model,
measurement SD is zero and there is no additional sampled measurement coordinate.

Collapsed inference integrates quality and period demand. Posterior policy
reconstruction first realizes any discrete role state, then recovers quality
using the total conditional covariance and recovers period demand conditional
on that same quality and role draw. Both recoveries use the projected relative
contrasts, observed margins when included in fitting, and the selected row's
measurement precision. Sampled-quality fits condition period demand on their
sampled quality. Outcomes excluded from the fitted likelihood do not become
conditioning data during recovery.

Version-4 policy plans store persistent quality, observation-level period
demand, measurement/demand scales, outside shares and role states. These fields
are bound to their posterior indices and included in integrity fingerprints.
The structural demand state is quality plus period demand; it remains fixed
through baseline, merger and sequential policies. Measurement noise remains
in replicated observations. Derived outside-share normalization uses the full
structural demand state. Cost/value inversion and structural solver conventions
are unchanged. Explicit `demand_anchoring = "observed"` retains its calibration
semantics; no numerical failure silently switches anchoring or statistical model.

`posterior_quality()` continues to return persistent quality alone.
`posterior_policy_draws()` adds period demand and recomputes the complete-state
normalization. Posterior predictive replication draws fresh period demand plus
measurement error through their total covariance; it does not add fresh noise
on top of a previously recovered period-demand realization.

This follows the structural interpretation of unobserved demand in
[Berry (1994)](https://www.its.caltech.edu/~mshum/gradio/papers/Berry1994.pdf):
mean utility contains an unobserved demand characteristic. Our persistent plus
period demand is its Bayesian structural analogue. An explicit measurement
equation permits imperfect reconstruction; the package retains its Bayesian
likelihood and identification assumptions rather than adopting IV/GMM estimation.

The independent transitory variances are not generally separately identified.
With free rho_d, different values of (omega, m, rho_d) can leave both omega^2+m^2
and rho_d*omega unchanged. Posterior updating therefore remains conditional on
the proper variance priors and correlation support. A narrower interval or an
improved conditional baseline alone does not establish separate identification
or out-of-sample predictive recovery.

Specification-14 saved fits remain the original combined model and use their
verified original kernels, source seals, data and priors. They are never
relabeled as new inference. Historical version-1–3 draw plans retain their
original semantics. Supplied `BayesSpecifiedFit` primitives retain supplied-state
calibration and do not acquire invented posterior quality or shocks.

The new default changes the statistical model and requires fresh inference and
appropriate evidence before adopting new candidate-model probabilities. No old
model weights are transferred or renormalized over successful candidates.
The completed [banking pilot](validation/pnb_measurement_prior_diagnostic_20261109.md)
and [variance experiments](validation/share_residual_experiment_record_20261109.md)
remain dated evidence; current integration validation is recorded separately.
