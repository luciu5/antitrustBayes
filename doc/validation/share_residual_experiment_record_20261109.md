# PNB share-residual experiments: recorded evidence

Status at 2026-11-09: these are completed diagnostic experiments, not instructions
to resume them. Raw fits, policies, seals and original reports are preserved.
No experimental variance or measurement model has become a production default;
the paper's estimates and model probabilities have not been replaced by them.

The earlier total/fraction alias repair fixes extraction of the actual quality
SD during recovery. It changes neither the statistical model nor its posterior,
and has not been shown to explain the share shortfalls. See the separate
[demand-state report](pnb_demand_shock_diagnostic_20261109.md) for its ten-draw,
no-refit verification.

The seven-model snapshot residual-prior experiment in mergerBayes is recorded in
`results/banking/share_residual_pilot_20261108/README.md` and
`doc/share_residual_model_fix.md`. Loosening that prior produced mixed share
improvements and worse margin predictions; its short fits do not establish a
preferred calibration. The history prior comparison was not executed.

| Completed experiment | What was tested | Recorded finding |
| --- | --- | --- |
| [Independent period demand, measurement SD fixed at .10](pnb_demand_shock_diagnostic_20261109.md) | Persistent quality plus jointly recovered period demand; measurement excluded from both equilibria | Strong snapshot and history share improvements, little margin change, zero divergences in new fits. The original reused snapshot reference was short; subsequent variance experiments supplied a stronger 4-chain reference. The fixed .10 is a diagnostic assumption. |
| [Unrestricted institution-assets variance slope](pnb_asset_variance_diagnostic_20261109.md) | One regularized coefficient on globally scaled full-institution assets | Positive fitted slope; snapshot baseline shares deteriorate severely. Not selected. |
| [Fixed inverse institution-assets variance](pnb_inverse_asset_variance_diagnostic_20261109.md) | Residual variance proportional to inverse full-institution assets | Snapshot large-bank errors fall below one percentage point. History remains less satisfactory and normalized likelihood worsens. Zero divergences in both new fits. |
| [Inverse assets plus linear time trend](pnb_variance_time_trend_diagnostic_20261109.md) | Regularized linear log-variance trend in history | History large-bank errors within 1.24 percentage points; little margin change; partial likelihood recovery; zero divergences. |
| [Inverse assets plus hierarchical year effects](pnb_year_variance_hierarchy_diagnostic_20261109.md) | Noncentered year log-variance contrasts with a proper hierarchy prior | History large-bank errors within 1.82 percentage points; improved year-specific residual calibration and likelihood relative to linear trend; zero divergences. The recorded hierarchy prior has an infinite prior mean variance multiplier and requires revision/sensitivity before production adoption. |

The inverse-size correction materially improves snapshot share predictions.
Hierarchical year effects improve historical temporal residual fit. Both changes
were propagated through the fitted covariance, conditional quality draws and
pre/post equilibrium reconstruction in those experiments. The reported policies
already use that correction; they are not observed-share anchored policies.
Full-institution assets are bank-year totals across markets, never summed
repeatedly over cases; history retains distinct market-year observations.

For scale, the observed 1960 shares of First Pennsylvania / PNB / Girard are
22.957 / 21.052 / 14.230 percent. The stronger original snapshot reference gives
18.675 / 17.373 / 9.121; inverse assets gives 22.346 / 20.457 / 13.274.
The historical original gives 22.021 / 19.101 / 14.110; inverse assets with
hierarchical years gives 23.101 / 19.235 / 14.206. These are 100-draw posterior
means after solving each draw's baseline equilibrium, not just intervals that
contain observations. All completed experiments retain their invalid-draw counts,
independent FOCs, no-op checks and solved-pre Paasche denominators.

The current authorized next diagnostic replaces the fixed .10 measurement SD
with a prior and **jointly estimates it**, using at most two new Bertrand fits
(snapshot first, then history if sampling is adequate). It is separate from the
variance experiments. A saved-fit conditional split is not a substitute for
that fitted model. The inverse-assets/year-effects route remains the fallback
if sampling or economic diagnostics are unsatisfactory. Package promotion and
other conduct models remain later work; no full rerun, bridge sampling, model
averaging, cluster job or push is authorized by this diagnostic.
