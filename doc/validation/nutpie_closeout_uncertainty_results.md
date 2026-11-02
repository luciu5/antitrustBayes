# Conditional-policy uncertainty audit

Collected 2026-10-21 02:07:50 UTC.

This saved-policy-only audit preserves every frozen point estimate, k weight, draw plan, failure and original approximate MCSE.
The frozen helper shares its latent-quality seed across components/backends, so independence of policy Monte Carlo errors is not assumed here.

For each route/regime, one k component at a time is resampled within each retained chain using circular blocks32 and200 fixed-seed repetitions.
The sum of the five component statistic SEs is a conservative first-order unknown-covariance bound; the retained evidence-weight SE is then added.
Comparison bounds add both routes' total SE bounds, including homogeneous comparisons that share latent-quality noise.
The bootstrap and normal bridge-weight error estimates remain approximate. Minkowski's inequality bounds the linearized sum without assuming independent component errors; it does not turn finite-sample nonlinear quantile bootstraps into exact confidence intervals.
In particular, atoms, low tail ESS, chain nonstationarity and insufficient block length can limit the approximation. Actual divergences/ESS/attrition remain in the original reports.

Matching exact constants retain unavailable MCSE and are exempt identities; different constants fail agreement. Existing .25price-point/.10relative-SD/.03probability precision caps are unchanged.
Within-component policy omissions retain original k weights; failure mass remains explicit. No claim of an unconditional policy distribution is made when any component has attrition.

Revised comparison rows: 28 ; agreement failures: 0 ; unavailable/imprecise: 13 .

Tracked tables: nutpie_closeout_uncertainty_components.csv, nutpie_closeout_uncertainty_bounds.csv and nutpie_closeout_uncertainty_comparisons.csv.
Raw component bootstraps, exact source/input seals and the audit result remain under .scratch/nutpie-closeout/uncertainty-audit.
No HMC, native target evaluation, policy reconstruction or new equilibrium solve was performed for this audit.
