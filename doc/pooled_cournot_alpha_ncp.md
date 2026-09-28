# Non-centering the pooled Cournot alpha hierarchy

Status: implemented, validated, not yet the subject of a fresh full-ladder
production run (that is `mergerBayes`'s responsibility after re-pinning).

## Motivation

The `mergerBayes` Supreme Court banking-ladder Rung 3 fit (`year_cross_section`,
6 markets, 176 products, pooled alpha, `s0_pooling="auto"`) showed persistent
HMC divergences on the pooled Cournot conduct model that were substantially
worse than the other three conduct models at every tested `adapt_delta`
(`.99`: MonCom 21 / Auction 14 / Bertrand 12 / **Cournot 37**; `.999`: 11 / 0 /
1 / **19**). Locating the divergent draws in the saved fits showed 46-84% of
Cournot's divergent draws sit in the bottom decile of `sigma_log_a_vec`
(posterior q05 ~ 0.02), while the other three conduct models — which use a
non-centered alpha hierarchy — have essentially none of their divergences
there. Inspecting `inst/stan/logit_supply.stan` confirmed the cause: pooled
Cournot used a *centered* alpha coordinate (`log_kappa_market`, with
`target += normal_lpdf(log_a_market_derived | mu_log_a, sigma_log_a)`) while
Bertrand/MonCom/Auction use the standard *non-centered* `r_market_a_raw`. This
is the textbook Neal's-funnel geometry: a centered hierarchical parameter
whose conditional scale shrinks toward zero.

`log_kappa_market` (kappa ~ log(alpha) + mean_log_s0) was originally designed
for the tiny-`s0` single-market PNB regime, where margins identify the product
alpha·s0 rather than alpha alone. On the Rung 3 data `s0` is order 0.2-0.3 (not
tiny) and alpha is well identified per market (posterior sd 0.10-0.22), so the
kappa coordinate buys nothing there and instead puts the sampler on a
0.83-0.96 correlation ridge between `log_kappa_market` and `mean_log_s0_market`
that the alpha coordinate itself does not have.

## The change

Both Logit kernels (`inst/stan/logit_supply.stan`,
`inst/stan/logit_stackelberg.stan`) had the identical construction; both were
edited identically. The **unpooled** Cournot path (`log_kappa_vec`, single
market-wide alpha, used when `pool_alpha_by_market == 0`) is untouched — it
remains the correct coordinate for the tiny-`s0` regime it was designed for.
Only the **pooled** path changed.

**parameters** — `r_market_a_raw`/`log_a_market_pool` are no longer restricted
to `is_cournot == 0`; `log_kappa_market` is deleted:
```stan
vector[(pool_alpha_by_market == 1 && use_centered_alpha == 0) ? N_market : 0] r_market_a_raw;
vector[(pool_alpha_by_market == 1 && use_centered_alpha == 1) ? N_market : 0] log_a_market_pool;
vector[(is_cournot == 1 && pool_alpha_by_market == 0) ? 1 : 0] log_kappa_vec;
```

**transformed parameters** — the pooled-Cournot branch now uses the same
coordinates as every other pooled conduct model:
```stan
} else {
  mu_log_a = mu_log_a_free[1];
  if (use_centered_alpha == 1) log_a_market_derived = log_a_market_pool;
  else log_a_market_derived = mu_log_a + sigma_log_a * r_market_a_raw;
  a_market = exp(log_a_market_derived);
}
```
`log_a_market_derived` and `mean_log_s0_market` stay declared and exported
(the former stays 0 for non-Cournot, as before), so no other conduct model's
draw schema or values change.

**model** — the pooled-alpha prior no longer excludes Cournot; the separate
centered-density term is gone:
```stan
if (pool_alpha_by_market == 1) {
  if (use_centered_alpha == 1) target += normal_lpdf(log_a_market_pool | mu_log_a, sigma_log_a);
  else target += std_normal_lpdf(r_market_a_raw);
}
```

## Exact equivalence

This is a reparameterization, not a model change. The map from the old
centered coordinate to the new noncentered one is
`kappa_m = mu_log_a + sigma_log_a * r_m + mean_log_s0_m`, a unit-Jacobian shift
of `r_m` (for fixed `mu_log_a`, `sigma_log_a`, `mean_log_s0_m`) composed with
the standard noncentering identity. The exact algebraic identity is

```
log_prob_new(r, ...) - log_prob_old(kappa, ...) = N_market * log(sigma_log_a)
```

at any matched economic state (same `a_market` for both), where the
right-hand side is the standard normal-vs-shift-and-scale Jacobian term. Under
`use_centered_alpha == 1` both kernels use the same centered density on
`log_a_market_pool`/`log_a_market_derived` with only a fixed unit-Jacobian
shift, so the gap is exactly 0.

Because the normalized target on every economic quantity (`a_market`,
`s0_prob`, `alpha_structural_market`, `sigma_log_a`, `mu_log_a`, ...) is
unchanged, `.antitrust_bayes_model_spec_version` (`R/helpers.R:2114`) is **not**
bumped — the same precedent as the earlier joint-s0 rotation
(`doc/s0_joint_rotation_diagnosis.md`: "Statistical model specification
version 10 remains appropriate because this adds sampler coordinates without
changing the normalized target or economic model").

## Validation

### Part A — analytic identity (pure R, `doc/pooled_cournot_alpha_ncp_native_validation.R`)

Both Stan formulas (OLD centered, NEW noncentered, NEW centered) were
transcribed exactly in R and checked against the identity above across a
5x3 grid of `sigma_log_a in {0.6, 0.15, 0.02, 0.005, 0.0005}` x
`mu_log_a in {log(0.84), -0.2, 0.35}`, at a matched economic state
(`log_a = mu + sigma*r` for a fixed `r` vector). All 45 checks passed,
including the `matched-state a_market` identity (<=1e-12), the exact log-prob
gap (<=1e-9 relative), and the `use_centered_alpha == 1` exact-match check
(<=1e-9), down to `sigma_log_a = 0.0005` — well past the funnel boundary that
caused the original divergences.

### Part B — compiled kernels (same script, real CmdStanR)

Both working-tree kernels (`logit_supply`, `logit_stackelberg`) were compiled
under CmdStan 2.39.0 and run with `fixed_param = TRUE` at an explicit
constrained state. All 8 checks passed (<=1e-9 tolerance, `sig_figs = 18`):
the compiled kernel's `log_a_market_derived`/`a_market` match the R formula
under both `use_centered_alpha` values, and a non-Cournot pooled state (the
untouched `is_cournot == 0` branch) is unaffected by this change, for both
kernels.

### Posterior-equivalence HMC smoke (`doc/pooled_cournot_alpha_ncp_hmc_smoke.R`)

Six-market, one-year, pooled-Cournot fixture (`logit_supply`), 4 chains x
500/500, `adapt_delta = 0.95`, same seed for both kernels:

| Kernel | Divergences | Treedepth hits |
|---|---:|---:|
| OLD (centered `log_kappa_market`) | 67 / 2000 (3.0%) | 0 |
| NEW (noncentered `r_market_a_raw`) | 1 / 2000 (0.05%) | 0 |

67x fewer divergences at an equal HMC budget. Economic posterior means agree
within combined MCSE for every quantity except `sigma_log_a_vec[1]` itself
(z = 2.85, still under the conventional z = 3 flag; expected, since that is
exactly the funnel-boundary parameter and the OLD run's divergent draws bias
its own posterior). `a_market`, `s0_prob`, `mu_log_a`: all z < 1.9. Full table
in `pooled_cournot_alpha_ncp_hmc_smoke.csv`.

## Side effect

Pooled Cournot now honors `use_centered_alpha` like every other pooled
conduct model (previously that flag had no effect on Cournot). The
mergerBayes production runner passes `use_centered_alpha = FALSE`
(non-centered), which is the branch validated above.

## Full test suite

`devtools::test()` clean (see commit); `tests/testthat/test-market-s0-rotation-api.R`
updated to assert `log_kappa_market` is gone from both Logit kernels and that
`log_a_market_derived = mu_log_a + sigma_log_a * r_market_a_raw` is present in
both.
