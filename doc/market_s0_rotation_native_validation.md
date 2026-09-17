# Independent market-s0 rotation gate

Date: 2026-09-16. All native calls use `rstan::sampling(..., algorithm="Fixed_param", iter=1, warmup=0)` followed by `rstan::log_prob`/`grad_log_prob`; no HMC or refitting was run. Original control is `/tmp/rotation-baseline-Rlib`; fresh rotation build is `/tmp/rotation-fresh-Rlib`. Each namespace ran in a separate R process.

## Results

- M=2,T=2, 176 supported cases per namespace: original-vs-rotated normalized `adjust_transform=TRUE` target max `2.842171e-14`; rotated AD-vs-central-FD max `6.864215e-09`; legacy-to-rotated coordinate chain max `4.156701e-09`; baseline-vs-rotated GQ max `1.776357e-15`; zero setup failures. Covers Logit/CES, Bertrand/Cournot, simultaneous/fringe/Stackelberg, none/market-year/market s0, centered/NC, sampled/collapsed quality, alpha common/pool, supported rho, missing margins, and exact BF/CF/BS/CS k cases.
- M=1,T=2 no-op, 256 cases: zero setup failures; baseline-vs-experimental target/GQ max `0`; experimental FD max `4.326377e-09`; identity legacy/rotated chain max `0`. This includes market-pool NC zero-width declarations.
- M=1,T=1 endpoint, 8 all-kernel cases: target/GQ max `0`, FD max `1.349636e-10`, chain max `0`, zero failures.
- M=2,T=1 endpoint, 8 active-rotation cases: target/GQ max `0`, FD max `9.919405e-10`, chain max `3.553685e-10`, zero failures.
- Student-t sampled-quality Logit B/C subset, 4 cases: target max `0`, GQ max `4.440892e-16`, FD max `1.439506e-09`, chain max `2.910929e-10`, zero failures.
- Centered alpha market-pool subset, 4 cases: target/GQ max `0`, FD max `1.661506e-09`, chain max `6.280652e-10`, zero failures.
- Input-side Logit market-pool NC B/C subset, 12 cases: target max `0`, common GQ max `0`, FD max `5.413214e-09`, chain max `2.830601e-10`, zero failures. Covers sampled/collapsed quality, alpha common/pool, and supported rho.

## Archived fit bridge reload

Fixture: `doc/statistical_toolkit_runs/20260916/hybrid/fits/logit_BF_noisy_271828__D.rds`; its archived native model name is `logit_stackelberg` (exact BF). Ordinary `readRDS`, class/metadata, and `get_inits` work. Direct `rstan::log_prob` and direct `bridge_sampler` on the deserialized `stanfit` fail in both namespaces with `the model object is not created or not valid`, demonstrating the serialized native pointer is not reusable.

A live Fixed_param context made from the archived data works in both namespaces (`log_prob=-203.8624`). The original control context plus archived samples bridges successfully (`logml=23.6902254754404`). The fresh rotated context with private `market_s0_parameterization=1L` fails when passed archived samples because current RStan exposes zero-width `c_s0`/`c_s0_raw` declarations absent from the archived sample schema. Adding those empty fields manually lets `unconstrain_pars` proceed, so a dedicated adapter could be implemented, but the current public bridge route does not provide it. The package must require a matching native context/schema adapter or refit/reconstruct the native instance for old serialized bridge work; it should not claim transparent old-fit bridge compatibility from the model-spec version alone.
