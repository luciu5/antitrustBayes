# Pilot decision

The bounded Logit Bertrand run supports a staged CmdStanR migration. On the
same 36-observation synthetic input and four-chain, 300-warmup/300-sampling
budget, the public RStan `bayes_est()` call took 48.3 seconds; CmdStanR sampling
took 18.1 seconds, with a separate cold model compile of 85.7 seconds. The
timers have different scopes: RStan includes package data preparation, whereas
CmdStanR's figure is sampling only. This is a promising speed signal, not a
production benchmark or a claim about all models.

The two compiled programs agreed on the normalized log target at three fixed
unconstrained parameter points to a maximum absolute difference of
`8.5e-14`. Selected posterior summaries and observation-level LOO/WAIC were
close. The short chains had 4 RStan and 10 CmdStanR divergences, and minimum
bulk ESS of 367 and 224, respectively; these fits are too short to establish
inferential quality or a stable ESS-per-second advantage. The saved CmdStanR
fit reloaded and produced posterior summaries in a fresh R session. CmdStanR's
bridge check was unavailable because the local CmdStan build lacks SUNDIALS
`-fPIC` flags and this sandbox cannot write to its installation for the
automatic rebuild of in-process model methods. `bridgesampling` provides a
`CmdStanMCMC` method, so this is an environment-specific check rather than
evidence of backend incompatibility. RStan bridge sampling returned a log
marginal likelihood of 53.08.

The CES model compiled, but the comparative run was stopped before sampling
finished. No CES speed or parity claim follows from this pilot.

Recommendation: make CmdStanR an opt-in backend behind the existing package
fit contract, starting with simultaneous Bertrand Logit and CES. Keep RStan
available for existing fits. During implementation, run one bridge check in a
writable CmdStan installation configured for model methods, then check the
remaining kernels and package distribution path before changing the default.
The current pilot is enough to justify that focused implementation; repeating
broad benchmarks now has low value.
# Historical pilot (superseded implementation recommendation)

This document records the pre-migration pilot and its original recommendation.
The current package implementation uses CmdStanR for HMC, with Pathfinder as
an explicit exploratory option. The pilot's short chains and unequal timing
scopes do not certify a speed advantage or sampler quality. Its raw outputs
are retained as engineering evidence.
