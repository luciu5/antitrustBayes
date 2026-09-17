# Statistical toolkit follow-up: bounded priorities

Status: requested work, specifications and experiments in progress. This note
records separation between an exact coordinate experiment and substantive
extensions. No likelihood or simulation default is changed by this plan.

1. Complete and validate the market/outside-share rotation without changing
   priors, pooling, likelihoods, conduct, or counterfactual definitions.
2. Derive and implement an explicit latent-quality counterfactual mode. Audit
   utility/outside normalization and cost recovery first. Draw quality once
   per common posterior draw and reuse it across pre/post interventions and
   panel units. Compare against the observed-share anchor with paired draws.
3. Execute the missing collapsed Gaussian-share/fixed-t5-margin benchmark on
   the seven frozen original fixtures. Keep rho off and all original priors
   and HMC settings. This is a separately labeled margin-likelihood experiment
   on the legacy-coordinate baseline; it does not establish rotation effects.
4. Formalize a model-universe layer around `bayes_compare()`. Represent shared
   endpoints once and interior conduct families explicitly. Restricted-family
   evidence must use normalized prior/posterior restriction masses, rather
   than dropping endpoint rows from a family comparison after fitting.
5. Add economically interpretable prior-predictive diagnostics, then resume
   exact-k recovery and PNB using configurations supported by executed evidence.

Hierarchical k, market-specific conduct, firm-market quality, nested outside
hierarchies, explicit cost equations, compositional share robustness, and
integrated out-of-market prediction remain separate research extensions.
The user explicitly deferred building market-specific conduct now. These
extensions require individual mathematical specifications and validation;
they are not silently added during the rotation experiment.

An economic nesting statement alone does not prove equality of marginal
likelihoods. Endpoint consolidation also requires the same observed outcomes,
parameter priors, likelihoods, normalization, and nuisance integration.
Different demand families or statistical specifications are not automatically
the same Bayesian hypothesis merely because their conduct labels agree.
