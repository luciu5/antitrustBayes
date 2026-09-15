# Inherited correlation-prior constant

This code audit precedes all new-model fits. The unchanged Logit/CES Stan
models set rho=.95*tanh(x), x=rho_raw, and add the prior kernel

    (eta-1)*log(1-rho^2) + log(1-tanh(x)^2),

where eta=prior_lkj. This is integrable and defines a proper prior after
normalization, but the manually written target omits its normalization
constant. Its integral on the raw parameter scale is

    I_eta = (1/.95) * integral[-.95,.95] (1-rho^2)^(eta-1) d rho.

For eta=2 this equals 2*(1-.95^2/3)=1.3983333333, so log(I_eta)=.3352811.
Direct numerical integration in both coordinates agrees within 1e-11 for
eta=1,2,4. This is an algebraic and numerical finding about inherited code,
not a posterior-estimation result.

This note records the initial version-6 correlation-only audit. The subsequent
full audit identifies other omitted constants from distribution statements and
positive-support priors. Matching correlation priors alone cannot justify raw
cross-structure model odds.

The final implementation follows the user's active-development instruction:
all four Stan programs retain normalized prior PDFs directly, subtract the
bounded-correlation normalizer above, and normalize each positive Normal's
support. Specification 7 requires earlier fits to be refitted. There is no
external bridge-offset path. The derivation and executed historical density
checks are preserved in
[evidence_normalization_contract.md](evidence_normalization_contract.md).
New compiled targets must equal independent full proper priors plus likelihood
and automatic constraint Jacobians; this is a separate required numerical gate.
