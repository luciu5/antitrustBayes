# Option B research target, 2026-10-22

This is an isolated prototype. No maintained source, priors, package default,
production experiment or saved fit is changed. Experimental fits do not carry a
production specification number. Scope is one market, one period, MC/B/C, with
Gaussian errors, sampled symmetric quality and rho fixed to zero. Eight fits (four matched A/B pairs) compare this target against the current free-intercept target under exactly
the same rho/quality choices; historical rho-active/collapsed results are not a
one-change control.

Let Q be J by (J-1), Q'Q=I and QQ'=P=I-11'/J. Quality is b=Q theta,
theta|tau ~ N(0,tau² I). Because Q is an isometry onto the zero-sum subspace,
its volume Jacobian is one; the normalized quality density is the J-1 dimensional
normal on theta. No extra sqrt(J) factor or redundant mean coordinate is added.
This matches the current package's Helmert coordinates and normalized prior.
The positive quality SD retains its current truncated N(1,1) prior; all active
PNB alpha, outside-share, residual-SD and conduct priors are copied unchanged.

The demand observation is z=Q'log(r). With x standardized price and market_sign
+1 for banks/input markets and -1 for output markets,

    z | theta,alpha,sigma_s ~ N(theta + market_sign*alpha_std*Q'x,
                               sigma_s² I_(J-1)).

The target retains the full -((J-1)/2)log(2*pi) -(J-1)log(sigma_s)
normalization. It does not insert the J centered components into an independent
J-dimensional normal. Rotating Q by an orthogonal matrix changes neither this
density nor the isotropic quality prior. The data density uses Euclidean dz.
For a density with respect to simplex coordinates dr_1 ... dr_(J-1), multiply
by 1/(sqrt(J)*prod(r)); that observed-data-only factor is common across conduct
models within Option B. Raw marginal-likelihood levels across A and B do not
share an observation reference measure and are not treated as Bayes factors.

There is no sampled inside-market intercept in B. Define v=b+market_sign*alpha_std*x.
The generated a=-logit(s0)-logsumexp(v) gives absolute utilities a+v and
outside utility zero. Then outside share is exactly s0 and conditional inside
shares are softmax(v). This is a calibration identity for the latent composition.
With measurement error the latent composition need not equal observed r. An
observed-share-anchored policy calculation uses observed r and recalibrates
utilities, matching the package's existing observed anchoring contract; this is
reported distinctly from latent-quality predictions.

Supply is unchanged: the exact maintained logit_inv_markup() function is copied
verbatim, including input price sign, physical price units, inverse-normalized
margin response, observed-margin mask and the existing 1e-9 stabilization. PNB
has 39 products and 32 observed margins. Missing margins are retained as missing,
not imputed. MC inverse markup is alpha, independent of s0. Thus at rho=0 all
other parameters/observations are independent of the outside-share logit in B:
its posterior exactly equals its proper prior. Cournot's existing log-kappa
coordinate is preserved: mu_log_alpha=log_kappa-log(s0), with conditional
Jacobian d(mu_log_alpha)/d(log_kappa)=1. Its normal prior remains a prior on
log-alpha, not on kappa. Bertrand/Cournot margins can identify s0.

For a future correlated implementation, let S select observed-margin rows.
If raw share and matched margin errors have covariance rho*sigma_s*sigma_m*S,
the transformed joint covariance is

    V = [ sigma_s² I                rho*sigma_s*sigma_m*Q'S ]
        [ rho*sigma_s*sigma_m*S'Q   sigma_m² I              ].

Conditioning margins on contrasts gives covariance
sigma_m²[I-rho² S' P S], generally dense. A rowwise conditional normal with
centered residuals and old variance sigma_m²(1-rho²) is incorrect. When theta
is analytically integrated, replace the upper-left block by
(sigma_s²+tau²)I; the cross block is unchanged. The full normalized multivariate
normal and its Schur factorization agree in deterministic tests with missing
margins at J=3/7/39; scalar numerical convolution separately checks the
rho=0 quality integral. These mathematical tests are not a production native
joint/collapsed implementation or validation for panels/repeated products.

Normalized bridge evidence uses the actual compiled experimental Stan target
through bridgesampling's CmdStanMCMC route: unconstrained parameters,
log_prob(jacobian=TRUE), normalized prior/likelihood constants and the standard
normal bridge proposal. Positive parameters use Stan's transform Jacobians;
logit outside share is already a real parameter with a normal prior. Every
requested repetition is retained. Only within-B conduct evidence is meaningful;
no seven-model BMA is claimed from a three-model prototype. Any three-model probabilities explicitly condition the original PNB prior on MC/B/C, yielding .6/.2/.2.

Runtime controls are sealed in manifest.rds: four chains, 1000 warmup plus
1000 retained per chain, .999 target acceptance, diagonal metric, depth 13,
fixed seed and at most four sampling cores. A single cumulative <=4h deadline
includes execution/compilation/evidence; controller kills the process tree at
the deadline. Completed and failed attempts are retained without retuning.
