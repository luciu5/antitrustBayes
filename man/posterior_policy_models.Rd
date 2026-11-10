\name{posterior_policy_models}
\alias{posterior_policy_models}
\title{Realize reusable structural models from a posterior policy draw plan}
\usage{
posterior_policy_models(
  fit, draw_plan = NULL, n_draws = NULL, seed = NULL,
  market = NULL, year = NULL,
  demand_anchoring = NULL,
  insideSize = NULL, priceOutside = NULL, labels = NULL,
  allow_experimental = FALSE
)
}
\arguments{
\item{fit}{A \code{BayesFit} for a supported output Logit/CES or input Logit model.}
\item{draw_plan}{A plan returned by \code{posterior_policy_draws()}. If omitted,
a new plan is constructed once.}
\item{n_draws}{Number of draws for a new plan; must agree with a supplied plan.}
\item{seed}{Optional seed for a new plan. Cannot accompany a supplied plan.}
\item{market, year}{Simulation-unit selectors. Panel fits require an explicit year selection.}
\item{demand_anchoring}{If omitted (\code{NULL}), estimated fits use
\code{"latent_quality"}; supplied \code{BayesSpecifiedFit} primitives retain
\code{"observed"} calibration. Explicit \code{"observed"} reproduces the
previous conditional observed-share reconstruction. For demand-shock fits,
latent anchoring carries jointly recovered persistent quality plus period
demand through the structural baseline; independent measurement noise remains
in the observation model.}
\item{insideSize, priceOutside, labels}{Optional unit overrides using \code{simulate()} conventions.}
\item{allow_experimental}{Explicitly allow exploratory Pathfinder draws.}
}
\value{
A \code{bayes_policy_models} list containing \code{models} nested by unit and
selected draw position, \code{draw_indices}, \code{units}, \code{draw_plan},
\code{demand_anchoring}, and a \code{failures} data frame. Failed models retain
their positions as \code{NULL}.
}
\description{
Constructs baseline structural models from reusable posterior draw plans for
deterministic policy simulation without fitting again or simulating an
unchanged counterfactual.
}
\details{
Support includes simultaneous Bertrand, Cournot, MonCom, second-score auction
and bargaining, and native exact fringe/Stackelberg output Logit/CES and input
Logit models. Ordinary models are
\code{AntitrustFit} objects; exact role models are native \pkg{coordination}
objects. Quality and roles are not resampled during realization. Latent
anchoring preserves the original structural cost primitive and may require
solving the genuine latent baseline. No tariff or policy equation is evaluated.
Attributes retain posterior index, local quality units, realized k and plan
fingerprint. Failures identify unit, posterior index and construction stage.
Second-score Logit auctions use their cost-based value
index follows \pkg{antitrust}: winning probabilities respond to effective
costs rather than bid prices. Latent reconstruction preserves costs, sets the
auction-normalized log share odds, and holds the implied value index fixed
through later policies, with the correct input/output coefficient sign.
Output CES auctions and supported bargaining retain their native equations.
Smooth-fringe models remain unsupported by this helper. In particular,
smooth-cutoff Logit auctions do not have an implemented latent cost-based
value-index reconstruction; no observed-share fallback is used.
}
\seealso{\code{\link{posterior_policy_draws}}, \code{\link{posterior_quality}}}
