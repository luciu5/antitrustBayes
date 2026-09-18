\name{posterior_policy_models}
\alias{posterior_policy_models}
\title{Realize reusable structural models from a posterior policy draw plan}
\usage{
posterior_policy_models(
  fit, draw_plan = NULL, n_draws = NULL, seed = NULL,
  market = NULL, year = NULL,
  demand_anchoring = c("observed", "latent_quality"),
  insideSize = NULL, priceOutside = NULL, labels = NULL
)
}
\arguments{
\item{fit}{A \code{BayesFit} for a supported output-market model.}
\item{draw_plan}{A plan returned by \code{posterior_policy_draws()}. If omitted,
a new plan is constructed once.}
\item{n_draws}{Number of draws for a new plan; must agree with a supplied plan.}
\item{seed}{Optional seed for a new plan. Cannot accompany a supplied plan.}
\item{market, year}{Simulation-unit selectors. Panel fits require an explicit year selection.}
\item{demand_anchoring}{Observed-share or conditional latent-quality anchoring.}
\item{insideSize, priceOutside, labels}{Optional unit overrides using \code{simulate()} conventions.}
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
Initial support includes simultaneous Bertrand, Cournot and MonCom and native
exact fringe/Stackelberg output Logit/CES models. Ordinary models are
\code{AntitrustFit} objects; exact role models are native \pkg{coordination}
objects. Quality and roles are not resampled during realization. Latent
anchoring preserves the original structural cost primitive and may require
solving the genuine latent baseline. No tariff or policy equation is evaluated.
Attributes retain posterior index, local quality units, realized k and plan
fingerprint. Failures identify unit, posterior index and construction stage.
Second-score output Logit auctions are also supported. Their cost-based value
index follows \pkg{antitrust}: winning probabilities respond to effective
costs rather than bid prices. Latent reconstruction preserves costs, sets the
auction-normalized log share odds, and holds the implied value index fixed
through later policies. CES auctions remain unsupported by this helper.
}
\seealso{\code{\link{posterior_policy_draws}}, \code{\link{posterior_quality}}}
