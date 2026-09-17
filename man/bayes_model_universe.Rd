% Generated manually for the model-universe API.
\name{bayes_model_universe}
\alias{bayes_model_universe}
\title{Declare a non-overlapping Bayesian model universe}
\usage{
bayes_model_universe(hypotheses, prior_weights)
}
\arguments{
  \item{hypotheses}{A named list. A character value names a full source fit;
  a list with \code{fit} and \code{k = "interior"} names a strict interior
  exact ordered-k family. Explicit integer k values or named market-to-k lists
  are also accepted for exact families.}
  \item{prior_weights}{Finite non-negative probabilities, one per named
  hypothesis, summing to one. These are the scientific prior weights and are
  never silently redistributed.}
}
\value{An object of class \code{bayes_model_universe}.}
\description{
Creates the explicit model-universe descriptor consumed by
\code{bayes_compare(..., universe = descriptor)}.
}
\details{
The descriptor represents shared economic endpoints exactly once. Endpoint
states are represented by the canonical simultaneous MC, B, and C hypotheses;
exact fringe and Stackelberg families must retain strictly interior states.
Evidence for a restricted family is corrected by its posterior restriction
mass divided by its recorded prior restriction mass.
}
\examples{
\donttest{
u <- bayes_model_universe(
  list(MC = "mc_fit", B = "b_fit",
       BF = list(fit = "bf_fit", k = "interior")),
  prior_weights = c(MC = .3, B = .4, BF = .3)
)
}
}
\keyword{models}
