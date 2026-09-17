% Generated manually for the model-universe API.
\name{bayes_dimensions}
\alias{bayes_dimensions}
\title{Aggregate model-universe posterior probabilities}
\usage{
bayes_dimensions(comparison,
  dimension = c("conduct", "structure", "strategic_core"))
}
\arguments{
  \item{comparison}{A \code{bayes_universe_comparison} returned by
  \code{bayes_compare} with a model-universe descriptor.}
  \item{dimension}{One or more dimensions: \code{"conduct"},
  \code{"structure"}, or \code{"strategic_core"}.}
}
\value{A data frame containing exact sums of atomic posterior probabilities
and the corresponding scientific prior probabilities.}
\description{Aggregates atomic model-universe probabilities without changing
the prior or renormalizing over an incomplete evidence set.}
\details{The strategic-core value \code{interior} denotes a strict k less
than n state in every market of the exact family, a global all-market event.
Auction remains its own conduct category and MC remains its own category.}
\keyword{models}
