% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/abundance.R
\name{jensenShannon}
\alias{jensenShannon}
\title{Calculate Jensen Shannon divergence between two probability distributions}
\usage{
jensenShannon(x, y, base = 2)
}
\arguments{
\item{x}{Vector of counts or proportions of N elements}

\item{y}{Second vector of counts or proportions of N elements}

\item{base}{Base of logarithm}
}
\value{
Jensen-Shannon divergence
}
\description{
Calculate Jensen Shannon divergence between two probability distributions
}
\examples{
jensenShannon(1:4,4:1)
}
