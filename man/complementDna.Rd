% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dna.R
\name{complementDna}
\alias{complementDna}
\title{Complement DNA}
\usage{
complementDna(dnas, ambigs = TRUE)
}
\arguments{
\item{dnas}{vector of sequences}

\item{ambigs}{if TRUE complement ambiguous bases}
}
\value{
vector with the DNA sequences complemented
}
\description{
Complement DNA
}
\examples{
complementDna(c('CTAG','ATCCAC'))
complementDna(c('CT[AC]G','ATNRY'))
}
