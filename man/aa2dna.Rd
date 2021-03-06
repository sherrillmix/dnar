% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dna.R
\name{aa2dna}
\alias{aa2dna}
\title{Convert amino acids to dna regex}
\usage{
aa2dna(aas)
}
\arguments{
\item{aas}{A vector of amino acids to turn to dna regex}
}
\value{
A vector with a regular expression for each element of aas
}
\description{
Convert amino acids to dna regex
}
\examples{
aa2dna(c('MATX','MGCATKRVX'))
}
