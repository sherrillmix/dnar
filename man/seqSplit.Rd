% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dna.R
\name{seqSplit}
\alias{seqSplit}
\title{Convenience function for splitting a bunch of sequences into a matrix}
\usage{
seqSplit(..., fill = ".")
}
\arguments{
\item{...}{Various sequences to split into a matrix}

\item{fill}{A character to pad ends of sequences. An error is generated if any sequence differs in length if NULL}
}
\value{
A matrix of single characters each row corresponding to a read
}
\description{
Convenience function for splitting a bunch of sequences into a matrix
}
\examples{
seqSplit('ACACA','ACA')
seqSplit('ACACA','ACA',fill='-')
seqSplit(c('ACACA','ACA'),'TCACA')
}
