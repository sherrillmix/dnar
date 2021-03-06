% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dna.R
\name{trimNs}
\alias{trimNs}
\title{Trim ambiguous sequence from ends of reads}
\usage{
trimNs(seqs, nonNStretch = c(10, 10), nChars = c("N"))
}
\arguments{
\item{seqs}{sequences to be trimmed}

\item{nonNStretch}{a one or two element (specifying start, end) to delete until finding a stretch of nonNs greater this length (0 means no trimming)}

\item{nChars}{a vector of characters to count as Ns}
}
\value{
Vector sequences with ends trimmed
}
\description{
Trim sequences to the first occurrence of a given length of unambiguous bases, e.g. to trim noisy starts and ends from Sanger sequences. Note that sequences shorter than nonNStretch will be completely deleted since they do not contain any stretch of unambiguous bases to pass the threshold.
}
\examples{
trimNs('NNNANNACTNGAGANNNNAAN',3)
}
