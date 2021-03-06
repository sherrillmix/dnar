% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dna.R
\name{pwm}
\alias{pwm}
\title{Generate a PWM for a set of sequences}
\usage{
pwm(seqs, chars = c("C", "G", "T", "A"), priors = table(chars) - 1)
}
\arguments{
\item{seqs}{vector of strings, sequences to form a position weight matrix from (all same length)}

\item{chars}{allowed characters}

\item{priors}{additional counts to add to each column with (default is a lazy way to get a named 0 vector)}
}
\value{
position weight matrix with length(chars) rows and nchar(seqs) columns
}
\description{
Generate a PWM for a set of sequences
}
\examples{
pwm(c('ACTG','ACTT','ACTT','ACGT'))
pwm(c('ACTG','ACTT','ACTT','ACGT'),priors=c('C'=1,'G'=1,'T'=1,'A'=1))
}
