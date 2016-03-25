% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/abundance.R
\name{kullback}
\alias{kullback}
\title{Calculate Kullback–Leibler divergence between two probability distributions}
\usage{
kullback(x, y, base = 2, standardize = TRUE)
}
\arguments{
\item{x}{Vector of counts or proportions of N elements}

\item{y}{Second vector of counts or proportions of N elements}

\item{base}{Base of logarithm}

\item{standardize}{If TRUE divide x and y by sum(x) and sum(y)}
}
\value{
Kullback-Leibler divergence
}
\description{
Calculate Kullback–Leibler divergence between two probability distributions
}
