% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dna.R
\name{index2range}
\alias{index2range}
\title{Find covered ranges with numerical index data}
\usage{
index2range(index)
}
\arguments{
\item{index}{numeric indices indicating positions of interest}
}
\value{
data.frame with rows for each contiguous region of interest with columns start and end of regions
}
\description{
Find covered ranges with numerical index data
}
\examples{
index2range(c(1:100,300:450))
index2range(c(1:100,2:110))
}
