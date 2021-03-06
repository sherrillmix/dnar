% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dna.R
\name{cigarToBlock}
\alias{cigarToBlock}
\title{Convert cigar and starts to qStarts, tStarts, blockSizes as in blat}
\usage{
cigarToBlock(cigars, starts, startEnds = FALSE)
}
\arguments{
\item{cigars}{vector of SAM cigar strings}

\item{starts}{vector of starting positions in target}

\item{startEnds}{if FALSE single line with comma separated starts, if TRUE data.frame with single start, end and id column}
}
\value{
dataframe with qStarts,tStarts,sizes if !startEnds or dataframe with starts, ends and ids if startEnds
}
\description{
Convert cigar and starts to qStarts, tStarts, blockSizes as in blat
}
\examples{
cigars<-c('100M','10H10M1000H','10M10I20M10D2M','10S10M10H','100M10000N100M')
starts<-1:5
cigarToBlock(cigars,starts)
cigarToBlock(cigars,starts,startEnds=TRUE)
}
