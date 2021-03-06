% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{insetScale}
\alias{insetScale}
\title{Add an inset color scale to a plot}
\usage{
insetScale(
  breaks,
  col,
  insetPos = c(0.025, 0.015, 0.04, 0.25),
  main = "",
  offset = 0.001,
  at = NULL,
  labels = NULL,
  cex = 1,
  labXOffset = 0,
  labYOffset = 0
)
}
\arguments{
\item{breaks}{a set of finite numeric breakpoints for the colours: must have one more breakpoint than colour and be in increasing order.}

\item{col}{a list of colors such as that generated by 'rainbow', 'heat.colors', 'topo.colors', 'terrain.colors' or similar functions.}

\item{insetPos}{a four element numeric vector giving the position of the bottom, left, top and right edges of the scale in nfc coordinates (0 = bottom/left of figure, 1 = top/right of figure)}

\item{main}{a character string giving label for color scale}

\item{offset}{add a small offset to rectangles to suppress pdf viewers showing a gap at exact border between rectangles}

\item{at}{positions for labels on the axis}

\item{labels}{labels for the axis}

\item{cex}{character size for axis and labels}

\item{labXOffset}{a vector of the same length as labels to be used to fine tune label positioning. Values given proportional to plot size}

\item{labYOffset}{a vector of the same length as labels to be used to fine tune label positioning. Values given proportional to plot size}
}
\description{
Adds an inset color scale to a plot without disrupting future plotting. Especially useful with \code{heatmap} or \code{image} plots.
}
\examples{
dists<-as.matrix(dist(sort(runif(40))))
breaks<-seq(0,1,.01)
cols<-rev(heat.colors(length(breaks)-1))
heatmap(dists,col=cols)
insetScale(breaks,cols,main='Distance')
}
\seealso{
\code{\link[stats]{heatmap}}, \code{\link[graphics]{image}}
}
