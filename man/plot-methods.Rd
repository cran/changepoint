\name{plot-methods}
\docType{methods}
\alias{plot-methods}
\alias{plot,ANY-method}
\alias{plot,cpt-method}
\alias{plot,cpt.reg-method}
\title{ ~~ Methods for Function plot in Package `graphics' ~~}
\description{
 ~~ Methods for function \code{plot} in Package `graphics' ~~
}
\section{Methods}{
\describe{

\item{\code{signature(x = "ANY")}}{
	Generic plot function, see graphics package description using ?plot
}

\item{\code{signature(x = "cpt")}}{
	Plots the data and identifies the changepoints using vertical lines (change in variance), horizontal lines (change in mean).  Optional arguments to control the lines: \code{cpt.col} equivilent to \code{col} to change the colour of the changepoint line; \code{cpt.width} equivilent to \code{lwd} to change the width of the changepoint line; \code{cpt.style} equivilent to \code{lty} to change the style of the line.
}
\item{\code{signature(x = "cpt.reg")}}{
	Plotting is only valid for one regressor.  Plots the regressor against the response and identifies the changepoints using horizontal lines.  Optional arguments to control the lines: \code{cpt.col} equivilent to \code{col} to change the colour of the changepoint line; \code{cpt.width} equivilent to \code{lwd} to change the width of the changepoint line; \code{cpt.style} equivilent to \code{lty} to change the style of the line.
}

}}
\keyword{methods}
\keyword{plot}
\keyword{cpt}
