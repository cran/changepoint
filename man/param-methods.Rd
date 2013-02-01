\name{param-methods}
\docType{methods}
\alias{param-methods}
\alias{param,cpt-method}
\alias{param,cpt.reg-method}
\title{ ~~ Methods for Function param  ~~}
\description{
 ~~ Methods for function \code{param}  ~~
}
\section{Methods}{
\describe{

\item{\code{signature(object = "cpt")}}{
	Replaces the \code{param.value} slot in object with the parameter estimates that are appropriate for the changepoint type (cpttype slot).  If the Gamma test statistic is used then the shape parameter is required as a variable.
}
\item{\code{signature(object = "cpt.reg")}}{
	Replaces the \code{param.value} slot in object with the parameter estimates that are appropriate for the changepoint test statistic (test.stat slot).
}
}}
\keyword{methods}
\keyword{cpt}
