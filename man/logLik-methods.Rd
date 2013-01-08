\name{logLik-methods}
\docType{methods}
\alias{logLik-methods}
\alias{logLik,cpt-method}
\title{ ~~ Methods for Function logLik  ~~}
\description{
 ~~ Methods for function \code{logLik}  ~~
}
\section{Methods}{
\describe{

\item{\code{signature(object = "cpt")}}{
	Returns the log-likelihood of the data with the fitted changepoints, two values are returned, the raw log-likelihood and the log-likelihood + penalty.  Only valid for cpttype="mean","variance" or "mean and variance".
}
}}
\keyword{methods}
\keyword{cpt}
