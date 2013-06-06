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
	Returns the log-likelihood of the data with the fitted changepoints, two values are returned, the raw scaled negative log-likelihood and the scaled negative log-likelihood + penalty.  Only valid for cpttype="mean","variance" or "mean and variance".
}
}}
\examples{
set.seed(1)
x=c(rnorm(50,0,1),rnorm(50,0,10),rnorm(50,0,5),rnorm(50,0,1))
out=cpt.var(x,penalty="Manual",pen.value="2*log(n)",method="BinSeg",Q=5)
likelihood(out) # returns the raw scaled negative likelihood (783.9144) and the scaled negative
#likelihood + penalty (805.1076)
}

\keyword{methods}
\keyword{cpt}
