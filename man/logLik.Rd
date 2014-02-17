\name{logLik}
\alias{logLik}
\title{
Generic Function - logLik
}
\description{
Generic function to calculate the log likelihood
}
\usage{
logLik(object)
}
\arguments{
  \item{object}{
	Depending on the class of \code{object} depends on the method used to calculate the log likelihood (and if one exists)
}
}
\details{
	Generic Function to calculate the log likelihood.
}
\value{
	Depends on the class of \code{object}, see individual methods
}
\author{
Rebecca Killick
}


\seealso{
\code{\link{logLik-methods}}, \code{\link{cpt.mean}},\code{\link{cpt.var}},\code{\link{cpt.meanvar}}
}
\examples{
set.seed(1)
x=c(rnorm(50,0,1),rnorm(50,0,10),rnorm(50,0,5),rnorm(50,0,1))
out=cpt.var(x,penalty="Manual",pen.value="2*log(n)",method="BinSeg",Q=5)
logLik(out) # returns the raw scaled negative likelihood and the scaled negative
#likelihood + penalty
}

\keyword{methods}
\keyword{cpt}
