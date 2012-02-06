\name{PELT.var.norm}
\alias{PELT.var.norm}
\title{
Multiple Changes in Variance using PELT pruned method - Normal Data
}
\description{
Calculates the optimal positioning and number of changepoints for Normal data using PELT pruned method.
}
\usage{
PELT.var.norm(data, pen=0, know.mean=FALSE, mu=NA)
}
\arguments{
  \item{data}{
	A vector containing the data within which you wish to find changepoints.
}
  \item{pen}{
	Numeric value of the linear penalty function.  This value is used in the decision for each individual changepoint so that in total the penalty is k*pen where k is the optimal number of changepoints detected.
}
  \item{know.mean}{
	Logical, if TRUE then the mean is assumed known and mu is taken as its value.  If FALSE, and mu=-1000 (default value) then the mean is estimated via maximum likelihood.  If FALSE and the value of mu is supplied, mu is not estimated but is counted as an estimated parameter for decisions.
}
  \item{mu}{
	Numerical value of the true mean of the data.  Either single value or vector of length nrow(data).  If data is a matrix and mu is a single value, the same mean is used for each row.
}
}
\details{
	This function is used to find a multiple changes in variance for data that is assumed to be normally distributed.  The value returned is the result of testing H0:existing number of changepoints against H1: one extra changepoint using the log of the likelihood ratio statistic coupled with the penalty supplied.  The PELT method keeps track of the optimal number and location of changepoints as it passes through the data.
}
\value{
	A vector of the changepoint locations is returned:
	\item{cpt}{Vector containing the changepoint locations for the penalty supplied.  This always ends with n.}
}
\references{
Change in variance: Chen, J. and Gupta, A. K. (2000) \emph{Parametric statistical change point analysis}, Birkhauser

PELT Algorithm: Killick, R. and Fearnhead, P. and Eckley, I.A. (2011) An exact linear time search algorithm for multiple changepoint detection, \emph{Submitted}
}
\author{
Rebecca Killick
}


\seealso{
\code{\link{PELT.mean.norm}},\code{\link{PELT.meanvar.norm}},\code{\link{cpt.var}},\code{\link{binseg.var.norm}},\code{\link{multiple.var.norm}},\code{\link{single.var.norm}},\code{\link{segneigh.var.norm}}
}
\examples{
# Example of multiple changes in variance at 50,100,150 in simulated normal data
set.seed(1)
x=c(rnorm(50,0,1),rnorm(50,0,10),rnorm(50,0,5),rnorm(50,0,1))
PELT.var.norm(x,pen=2*log(200)) # returns c(50,99,150,200)

# Example no change in variance
set.seed(10)
x=rnorm(200,0,1)
PELT.var.norm(x,pen=2*log(200)) # returns 200 to show no change in variance has been found
}

\keyword{methods}
\keyword{univar}
\keyword{models}
\keyword{ts}
