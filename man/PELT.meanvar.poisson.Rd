\name{PELT.meanvar.poisson}
\alias{PELT.meanvar.poisson}
\title{
Multiple Changes in Mean and Variance using PELT pruned method - Poisson Data
}
\description{
Calculates the optimal positioning and number of changepoints for Poisson data using PELT pruned method.
}
\usage{
PELT.meanvar.poisson(data, pen=0)
}
\arguments{
  \item{data}{
	A vector containing the data within which you wish to find changepoints.
}
  \item{pen}{
	Numeric value of the linear penalty function.  This value is used in the decision for each individual changepoint so that in total the penalty is k*pen where k is the optimal number of changepoints detected.
}
}
\details{
	This function is used to find a multiple changes in mean and variance for data that is assumed to be Poisson distributed.  The value returned is the result of testing H0:existing number of changepoints against H1: one extra changepoint using the log of the likelihood ratio statistic coupled with the penalty supplied.  The PELT method keeps track of the optimal number and location of changepoints as it passes through the data.
}
\value{
	A vector of the changepoint locations is returned:
	\item{cpt}{Vector containing the changepoint locations for the penalty supplied.  This always ends with n.}
}
\references{
Change in Poisson model: Chen, J. and Gupta, A. K. (2000) \emph{Parametric statistical change point analysis}, Birkhauser

PELT Algorithm: Killick, R. and Fearnhead, P. and Eckley, I.A. (2011) An exact linear time search algorithm for multiple changepoint detection, \emph{Submitted}
}
\author{
Rebecca Killick
}


\seealso{
\code{\link{cpt.meanvar}},\code{\link{binseg.meanvar.poisson}},\code{\link{multiple.meanvar.poisson}},\code{\link{single.meanvar.poisson}},\code{\link{segneigh.meanvar.poisson}}
}
\examples{
# Example of multiple changes in mean and variance at 50,100,150 in simulated Poisson data
set.seed(1)
x=c(rpois(50,lambda=1),rpois(50,lambda=3),rpois(50,lambda=1),rpois(50,lambda=10))
PELT.meanvar.poisson(x,pen=2*log(200)) # returns c(50,100,150,200)

# Example no change in parameter
set.seed(1)
x=rpois(200,lambda=1)
PELT.meanvar.poisson(x,pen=2*log(200)) # returns 200 to show no change in mean or variance has been found
}

\keyword{methods}
\keyword{univar}
\keyword{models}
\keyword{ts}
