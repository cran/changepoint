\name{PELT.meanvar.exp}
\alias{PELT.meanvar.exp}
\title{
Multiple Changes in Mean and Variance using PELT pruned method - Exponential Data
}
\description{
Calculates the optimal positioning and number of changepoints for Exponential data using PELT pruned method.
}
\usage{
PELT.meanvar.exp(data, pen=0)
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
	This function is used to find a multiple changes in mean and variance for data that is assumed to be Exponential distributed.  The value returned is the result of testing H0:existing number of changepoints against H1: one extra changepoint using the log of the likelihood ratio statistic coupled with the penalty supplied.  The PELT method keeps track of the optimal number and location of changepoints as it passes through the data.
}
\value{
	A vector of the changepoint locations is returned:
	\item{cpt}{Vector containing the changepoint locations for the penalty supplied.  This always ends with n.}
}
\references{
Change in Exponential model: Chen, J. and Gupta, A. K. (2000) \emph{Parametric statistical change point analysis}, Birkhauser

PELT Algorithm: Killick, R. and Fearnhead, P. and Eckley, I.A. (2011) An exact linear time search algorithm for multiple changepoint detection, \emph{Submitted}
}
\author{
Rebecca Killick
}


\seealso{
\code{\link{PELT.meanvar.norm}},\code{\link{PELT.meanvar.gamma}},\code{\link{cpt.meanvar}},\code{\link{binseg.meanvar.exp}},\code{\link{multiple.meanvar.exp}},\code{\link{single.meanvar.exp}},\code{\link{segneigh.meanvar.exp}}
}
\examples{
# Example of multiple changes in mean and variance at 50,100,150 in simulated Exponential data
set.seed(1)
x=c(rexp(50,rate=1),rexp(50,rate=3),rexp(50,rate=1),rexp(50,rate=10))
PELT.meanvar.exp(x,pen=2*log(200)) # returns c(53,100,150,200)

# Example no change in scale parameter
set.seed(1)
x=rexp(200,rate=1)
PELT.meanvar.exp(x,pen=2*log(200)) # returns 200 to show no change in mean or variance has been found
}

\keyword{methods}
\keyword{univar}
\keyword{models}
\keyword{ts}
