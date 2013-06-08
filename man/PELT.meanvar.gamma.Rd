\name{PELT.meanvar.gamma}
\alias{PELT.meanvar.gamma}
\title{
Multiple Changes in Mean and Variance using PELT pruned method - Gamma Data (i.e. change in scale parameter)
}
\description{
Calculates the optimal positioning and number of changepoints for Gamma data using PELT pruned method.
}
\usage{
PELT.meanvar.gamma(data, shape=1, pen=0)
}
\arguments{
  \item{data}{
	A vector containing the data within which you wish to find changepoints.
}
  \item{shape}{
	Numerical value of the true shape parameter for the data.  Either single value or vector of length nrow(data).  If data is a matrix and shape is a single value, the same shape parameter is used for each row.
}
  \item{pen}{
	Numeric value of the linear penalty function.  This value is used in the decision for each individual changepoint so that in total the penalty is k*pen where k is the optimal number of changepoints detected.
}
}
\details{
	This function is used to find a multiple changes in mean and variance for data that is assumed to be Gamma distributed.  The value returned is the result of testing H0:existing number of changepoints against H1: one extra changepoint using the log of the likelihood ratio statistic coupled with the penalty supplied.  The PELT method keeps track of the optimal number and location of changepoints as it passes through the data.
}
\value{
	A vector of the changepoint locations is returned:
	\item{cpt}{Vector containing the changepoint locations for the penalty supplied.  This always ends with n.}
}
\references{
Change in Gamma shape parameter: Chen, J. and Gupta, A. K. (2000) \emph{Parametric statistical change point analysis}, Birkhauser

PELT Algorithm: Killick R, Fearnhead P, Eckley IA (2012) Optimal detection of changepoints with a linear computational cost, \emph{JASA} \bold{107(500)}, 1590--1598
}
\author{
Rebecca Killick
}


\seealso{
\code{\link{PELT.meanvar.norm}},\code{\link{cpt.meanvar}},\code{\link{binseg.meanvar.gamma}},\code{\link{multiple.meanvar.gamma}},\code{\link{single.meanvar.gamma}},\code{\link{segneigh.meanvar.gamma}}
}
\examples{
# Example of multiple changes in scale parameter at 50,100,150 in simulated gamma data
set.seed(1)
x=c(rgamma(50,shape=1,rate=1),rgamma(50,shape=1,rate=3),rgamma(50,shape=1,rate=1),rgamma(50,shape=1,
rate=10))
PELT.meanvar.gamma(x,pen=2*log(200)) # returns c(47,102,150,200)

# Example no change in scale parameter
set.seed(1)
x=rgamma(200,shape=1,rate=1)
PELT.meanvar.gamma(x,pen=2*log(200)) # returns 200 to show no change in mean or variance has been
#found
}

\keyword{methods}
\keyword{univar}
\keyword{models}
\keyword{ts}
