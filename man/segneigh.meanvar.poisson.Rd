\name{segneigh.meanvar.poisson}
\alias{segneigh.meanvar.poisson}
\title{
Multiple Changes in Mean and Variance using Segment Neighbourhood method - Poisson Data
}
\description{
Calculates the optimal positioning and number of changepoints for Poisson data using Segment Neighbourhood method.  Note that this gives the same results as PELT method but takes more computational time.
}
\usage{
segneigh.meanvar.poisson(data, Q=5, pen=0)
}
\arguments{
  \item{data}{
	A vector containing the data within which you wish to find changepoints.
}
  \item{Q}{
	Numeric value of the maximum number of segments (number of changepoints +1) you wish to search for, default is 5.
}
  \item{pen}{
	Numeric value of the linear penalty function.  This value is used in the final decision as to the optimal number of changepoints, used as k*pen where k is the number of changepoints to be tested.
}
}
\details{
	This function is used to find a multiple changes in mean and variance for data that is assumed to be Poisson distributed.  The value returned is the result of finding the optimal location of up to Q changepoints using the log of the likelihood ratio statistic.  Once all changepoint locations have been calculated, the optimal number of changepoints is decided using k*pen as the penalty function where k is the number of changepoints tested (k in (1,Q)).
}
\value{
	A list is returned containing the following items
	\item{cps}{Matrix containing the changepoint positions for 1,...,Q changepoints.}
	\item{op.cpts}{The optimal changepoint locations for the penalty supplied.}
	\item{like}{Value of the -2*log(likelihood ratio) + penalty for the optimal number of changepoints selected.}
}
\references{
Change in Poisson model: Chen, J. and Gupta, A. K. (2000) \emph{Parametric statistical change point analysis}, Birkhauser

Segment Neighbourhoods: Auger, I. E. And Lawrence, C. E. (1989) Algorithms for the Optimal Identification of Segment Neighborhoods, \emph{Bulletin of Mathematical Biology} \bold{51(1)}, 39--54
}
\author{
Rebecca Killick
}


\seealso{
\code{\link{cpt.meanvar}},\code{\link{PELT.meanvar.poisson}},\code{\link{multiple.meanvar.poisson}},\code{\link{single.meanvar.poisson}},\code{\link{binseg.meanvar.poisson}}
}
\examples{
# Example of multiple changes in mean and variance at 50,100,150 in simulated Poisson data
set.seed(1)
x=c(rpois(50,lambda=1),rpois(50,lambda=5),rpois(50,lambda=1),rpois(50,lambda=10))
segneigh.meanvar.poisson(x,Q=5, pen=2*log(200)) # returns optimal number as 3 and the locations as
#c(50,100,150)
segneigh.meanvar.poisson(x,Q=3, pen=2*log(200)) # returns optimal number as 2 as this is the maximum
#number of changepoints it can find.  If you get the maximum number, you need to increase Q until
#this is not the case.

# Example no change in mean or variance
set.seed(1)
x=rpois(200,lambda=1)
segneigh.meanvar.poisson(x,pen=2*log(200)) # returns optimal number as 0
}

\keyword{methods}
\keyword{univar}
\keyword{models}
\keyword{ts}
