\name{segneigh.meanvar.gamma}
\alias{segneigh.meanvar.gamma}
\title{
Multiple Changes in Mean and Variance using Segment Neighbourhood method - Gamma Data (i.e. change in scale parameter)
}
\description{
Calculates the optimal positioning and number of changepoints for Gamma data using Segment Neighbourhood method.  Note that this gives the same results as PELT method but takes more computational time.
}
\usage{
segneigh.meanvar.gamma(data, shape=1, Q=5, pen=0)
}
\arguments{
  \item{data}{
	A vector containing the data within which you wish to find changepoints.
}
  \item{shape}{
	Numerical value of the true shape parameter for the data.  Either single value or vector of length nrow(data).  If data is a matrix and shape is a single value, the same shape parameter is used for each row.
}
  \item{Q}{
	Numeric value of the maximum number of segments (number of changepoints +1) you wish to search for, default is 5.
}
  \item{pen}{
	Numeric value of the linear penalty function.  This value is used in the final decision as to the optimal number of changepoints, used as k*pen where k is the number of changepoints to be tested.
}
}
\details{
	This function is used to find a multiple changes in mean and variance for data that is assumed to be Gamma distributed.  The value returned is the result of finding the optimal location of up to Q changepoints using the log of the likelihood ratio statistic.  Once all changepoint locations have been calculated, the optimal number of changepoints is decided using k*pen as the penalty function where k is the number of changepoints tested (k in (1,Q)).
}
\value{
	A list is returned containing the following items
	\item{cps}{Matrix containing the changepoint positions for 1,...,Q changepoints.}
	\item{op.cpts}{The optimal changepoint locations for the penalty supplied.}
	\item{like}{Value of the -2*log(likelihood ratio) + penalty for the optimal number of changepoints selected.}
}
\references{
Change in Gamma shape parameter: Chen, J. and Gupta, A. K. (2000) \emph{Parametric statistical change point analysis}, Birkhauser

Segment Neighbourhoods: Auger, I. E. And Lawrence, C. E. (1989) Algorithms for the Optimal Identification of Segment Neighborhoods, \emph{Bulletin of Mathematical Biology} \bold{51(1)}, 39--54
}
\author{
Rebecca Killick
}


\seealso{
\code{\link{segneigh.meanvar.norm}},\code{\link{cpt.meanvar}},\code{\link{PELT.meanvar.gamma}},\code{\link{multiple.meanvar.gamma}},\code{\link{single.meanvar.gamma}},\code{\link{binseg.meanvar.gamma}}
}
\examples{
# Example of multiple changes in mean and variance at 50,100,150 in simulated Gamma data
set.seed(1)
x=c(rgamma(50,shape=1,rate=1),rgamma(50,shape=1,rate=3),rgamma(50,shape=1,rate=1),rgamma(50,shape=1,
rate=10))
segneigh.meanvar.gamma(x,shape=1, Q=5, pen=2*log(200)) # returns optimal number as 3 and the
#locations as c(47,102,150)
segneigh.meanvar.gamma(x,shape=1, Q=3, pen=2*log(200)) # returns optimal number as 2 as this is the
#maximum number of changepoints it can find.  If you get the maximum number, you need to increase Q
#until this is not the case.

# Example no change in mean or variance
set.seed(1)
x=rgamma(200,shape=1,rate=1)
segneigh.meanvar.gamma(x,shape=1,pen=2*log(200)) # returns optimal number as 0
}

\keyword{methods}
\keyword{univar}
\keyword{models}
\keyword{ts}
