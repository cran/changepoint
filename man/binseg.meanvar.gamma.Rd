\name{binseg.meanvar.gamma}
\alias{binseg.meanvar.gamma}
\title{
Multiple Changes in Mean and Variance using Binary Segmentation method - Gamma Data (i.e. change in scale parameter)
}
\description{
Calculates the optimal positioning and number of changepoints for Gamma data using Binary Segmentation method.  Note that this is an approximate method.
}
\usage{
binseg.meanvar.gamma(data, shape=1, Q=5, pen=0)
}
\arguments{
  \item{data}{
	A vector containing the data within which you wish to find changepoints.
}
  \item{shape}{
	Numerical value of the true shape parameter for the data.  Either single value or vector of length nrow(data).  If data is a matrix and shape is a single value, the same shape parameter is used for each row.
}
  \item{Q}{
	Numeric value of the maximum number of changepoints you wish to search for, default is 5.
}
  \item{pen}{
	Numeric value of the linear penalty function.  This value is used in the decision as to the optimal number of changepoints.
}
}
\details{
	This function is used to find a multiple changes in mean and variance for data that is assumed to be Gamma distributed.  The value returned is the result of finding the optimal location of up to Q changepoints using the log of the likelihood ratio statistic.  Once all changepoint locations have been calculated, the optimal number of changepoints is decided using pen as the penalty function.
}
\value{
	A list is returned containing the following items
	\item{cps}{2xQ Matrix containing the changepoint positions on the first row and the test statistic on the second row.}
	\item{op.cpts}{The optimal changepoint locations for the penalty supplied.}
	\item{pen}{Penalty used to find the optimal number of changepoints.}
}
\references{
Binary Segmentation: Scott, A. J. and Knott, M. (1974) A Cluster Analysis Method for Grouping Means in the Analysis of Variance, \emph{Biometrics} \bold{30(3)}, 507--512

Change in Gamma scale parameter: Chen, J. and Gupta, A. K. (2000) \emph{Parametric statistical change point analysis}, Birkhauser
}
\author{
Rebecca Killick
}


\seealso{
\code{\link{binseg.meanvar.norm}},\code{\link{cpt.meanvar}},\code{\link{PELT.meanvar.gamma}},\code{\link{multiple.meanvar.gamma}},\code{\link{single.meanvar.gamma}},\code{\link{segneigh.meanvar.gamma}}
}
\examples{
# Example of multiple changes in mean and variance at 50,100,150 in simulated Gamma data
set.seed(1)
x=c(rgamma(50,shape=1,rate=1),rgamma(50,shape=1,rate=3),rgamma(50,shape=1,rate=1),rgamma(50,shape=1,
rate=10))
binseg.meanvar.gamma(x,shape=1,Q=5,pen=2*log(200)) # returns optimal number as 3 and the locations
#as c(47,104,150,152)
binseg.meanvar.gamma(x,shape=1,Q=2,pen=2*log(200)) # returns optimal number as 2 as this is the
#maximum number of changepoints it can find.  If you get the maximum number, you need to increase Q
#until this is not the case.

# Example no change in mean or variance
set.seed(1)
x=rgamma(200,shape=1,rate=1)
binseg.meanvar.gamma(x,shape=1,pen=2*log(200)) # returns optimal number as 0
}

\keyword{methods}
\keyword{univar}
\keyword{models}
\keyword{ts}
