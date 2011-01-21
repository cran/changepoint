\name{binseg.meanvar.norm}
\alias{binseg.meanvar.norm}
\title{
Multiple Changes in Mean and Variance using Binary Segmentation method - Normal Data
}
\description{
Calculates the optimal positioning and number of changepoints for Normal data using Binary Segmentation method.  Note that this is an approximate method.
}
\usage{
binseg.meanvar.norm(data, Q=5, pen=0)
}
\arguments{
  \item{data}{
	A vector containing the data within which you wish to find changepoints.
}
  \item{Q}{
	Numeric value of the maximum number of changepoints you wish to search for, default is 5.
}
  \item{pen}{
	Numeric value of the linear penalty function.  This value is used in the decision as to the optimal number of changepoints.
}
}
\details{
	This function is used to find a multiple changes in mean and variance for data that is assumed to be normally distributed.  The value returned is the result of finding the optimal location of up to Q changepoints using the log of the likelihood ratio statistic.  Once all changepoint locations have been calculated, the optimal number of changepoints is decided using pen as the penalty function.
}
\value{
	A list is returned containing the following items
	\item{cps}{2xQ Matrix containing the changepoint positions on the first row and the test statistic on the second row.}
	\item{op.cpts}{The optimal changepoint locations for the penalty supplied.}
	\item{pen}{Penalty used to find the optimal number of changepoints.}
}
\references{
Binary Segmentation: Scott, A. J. and Knott, M. (1974) A Cluster Analysis Method for Grouping Means in the Analysis of Variance, \emph{Biometrics} \bold{30(3)}, 507--512

Change in mean and variance: Chen, J. and Gupta, A. K. (2000) \emph{Parametric statistical change point analysis}, Birkhauser
}
\author{
Rebecca Killick
}


\seealso{
\code{\link{binseg.var.norm}},\code{\link{binseg.mean.norm}},\code{\link{cpt.meanvar}},\code{\link{PELT.meanvar.norm}},\code{\link{multiple.meanvar.norm}},\code{\link{single.meanvar.norm}},\code{\link{segneigh.meanvar.norm}}
}
\examples{
# Example of multiple changes in mean and variance at 50,100,150 in simulated normal data
set.seed(1)
x=c(rnorm(50,0,1),rnorm(50,5,3),rnorm(50,10,1),rnorm(50,3,10))
binseg.meanvar.norm(x,Q=5, pen=4*log(200)) # returns optimal number as 4 and the locations as c(50,100,150,152)
binseg.meanvar.norm(x,Q=3, pen=4*log(200)) # returns optimal number as 2 as this is the maximum number of changepoints it can find.  If you get the maximum number, you need to increase Q until this is not the case.

# Example no change in mean or variance
set.seed(1)
x=rnorm(200,0,1)
binseg.meanvar.norm(x,pen=4*log(200)) # returns optimal number as 0
}

\keyword{methods}
\keyword{univar}
\keyword{models}
\keyword{ts}
