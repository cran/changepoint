\name{binseg.mean.cusum}
\alias{binseg.mean.cusum}
\title{
Multiple Changes in Mean using Binary Segmentation method - Cumulative Sums
}
\description{
Calculates the optimal positioning and number of changepoints for the cumulative sums test statistic using Binary Segmentation method.  Note that this is an approximate method.
}
\usage{
binseg.mean.cusum(data, Q=5, pen=0)
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
	This function is used to find a multiple changes in mean for data where no assumption about the distribution is made.  The value returned is the result of finding the optimal location of up to Q changepoints using the cumulative sums test statistic.  Once all changepoint locations have been calculated, the optimal number of changepoints is decided using pen as the penalty function.
}
\value{
	A list is returned containing the following items
	\item{cps}{2xQ Matrix containing the changepoint positions on the first row and the test statistic on the second row.}
	\item{op.cpts}{The optimal changepoint locations for the penalty supplied.}
	\item{pen}{Penalty used to find the optimal number of changepoints.}
}
\references{
Binary Segmentation: Scott, A. J. and Knott, M. (1974) A Cluster Analysis Method for Grouping Means in the Analysis of Variance, \emph{Biometrics} \bold{30(3)}, 507--512

M. Csorgo, L. Horvath (1997) Limit Theorems in Change-Point Analysis, \emph{Wiley}

E. S. Page (1954) Continuous Inspection Schemes, \emph{Biometrika} \bold{41(1/2)}, 100--115
}
\author{
Rebecca Killick
}


\seealso{
\code{\link{binseg.mean.norm}},\code{\link{cpt.mean}},\code{\link{multiple.mean.cusum}},\code{\link{single.mean.cusum}},\code{\link{segneigh.mean.cusum}}
}
\examples{
# Example of multiple changes in mean at 50,100,150 in simulated data
set.seed(1)
x=c(rnorm(50,0,1),rnorm(50,5,1),rnorm(50,10,1),rnorm(50,3,1))
binseg.mean.cusum(x,Q=5,pen=0.8) # returns optimal number as 3 and the locations as c(50,100,150)
binseg.mean.cusum(x,Q=2,pen=0.8) # returns optimal number as 2 as this is the maximum number of
#changepoints it can find.  If you get the maximum number, you need to increase Q until this is not
#the case.

# Example no change in mean
set.seed(10)
x=rnorm(200,0,1)
binseg.mean.cusum(x,Q=5,pen=0.8) # returns optimal number as 0
}

\keyword{methods}
\keyword{univar}
\keyword{models}
\keyword{ts}
