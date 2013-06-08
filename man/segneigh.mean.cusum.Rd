\name{segneigh.mean.cusum}
\alias{segneigh.mean.cusum}
\title{
Multiple Changes in Mean using Segment Neighbourhood method - Cumulative Sums
}
\description{
Calculates the optimal positioning and number of changepoints for Cumulative Sums test statistic using Segment Neighbourhood method.
}
\usage{
segneigh.mean.cusum(data, Q=5, pen=0)
}
\arguments{
  \item{data}{
	A vector containing the data within which you wish to find changepoints.
}
  \item{Q}{
	Numeric value of the maximum number of segments (number of changepoints +1) you wish to search for, default is 5.
}
  \item{pen}{
	Numeric value of the linear penalty function.  This value is used in the final decision as to the optimal number of changepoints.
}
}
\details{
	This function is used to find a multiple changes in mean for data that is not assumed to have a particular distribution.  The value returned is the result of finding the optimal location of up to Q changepoints using the cumulative sums test statistic.  Once all changepoint locations have been calculated, the optimal number of changepoints is decided using pen as the penalty function.
}
\value{
	A list is returned containing the following items
	\item{cps}{Matrix containing the changepoint positions for 1,...,Q changepoints.}
	\item{op.cpts}{The optimal changepoint locations for the penalty supplied.}
}
\references{
M. Csorgo, L. Horvath (1997) Limit Theorems in Change-Point Analysis, \emph{Wiley}

E. S. Page (1954) Continuous Inspection Schemes, \emph{Biometrika} \bold{41(1/2)}, 100--115

Segment Neighbourhoods: Auger, I. E. And Lawrence, C. E. (1989) Algorithms for the Optimal Identification of Segment Neighborhoods, \emph{Bulletin of Mathematical Biology} \bold{51(1)}, 39--54
}
\author{
Rebecca Killick
}


\seealso{
\code{\link{segneigh.mean.cusum}},\code{\link{cpt.mean}},\code{\link{multiple.mean.cusum}},\code{\link{single.mean.cusum}},\code{\link{binseg.mean.cusum}}
}
\examples{
# Example of multiple changes in mean at 50,100,150 in simulated normal data
set.seed(1)
x=c(rnorm(50,0,1),rnorm(50,5,1),rnorm(50,10,1),rnorm(50,3,1))
segneigh.mean.cusum(x,Q=5,pen=1) # returns optimal number as 3 and the locations as c(50,101,150)
segneigh.mean.cusum(x,Q=3,pen=1) # returns optimal number as 2 as this is the maximum number of
#changepoints it can find.  If you get the maximum number, you need to increase Q until this is not
#the case.

# Example no change in mean
set.seed(10)
x=rnorm(200,0,1)
segneigh.mean.cusum(x,Q=5,pen=1) # returns optimal number as 0
}

\keyword{methods}
\keyword{univar}
\keyword{models}
\keyword{ts}
