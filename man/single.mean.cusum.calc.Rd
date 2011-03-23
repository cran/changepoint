\name{single.mean.cusum.calc}
\alias{single.mean.cusum.calc}
\title{
Single Change in Mean - Cumulative Sums
}
\description{
Calculates the cumulative sums (cusum) test statistic for all possible changepoint locations and returns the single most probable (max).
}
\usage{
single.mean.cusum.calc(data, extrainf = TRUE)
}
\arguments{
  \item{data}{
	A vector or matrix containing the data within which you wish to find a changepoint.  If data is a matrix, each row is considered a separate dataset.
}
  \item{extrainf}{
	Logical, if TRUE the test statistic is returned along with the changepoint location.  If FALSE, only the changepoint location is returned.
}
}
\details{
	This function is used to find a single change in mean for data where no distributional assumption is made.  The changepoint returned is simply the location where the test statistic is maximised, there is no test performed as to whether this location is a true changepoint or not.

In reality this function should not be used unless you are performing a changepoint test using the output supplied.  This function is used in the "see also" functions that perform various changepoint tests, ideally these should be used.
}
\value{
	If data is a vector (single dataset) and extrainf=FALSE then a single value is returned:
	\item{cpt}{The most probable location of a changepoint}
	If data is a vector (single dataset) and extrainf=TRUE then a vector with two elements is returned:
	\item{test statistic}{The cumulative sums test statistic}
	If data is an mxn matrix (multiple datasets) and extrainf=FALSE then a vector is returned:
	\item{cpt}{Vector of length m containing the most probable location of a changepoint for each row in data.  cpt[1] is the most probable changepoint of the first row in data and cpt[m] is the most probable changepoint for the final row in data.}
	If data is a matrix (multiple datasets) and extrainf=TRUE then a matrix is returned where the first column is the changepoint location for each row in data, the second column is the test statistic for each row in data.
}
\references{
M. Csorgo, L. Horvath (1997) Limit Theorems in Change-Point Analysis, \emph{Wiley}

E. S. Page (1954) Continuous Inspection Schemes, \emph{Biometrika} \bold{41(1/2)}, 100--115
}
\author{
Rebecca Killick
}


\seealso{
\code{\link{single.mean.cusum}}, \code{\link{cpt.mean}}
}
\examples{
# Example of a change in mean at 100 in simulated data
set.seed(1)
x=c(rnorm(100,0,1),rnorm(100,10,1))
single.mean.cusum.calc(x,extrainf=FALSE) # finds change at 101
single.mean.cusum.calc(x) # finds change at 101 and gives the test statistic as 2.463326


# Example of no change in mean in simulated data
set.seed(1)
x=rnorm(100,0,1)
single.mean.cusum.calc(x,extrainf=FALSE) # finds change at 97, this is the most probable point of change but if a changepoint test is performed then no change will be found.
single.mean.cusum.calc(x)# change at 97, test statistic is 0.0398342
}

\keyword{methods}
\keyword{univar}
\keyword{models}
\keyword{ts}
