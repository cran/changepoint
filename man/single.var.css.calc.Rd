\name{single.var.css.calc}
\alias{single.var.css.calc}
\title{
Single Change in Variance - Cumulative Sums of Squares
}
\description{
Calculates the cumulative sums of squares (css) test statistic for all possible changepoint locations and returns the single most probable (max).
}
\usage{
single.var.css.calc(data, extrainf = TRUE)
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
	This function is used to find a single change in variance for data where no distributional assumption is made.  The changepoint returned is simply the location where the test statistic is maximised, there is no test performed as to whether this location is a true changepoint or not.

In reality this function should not be used unless you are performing a changepoint test using the output supplied.  This function is used in the "see also" functions that perform various changepoint tests, ideally these should be used.
}
\value{
	If data is a vector (single dataset) and extrainf=FALSE then a single value is returned:
	\item{cpt}{The most probable location of a changepoint}
	If data is a vector (single dataset) and extrainf=TRUE then a vector with two elements is returned:
	\item{test statistic}{The cumulative sums of squares test statistic}
	If data is an mxn matrix (multiple datasets) and extrainf=FALSE then a vector is returned:
	\item{cpt}{Vector of length m containing the most probable location of a changepoint for each row in data.  cpt[1] is the most probable changepoint of the first row in data and cpt[m] is the most probable changepoint for the final row in data.}
	If data is a matrix (multiple datasets) and extrainf=TRUE then a matrix is returned where the first column is the changepoint location for each row in data, the second column is the test statistic for each row in data.
}
\references{
C. Inclan, G. C. Tiao (1994) Use of Cumulative Sums of Squares for Retrospective Detection of Changes of Variance, \emph{Journal of the American Statistical Association} \bold{89(427)}, 913--923

R. L. Brown, J. Durbin, J. M. Evans (1975) Techniques for Testing the Constancy of Regression Relationships over Time, \emph{Journal of the Royal Statistical Society B} \bold{32(2)}, 149--192
}
\author{
Rebecca Killick
}


\seealso{
\code{\link{single.var.css}}, \code{\link{cpt.var}}
}
\examples{
# Example of a change in variance at 100 in simulated normal data
set.seed(1)
x=c(rnorm(100,0,1),rnorm(100,0,10))
single.var.css.calc(x,extrainf=FALSE) # finds change at 105
single.var.css.calc(x) # finds change at 105 and gives test statistic as 4.979771


# Example of no change in variance in simulated normal data
set.seed(1)
x=rnorm(100,0,1)
single.var.css.calc(x,extrainf=FALSE) # finds change at 53, this is the most probable point of
#change but if a changepoint test is performed then no change will be found.
single.var.css.calc(x)# change at 53, test statistic is 0.6922931
}

\keyword{methods}
\keyword{univar}
\keyword{models}
\keyword{ts}
