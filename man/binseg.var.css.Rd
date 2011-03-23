\name{binseg.var.css}
\alias{binseg.var.css}
\title{
Multiple Changes in Variance using Binary Segmentation method - Cumulative Sums of Squares
}
\description{
Calculates the optimal positioning and number of changepoints for the cumulative sums of squares test statistic using Binary Segmentation method.  Note that this is an approximate method.
}
\usage{
binseg.var.css(data, Q=5, pen=0)
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
	This function is used to find a multiple changes in variance for data where no assumption about the distribution is made.  The value returned is the result of finding the optimal location of up to Q changepoints using the cumulative sums of squares test statistic.  Once all changepoint locations have been calculated, the optimal number of changepoints is decided using pen as the penalty function.
}
\value{
	A list is returned containing the following items
	\item{cps}{2xQ Matrix containing the changepoint positions on the first row and the test statistic on the second row.}
	\item{op.cpts}{The optimal changepoint locations for the penalty supplied.}
	\item{pen}{Penalty used to find the optimal number of changepoints.}
}
\references{
Binary Segmentation: Scott, A. J. and Knott, M. (1974) A Cluster Analysis Method for Grouping Means in the Analysis of Variance, \emph{Biometrics} \bold{30(3)}, 507--512

C. Inclan, G. C. Tiao (1994) Use of Cumulative Sums of Squares for Retrospective Detection of Changes of Variance, \emph{Journal of the American Statistical Association} \bold{89(427)}, 913--923

R. L. Brown, J. Durbin, J. M. Evans (1975) Techniques for Testing the Constancy of Regression Relationships over Time, \emph{Journal of the Royal Statistical Society B} \bold{32(2)}, 149--192
}
\author{
Rebecca Killick
}


\seealso{
\code{\link{binseg.var.norm}},\code{\link{cpt.var}},\code{\link{multiple.var.css}},\code{\link{single.var.css}},\code{\link{segneigh.var.css}}
}
\examples{
# Example of multiple changes in variance at 50,100,150 in simulated normal data
set.seed(10)
x=c(rnorm(50,0,1),rnorm(50,0,10),rnorm(50,0,5),rnorm(50,0,1))
binseg.var.css(x,Q=5, pen=1.358) # returns optimal number as 4 and the locations as c(50,52,100,149)
binseg.var.css(x,Q=2, pen=1.358) # returns optimal number as 2 as this is the maximum number of changepoints it can find.  If you get the maximum number, you need to increase Q until this is not the case.
# 1.358 is the asymptotic value of the penalty for 95% confidence

# Example no change in variance
set.seed(1)
x=rnorm(200,0,1)
binseg.var.css(x,Q=5, pen=1.358) # returns optimal number as 0
}

\keyword{methods}
\keyword{univar}
\keyword{models}
\keyword{ts}
