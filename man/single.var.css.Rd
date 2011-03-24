\name{single.var.css}
\alias{single.var.css}
\title{
Single Change in Variance - Cumulative Sums of Squares
}
\description{
Calculates the cumulative sums of squares (css) test statistic for all possible changepoint locations and returns the single most probable (max).
}
\usage{
single.var.css(data,penalty="SIC",value=0,class=TRUE,param.estimates=TRUE)
}
\arguments{
  \item{data}{
	A vector or matrix containing the data within which you wish to find a changepoint.  If data is a matrix, each row is considered a separate dataset.
}
  \item{penalty}{
	Choice of "None", "SIC", "BIC", "AIC", "Hannan-Quinn", "Asymptotic" and "Manual" penalties.  If Manual is specified, the manual penalty is contained in the value parameter. If Asymptotic is specified, the theoretical type I error is contained in the value parameter.
}
  \item{value}{
	The theoretical type I error e.g.0.05 when using the Asymptotic penalty (options are 0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95).  The value of the penalty when using the Manual penalty option.  This can be a numeric value or text giving the formula to use.  Available variables are, n=length of original data, null=test statistic, tau=proposed changepoint, diffparam=difference in number of alternatve and null parameters.
}
  \item{class}{
	Logical.  If TRUE then an object of class \code{cpt} is returned. If FALSE a vector is returned.
}
  \item{param.estimates}{
	Logical.  If TRUE and class=TRUE then parameter estimates are returned. If FALSE or class=FALSE no parameter estimates are returned.
}
}
\details{
	This function is used to find a single change in variance for data that is is not assumed to follow a specific distribtuion.  The value returned is the result of testing H0:no change in variance against H1: single change in variance using the cumulative sums of squares test statistic coupled with the penalty supplied.
}
\value{
	If class=TRUE then an object of S4 class "cpt" is returned.  The slot \code{cpts} contains the changepoints that are solely returned if class=FALSE.  The structure of \code{cpts} is as follows.

	If data is a vector (single dataset) then a single value is returned:
	\item{cpt}{The most probable location of a changepoint if H0 was rejected or NA if H1 was rejected.}
	If data is an mxn matrix (multiple datasets) then a vector is returned:
	\item{cpt}{Vector of length m containing where each element is the result of the test for data[m,].  If cpt[m] is a number then it is the most probable location of a changepoint under H1.  Otherwise cpt[m] has the value NA and indicates that H1 was rejected.}
}
\references{
C. Inclan, G. C. Tiao (1994) Use of Cumulative Sums of Squares for Retrospective Detection of Changes of Variance, \emph{Journal of the American Statistical Association} \bold{89(427)}, 913--923

R. L. Brown, J. Durbin, J. M. Evans (1975) Techniques for Testing the Constancy of Regression Relationships over Time, \emph{Journal of the Royal Statistical Society B} \bold{32(2)}, 149--192
}
\author{
Rebecca Killick
}


\seealso{
\code{\link{cpt.var}},\code{\linkS4class{cpt}}
}
\examples{
# Example of a change in variance at 100 in simulated normal data
set.seed(1)
x=c(rnorm(100,0,1),rnorm(100,0,10))
single.var.css(x,penalty="Asymptotic",value=0.05,class=FALSE) # returns 105 to show that the null hypothesis was rejected and the change in variance is at 105
ans=single.var.css(x,penalty="Asymptotic",value=0.01) 
cpts(ans) # returns 105 to show that the null hypothesis was rejected, the change in variance is at 105 and we are 99% confident of this result

# Example of a data matrix containing 2 rows, row 1 has a change in variance and row 2 had no change in variance
set.seed(10)
x=c(rnorm(100,0,1),rnorm(100,0,10))
y=rnorm(200,0,1)
z=rbind(x,y)
single.var.css(z,penalty="Asymptotic",value=0.05,class=FALSE) # returns vector c(102,200) which shows that the first dataset has a change in variance at 102 and the second dataset rejected H1 and has no change in variance
ans=single.var.css(z,penalty="Manual",value=2) 
cpts(ans[[1]]) # test using a manual penalty which is the same as the AIC penalty for this example
cpts(ans[[2]]) # result is the same as above, c(102, 200)
}

\keyword{methods}
\keyword{univar}
\keyword{models}
\keyword{ts}
