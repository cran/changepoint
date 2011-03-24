\name{single.mean.cusum}
\alias{single.mean.cusum}
\title{
Single Change in Mean - Cumulative Sums
}
\description{
Calculates the cumulative sums (cusum) test statistic for all possible changepoint locations and returns the single most probable (max).
}
\usage{
single.mean.cusum(data,penalty="Asymptotic",value=0.05,class=TRUE,param.estimates=TRUE)
}
\arguments{
  \item{data}{
	A vector or matrix containing the data within which you wish to find a changepoint.  If data is a matrix, each row is considered a separate dataset.
}
  \item{penalty}{
	Choice of "None", "SIC", "BIC", "AIC", "Hannan-Quinn" and "Manual" penalties.  If Manual is specified, the manual penalty is contained in the value parameter.
}
  \item{value}{
	The value of the penalty when using the Manual penalty option.  This can be a numeric value or text giving the formula to use.  Available variables are, n=length of original data, null=test statistic, tau=proposed changepoint, diffparam=difference in number of alternatve and null parameters.
}
  \item{class}{
	Logical.  If TRUE then an object of class \code{cpt} is returned. If FALSE a vector is returned.
}
  \item{param.estimates}{
	Logical.  If TRUE and class=TRUE then parameter estimates are returned. If FALSE or class=FALSE no parameter estimates are returned.
}
}
\details{
	This function is used to find a single change in mean for data that is is not assumed to follow a specific distribtuion.  The value returned is the result of testing H0:no change in mean against H1: single change in mean using the cumulative sums test statistic coupled with the penalty supplied.

Warning: The prescribed penalty values are not defined for use on CUSUM tests.  The values tend to be too large and thus manual penalties are preferred.
}
\value{
	If class=TRUE then an object of S4 class "cpt" is returned.  The slot \code{cpts} contains the changepoints that are solely returned if class=FALSE.  The structure of \code{cpts} is as follows.

	If data is a vector (single dataset) then a single value is returned:
	\item{cpt}{The most probable location of a changepoint if H0 was rejected or NA if H1 was rejected.}
	If data is an mxn matrix (multiple datasets) then a vector is returned:
	\item{cpt}{Vector of length m containing where each element is the result of the test for data[m,].  If cpt[m] is a number then it is the most probable location of a changepoint under H1.  Otherwise cpt[m] has the value NA and indicates that H1 was rejected.}
}
\references{
M. Csorgo, L. Horvath (1997) Limit Theorems in Change-Point Analysis, \emph{Wiley}

E. S. Page (1954) Continuous Inspection Schemes, \emph{Biometrika} \bold{41(1/2)}, 100--115
}
\author{
Rebecca Killick
}


\seealso{
\code{\link{cpt.mean}},\code{\linkS4class{cpt}}
}
\examples{
# Example of a change in mean at 100 in simulated normal data
set.seed(1)
x=c(rnorm(100,0,1),rnorm(100,10,1))
single.mean.cusum(x,penalty="Manual",value=1,class=FALSE) # returns 101 to show that the null hypothesis was rejected and the change in mean is at 101
ans=single.mean.cusum(x,penalty="Manual",value=1) 
cpts(ans) # returns 101 to show that the null hypothesis was rejected, the change in mean is at 101


# Example of a data matrix containing 2 rows, row 1 has a change in mean and row 2 had no change in mean
set.seed(1)
x=c(rnorm(100,0,1),rnorm(100,10,1))
y=rnorm(200,0,1)
z=rbind(x,y)
single.mean.cusum(z,penalty="Manual",value=1,class=FALSE) # returns vector c(101,200) which shows that the first dataset has a change in mean at 101 and the second dataset rejected H1 and has no change in mean
ans=single.mean.cusum(z,penalty="Manual",value="log(log(n))") 
cpts(ans[[1]]) # test using a manual penalty which is the log(SIC) penalty and gives the same values for this example 
cpts(ans[[2]]) # result is the same as above, c(101, 200)
}

\keyword{methods}
\keyword{univar}
\keyword{models}
\keyword{ts}
