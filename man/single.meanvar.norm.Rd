\name{single.meanvar.norm}
\alias{single.meanvar.norm}
\title{
Single Change in Mean and Variance - Normal Data
}
\description{
Calculates the scaled log-likelihood (assuming the data is normally distributed) for all possible changepoint locations and returns the single most probable (max).
}
\usage{
single.meanvar.norm(data, penalty="SIC",value=0,class=TRUE,param.estimates=TRUE)
}
\arguments{
  \item{data}{
	A vector or matrix containing the data within which you wish to find a changepoint.  If data is a matrix, each row is considered a separate dataset.
}
  \item{penalty}{
	Choice of "None", "SIC", "BIC", "AIC", "Hannan-Quinn", "Asymptotic" and "Manual" penalties.  If Manual is specified, the manual penalty is contained in the value parameter. If Asymptotic is specified, the theoretical type I error is contained in the value parameter.  The predefined penalties listed do NOT count the changepoint as a parameter, postfix a 1 e.g."SIC1" to count the changepoint as a parameter.
}
  \item{value}{
	The theoretical type I error e.g.0.05 when using the Asymptotic penalty.  The value of the penalty when using the Manual penalty option.  This can be a numeric value or text giving the formula to use.  Available variables are, n=length of original data, null=null likelihood, alt=alternative likelihood, tau=proposed changepoint, diffparam=difference in number of alternatve and null parameters.
}
  \item{class}{
	Logical.  If TRUE then an object of class \code{cpt} is returned. If FALSE a vector is returned.
}
  \item{param.estimates}{
	Logical.  If TRUE and class=TRUE then parameter estimates are returned. If FALSE or class=FALSE no parameter estimates are returned.
}
}
\details{
	This function is used to find a single change in mean and variance for data that is assumed to be normally distributed.  The value returned is the result of testing H0:no change in mean or variance against H1: single change in mean and/or variance using the log of the likelihood ratio statistic coupled with the penalty supplied.
}
\value{
	If class=TRUE then an object of S4 class "cpt" is returned.  The slot \code{cpts} contains the changepoints that are solely returned if class=FALSE.  The structure of \code{cpts} is as follows.

	If data is a vector (single dataset) then a single value is returned:
	\item{cpt}{The most probable location of a changepoint if H0 was rejected or NA if H1 was rejected.}
	If data is an mxn matrix (multiple datasets) then a vector is returned:
	\item{cpt}{Vector of length m containing where each element is the result of the test for data[m,].  If cpt[m] is a number then it is the most probable location of a changepoint under H1.  Otherwise cpt[m] has the value NA and indicates that H1 was rejected.}
}
\references{
Change in Normal mean and variance: Chen, J. and Gupta, A. K. (2000) \emph{Parametric statistical change point analysis}, Birkhauser
}
\author{
Rebecca Killick
}


\seealso{
\code{\link{cpt.meanvar}},\code{\link{single.mean.norm}},\code{\link{single.var.norm}},\code{\linkS4class{cpt}}
}
\examples{
# Example of a change in mean and variance at 100 in simulated normal data
set.seed(1)
x=c(rnorm(100,0,1),rnorm(100,1,10))
single.meanvar.norm(x,penalty="SIC",class=FALSE) # returns 99 to show that the null hypothesis was rejected and the change in mean and variance is at 99
ans=single.meanvar.norm(x,penalty="Asymptotic",value=0.01) 
cpts(ans) # returns 99 to show that the null hypothesis was rejected, the change in mean and variance is at 99 and we are 99% confident of this result

# Example of a data matrix containing 2 rows, row 1 has a change in mean and variance and row 2 had no change in mean or variance
set.seed(10)
x=c(rnorm(100,0,1),rnorm(100,1,10))
y=rnorm(200,0,1)
z=rbind(x,y)
single.meanvar.norm(z,penalty="Asymptotic",value=0.01,class=FALSE) # returns vector c(99,200) which shows that the first dataset has a change in mean and variance at 99 and the second dataset rejected H1 and has no change in mean or variance
ans=single.meanvar.norm(z,penalty="Manual",value="diffparam*log(n)") # list returned
cpts(ans[[1]]) # test using a manual penalty which is the same as the SIC penalty for this example, result shows that the penalty is too small.  The same changepoint is detected for the first dataset 
cpts(ans[[2]]) # but the second dataset returns a change in mean and variance at 198 which was rejected under the asymptotic penalty above.
}

\keyword{methods}
\keyword{univar}
\keyword{models}
\keyword{ts}
