\name{single.var.norm}
\alias{single.var.norm}
\title{
Single Change in Variance - Normal Data
}
\description{
Calculates the scaled log-likelihood (assuming the data is normally distributed) for all possible changepoint locations and returns the single most probable (max).
}
\usage{
single.var.norm(data,penalty="SIC",pen.value=0,know.mean=FALSE,mu=NA,class=TRUE,
param.estimates=TRUE)
}
\arguments{
  \item{data}{
	A vector, ts object or matrix containing the data within which you wish to find a changepoint.  If data is a matrix, each row is considered a separate dataset.
}
  \item{penalty}{
	Choice of "None", "SIC", "BIC", "AIC", "Hannan-Quinn", "Asymptotic" and "Manual" penalties.  If Manual is specified, the manual penalty is contained in the pen.value parameter. If Asymptotic is specified, the theoretical type I error is contained in the pen.value parameter.  The predefined penalties listed do NOT count the changepoint as a parameter, postfix a 1 e.g."SIC1" to count the changepoint as a parameter.
}
  \item{pen.value}{
	The theoretical type I error e.g.0.05 when using the Asymptotic penalty.  The value of the penalty when using the Manual penalty option.  This can be a numeric value or text giving the formula to use.  Available variables are, n=length of original data, null=null likelihood, alt=alternative likelihood, tau=proposed changepoint, diffparam=difference in number of alternatve and null parameters.
}
  \item{know.mean}{
	Logical, if TRUE then the mean is assumed known and mu is taken as its value.  If FALSE, and mu=-1000 (default value) then the mean is estimated via maximum likelihood.  If FALSE and the value of mu is supplied, mu is not estimated but is counted as an estimated parameter for decisions.
}
  \item{mu}{
	Numerical value of the true mean of the data.  Either single value or vector of length nrow(data).  If data is a matrix and mu is a single value, the same mean is used for each row.
}
  \item{class}{
	Logical.  If TRUE then an object of class \code{cpt} is returned. If FALSE a vector is returned.
}
  \item{param.estimates}{
	Logical.  If TRUE and class=TRUE then parameter estimates are returned. If FALSE or class=FALSE no parameter estimates are returned.
}
}
\details{
	This function is used to find a single change in variance for data that is assumed to be normally distributed.  The value returned is the result of testing H0:no change in variance against H1: single change in variance using the log of the likelihood ratio statistic coupled with the penalty supplied.
}
\value{
	If class=TRUE then an object of S4 class "cpt" is returned.  The slot \code{cpts} contains the changepoints that are returned (with p-value) if class=FALSE.  The structure of \code{cpts} is as follows.

	If data is a vector (single dataset) then a single value is returned:
	\item{cpt}{The most probable location of a changepoint if H0 was rejected or NA if H1 was rejected.}
	If data is an mxn matrix (multiple datasets) then a vector is returned:
	\item{cpt}{Vector of length m containing where each element is the result of the test for data[m,].  If cpt[m] is a number then it is the most probable location of a changepoint under H1.  Otherwise cpt[m] has the value NA and indicates that H1 was rejected.}
}
\references{
Change in variance: Chen, J. and Gupta, A. K. (2000) \emph{Parametric statistical change point analysis}, Birkhauser
}
\author{
Rebecca Killick
}


\seealso{
\code{\link{cpt.var}},\code{\link{single.meanvar.norm}},\code{\link{single.mean.norm}},\code{\linkS4class{cpt}}
}
\examples{
# Example of a change in variance at 100 in simulated normal data
set.seed(1)
x=c(rnorm(100,0,1),rnorm(100,0,10))
single.var.norm(x,penalty="SIC",class=FALSE) # returns 100 to show that the null hypothesis was
#rejected and the change in variance is at 100
ans=single.var.norm(x,penalty="Asymptotic",pen.value=0.01) 
cpts(ans) # returns 100 to show that the null hypothesis was rejected, the change in variance is at
#100 and we are 99% confident of this result

# Example of a data matrix containing 2 rows, row 1 has a change in variance and row 2 had no change
#in variance
set.seed(10)
x=c(rnorm(100,0,1),rnorm(100,0,10))
y=rnorm(200,0,1)
z=rbind(x,y)
single.var.norm(z,penalty="SIC",class=FALSE) # returns vector c(100,200) which shows that the first
#dataset has a change in variance at 100 and the second dataset rejected H1 and has no change in
#variance
ans=single.var.norm(z,penalty="Manual",pen.value="diffparam*log(n)") 
cpts(ans[[1]]) # test using a manual penalty which is the same as the SIC penalty for this example
cpts(ans[[2]]) # result is the same as above, c(100, 200)
}

\keyword{methods}
\keyword{univar}
\keyword{models}
\keyword{ts}
