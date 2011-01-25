\name{single.reg.norm}
\alias{single.reg.norm}
\title{
Single Change in Regression - Normal Data
}
\description{
Calculates the scaled log-likelihood (assuming the errors are zero-mean normally distributed) for all possible changepoint locations and returns the single most probable (max).
}
\usage{
single.reg.norm(data, penalty="SIC",value=0,class=TRUE,param.estimates=TRUE)
}
\arguments{
  \item{data}{
	A matrix or 3-d array containing the data within which you wish to find a changepoint.  If data is a 3-d array, each first dimension is considered a separate dataset.  Within each dataset the first column is considered the response and the further columns are the covariates.  This function does not create a constant term so if one is required a column of ones should appear in the covariates.

NOTE: each model must have the same number of parameters otherwise singularities appear
}
  \item{penalty}{
	Choice of "None", "SIC", "BIC", "AIC", "Hannan-Quinn", "Asymptotic" and "Manual" penalties.  If Manual is specified, the manual penalty is contained in the value parameter. If Asymptotic is specified, the theoretical type I error is contained in the value parameter.
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
	This function is used to find a single change in regression for errors that are assumed to be normally distributed with zero mean.  The value returned is the result of testing H0:no change in regression against H1: single change in regression using the log of the likelihood ratio statistic coupled with the penalty supplied.
}
\value{
	If class=TRUE then an object of S4 class "cpt.reg" is returned.  The slot \code{cpts} contains the changepoints that are solely returned if class=FALSE.  If class=FALSE, the structure of \code{cpts} is as follows.

	If data is an nxp matrix (single dataset) then a single value is returned:
	\item{cpt}{The most probable location of a changepoint if H0 was rejected or n if H1 was rejected.}
	If data is an mxnxp matrix (multiple datasets) then a vector is returned:
	\item{cpt}{Vector of length m containing where each element is the result of the test for data[m,].  If cpt[m] is a number then it is the most probable location of a changepoint under H1.  Otherwise cpt[m] has the value n and indicates that H1 was rejected.}

}
\references{
Change in regression: Chen, J. and Gupta, A. K. (2000) \emph{Parametric statistical change point analysis}, Birkhauser
}
\author{
Rebecca Killick
}


\seealso{
\code{\link{cpt.reg}},\code{\link{single.mean.norm}},\code{\link{single.meanvar.norm}},\code{\link{single.var.norm}},\code{\linkS4class{cpt.reg}}
}
\examples{
# Example of a change in regression at 100 in simulated normal data
set.seed(1)
x=1:250
y=c(0.01*x[1:100],1.5-0.02*(x[101:250]-101))
ynoise=y+rnorm(250,0,0.2)
data=cbind(ynoise,1,x)
single.reg.norm(data,penalty="SIC",class=FALSE) # returns 100 to show that the null hypothesis was rejected and the change in regression is at 100
ans=single.reg.norm(data,penalty="Asymptotic",value=0.01) 
cpts(ans) # returns 100 to show that the null hypothesis was rejected, the change in regression is at 100 and we are 99% confident of this result


# Example of a data matrix containing 2 rows, row 1 has a change in regression and row 2 had no change in regression
set.seed(1)
x=1:250
y=c(0.01*x[1:100],1.5-0.02*(x[101:250]-101))
ynoise=y+rnorm(250,0,0.2)
data1=cbind(ynoise,1,x)
x=1:250
y=0.01*x+1
ynoise=y+rnorm(250,0,0.2)
data2=cbind(ynoise,1,x)

z=array(0,dim=c(2,250,3)) # note that the same model is fitted to all datasets
z[1,,]=data1
z[2,,]=data2
single.reg.norm(z,penalty="SIC",class=FALSE) # returns vector c(100,250) which shows that the first dataset has a change in regression at 100 and the second dataset rejected H1 and has no change in regression
ans=single.reg.norm(z,penalty="Manual",value="diffparam*log(n)") 
cpts(ans[[1]]) # test using a manual penalty which is the same as the SIC penalty for this example
cpts(ans[[2]]) # result is the same as above, c(100, 250) and 250.
}

\keyword{methods}
\keyword{univar}
\keyword{models}
\keyword{ts}
