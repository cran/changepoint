\name{cpt.reg}
\alias{cpt.reg}
\title{
Identifying Changes in Regression
}
\description{
Calculates the optimal positioning and (potentially) number of changepoints for data using the user specified method.
}
\usage{
cpt.reg(data,penalty="SIC",value=0,method="AMOC",Q=5,dist="Normal",class=TRUE,param.estimates=TRUE)
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
  \item{method}{
	Choice of "AMOC", "PELT", "SegNeigh" or "BinSeg".
}
  \item{Q}{
	The maximum number of changepoints to search for using the "BinSeg" method.  The maximum number of segments (number of changepoints + 1) to search for using the "SegNeigh" method.
}
  \item{dist}{
	The assumed distribution of the data.  Currently only "Normal" is supported.
}
  \item{class}{
	Logical.  If TRUE then an object of class \code{cpt} is returned.
}
  \item{param.estimates}{
	Logical.  If TRUE and class=TRUE then parameter estimates are returned. If FALSE or class=FALSE no parameter estimates are returned.
}
}
\details{
	This function is used to find changes in regression for data that is assumed to be distributed as the dist parameter.  The changes are found using the method supplied which can be single changepoint (AMOC) or multiple changepoints using exact (PELT or SegNeigh) or approximate (BinSeg) methods.
}
\value{
	If class=TRUE then an object of S4 class "cpt" is returned.  The slot \code{cpts} contains the changepoints that are solely returned if class=FALSE.  The structure of \code{cpts} is as follows.

	If data is an nxp matrix (single dataset) then a single value is returned:
	\item{cpt}{The most probable location of a changepoint if H0 was rejected or n if H1 was rejected.}
	If data is an mxnxp matrix (multiple datasets) then a vector is returned:
	\item{cpt}{Vector of length m containing where each element is the result of the test for data[m,,].  If cpt[m] is less than n then it is the most probable location of a changepoint under H1.  Otherwise cpt[m] has the value n and indicates that H1 was rejected.}
}
\references{
Change in regression: Chen, J. and Gupta, A. K. (2000) \emph{Parametric statistical change point analysis}, Birkhauser
}
\author{
Rebecca Killick
}


\seealso{
\code{\link{cpt.var}},\code{\link{cpt.meanvar}},\code{\link{cpt.mean}},\code{\linkS4class{cpt.reg}}
}
\examples{
# Example of a change in regression at 100 in simulated data with zero-mean normal errors
set.seed(1)
x=1:250
y=c(0.01*x[1:100],1.5-0.02*(x[101:250]-101))
ynoise=y+rnorm(250,0,0.2)
data=cbind(ynoise,1,x)
cpt.reg(data,penalty="SIC",class=FALSE) # returns 100 to show that the null hypothesis was rejected and the change in regression is at 100
ans=cpt.reg(data,penalty="Asymptotic",value=0.01) 
cpts(ans) # returns 100 to show that the null hypothesis was rejected, the change in regression is at 100 and we are 99% confident of this result

# Example of multiple changes in regression at 100,250 in simulated data with zero-mean normal errors
set.seed(1)
x=1:400
y=c(0.01*x[1:100],3.5-0.02*x[101:250],-15+0.05*x[251:400])
ynoise=y+rnorm(400,0,0.2)
yx=cbind(ynoise,1,x)
cpt.reg(yx,method="BinSeg",penalty="Manual",value="4*log(n)",Q=5,class=FALSE) # returns optimal number of changepoints is 2, locations are 100,250.

# Example multiple datasets where the first has multiple changes in regression and the second has no change in regression
set.seed(1)
x1=1:400
y1=c(0.01*x1[1:100],3.5-0.02*x1[101:250],-15+0.05*x1[251:400])
ynoise1=y1+rnorm(400,0,0.2)
yx1=cbind(ynoise1,1,x1)

x2=1:400
y2=0.01*x2
ynoise2=y2+rnorm(400,0,0.2)
yx2=cbind(ynoise2,1,x2)

data=array(0,dim=c(2,400,3))
data[1,,]=yx1; data[2,,]=yx2

cpt.reg(data,method="SegNeigh",penalty="Asymptotic",value=0.01,Q=5,class=FALSE) # returns list that has two elements, the first has 2 changes in regression at 100,250 and the second has no changes in regression
ans=cpt.reg(data,method="PELT",penalty="Asymptotic",value=0.01) 
cpts(ans[[1]]) # same results as for the SegNeigh method.
cpts(ans[[2]]) # same results as for the SegNeigh method.
}

\keyword{methods}
\keyword{univar}
\keyword{models}
\keyword{ts}
