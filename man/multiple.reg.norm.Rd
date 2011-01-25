\name{multiple.reg.norm}
\alias{multiple.reg.norm}
\title{
Multiple Changes in Regression - Normal Errors
}
\description{
Calculates the optimal positioning and number of changepoints for regression data with zero-mean Normal errors using the user specified method.
}
\usage{
multiple.reg.norm(data,mul.method="PELT",penalty="SIC",value=0,Q=5,class=TRUE,param.estimates=TRUE)
}
\arguments{
  \item{data}{
	A matrix or 3-d array containing the data within which you wish to find a changepoint.  If data is a 3-d array, each first dimension is considered a separate dataset.  Within each dataset the first column is considered the response and the further columns are the covariates.  This function does not create a constant term so if one is required a column of ones should appear in the covariates.

NOTE: each model must have the same number of parameters otherwise singularities appear
}
  \item{mul.method}{
	Choice of "PELT", "SegNeigh" or "BinSeg".
}
  \item{penalty}{
	Choice of "None", "SIC", "BIC", "AIC", "Hannan-Quinn", "Asymptotic" and "Manual" penalties.  If Manual is specified, the manual penalty is contained in the value parameter. If Asymptotic is specified, the theoretical type I error is contained in the value parameter.
}
  \item{value}{
	The theoretical type I error e.g.0.05 when using the Asymptotic penalty.  The value of the penalty when using the Manual penalty option.  This can be a numeric value or text giving the formula to use.  Available variables are, n=length of original data, null=null likelihood, alt=alternative likelihood, tau=proposed changepoint, diffparam=difference in number of alternatve and null parameters.
}
  \item{Q}{
	The maximum number of changepoints to search for using the "BinSeg" method.  The maximum number of segments (number of changepoints + 1) to search for using the "SegNeigh" method.
}
  \item{class}{
	Logical.  If TRUE then an object of class \code{cpt} is returned.
}
  \item{param.estimates}{
	Logical.  If TRUE and class=TRUE then parameter estimates are returned. If FALSE or class=FALSE no parameter estimates are returned.
}
}
\details{
	This function is used to find multiple changes in regression for data that is assumed to have zero-mean normally distributed errors.  The changes are found using the method supplied which can be exact (PELT or SegNeigh) or approximate (BinSeg).
}
\value{
	If class=TRUE then an object of S4 class "cpt" is returned.  The slot \code{cpts} contains the changepoints that are solely returned if class=FALSE.  The structure of \code{cpts} is as follows.
	
	If data is a matrix (single dataset) then a vector/list is returned depending on the value of mul.method.  If data is an array (multiple datasets) then a list is returned where each element in the list is either a vector or list depending on the value of mul.method.

	If mul.method is PELT then a vector is returned:
	\item{cpt}{Vector containing the changepoint locations for the penalty supplied.  This always ends with n.}
	If mul.method is SegNeigh then a list is returned with elements:
	\item{cps}{Matrix containing the changepoint positions for 1,...,Q changepoints.}
	\item{op.cpts}{The optimal changepoint locations for the penalty supplied.}
	\item{like}{Value of the -2*log(likelihood ratio) + penalty for the optimal number of changepoints selected.}
	If mul.method is BinSeg then a list is returned with elements:
	\item{cps}{2xQ Matrix containing the changepoint positions on the first row and the test statistic on the second row.}
	\item{op.cpts}{The optimal changepoint locations for the penalty supplied.}
	\item{pen}{Penalty used to find the optimal number of changepoints.}
}
\references{
Change in regression: Chen, J. and Gupta, A. K. (2000) \emph{Parametric statistical change point analysis}, Birkhauser

PELT Algorithm: Killick, R. and Fearnhead, P. and Eckley, I.A. (2011) An exact linear time search algorithm for multiple changepoint detection, \emph{Submitted}

Binary Segmentation: Scott, A. J. and Knott, M. (1974) A Cluster Analysis Method for Grouping Means in the Analysis of Variance, \emph{Biometrics} \bold{30(3)}, 507--512

Segment Neighbourhoods: Auger, I. E. And Lawrence, C. E. (1989) Algorithms for the Optimal Identification of Segment Neighborhoods, \emph{Bulletin of Mathematical Biology} \bold{51(1)}, 39--54
}
\author{
Rebecca Killick
}


\seealso{
\code{\link{multiple.mean.norm}},\code{\link{multiple.var.norm}},\code{\link{multiple.meanvar.norm}},\code{\link{cpt.reg}},\code{\link{PELT.reg.norm}},\code{\link{binseg.reg.norm}},\code{\link{single.reg.norm}},\code{\link{segneigh.reg.norm}},\code{\linkS4class{cpt}}
}
\examples{
# Example of multiple changes in regression at 100,250 in simulated data with zero-mean normal errors
set.seed(1)
x=1:400
y=c(0.01*x[1:100],3.5-0.02*x[101:250],-15+0.05*x[251:400])
ynoise=y+rnorm(400,0,0.2)
yx=cbind(ynoise,1,x)
multiple.reg.norm(yx,mul.method="BinSeg",penalty="Manual",value="4*log(n)",Q=5,class=FALSE) # returns optimal number of changepoints is 2, locations are 100,250.

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

multiple.reg.norm(data,mul.method="SegNeigh",penalty="Asymptotic",value=0.01,Q=5,class=FALSE) # returns list that has two elements, the first has 2 changes in regression at 100,250 and the second has no changes in regression
ans=multiple.reg.norm(data,mul.method="PELT",penalty="Asymptotic",value=0.01) 
cpts(ans[[1]]) # same results as for the SegNeigh method.
cpts(ans[[2]]) # same results as for the SegNeigh method.
}

\keyword{methods}
\keyword{univar}
\keyword{models}
\keyword{ts}
