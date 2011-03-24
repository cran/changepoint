\name{cpt.meanvar}
\alias{cpt.meanvar}
\title{
Identifying Changes in Mean and Variance
}
\description{
Calculates the optimal positioning and (potentially) number of changepoints for data using the user specified method.
}
\usage{
cpt.meanvar(data,penalty="SIC",value=0,method="AMOC",Q=5,dist="Normal",class=TRUE,param.estimates=TRUE,shape=0)
}
\arguments{
  \item{data}{
	A vector or matrix containing the data within which you wish to find a changepoint.  If data is a matrix, each row is considered a separate dataset.
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
	The assumed distribution of the data.  Currently only "Normal", "Gamma" and "Exponential" are supported.
}
  \item{class}{
	Logical.  If TRUE then an object of class \code{cpt} is returned.
}
  \item{param.estimates}{
	Logical.  If TRUE and class=TRUE then parameter estimates are returned. If FALSE or class=FALSE no parameter estimates are returned.
}
  \item{shape}{
	Value of the assumed known shape parameter required when dist="Gamma".
}
}
\details{
	This function is used to find changes in mean and variance for data that is assumed to be distributed as the dist parameter.  The changes are found using the method supplied which can be single changepoint (AMOC) or multiple changepoints using exact (PELT or SegNeigh) or approximate (BinSeg) methods.
}
\value{
	If class=TRUE then an object of S4 class "cpt" is returned.  The slot \code{cpts} contains the changepoints that are solely returned if class=FALSE.  The structure of \code{cpts} is as follows.

	If data is a vector (single dataset) then a vector/list is returned depending on the value of method.  If data is a matrix (multiple datasets) then a list is returned where each element in the list is either a vector or list depending on the value of method.

	If method is AMOC then a single value (one dataset) or vector (multiple datasets) is returned:
	\item{cpt}{The most probable location of a changepoint if a change was identified or NA if no changepoint.}
	If method is PELT then a vector is returned:
	\item{cpt}{Vector containing the changepoint locations for the penalty supplied.  This always ends with n.}
	If method is SegNeigh then a list is returned with elements:
	\item{cps}{Matrix containing the changepoint positions for 1,...,Q changepoints.}
	\item{op.cpts}{The optimal changepoint locations for the penalty supplied.}
	\item{like}{Value of the -2*log(likelihood ratio) + penalty for the optimal number of changepoints selected.}
	If method is BinSeg then a list is returned with elements:
	\item{cps}{2xQ Matrix containing the changepoint positions on the first row and the test statistic on the second row.}
	\item{op.cpts}{The optimal changepoint locations for the penalty supplied.}
	\item{pen}{Penalty used to find the optimal number of changepoints.}
}
\references{
Change in Normal mean and variance: Chen, J. and Gupta, A. K. (2000) \emph{Parametric statistical change point analysis}, Birkhauser

Change in Gamma shape parameter: Chen, J. and Gupta, A. K. (2000) \emph{Parametric statistical change point analysis}, Birkhauser

Change in Exponential model: Chen, J. and Gupta, A. K. (2000) \emph{Parametric statistical change point analysis}, Birkhauser

PELT Algorithm: Killick, R. and Fearnhead, P. and Eckley, I.A. (2011) An exact linear time search algorithm for multiple changepoint detection, \emph{Submitted}

Binary Segmentation: Scott, A. J. and Knott, M. (1974) A Cluster Analysis Method for Grouping Means in the Analysis of Variance, \emph{Biometrics} \bold{30(3)}, 507--512

Segment Neighbourhoods: Auger, I. E. And Lawrence, C. E. (1989) Algorithms for the Optimal Identification of Segment Neighborhoods, \emph{Bulletin of Mathematical Biology} \bold{51(1)}, 39--54
}
\author{
Rebecca Killick
}


\seealso{
\code{\link{cpt.var}},\code{\link{cpt.mean}},\code{\link{plot-methods}},\code{\linkS4class{cpt}}
}
\examples{
# Example of a change in scale parameter (mean and variance) at 100 in simulated gamma data
set.seed(1)
x=c(rgamma(100,shape=1,rate=1),rgamma(100,shape=1,rate=5))
cpt.meanvar(x,penalty="SIC",method="AMOC",dist="Gamma",class=FALSE,shape=1) # returns 97 to show that the null hypothesis was rejected and the change in scale parameter is at 97
ans=cpt.meanvar(x,penalty="AIC",method="AMOC",dist="Gamma",shape=1)
cpts(ans) # returns 97 to show that the null hypothesis was rejected, the change in scale parameter is at 97

# Example of multiple changes in mean and variance at 50,100,150 in simulated normal data
set.seed(1)
x=c(rnorm(50,0,1),rnorm(50,5,3),rnorm(50,10,1),rnorm(50,3,10))
cpt.meanvar(x,penalty="Manual",value="4*log(n)",method="BinSeg",Q=5,class=FALSE) # returns optimal number of changepoints is 4, locations are 50,100,150,152.

# Example multiple datasets where the first row has multiple changes in mean and variance and the second row has no change in mean or variance
set.seed(1)
x=c(rnorm(50,0,1),rnorm(50,5,3),rnorm(50,10,1),rnorm(50,3,10))
y=rnorm(200,0,1)
z=rbind(x,y)
cpt.meanvar(z,penalty="Asymptotic",value=0.01,method="SegNeigh",Q=5,class=FALSE) # returns list that has two elements, the first has 3 changes in mean and variance at 50,100,150 and the second has no changes in mean or variance
ans=cpt.meanvar(z,penalty="Asymptotic",value=0.01,method="PELT") 
cpts(ans[[1]]) # same results as for the SegNeigh method.
cpts(ans[[2]]) # same results as for the SegNeigh method.
}

\keyword{methods}
\keyword{univar}
\keyword{models}
\keyword{ts}
