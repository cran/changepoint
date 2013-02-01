\name{multiple.var.norm}
\alias{multiple.var.norm}
\title{
Multiple Changes in Variance - Normal Data
}
\description{
Calculates the optimal positioning and number of changepoints for Normal data using the user specified method.
}
\usage{
multiple.var.norm(data,mul.method="PELT",penalty="SIC",pen.value=0,Q=5,know.mean=FALSE,mu=NA,class=TRUE,param.estimates=TRUE)
}
\arguments{
  \item{data}{
	A vector, ts object or matrix containing the data within which you wish to find a changepoint.  If data is a matrix, each row is considered a separate dataset.
}
  \item{mul.method}{
	Choice of "PELT", "SegNeigh" or "BinSeg".
}
  \item{penalty}{
	Choice of "None", "SIC", "BIC", "AIC", "Hannan-Quinn", "Asymptotic" and "Manual" penalties.  If Manual is specified, the manual penalty is contained in the pen.value parameter. If Asymptotic is specified, the theoretical type I error is contained in the pen.value parameter.  The predefined penalties listed do NOT count the changepoint as a parameter, postfix a 1 e.g."SIC1" to count the changepoint as a parameter.
}
  \item{pen.value}{
	The theoretical type I error e.g.0.05 when using the Asymptotic penalty.  The value of the penalty when using the Manual penalty option.  This can be a numeric value or text giving the formula to use.  Available variables are, n=length of original data, null=null likelihood, alt=alternative likelihood, tau=proposed changepoint, diffparam=difference in number of alternatve and null parameters.
}
  \item{Q}{
	The maximum number of changepoints to search for using the "BinSeg" method.  The maximum number of segments (number of changepoints + 1) to search for using the "SegNeigh" method.
}
  \item{know.mean}{
	Logical, if TRUE then the mean is assumed known and mu is taken as its value.  If FALSE, and mu=-1000 (default value) then the mean is estimated via maximum likelihood.  If FALSE and the value of mu is supplied, mu is not estimated but is counted as an estimated parameter for decisions.
}
  \item{mu}{
	Numerical value of the true mean of the data.  Either single value or vector of length nrow(data).  If data is a matrix and mu is a single value, the same mean is used for each row.
}
  \item{class}{
	Logical.  If TRUE then an object of class \code{cpt} is returned.
}
  \item{param.estimates}{
	Logical.  If TRUE and class=TRUE then parameter estimates are returned. If FALSE or class=FALSE no parameter estimates are returned.
}
}
\details{
	This function is used to find multiple changes in variance for data that is assumed to be normally distributed.  The changes are found using the method supplied which can be exact (PELT or SegNeigh) or approximate (BinSeg).
}
\value{
	If class=TRUE then an object of S4 class "cpt" is returned.  The slot \code{cpts} contains the changepoints that are solely returned if class=FALSE.  The structure of \code{cpts} is as follows.

	If data is a vector (single dataset) then a vector/list is returned depending on the value of mul.method.  If data is a matrix (multiple datasets) then a list is returned where each element in the list is either a vector or list depending on the value of mul.method.

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
Change in variance: Chen, J. and Gupta, A. K. (2000) \emph{Parametric statistical change point analysis}, Birkhauser

PELT Algorithm: Killick, R. and Fearnhead, P. and Eckley, I.A. (2011) An exact linear time search algorithm for multiple changepoint detection, \emph{Submitted}

Binary Segmentation: Scott, A. J. and Knott, M. (1974) A Cluster Analysis Method for Grouping Means in the Analysis of Variance, \emph{Biometrics} \bold{30(3)}, 507--512

Segment Neighbourhoods: Auger, I. E. And Lawrence, C. E. (1989) Algorithms for the Optimal Identification of Segment Neighborhoods, \emph{Bulletin of Mathematical Biology} \bold{51(1)}, 39--54

}
\author{
Rebecca Killick
}


\seealso{
\code{\link{multiple.mean.norm}},\code{\link{multiple.meanvar.norm}},\code{\link{cpt.var}},\code{\link{PELT.var.norm}},\code{\link{binseg.var.norm}},\code{\link{single.var.norm}},\code{\link{segneigh.var.norm}},\code{\linkS4class{cpt}}
}
\examples{
# Example of multiple changes in variance at 50,100,150 in simulated normal data
set.seed(1)
x=c(rnorm(50,0,1),rnorm(50,0,10),rnorm(50,0,5),rnorm(50,0,1))
multiple.var.norm(x,mul.method="BinSeg",penalty="Manual",pen.value="2*log(n)",Q=5,class=FALSE) # returns optimal number of changepoints is 3, locations are 50,99,150.

# Example multiple datasets where the first row has multiple changes in variance and the second row has no change in variance
set.seed(10)
x=c(rnorm(50,0,1),rnorm(50,0,10),rnorm(50,0,5),rnorm(50,0,1))
y=rnorm(200,0,1)
z=rbind(x,y)
multiple.var.norm(z,mul.method="SegNeigh",penalty="Asymptotic",pen.value=0.01,Q=5,class=FALSE) # returns list that has two elements, the first has 3 changes in variance at 50,100,149 and the second has no changes in variance
ans=multiple.var.norm(z,mul.method="PELT",penalty="Asymptotic",pen.value=0.01)
cpts(ans[[1]]) # same results as for the SegNeigh method.
cpts(ans[[2]]) # same results as for the SegNeigh method.
}

\keyword{methods}
\keyword{univar}
\keyword{models}
\keyword{ts}
