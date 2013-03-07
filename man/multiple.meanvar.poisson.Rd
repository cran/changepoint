\name{multiple.meanvar.poisson}
\alias{multiple.meanvar.poisson}
\title{
Multiple Changes in Mean and Variance - Poisson Data
}
\description{
Calculates the optimal positioning and number of changepoints for Poisson data using the user specified method.
}
\usage{
multiple.meanvar.poisson(data,mul.method="PELT",penalty="SIC",pen.value=0,Q=5,class=TRUE,param.estimates=TRUE)
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
  \item{class}{
	Logical.  If TRUE then an object of class \code{cpt} is returned.
}
  \item{param.estimates}{
	Logical.  If TRUE and class=TRUE then parameter estimates are returned. If FALSE or class=FALSE no parameter estimates are returned.
}
}
\details{
	This function is used to find multiple changes in mean and variance for data that is assumed to be Poisson distributed.  The changes are found using the method supplied which can be exact (PELT or SegNeigh) or approximate (BinSeg).
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
Change in Poisson model: Chen, J. and Gupta, A. K. (2000) \emph{Parametric statistical change point analysis}, Birkhauser

PELT Algorithm: Killick R, Fearnhead P, Eckley IA (2012) Optimal detection of changepoints with a linear computational cost, \emph{JASA} \bold{107(500)}, 1590--1598

Binary Segmentation: Scott, A. J. and Knott, M. (1974) A Cluster Analysis Method for Grouping Means in the Analysis of Variance, \emph{Biometrics} \bold{30(3)}, 507--512

Segment Neighbourhoods: Auger, I. E. And Lawrence, C. E. (1989) Algorithms for the Optimal Identification of Segment Neighborhoods, \emph{Bulletin of Mathematical Biology} \bold{51(1)}, 39--54
}
\author{
Rebecca Killick
}


\seealso{
\code{\link{cpt.meanvar}},\code{\link{PELT.meanvar.poisson}},\code{\link{binseg.meanvar.poisson}},\code{\link{single.meanvar.poisson}},\code{\link{segneigh.meanvar.poisson}},\code{\linkS4class{cpt}}
}
\examples{
# Example of multiple changes in mean and variance at 50,100,150 in simulated Poisson data
set.seed(1)
x=c(rpois(50,lambda=1),rpois(50,lambda=5),rpois(50,lambda=1),rpois(50,lambda=10))
multiple.meanvar.poisson(x,mul.method="BinSeg",penalty="Manual",pen.value="2*log(n)",Q=5,class=FALSE) # returns optimal number of changepoints is 3, locations are 50,100,150.

# Example multiple datasets where the first row has multiple changes in mean and variance and the second row has no change in mean or variance
set.seed(1)
x=c(rpois(50,lambda=1),rpois(50,lambda=5),rpois(50,lambda=1),rpois(50,lambda=10))
y=rpois(200,lambda=1)
z=rbind(x,y)
multiple.meanvar.poisson(z,mul.method="SegNeigh",penalty="Manual",pen.value=2*log(200),Q=5,class=FALSE) # returns list that has two elements, the first has 3 changes in mean and variance at 50,100,150 and the second has no changes in mean or variance
ans=multiple.meanvar.poisson(z,mul.method="PELT",penalty="Manual",pen.value=2*log(200)) 
cpts(ans[[1]]) # same results as for the SegNeigh method.
cpts(ans[[2]]) # same results as for the SegNeigh method.
}

\keyword{methods}
\keyword{univar}
\keyword{models}
\keyword{ts}
