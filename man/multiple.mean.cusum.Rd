\name{multiple.mean.cusum}
\alias{multiple.mean.cusum}
\title{
Multiple Changes in Mean - Cumulative Sums
}
\description{
Calculates the optimal positioning and number of changepoints for the cumulative sums test statistic using the user specified method.
}
\usage{
multiple.mean.cusum(data,mul.method="BinSeg",penalty="Asymptotic",value=0.05,Q=5,class=TRUE,param.estimates=TRUE)
}
\arguments{
  \item{data}{
	A vector or matrix containing the data within which you wish to find a changepoint.  If data is a matrix, each row is considered a separate dataset.
}
  \item{mul.method}{
	Choice of "SegNeigh" or "BinSeg".
}
  \item{penalty}{
	Choice of "None", "SIC", "BIC", "AIC", "Hannan-Quinn" and "Manual" penalties.  If Manual is specified, the manual penalty is contained in the value parameter.  The predefined penalties listed do NOT count the changepoint as a parameter, postfix a 1 e.g."SIC1" to count the changepoint as a parameter.
}
  \item{value}{
	The value of the penalty when using the Manual penalty option.  This can be a numeric value or text giving the formula to use.  Available variables are, n=length of original data, null=null likelihood, alt=alternative likelihood, tau=proposed changepoint, diffparam=difference in number of alternatve and null parameters.
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
	This function is used to find multiple changes in mean for data where no assumption about the distribution is made.  The changes are found using the method supplied which can be exact (SegNeigh) or approximate (BinSeg).  Note that the programmed penalty values are not designed to be used with the CUSUM method, it is advised to use Asymptotic or Manual penalties.
}
\value{
	If class=TRUE then an object of S4 class "cpt" is returned.  The slot \code{cpts} contains the changepoints that are solely returned if class=FALSE.  The structure of \code{cpts} is as follows.

	If data is a vector (single dataset) then a vector/list is returned depending on the value of mul.method.  If data is a matrix (multiple datasets) then a list is returned where each element in the list is either a vector or list depending on the value of mul.method.

	If mul.method is SegNeigh then a list is returned with elements:
	\item{cps}{Matrix containing the changepoint positions for 1,...,Q changepoints.}
	\item{op.cpts}{The optimal changepoint locations for the penalty supplied.}
	\item{pen}{Penalty used to find the optimal number of changepoints.}
	If mul.method is BinSeg then a list is returned with elements:
	\item{cps}{2xQ Matrix containing the changepoint positions on the first row and the test statistic on the second row.}
	\item{op.cpts}{The optimal changepoint locations for the penalty supplied.}
	\item{pen}{Penalty used to find the optimal number of changepoints.}
}
\references{
Binary Segmentation: Scott, A. J. and Knott, M. (1974) A Cluster Analysis Method for Grouping Means in the Analysis of Variance, \emph{Biometrics} \bold{30(3)}, 507--512

Segment Neighbourhoods: Auger, I. E. And Lawrence, C. E. (1989) Algorithms for the Optimal Identification of Segment Neighborhoods, \emph{Bulletin of Mathematical Biology} \bold{51(1)}, 39--54

M. Csorgo, L. Horvath (1997) Limit Theorems in Change-Point Analysis, \emph{Wiley}

E. S. Page (1954) Continuous Inspection Schemes, \emph{Biometrika} \bold{41(1/2)}, 100--115
}
\author{
Rebecca Killick
}


\seealso{
\code{\link{multiple.var.css}},\code{\link{cpt.mean}},\code{\link{binseg.mean.cusum}},\code{\link{single.mean.cusum}},\code{\link{segneigh.mean.cusum}},\code{\linkS4class{cpt}}
}
\examples{
# Example of multiple changes in mean at 50,100,150 in simulated data
set.seed(1)
x=c(rnorm(50,0,1),rnorm(50,5,1),rnorm(50,10,1),rnorm(50,3,1))
multiple.mean.cusum(x,mul.method="BinSeg",penalty="Manual",value=0.8,Q=5,class=FALSE) # returns optimal number of changepoints is 3, locations are 50,100,150.

# Example multiple datasets where the first row has multiple changes in mean and the second row has no change in mean
set.seed(1)
x=c(rnorm(50,0,1),rnorm(50,5,1),rnorm(50,10,1),rnorm(50,3,1))
y=rnorm(200,0,1)
z=rbind(x,y)
multiple.mean.cusum(z,mul.method="SegNeigh",penalty="Manual",value=1,Q=5,class=FALSE) # returns list that has two elements, the first has 3 changes in mean at 50,101,150 and the second has no changes in mean
ans=multiple.mean.cusum(z,mul.method="BinSeg",penalty="Manual",value=0.8) 
cpts(ans[[1]]) # same results as for the SegNeigh method.
cpts(ans[[2]]) # same results as for the SegNeigh method.
}

\keyword{methods}
\keyword{univar}
\keyword{models}
\keyword{ts}
