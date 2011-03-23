\name{multiple.var.css}
\alias{multiple.var.css}
\title{
Multiple Changes in Variance - Cumulative Sums of Squares
}
\description{
Calculates the optimal positioning and number of changepoints for the cumulative sums of squares test statistic using the user specified method.
}
\usage{
multiple.var.css(data,mul.method="BinSeg",penalty="SIC",value=0,Q=5,class=TRUE,param.estimates=TRUE)
}
\arguments{
  \item{data}{
	A vector or matrix containing the data within which you wish to find a changepoint.  If data is a matrix, each row is considered a separate dataset.
}
  \item{mul.method}{
	Choice of "SegNeigh" or "BinSeg".
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
	This function is used to find multiple changes in variance for data where no assumption about the distribution is made.  The changes are found using the method supplied which can be exact (SegNeigh) or approximate (BinSeg).  Note that the penalty values are log(.) to be comparable with the distributional penalties.
}
\value{
	If class=TRUE then an object of S4 class "cpt" is returned.  The slot \code{cpts} contains the changepoints that are solely returned if class=FALSE.  The structure of \code{cpts} is as follows.

	If data is a vector (single dataset) then a vector/list is returned depending on the value of mul.method.  If data is a matrix (multiple datasets) then a list is returned where each element in the list is either a vector or list depending on the value of mul.method.

	If mul.method is SegNeigh then a list is returned with elements:
	\item{cps}{Matrix containing the changepoint positions for 1,...,Q changepoints.}
	\item{op.cpts}{The optimal changepoint locations for the penalty supplied.}
	\item{pen}{Penalty used to find the optimal number of changepoints.}	If mul.method is BinSeg then a list is returned with elements:
	\item{cps}{2xQ Matrix containing the changepoint positions on the first row and the test statistic on the second row.}
	\item{op.cpts}{The optimal changepoint locations for the penalty supplied.}
	\item{pen}{Penalty used to find the optimal number of changepoints.}
}
\references{
Binary Segmentation: Scott, A. J. and Knott, M. (1974) A Cluster Analysis Method for Grouping Means in the Analysis of Variance, \emph{Biometrics} \bold{30(3)}, 507--512

Segment Neighbourhoods: Auger, I. E. And Lawrence, C. E. (1989) Algorithms for the Optimal Identification of Segment Neighborhoods, \emph{Bulletin of Mathematical Biology} \bold{51(1)}, 39--54

C. Inclan, G. C. Tiao (1994) Use of Cumulative Sums of Squares for Retrospective Detection of Changes of Variance, \emph{Journal of the American Statistical Association} \bold{89(427)}, 913--923

R. L. Brown, J. Durbin, J. M. Evans (1975) Techniques for Testing the Constancy of Regression Relationships over Time, \emph{Journal of the Royal Statistical Society B} \bold{32(2)}, 149--192
}
\author{
Rebecca Killick
}


\seealso{
\code{\link{multiple.mean.cusum}},\code{\link{cpt.var}},\code{\link{binseg.var.css}},\code{\link{single.var.css}},\code{\link{segneigh.var.css}},\code{\linkS4class{cpt}}
}
\examples{
# Example of multiple changes in variance at 50,100,150 in simulated data
set.seed(1)
x=c(rnorm(50,0,1),rnorm(50,0,10),rnorm(50,0,5),rnorm(50,0,1))
multiple.var.css(x,mul.method="BinSeg",penalty="Manual",value="log(2*log(n))",Q=5,class=FALSE) # returns optimal number of changepoints is 4, locations are 50,53,99,150.

# Example multiple datasets where the first row has multiple changes in variance and the second row has no change in variance
set.seed(10)
x=c(rnorm(50,0,1),rnorm(50,0,10),rnorm(50,0,5),rnorm(50,0,1))
y=rnorm(200,0,1)
z=rbind(x,y)
multiple.var.css(z,mul.method="SegNeigh",penalty="Asymptotic",value=0.01,Q=5,class=FALSE) # returns list that has two elements, the first has 3 changes in variance at 52,100,149 and the second has no changes in variance
ans=multiple.var.css(z,mul.method="SegNeigh",penalty="Asymptotic",value=0.01)
cpts(ans[[1]]) # same results as for class=FALSE.
cpts(ans[[2]]) # same results as for class=FALSE.
}

\keyword{methods}
\keyword{univar}
\keyword{models}
\keyword{ts}
