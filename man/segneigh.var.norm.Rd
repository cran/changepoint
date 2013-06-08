\name{segneigh.var.norm}
\alias{segneigh.var.norm}
\title{
Multiple Changes in Variance using Segment Neighbourhood method - Normal Data
}
\description{
Calculates the optimal positioning and number of changepoints for Normal data using Segment Neighbourhood method.  Note that this gives the same results as PELT method but takes more computational time.
}
\usage{
segneigh.var.norm(data, Q=5, pen=0, know.mean=FALSE, mu=NA)
}
\arguments{
  \item{data}{
	A vector containing the data within which you wish to find changepoints.
}
  \item{Q}{
	Numeric value of the maximum number of segments (number of changepoints +1) you wish to search for, default is 5.
}
  \item{pen}{
	Numeric value of the linear penalty function.  This value is used in the final decision as to the optimal number of changepoints, used as k*pen where k is the number of changepoints to be tested.
}
  \item{know.mean}{
	Logical, if TRUE then the mean is assumed known and mu is taken as its value.  If FALSE, and mu=-1000 (default value) then the mean is estimated via maximum likelihood.  If FALSE and the value of mu is supplied, mu is not estimated but is counted as an estimated parameter for decisions.
}
  \item{mu}{
	Numerical value of the true mean of the data.  Either single value or vector of length nrow(data).  If data is a matrix and mu is a single value, the same mean is used for each row.
}
}
\details{
	This function is used to find a multiple changes in variance for data that is assumed to be normally distributed.  The value returned is the result of finding the optimal location of up to Q changepoints using the log of the likelihood ratio statistic.  Once all changepoint locations have been calculated, the optimal number of changepoints is decided using k*pen as the penalty function where k is the number of changepoints tested (k in (1,Q)).
}
\value{
	A list is returned containing the following items
	\item{cps}{Matrix containing the changepoint positions for 1,...,Q changepoints.}
	\item{op.cpts}{The optimal changepoint locations for the penalty supplied.}
	\item{like}{Value of the -2*log(likelihood ratio) + penalty for the optimal number of changepoints selected.}
}
\references{
Change in variance: Chen, J. and Gupta, A. K. (2000) \emph{Parametric statistical change point analysis}, Birkhauser

Segment Neighbourhoods: Auger, I. E. And Lawrence, C. E. (1989) Algorithms for the Optimal Identification of Segment Neighborhoods, \emph{Bulletin of Mathematical Biology} \bold{51(1)}, 39--54
}
\author{
Rebecca Killick
}


\seealso{
\code{\link{segneigh.mean.norm}},\code{\link{segneigh.meanvar.norm}},\code{\link{cpt.var}},\code{\link{PELT.var.norm}},\code{\link{multiple.var.norm}},\code{\link{single.var.norm}},\code{\link{binseg.var.norm}}
}
\examples{
# Example of multiple changes in variance at 50,100,150 in simulated normal data
set.seed(1)
x=c(rnorm(50,0,1),rnorm(50,0,10),rnorm(50,0,5),rnorm(50,0,1))
segneigh.var.norm(x,Q=5, pen=2*log(200)) # returns optimal number as 3 and the locations as
#c(50,99,150)
segneigh.var.norm(x,Q=3, pen=2*log(200)) # returns optimal number as 2 as this is the maximum number
#of changepoints it can find.  If you get the maximum number, you need to increase Q until this is
#not the case.

# Example no change in variance
set.seed(10)
x=rnorm(200,0,1)
segneigh.var.norm(x,Q=5, pen=2*log(200)) # returns optimal number as 0
}

\keyword{methods}
\keyword{univar}
\keyword{models}
\keyword{ts}
