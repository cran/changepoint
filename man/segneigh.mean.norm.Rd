\name{segneigh.mean.norm}
\alias{segneigh.mean.norm}
\title{
Multiple Changes in Mean using Segment Neighbourhood method - Normal Data
}
\description{
Calculates the optimal positioning and number of changepoints for Normal data using Segment Neighbourhood method.  Note that this gives the same results as PELT method but takes more computational time.
}
\usage{
segneigh.mean.norm(data, Q=5, pen=0)
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
}
\details{
	This function is used to find a multiple changes in mean for data that is assumed to be normally distributed.  The value returned is the result of finding the optimal location of up to Q changepoints using the log of the likelihood ratio statistic.  Once all changepoint locations have been calculated, the optimal number of changepoints is decided using k*pen as the penalty function where k is the number of changepoints tested (k in (1,Q)).
}
\value{
	A list is returned containing the following items
	\item{cps}{Matrix containing the changepoint positions for 1,...,Q changepoints.}
	\item{op.cpts}{The optimal changepoint locations for the penalty supplied.}
	\item{like}{Value of the -2*log(likelihood ratio) + penalty for the optimal number of changepoints selected.}
}
\references{
Change in Normal mean: Hinkley, D. V. (1970) Inference About the Change-Point in a Sequence of Random Variables, \emph{Biometrika} \bold{57}, 1--17

Segment Neighbourhoods: Auger, I. E. And Lawrence, C. E. (1989) Algorithms for the Optimal Identification of Segment Neighborhoods, \emph{Bulletin of Mathematical Biology} \bold{51(1)}, 39--54
}
\author{
Rebecca Killick
}


\seealso{
\code{\link{segneigh.var.norm}},\code{\link{segneigh.meanvar.norm}},\code{\link{cpt.mean}},\code{\link{PELT.mean.norm}},\code{\link{multiple.mean.norm}},\code{\link{single.mean.norm}},\code{\link{binseg.mean.norm}}
}
\examples{
# Example of multiple changes in mean at 50,100,150 in simulated normal data
set.seed(1)
x=c(rnorm(50,0,1),rnorm(50,5,1),rnorm(50,10,1),rnorm(50,3,1))
segneigh.mean.norm(x,Q=5,pen=2*log(200)) # returns optimal number as 3 and the locations as c(50,100,150)
segneigh.mean.norm(x,Q=3,pen=2*log(200)) # returns optimal number as 2 as this is the maximum number of changepoints it can find.  If you get the maximum number, you need to increase Q until this is not the case.

# Example no change in mean
set.seed(10)
x=rnorm(200,0,1)
segneigh.mean.norm(x,Q=5,pen=2*log(200)) # returns optimal number as 0
}

\keyword{methods}
\keyword{univar}
\keyword{models}
\keyword{ts}
