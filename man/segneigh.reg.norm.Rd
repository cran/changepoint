\name{segneigh.reg.norm}
\alias{segneigh.reg.norm}
\title{
Multiple Changes in Regression using Segment Neighbourhood method - Normal Errors
}
\description{
Calculates the optimal positioning and number of changepoints for regression data with zero-mean Normal errors using Segment Neighbourhood method.  Note that this gives the same results as PELT method but takes more computational time.
}
\usage{
segneigh.reg.norm(data, Q=5, pen=0)
}
\arguments{
  \item{data}{
	A matrix or 3-d array containing the data within which you wish to find a changepoint.  If data is a 3-d array, each first dimension is considered a separate dataset.  Within each dataset the first column is considered the response and the further columns are the covariates.  This function does not create a constant term so if one is required a column of ones should appear in the covariates.

NOTE: each model must have the same number of parameters otherwise singularities appear
}
  \item{Q}{
	Numeric value of the maximum number of segments (number of changepoints +1) you wish to search for, default is 5.
}
  \item{pen}{
	Numeric value of the linear penalty function.  This value is used in the final decision as to the optimal number of changepoints, used as k*pen where k is the number of changepoints to be tested.
}
}
\details{
	This function is used to find a multiple changes in regression for data that is assumed to have zero-mean normally distributed errors.  The value returned is the result of finding the optimal location of up to Q changepoints using the log of the likelihood ratio statistic.  Once all changepoint locations have been calculated, the optimal number of changepoints is decided using k*pen as the penalty function where k is the number of changepoints tested (k in (1,Q)).
}
\value{
	A list is returned containing the following items
	\item{cps}{Matrix containing the changepoint positions for 1,...,Q changepoints.}
	\item{op.cpts}{The optimal changepoint locations for the penalty supplied.}
	\item{like}{Value of the -2*log(likelihood ratio) + penalty for the optimal number of changepoints selected.}
}
\references{
Change in regression: Chen, J. and Gupta, A. K. (2000) \emph{Parametric statistical change point analysis}, Birkhauser

Segment Neighbourhoods: Auger, I. E. And Lawrence, C. E. (1989) Algorithms for the Optimal Identification of Segment Neighborhoods, \emph{Bulletin of Mathematical Biology} \bold{51(1)}, 39--54
}
\author{
Rebecca Killick
}


\seealso{
\code{\link{segneigh.mean.norm}},\code{\link{segneigh.var.norm}},\code{\link{segneigh.meanvar.norm}},\code{\link{cpt.reg}},\code{\link{PELT.reg.norm}},\code{\link{multiple.reg.norm}},\code{\link{single.reg.norm}},\code{\link{binseg.reg.norm}}
}
\examples{
# Example of multiple changes in regression at 100,250 in simulated data with zero-mean normal errors
set.seed(1)
x=1:400
y=c(0.01*x[1:100],3.5-0.02*x[101:250],-15+0.05*x[251:400])
ynoise=y+rnorm(400,0,0.2)
yx=cbind(ynoise,1,x)
segneigh.reg.norm(yx,Q=5,pen=4*log(400)) # returns optimal number as 2 and the locations as c(100,250)
segneigh.reg.norm(yx,Q=2,pen=4*log(400)) # returns optimal number as 1 as this is the maximum number of changepoints it can find.  If you get the maximum number, you need to increase Q until this is not the case.

# Example no change in regression
set.seed(10)
x=1:400
y=0.01*x
ynoise=y+rnorm(400,0,0.2)
yx=cbind(ynoise,1,x)
segneigh.reg.norm(yx,Q=5,pen=2*log(400)) # returns optimal number as 0
}

\keyword{methods}
\keyword{univar}
\keyword{models}
\keyword{ts}
