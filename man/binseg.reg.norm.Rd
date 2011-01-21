\name{binseg.reg.norm}
\alias{binseg.reg.norm}
\title{
Multiple Changes in Regression using Binary Segmentation method - Normal Errors
}
\description{
Calculates the optimal positioning and number of changepoints for regression data with zero-mean Normal errors using Binary Segmentation method.  Note that this is an approximate method.
}
\usage{
binseg.reg.norm(data, Q=5, pen=0)
}
\arguments{
  \item{data}{
	A matrix or 3-d array containing the data within which you wish to find a changepoint.  If data is a 3-d array, each first dimension is considered a separate dataset.  Within each dataset the first column is considered the response and the further columns are the covariates.  This function does not create a constant term so if one is required a column of ones should appear in the covariates.

NOTE: each model must have the same number of parameters otherwise singularities appear
}
  \item{Q}{
	Numeric value of the maximum number of changepoints you wish to search for, default is 5.
}
  \item{pen}{
	Numeric value of the linear penalty function.  This value is used in the decision as to the optimal number of changepoints.
}
}
\details{
	This function is used to find a multiple changes in regression for data that is assumed to have zero-mean normally distributed errors.  The value returned is the result of finding the optimal location of up to Q changepoints using the log of the likelihood ratio statistic.  Once all changepoint locations have been calculated, the optimal number of changepoints is decided using pen as the penalty function.
}
\value{
	A list is returned containing the following items
	\item{cps}{2xQ Matrix containing the changepoint positions on the first row and the test statistic on the second row.}
	\item{op.cpts}{The optimal changepoint locations for the penalty supplied.}
	\item{pen}{Penalty used to find the optimal number of changepoints.}
}
\references{
Binary Segmentation: Scott, A. J. and Knott, M. (1974) A Cluster Analysis Method for Grouping Means in the Analysis of Variance, \emph{Biometrics} \bold{30(3)}, 507--512

Change in regression: Chen, J. and Gupta, A. K. (2000) \emph{Parametric statistical change point analysis}, Birkhauser
}
\author{
Rebecca Killick
}


\seealso{
\code{\link{binseg.mean.norm}},\code{\link{binseg.var.norm}},\code{\link{binseg.meanvar.norm}},\code{\link{cpt.reg}},\code{\link{PELT.reg.norm}},\code{\link{multiple.reg.norm}},\code{\link{single.reg.norm}},\code{\link{segneigh.reg.norm}}
}
\examples{
# Example of multiple changes in regression at 100,250 in simulated data with zero-mean normal errors
set.seed(1)
x=1:400
y=c(0.01*x[1:100],3.5-0.02*x[101:250],-15+0.05*x[251:400])
ynoise=y+rnorm(400,0,0.2)
yx=cbind(ynoise,1,x)
binseg.reg.norm(yx,Q=5,pen=4*log(400)) # returns optimal number as 2 and the locations as c(100,250)
binseg.reg.norm(yx,Q=1,pen=4*log(400)) # returns optimal number as 1 as this is the maximum number of changepoints it can find.  If you get the maximum number, you need to increase Q until this is not the case.

# Example no change in regression
set.seed(10)
x=1:400
y=0.01*x
ynoise=y+rnorm(400,0,0.2)
yx=cbind(ynoise,1,x)
binseg.reg.norm(yx,Q=5,pen=4*log(400)) # returns optimal number as 0
}

\keyword{methods}
\keyword{univar}
\keyword{models}
\keyword{ts}
