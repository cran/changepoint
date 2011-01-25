\name{PELT.reg.norm}
\alias{PELT.reg.norm}
\title{
Multiple Changes in Regression using PELT pruned method - Normal Errors
}
\description{
Calculates the optimal positioning and number of changepoints for regression data with zero-mean Normal errors using PELT pruned method.
}
\usage{
PELT.reg.norm(data, pen=0, nprune=FALSE)
}
\arguments{
  \item{data}{
	A matrix or 3-d array containing the data within which you wish to find a changepoint.  If data is a 3-d array, each first dimension is considered a separate dataset.  Within each dataset the first column is considered the response and the further columns are the covariates.  This function does not create a constant term so if one is required a column of ones should appear in the covariates.

NOTE: each model must have the same number of parameters otherwise singularities appear
}
  \item{pen}{
	Numeric value of the linear penalty function.  This value is used in the decision for each individual changepoint so that in total the penalty is k*pen where k is the optimal number of changepoints detected.
}
  \item{nprune}{
	Logical.  If nprune=FALSE, the locations of the changepoints are returned.  If nprune=TRUE then the number of points kept at each stage of pruning is returned instead.
}
}
\details{
	This function is used to find a multiple changes in regression for data that is assumed to have zero-mean normally distributed errors.  The value returned is the result of testing H0:existing number of changepoints against H1: one extra changepoint using the log of the likelihood ratio statistic coupled with the penalty supplied.  The PELT method keeps track of the optimal number and location of changepoints as it passes through the data.
}
\value{
	If nprune=FALSE then a vector of the changepoint locations is returned:
	\item{cpt}{Vector containing the changepoint locations for the penalty supplied.  This always ends with n.}
	If nprune=TRUE then a vector of the number of points kept during pruning at each iteration is returned:
	\item{nprune}{Vector containing the number of points kept at each stage of pruning.}
}
\references{
Change in regression: Chen, J. and Gupta, A. K. (2000) \emph{Parametric statistical change point analysis}, Birkhauser

PELT Algorithm: Killick, R. and Fearnhead, P. and Eckley, I.A. (2011) An exact linear time search algorithm for multiple changepoint detection, \emph{Submitted}

}
\author{
Rebecca Killick
}


\seealso{
\code{\link{PELT.mean.norm}},\code{\link{PELT.var.norm}},\code{\link{PELT.meanvar.norm}},\code{\link{cpt.reg}},\code{\link{binseg.reg.norm}},\code{\link{multiple.reg.norm}},\code{\link{single.reg.norm}},\code{\link{segneigh.reg.norm}}
}
\examples{
# Example of multiple changes in regression at 100,250 in simulated data with zero-mean normal errors
set.seed(1)
x=1:400
y=c(0.01*x[1:100],3.5-0.02*x[101:250],-15+0.05*x[251:400])
ynoise=y+rnorm(400,0,0.2)
yx=cbind(ynoise,1,x)
PELT.reg.norm(yx,pen=4*log(400)) # returns optimal number as 2 and the locations as c(100,249)
PELT.reg.norm(yx,pen=4*log(400),nprune=TRUE) # returns vector showing how many points are kept during pruing at each iteration of the method

# Example no change in regression
set.seed(10)
x=1:400
y=0.01*x
ynoise=y+rnorm(400,0,0.2)
yx=cbind(ynoise,1,x)
PELT.reg.norm(yx,pen=4*log(400)) # returns 400 to show no changepoint were found
PELT.reg.norm(yx,pen=4*log(200),nprune=TRUE) # notice how the number of points kept it steadily increasing compared to the previous example where it almost resets when a true change has been found
}

\keyword{methods}
\keyword{univar}
\keyword{models}
\keyword{ts}
