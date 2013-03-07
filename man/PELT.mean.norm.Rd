\name{PELT.mean.norm}
\alias{PELT.mean.norm}
\title{
Multiple Changes in Mean using PELT pruned method - Normal Data
}
\description{
Calculates the optimal positioning and number of changepoints for Normal data using PELT pruned method.
}
\usage{
PELT.mean.norm(data, pen=0)
}
\arguments{
  \item{data}{
	A vector containing the data within which you wish to find changepoints.
}
  \item{pen}{
	Numeric value of the linear penalty function.  This value is used in the decision for each individual changepoint so that in total the penalty is k*pen where k is the optimal number of changepoints detected.
}
}
\details{
	This function is used to find a multiple changes in mean for data that is assumed to be normally distributed.  The value returned is the result of testing H0:existing number of changepoints against H1: one extra changepoint using the log of the likelihood ratio statistic coupled with the penalty supplied.  The PELT method keeps track of the optimal number and location of changepoints as it passes through the data.
}
\value{
	A vector of the changepoint locations is returned:
	\item{cpt}{Vector containing the changepoint locations for the penalty supplied.  This always ends with n.}
}
\references{
Change in Normal mean: Hinkley, D. V. (1970) Inference About the Change-Point in a Sequence of Random Variables, \emph{Biometrika} \bold{57}, 1--17

PELT Algorithm: Killick R, Fearnhead P, Eckley IA (2012) Optimal detection of changepoints with a linear computational cost, \emph{JASA} \bold{107(500)}, 1590--1598

}
\author{
Rebecca Killick
}


\seealso{
\code{\link{PELT.var.norm}},\code{\link{PELT.meanvar.norm}},\code{\link{cpt.mean}},\code{\link{binseg.mean.norm}},\code{\link{multiple.mean.norm}},\code{\link{single.mean.norm}},\code{\link{segneigh.mean.norm}}
}
\examples{
# Example of multiple changes in mean at 50,100,150 in simulated normal data
set.seed(1)
x=c(rnorm(50,0,1),rnorm(50,5,1),rnorm(50,10,1),rnorm(50,3,1))
PELT.mean.norm(x,pen=2*log(200)) # returns c(50,100,150,200)

# Example no change in mean
set.seed(10)
x=rnorm(200,0,1)
PELT.mean.norm(x,pen=2*log(200)) # returns 200 to show no change in mean has been found
}

\keyword{methods}
\keyword{univar}
\keyword{models}
\keyword{ts}
