\name{single.meanvar.exp.calc}
\alias{single.meanvar.exp.calc}
\title{
Single Change in Mean and Variance - Exponential Data
}
\description{
Calculates the scaled log-likelihood (assuming the data is Exponential distributed) for all possible changepoint locations and returns the single most probable (max).
}
\usage{
single.meanvar.exp.calc(data, extrainf = TRUE)
}
\arguments{
  \item{data}{
	A vector or matrix containing the data within which you wish to find a changepoint.  If data is a matrix, each row is considered a separate dataset.
}
  \item{extrainf}{
	Logical, if TRUE the scaled null and alternative likelihood values are returned along with the changepoint location.  If FALSE, only the changepoint location is returned.
}
}
\details{
	This function is used to find a single change in mean and variance for data that is assumed to be Exponential distributed.  The changepoint returned is simply the location where the log likelihood is maximised, there is no test performed as to whether this location is a true changepoint or not.

The returned likelihoods are scaled so that a test can be directly performed using the log of the likelihood ratio, \eqn{\lambda=\left\{null-alt\right\}}, which should be maximised.

In reality this function should not be used unless you are performing a changepoint test using the output supplied.  This function is used in the "see also" functions that perform various changepoint tests, ideally these should be used.
}
\value{
	If data is a vector (single dataset) and extrainf=FALSE then a single value is returned:
	\item{cpt}{The most probable location of a changepoint (scaled max log likelihood over all possible changepoint locations}
	If data is a vector (single dataset) and extrainf=TRUE then a vector with three elements is returned:
	\item{cpt}{The most probable location of a changepoint (scaled max log likelihood over all possible changepoint locations}
	\item{null}{The scaled null likelihood (log likelihood of entire data with no change)}
	\item{altlike}{The scaled alternative liklihood at cpt (log likelihood of entire data with a change at cpt)}
	If data is an mxn matrix (multiple datasets) and extrainf=FALSE then a vector is returned:
	\item{cpt}{Vector of length m containing the most probable location of a changepoint (scaled max log likelihood over all possible changepoint locations for each row in data.  cpt[1] is the most probable changepoint of the first row in data and cpt[m] is the most probable changepoint for the final row in data.}
	If data is a matrix (multiple datasets) and extrainf=TRUE then a matrix is returned where the first column is the changepoint location for each row in data, the second column is the scaled null likelihood for each row in data, the final column is the scaled maximum of the alternative likelihoods for each row in data.
}
\references{
Change in Exponential Model: Chen, J. and Gupta, A. K. (2000) \emph{Parametric statistical change point analysis}, Birkhauser
}
\author{
Rebecca Killick
}


\seealso{
\code{\link{single.meanvar.exp}}, \code{\link{cpt.meanvar}}
}
\examples{
# Example of a change in mean and variance at 100 in simulated exponential data
set.seed(10)
x=c(rexp(100,rate=1),rexp(100,rate=5))
single.meanvar.exp.calc(x,extrainf=FALSE) # finds change at 99
single.meanvar.exp.calc(x) # finds change at 99 and gives null likelihood as -192.5052 and alternative likelihood as -299.7263


# Example of no change in scale parameter (mean or variance) in simulated gamma data
set.seed(10)
x=rexp(100,rate=1)
single.meanvar.exp.calc(x,extrainf=FALSE) # finds change at 90, this is the most probable point of change but if a changepoint test is performed then no change will be found.
single.meanvar.exp.calc(x)# change at 90, null liklihood is 3.196690 and alternative is 1.060281
}

\keyword{methods}
\keyword{univar}
\keyword{models}
\keyword{ts}
