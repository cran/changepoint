\name{single.var.norm.calc}
\alias{single.var.norm.calc}
\title{
Single Change in Variance - Normal Data
}
\description{
Calculates the scaled negative log-likelihood (assuming the data is normally distributed) for all possible changepoint locations and returns the single most probable (min).
}
\usage{
single.var.norm.calc(data, know.mean=FALSE, mu=-1000, extrainf = TRUE)
}
\arguments{
  \item{data}{
	A vector or matrix containing the data within which you wish to find a changepoint.  If data is a matrix, each row is considered a separate dataset.
}
  \item{know.mean}{
	Logical, if TRUE then the mean is assumed known and mu is taken as its value.  If FALSE, and mu=-1000 (default value) then the mean is estimated via maximum likelihood.  If FALSE and the value of mu is supplied, mu is not estimated but is counted as an estimated parameter for decisions.
}
  \item{mu}{
	Numerical value of the true mean of the data.  Either single value or vector of length nrow(data).  If data is a matrix and mu is a single value, the same mean is used for each row.
}
  \item{extrainf}{
	Logical, if TRUE the scaled null and alternative negative likelihood values are returned along with the changepoint location.  If FALSE, only the changepoint location is returned.
}
}
\details{
	This function is used to find a single change in variance for data that is assumed to be normally distributed.  The changepoint returned is simply the location where the log likelihood ratio is maximised, there is no test performed as to whether this location is a true changepoint or not.

The returned negative log likelihoods are scaled so that a test can be directly performed using the log of the likelihood ratio, \eqn{\lambda=\left\{null-alt\right\}}, which should be maximised.

In reality this function should not be used unless you are performing a changepoint test using the output supplied.  This function is used in the "see also" functions that perform various changepoint tests, ideally these should be used.
}
\value{
	If data is a vector (single dataset) and extrainf=FALSE then a single value is returned:
	\item{cpt}{The most probable location of a changepoint (scaled max log likelihood over all possible changepoint locations}
	If data is a vector (single dataset) and extrainf=TRUE then a vector with three elements is returned:
	\item{null}{The scaled null likelihood (negative log likelihood of entire data with no change)}
	\item{alt}{The scaled alternative liklihood at cpt (negative log likelihood of entire data with a change at cpt)}
	If data is an mxn matrix (multiple datasets) and extrainf=FALSE then a vector is returned:
	\item{cpt}{Vector of length m containing the most probable location of a changepoint (scaled max log likelihood over all possible changepoint locations for each row in data.  cpt[1] is the most probable changepoint of the first row in data and cpt[m] is the most probable changepoint for the final row in data.}
	If data is a matrix (multiple datasets) and extrainf=TRUE then a matrix is returned where the first column is the changepoint location for each row in data, the second column is the scaled null likelihood for each row in data, the final column is the scaled maximum of the alternative likelihoods for each row in data.
}
\references{
Change in variance: Chen, J. and Gupta, A. K. (2000) \emph{Parametric statistical change point analysis}, Birkhauser
}
\author{
Rebecca Killick
}


\seealso{
\code{\link{single.var.norm}}, \code{\link{cpt.var}}
}
\examples{
# Example of a change in variance at 100 in simulated normal data
set.seed(1)
x=c(rnorm(100,0,1),rnorm(100,0,10))
single.var.norm.calc(x,extrainf=FALSE) # finds change at 100
single.var.norm.calc(x) # finds change at 100 and gives null likelihood as 765.1905 and alternative likelihood as 435.6529


# Example of no change in variance in simulated normal data
set.seed(1)
x=rnorm(100,0,1)
single.var.norm.calc(x,extrainf=FALSE) # finds change at 53, this is the most probable point of change but if a changepoint test is performed then no change will be found.
single.var.norm.calc(x)# change at 53, null liklihood is -22.47768 and alternative is -24.39894. 
}

\keyword{methods}
\keyword{univar}
\keyword{models}
\keyword{ts}
