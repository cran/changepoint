\name{single.mean.norm.calc}
\alias{single.mean.norm.calc}
\title{
Single Change in Mean - Normal Data
}
\description{
Calculates the scaled log-likelihood (assuming the data is normally distributed) for all possible changepoint locations and returns the single most probable (max).
}
\usage{
single.mean.norm.calc(data, extrainf = TRUE)
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
	This function is used to find a single change in mean for data that is assumed to be normally distributed.  The changepoint returned is simply the location where the log likelihood is maximised, there is no test performed as to whether this location is a true changepoint or not.

The returned likelihoods are scaled so that a test can be directly performed using the log of the likelihood ratio, \eqn{\lambda=\frac{1}{\sigma^2}\left\{null-alt\right\}}, which should be maximised.

In reality this function should not be used unless you are performing a changepoint test using the output supplied.  This function is used in the "see also" functions that perform various changepoint tests, ideally these should be used.
}
\value{
	If data is a vector (single dataset) and extrainf=FALSE then a single value is returned:
	\item{cpt}{The most probable location of a changepoint (scaled max log likelihood over all possible changepoint locations}
	If data is a vector (single dataset) and extrainf=TRUE then a vector with three elements is returned:
	\item{cpt}{The most probable location of a changepoint (scaled max log likelihood over all possible changepoint locations}
	\item{null}{The scaled null likelihood (log likelihood of entire data with no change)}
	\item{alt}{The scaled alternative liklihood at cpt (log likelihood of entire data with a change at cpt)}
	If data is an mxn matrix (multiple datasets) and extrainf=FALSE then a vector is returned:
	\item{cpt}{Vector of length m containing the most probable location of a changepoint (scaled max log likelihood over all possible changepoint locations for each row in data.  cpt[1] is the most probable changepoint of the first row in data and cpt[m] is the most probable changepoint for the final row in data.}
	If data is a matrix (multiple datasets) and extrainf=TRUE then a matrix is returned where the first column is the changepoint location for each row in data, the second column is the scaled null likelihood for each row in data, the final column is the scaled maximum of the alternative likelihoods for each row in data.
}
\references{
Change in Normal mean: Hinkley, D. V. (1970) Inference About the Change-Point in a Sequence of Random Variables, \emph{Biometrika} \bold{57}, 1--17
}
\author{
Rebecca Killick
}


\seealso{
\code{\link{single.mean.norm}}, \code{\link{cpt.mean}}
}
\examples{
# Example of a change in mean at 100 in simulated normal data
set.seed(1)
x=c(rnorm(100,0,1),rnorm(100,10,1))
single.mean.norm.calc(x,extrainf=FALSE) # finds change at 100
single.mean.norm.calc(x) # finds change at 100 and gives null likelihood as 5025.0857 and alternative likelihood as 170.7051


# Example of no change in mean in simulated normal data
set.seed(1)
x=rnorm(100,0,1)
single.mean.norm.calc(x,extrainf=FALSE) # finds change at 96, this is the most probable point of change but if a changepoint test is performed then no change will be found.
single.mean.norm.calc(x)# change at 96, null liklihood is 79.86945 and alternative is 75.73725
}

\keyword{methods}
\keyword{univar}
\keyword{models}
\keyword{ts}
