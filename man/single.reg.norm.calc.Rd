\name{single.reg.norm.calc}
\alias{single.reg.norm.calc}
\title{
Single Change in Regression - Normal Data
}
\description{
Calculates the scaled log-likelihood (assuming the errors are zero-mean normally distributed) for all possible changepoint locations and returns the single most probable (max).
}
\usage{
single.reg.norm.calc(data, extrainf = TRUE)
}
\arguments{
  \item{data}{
	A matrix or 3-d array containing the data within which you wish to find a changepoint.  If data is a 3-d array, each first dimension is considered a separate dataset.  Within each dataset the first column is considered the response and the further columns are the covariates.  This function does not create a constant term so if one is required a column of ones should appear in the covariates.

NOTE: each model must have the same number of parameters otherwise singularities appear
}
  \item{extrainf}{
	Logical, if TRUE the scaled null and alternative likelihood values are returned along with the changepoint location.  If FALSE, only the changepoint location is returned.
}
}
\details{
	This function is used to find a single change in regression for data that is assumed to have normally distributed errors with zero mean.  The changepoint returned is simply the location where the log likelihood is maximised, there is no test performed as to whether this location is a true changepoint or not.

The returned likelihoods are scaled so that a test can be directly performed using the log of the likelihood ratio, \eqn{\lambda=\left\{null-alt\right\}}, which should be maximised.

In reality this function should not be used unless you are performing a changepoint test using the output supplied.  This function is used in the "see also" functions that perform various changepoint tests, ideally these should be used.
}
\value{
	If data is a nxp matrix (single dataset) and extrainf=FALSE then a single value is returned:
	\item{cpt}{The most probable location of a changepoint (scaled max log likelihood over all possible changepoint locations}
	If data is a nxp matrix (single dataset) and extrainf=TRUE then a vector with three elements is returned:
	\item{cpt}{The most probable location of a changepoint (scaled max log likelihood over all possible changepoint locations}
	\item{null}{The scaled null likelihood (log likelihood of entire data with no change)}
	\item{alt}{The scaled alternative liklihood at cpt (log likelihood of entire data with a change at cpt)}
	If data is an mxnxp array (multiple datasets) and extrainf=FALSE then a vector is returned:
	\item{cpt}{Vector of length m containing the most probable location of a changepoint (scaled max log likelihood over all possible changepoint locations for each row in data.  cpt[1] is the most probable changepoint of data[1,,] and cpt[m] is the most probable changepoint for data[m,,].}
	If data is an mxnxp array (multiple datasets) and extrainf=TRUE then a matrix is returned where the first column is the changepoint location for each dataset, the second column is the scaled null likelihood for each dataset, the final column is the scaled maximum of the alternative likelihoods for each dataset.
}
\references{
Change in regression: Chen, J. and Gupta, A. K. (2000) \emph{Parametric statistical change point analysis}, Birkhauser
}
\author{
Rebecca Killick
}


\seealso{
\code{\link{single.reg.norm}}, \code{\link{cpt.reg}}
}
\examples{
# Example of a change in regression at 100 in simulated normal data
set.seed(1)
x=1:250
y=c(0.01*x[1:100],1.5-0.02*(x[101:250]-101))
ynoise=y+rnorm(250,0,0.2)
data=cbind(ynoise,1,x) # including an intercept
single.reg.norm.calc(data,extrainf=FALSE) # finds change at 100
single.reg.norm.calc(data) # finds change at 100 and gives null likelihood as 1119.3264 and alternative likelihood as 554.1691


# Example of no change in regression in simulated normal data
set.seed(1)
x=1:250
y=0.01*x
ynoise=y+rnorm(250,0,0.2)
data=cbind(ynoise,x) # no intercept included
single.reg.norm.calc(data,extrainf=FALSE) # finds change at 96, this is the most probable point of change but if a changepoint test is performed then no change will be found.
single.reg.norm.calc(data)# change at 96, null liklihood is 555.6434 and alternative is 552.6056
}

\keyword{methods}
\keyword{univar}
\keyword{models}
\keyword{ts}
\keyword{lm}
