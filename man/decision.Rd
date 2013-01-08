\name{decision}
\alias{decision}
\title{
Likelihood Ratio Decision Function
}
\description{
Uses the function parameters to decide if a proposed changepoint is a true changepoint or due to random variability.  Test is conducted using the user specified penalty.}
\usage{
decision(tau,null,alt=NA,penalty="SIC",n=0,diffparam=1,pen.value=0)
}
\arguments{
  \item{tau}{
	A numeric value or vector specifying the proposed changepoint location(s).
}
  \item{null}{
	The value of the null likelihood to be used in the liklihood ratio calculation.  If tau is a vector, so is null.  If the test statistic is already known, replace the null argument with the test statistic.
}
  \item{alt}{
	The value of the alternative likelihood (at tau) to be used in the likelihood ratio calculation.  If tau is a vector, so is alt.  If the test statistic is already known, then it is used in replacement of the null argument and the alternative should not be specified (default NA to account for this)
}
  \item{penalty}{
	Choice of "None", "SIC", "BIC", "AIC", "Hannan-Quinn", "Asymptotic" and "Manual" penalties.  If Manual is specified, the manual penalty is contained in the pen.value parameter.  The predefined penalties listed do NOT count the changepoint as a parameter, postfix a 1 e.g."SIC1" to count the changepoint as a parameter.
}
  \item{n}{
	The length of the original data, required to give sensible "no changepoint".
}
  \item{diffparam}{
	The difference in the number of parameters in the null and alternative hypotheses, required for the SIC, BIC, AIC, Hanna-Quinn and possibly Manual penalties.
}
  \item{pen.value}{
	The value of the penalty when using the Manual or Asymptotic option.  This can be a numeric value or text giving the formula to use.  Available variables are, n=length of original data, null=null likelihood, alt=alternative likelihood, tau=proposed changepoint, diffparam=difference in number of alternatve and null parameters.
}
}
\details{
	This function is used to test whether tau is a true changepoint or not.  This test uses the likelihood ratio as the test statistic and performs the test where the null hypothesis is no change point and the alternative hypothesis is a single changepoint at tau.  The test is (null-alt)>=penalty, if TRUE then the changepoint is deemed a true changepoint, if FALSE then n (length of data) is returned.

If the test statistic is already known then it replaces the null value and the alternative is not required (default NA).  In this case the test is null>=penalty, if TRUE then the changepoint is deemed a true changepoint, if FALSE then n (length of data) is returned.

In reality this function should not be used unless you are performing a changepoint test using output from other functions.  This function is used in the "see also" functions that perform various changepoint tests, ideally these should be used.
}
\value{
	A list is returned with two elements, cpt and pen.
	\item{cpt}{ If tau is a single value then a single value is returned: Either the value of the true changepoint location or n (length of data) if no changepoint is found.
	
	If tau is a vector of length m then a vector of length m is returned:Each element is either the value of the true changepoint location or n (length of data) if no changepoint is found.  The first element is for the first value of tau and the final element is for the final value of tau.
	}
	\item{pen}{The numeric value of the penalty used for the test(s).}
}
\references{
SIC/BIC: Schwarz, G. (1978) Estimating the Dimension of a Model, \emph{The Annals of Statistics} \bold{6(2)}, 461--464

AIC: Akaike, H. (1974) A new look at the statistical model identification, \emph{Automatic Control, IEEE Transactions on} \bold{19(6)}, 716--723

Hannan-Quinn: Hannan, E. J. and B. G. Quinn (1979) The Determination of the Order of an Autoregression, \emph{Journal of the Royal Statistical Society, B} \bold{41}, 190--195
}
\author{
Rebecca Killick
}


\seealso{
\code{\link{cpt.mean}},\code{\link{cpt.var}},\code{\link{cpt.meanvar}},\code{\link{single.mean.norm}},\code{\link{single.var.norm}},\code{\link{single.meanvar.norm}}
}
\examples{
# Example of finding a change
out=c(100,765.1905,435.6529) # tau, null, alt
decision(out[1],out[2],out[3],penalty="SIC",n=200,diffparam=1) # returns 100 as a true changepoint

# Example of no change found
out=c(53,-22.47768,-24.39894) # tau, null, alt
decision(out[1],out[2],out[3],penalty="Manual",n=200,diffparam=1,pen.value="2*log(n)")
}

\keyword{methods}
\keyword{univar}
\keyword{models}
\keyword{ts}
