\name{changepoint-package}
\alias{changepoint-package}
\alias{changepoint}
\docType{package}
\title{
Contains funcions that run various single and multiple changepoint methods
}
\description{
Implements various mainstream and specialised changepoint methods for finding single and multiple changepoints within data.  Many popular non-parametric and frequentist methods are included.  Users should start by looking at the documentation for cpt.mean, cpt.var and cpt.meanvar.
}
\details{
\tabular{ll}{
Package: \tab changepoint\cr
Type: \tab Package\cr
Version: \tab 0.7 \cr
Date: \tab 2012-05-13\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
}

}
\author{
Rebecca Killick <r.killick@lancs.ac.uk>, Idris A. Eckley <i.eckley@lancs.ac.uk> with contributions from Paul Fearnhead <p.fearnhead@lancs.ac.uk>.

Maintainer: Rebecca Killick <r.killick@lancs.ac.uk>
}

\references{
Chen, J. and Gupta, A. K. (2000) \emph{Parametric statistical change point analysis}, Birkhauser

PELT Algorithm: Killick, R. and Fearnhead, P. and Eckley, I.A. (2011) An exact linear time search algorithm for multiple changepoint detection, \emph{Submitted}

Binary Segmentation: Scott, A. J. and Knott, M. (1974) A Cluster Analysis Method for Grouping Means in the Analysis of Variance, \emph{Biometrics} \bold{30(3)}, 507--512

Segment Neighbourhoods: Auger, I. E. And Lawrence, C. E. (1989) Algorithms for the Optimal Identification of Segment Neighborhoods, \emph{Bulletin of Mathematical Biology} \bold{51(1)}, 39--54

}
\keyword{changepoint}
\keyword{segmentation}
\seealso{
\code{\link{cpt.mean}},\code{\link{cpt.var}},\code{\link{cpt.meanvar}}
}
\examples{
# change in variance
x=c(rnorm(100,0,1),rnorm(100,0,10))
ansvar=cpt.var(x)
plot(ansvar)
print(ansvar)

# change in mean
y=c(rnorm(100,0,1),rnorm(100,5,1))
ansmean=cpt.mean(y)
plot(ansmean,cpt.col='blue')
print(ansmean)

# change in mean and variance
z=c(rnorm(100,0,1),rnorm(100,2,10))
ansmeanvar=cpt.meanvar(z)
plot(ansmeanvar,cpt.width=3)
print(ansmeanvar)
}
