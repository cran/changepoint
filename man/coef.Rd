\name{coef}
\alias{coef}
\title{
Generic Function - coef
}
\description{
Generic function
}
\usage{
coef(object)
}
\arguments{
  \item{object}{
	Depending on the class of \code{object} depends on the method used (and if one exists)
}
}
\details{
	Generic Function
}
\value{
	Depends on the class of \code{object}, see individual methods
}
\author{
Rebecca Killick
}


\seealso{
\code{\link{coef-methods}}
}
\examples{
x=new("cpt") # new cpt object
coef(x) # retrieves the param.est slot from x
  # which shows the coefficients used in the model
}

\keyword{methods}
\keyword{cpt}
