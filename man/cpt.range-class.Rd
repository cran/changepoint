\name{cpt.range-class}
\Rdversion{1.1}
\docType{class}
\alias{cpt.range-class}

\title{Class "cpt.range"}
\description{
A class for changepoint objects that return more than 1 segmentation. Inherits from cpt class.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("cpt.range", ...)}.
  \describe{
    \item{\code{new("cpt.range", ...)}:}{creates a new object with class cpt.range }
  }
}
\section{Slots}{
  \describe{
    \item{\code{cpts.full}:}{Object of class \code{"matrix"}, each row of the matrix is a different segmentation of the data (different set of changepoints).}
    \item{\code{pen.value.full}:}{Object of class \code{"vector"}, each element is the penalty used to create the set of changepoints in the corresponding row of \code{cpts.full}.}
    The remaining slots are inherited from the \code{cpt} class.
    \item{\code{data.set}:}{Object of class \code{"ts"}, a coerced time series of the original data. Inherited from cpt class. }
    \item{\code{cpttype}:}{Object of class \code{"character"}, the type of changepoint that was identified. Inherited from cpt class. }
    \item{\code{method}:}{Object of class \code{"character"}, the method that was used to search for changepoints. Inherited from cpt class. }
    \item{\code{test.stat}:}{Object of class \code{"character"}, the test statistic for the analysis of the data. Inherited from cpt class. }
    \item{\code{pen.type}:}{Object of class \code{"character"}, the penalty type specified in the analysis. Inherited from cpt class.}
    \item{\code{pen.value}:}{Object of class \code{"numeric"}, the value of the penalty used in the analysis. Inherited from cpt class.}
    \item{\code{minseglen}:}{Object of class \code{"numeric"}, the minimum segment length (no. of observations between changepoints) used in the analysis. Inherited from cpt class.}
    \item{\code{cpts}:}{Object of class \code{"numeric"}, vector of optimal changepoints identified. Inherited from cpt class. }
    \item{\code{ncpts.max}:}{Object of class \code{"numeric"}, maximum number of changepoint that can be identified. Inherited from cpt class. }
    \item{\code{param.est}:}{Object of class \code{"list"}, list where each element is a vector of parameter estimates, if requested. Inherited from cpt class. }
    \item{\code{date}:}{Object of class \code{"character"}, date and time the changepoint analysis was run. Inherited from cpt class. }
  }
}
\section{Methods}{
  \describe{
    \item{cpts.full}{\code{signature(object = "cpt.range")}: retrieves cpts.full slot }
    \item{pen.value.full}{\code{signature(object = "cpt.range")}: retrieves pen.value.full slot }
    \item{cpts.full<-}{\code{signature(object = "cpt.range")}: replaces cpts.full slot }
    \item{pen.value.full<-}{\code{signature(object = "cpt.range")}: replaces pen.value.full slot }
    \item{print}{\code{signature(object = "cpt.range")}: prints details of the cpt.range object including summary}
    \item{summary}{\code{signature(object = "cpt.range")}: prints a summary of the cpt.range object }
  }
}

\author{
Rebecca Killick
}

\seealso{
\code{\link{cpts.full-methods}},\code{\linkS4class{cpt}},\code{\link{cpt.mean}},\code{\link{cpt.var}},\code{\link{cpt.meanvar}}
}
\examples{
showClass("cpt.range") # shows the structure of the cpt.range class

x=new("cpt.range") # creates a new object with the cpt.range class defaults
cpts(x) # retrieves the cpts slot from x
cpts(x)<-c(10,50,100) # replaces the cpts slot from x with c(10,50,100)

# Example of a change in variance at 100 in simulated normal data
set.seed(1)
x=c(rnorm(100,0,1),rnorm(100,0,10))
ans=cpt.var(x,method="PELT",penalty="CROPS",pen.value=c(3,10)) 
print(ans) # prints details of the analysis including a summary
summary(ans)
plot(ans) # plots the data with change (vertical line) at 100
logLik(ans) # raw likelihood of the data with changepoints, second value is likelihood + penalty
}
\keyword{classes}
\keyword{cpt}
