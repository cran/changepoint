\name{cpt.reg-class}
\Rdversion{1.1}
\docType{class}
\alias{cpt.reg-class}
\alias{cpt.reg-method,cpts}
\alias{cpt.reg-method,cpttype}
\alias{cpt.reg-method,data.set}
\alias{cpt.reg-method,distribution}
\alias{cpt.reg-method,method}
\alias{cpt.reg-method,ncpts.max}
\alias{cpt.reg-method,param.est}
\alias{cpt.reg-method,pen.type}
\alias{cpt.reg-method,pen.value}
\alias{cpt.reg-method,cpts<-}
\alias{cpt.reg-method,cpttype<-}
\alias{cpt.reg-method,data.set<-}
\alias{cpt.reg-method,distribution<-}
\alias{cpt.reg-method,method<-}
\alias{cpt.reg-method,ncpts.max<-}
\alias{cpt.reg-method,param.est<-}
\alias{cpt.reg-method,pen.type<-}
\alias{cpt.reg-method,pen.value<-}
\alias{cpt.reg-method,print}
\alias{cpt.reg-method,summary}
\alias{cpt.reg-method,param}

\title{Class "cpt.reg"}
\description{
A class for changepoint objects, specifically change in regression.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("cpt", ...)}.
  \describe{
    \item{\code{new("cpt", ...)}:}{creates a new object with class cpt }
  }
}
\section{Slots}{
  \describe{
    \item{\code{data.set}:}{Object of class \code{"numeric"}, the original vector of data }
    \item{\code{cpttype}:}{Object of class \code{"character"}, the type of changepoint that was identified }
    \item{\code{method}:}{Object of class \code{"character"}, the method that was used to search for changepoints, default change in regression }
    \item{\code{distribution}:}{Object of class \code{"character"}, the assumed distribution of the data }
    \item{\code{pen.type}:}{Object of class \code{"character"}, the penalty type specified in the analysis}
    \item{\code{pen.value}:}{Object of class \code{"numeric"}, the value of the penalty used in the analysis}
    \item{\code{cpts}:}{Object of class \code{"numeric"}, vector of changepoints identified }
    \item{\code{ncpts.max}:}{Object of class \code{"numeric"}, maximum number of changepoint that can be identified }
    \item{\code{param.est}:}{Object of class \code{"list"}, list where each element is a vector of parameter estimates, if requested }
    \item{\code{date}:}{Object of class \code{"character"}, date and time the changepoint analysis was run }
  }
}
\section{Methods}{
  \describe{
    \item{cpts}{\code{signature(object = "cpt.reg")}: retrieves cpts slot }
    \item{cpttype}{\code{signature(object = "cpt.reg")}: retrieves cpttype slot }
    \item{data.set}{\code{signature(object = "cpt.reg")}: retrieves data.set slot }
    \item{distribution}{\code{signature(object = "cpt.reg")}: retrieves distribution slot }
    \item{ncpts.max}{\code{signature(object = "cpt.reg")}: retrieves ncpts.max slot }
    \item{method}{\code{signature(object = "cpt.reg")}: retrieves method slot }
    \item{param.est}{\code{signature(object = "cpt.reg")}: retrieves param.est slot }
    \item{pen.type}{\code{signature(object = "cpt.reg")}: retrieves pen.type slot }
    \item{pen.value}{\code{signature(object = "cpt.reg")}: retrieves pen.value slot }
    \item{cpts<-}{\code{signature(object = "cpt.reg")}: replaces cpts slot }
    \item{cpttype<-}{\code{signature(object = "cpt.reg")}: replaces cpttype slot }
    \item{data.set<-}{\code{signature(object = "cpt.reg")}: replaces data.set slot }
    \item{distribution<-}{\code{signature(object = "cpt.reg")}: replaces distribution slot }
    \item{ncpts.max<-}{\code{signature(object = "cpt.reg")}: replaces ncpts.max slot }
    \item{method<-}{\code{signature(object = "cpt.reg")}: replaces method slot }
    \item{param.est<-}{\code{signature(object = "cpt.reg")}: replaces param.est slot }
    \item{pen.type<-}{\code{signature(object = "cpt.reg")}: replaces pen.type slot }
    \item{pen.value<-}{\code{signature(object = "cpt.reg")}: replaces pen.value slot }
    \item{print}{\code{signature(object = "cpt.reg")}: prints details of the cpt object including summary}
    \item{summary}{\code{signature(object = "cpt.reg")}: prints a summary of the cpt object }
    \item{param}{\code{signature(object = "cpt.reg")}: calculates the parameter estimates for the cpt object}
	 }
}
\author{
Rebecca Killick
}

\seealso{
\code{\link{plot-methods}},\code{\link{cpts-methods}},\code{\linkS4class{cpt}},\code{\link{cpt.reg}}
}
\examples{
showClass("cpt.reg")

x=new("cpt.reg") # creates a new object with the cpt.reg class defaults
data.set(x)  # retrieves the data.set slot from x
data.set(x)<-matrix(1:10,nrow=5,ncol=2) # replaces the data.set slot from x with a matrix

# single change in regression example
set.seed(1)
x=1:250
y=c(0.01*x[1:100],1.5-0.02*(x[101:250]-101))
ynoise=y+rnorm(250,0,0.2)
data=cbind(ynoise,x)
out=cpt.reg(data) # analyses x using default values
print(out) # gives details of analysis including default values used
summary(out) # shorter version of print(out)
}
\keyword{classes}
\keyword{cpt}
