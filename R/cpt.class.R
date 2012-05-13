	setClass("cpt",representation(data.set="numeric", cpttype="character", method="character", 	distribution="character",pen.type="character",pen.value="numeric",cpts="numeric",ncpts.max="numeric",param.est="list",date="character"),prototype(date=date()))

	setClass("cpt.reg",representation(data.set="matrix", cpttype="character", method="character", distribution="character",pen.type="character",pen.value="numeric",cpts="numeric",ncpts.max="numeric",param.est="list",date="character"),prototype(cpttype="regression",date=date()))

# retrival functions for slots
	if(!isGeneric("data.set")) {
		if (is.function("data.set")){
			fun <- data.set
		}
		else {fun <- function(object){
				standardGeneric("data.set")
			}
		}
		setGeneric("data.set", fun)
	}
	setMethod("data.set","cpt",function(object) object@data.set)
	setMethod("data.set","cpt.reg",function(object) object@data.set)

	if(!isGeneric("cpttype")) {
		if (is.function("cpttype")){
			fun <- cpttype
		}
		else {fun <- function(object){
				standardGeneric("cpttype")
			}
		}
		setGeneric("cpttype", fun)
	}
	setMethod("cpttype","cpt",function(object) object@cpttype)
	setMethod("cpttype","cpt.reg",function(object) object@cpttype)

	if(!isGeneric("method")) {
		if (is.function("method")){
			fun <- method
		}
		else {fun <- function(object){
				standardGeneric("method")
			}
		}
		setGeneric("method", fun)
	}
	setMethod("method","cpt",function(object) object@method)
	setMethod("method","cpt.reg",function(object) object@method)
	
	if(!isGeneric("distribution")) {
		if (is.function("distribution")){
			fun <- distribution
		}
		else {fun <- function(object){
				standardGeneric("distribution")
			}
		}
		setGeneric("distribution", fun)
	}
	setMethod("distribution","cpt",function(object) object@distribution)
	setMethod("distribution","cpt.reg",function(object) object@distribution)
	
	if(!isGeneric("pen.type")) {
		if (is.function("pen.type")){
			fun <- pen.type
		}
		else {fun <- function(object){
				standardGeneric("pen.type")
			}
		}
		setGeneric("pen.type", fun)
	}
	setMethod("pen.type","cpt",function(object) object@pen.type)
	setMethod("pen.type","cpt.reg",function(object) object@pen.type)
	
	if(!isGeneric("pen.value")) {
		if (is.function("pen.value")){
			fun <- pen.value
		}
		else {fun <- function(object){
				standardGeneric("pen.value")
			}
		}
		setGeneric("pen.value", fun)
	}
	setMethod("pen.value","cpt",function(object) object@pen.value)
	setMethod("pen.value","cpt.reg",function(object) object@pen.value)
	
	if(!isGeneric("cpts")) {
		if (is.function("cpts")){
			fun <- cpts
		}
		else {fun <- function(object){
				standardGeneric("cpts")
			}
		}
		setGeneric("cpts", fun)
	}
	setMethod("cpts","cpt",function(object) object@cpts)
	setMethod("cpts","cpt.reg",function(object) object@cpts)
	
	if(!isGeneric("ncpts.max")) {
		if (is.function("ncpts.max")){
			fun <- cpts
		}
		else {fun <- function(object){
				standardGeneric("ncpts.max")
			}
		}
		setGeneric("ncpts.max", fun)
	}
	setMethod("ncpts.max","cpt",function(object) object@ncpts.max)
	setMethod("ncpts.max","cpt.reg",function(object) object@ncpts.max)

	if(!isGeneric("param.est")) {
		if (is.function("param.est")){
			fun <- param.est
		}
		else {fun <- function(object){
				standardGeneric("param.est")
			}
		}
		setGeneric("param.est", fun)
	}
	setMethod("param.est","cpt",function(object) object@param.est)
	setMethod("param.est","cpt.reg",function(object) object@param.est)
	
# ncpts function
	if(!isGeneric("ncpts")) {
		if (is.function("ncpts")){
			fun <- ncpts
		}
		else {fun <- function(object){
				standardGeneric("ncpts")
			}
		}
		setGeneric("ncpts", fun)
	}
	setMethod("ncpts","cpt",function(object) length(cpts(object))-1)
	setMethod("ncpts","cpt.reg",function(object) length(cpts(object))-1)

# replacement functions for slots
	setGeneric("data.set<-", function(object, value) standardGeneric("data.set<-"))
	setReplaceMethod("data.set", "cpt", function(object, value) {
		object@data.set <- value
		return(object)
	})
	setReplaceMethod("data.set", "cpt.reg", function(object, value) {
		object@data.set <- value
		return(object)
	})
	
	setGeneric("cpttype<-", function(object, value) standardGeneric("cpttype<-"))
	setReplaceMethod("cpttype", "cpt", function(object, value) {
		object@cpttype <- value
		return(object)
	})
	setReplaceMethod("cpttype", "cpt.reg", function(object, value) {
		object@cpttype <- value
		return(object)
	})
	
	setGeneric("method<-", function(object, value) standardGeneric("method<-"))
	setReplaceMethod("method", "cpt", function(object, value) {
		object@method <- value
		return(object)
	})
	setReplaceMethod("method", "cpt.reg", function(object, value) {
		object@method <- value
		return(object)
	})
	
	setGeneric("distribution<-", function(object, value) standardGeneric("distribution<-"))
	setReplaceMethod("distribution", "cpt", function(object, value) {
		object@distribution <- value
		return(object)
	})
	setReplaceMethod("distribution", "cpt.reg", function(object, value) {
		object@distribution <- value
		return(object)
	})
	
	setGeneric("pen.type<-", function(object, value) standardGeneric("pen.type<-"))
	setReplaceMethod("pen.type", "cpt", function(object, value) {
		object@pen.type <- value
		return(object)
	})
	setReplaceMethod("pen.type", "cpt.reg", function(object, value) {
		object@pen.type <- value
		return(object)
	})
	
	setGeneric("pen.value<-", function(object, value) standardGeneric("pen.value<-"))
	setReplaceMethod("pen.value", "cpt", function(object, value) {
		object@pen.value <- value
		return(object)
	})
	setReplaceMethod("pen.value", "cpt.reg", function(object, value) {
		object@pen.value <- value
		return(object)
	})
	
	setGeneric("cpts<-", function(object, value) standardGeneric("cpts<-"))
	setReplaceMethod("cpts", "cpt", function(object, value) {
		object@cpts <- value
		return(object)
	})
	setReplaceMethod("cpts", "cpt.reg", function(object, value) {
		object@cpts <- value
		return(object)
	})

	setGeneric("ncpts.max<-", function(object, value) standardGeneric("ncpts.max<-"))
	setReplaceMethod("ncpts.max", "cpt", function(object, value) {
		object@ncpts.max <- value
		return(object)
	})
	setReplaceMethod("ncpts.max", "cpt.reg", function(object, value) {
		object@ncpts.max <- value
		return(object)
	})
	
	setGeneric("param.est<-", function(object, value) standardGeneric("param.est<-"))
	setReplaceMethod("param.est", "cpt", function(object, value) {
		object@param.est <- value
		return(object)
	})
	setReplaceMethod("param.est", "cpt.reg", function(object, value) {
		object@param.est <- value
		return(object)
	})

# parameter functions
	setGeneric("param", function(object,...) standardGeneric("param"))
	setMethod("param", "cpt", function(object,shape,...) {			
		param.mean=function(object){
			cpts=c(0,cpts(object))
			nseg=length(cpts)-1
			data=data.set(object)
			tmpmean=NULL
			for(j in 1:nseg){
				tmpmean[j]=mean(data[(cpts[j]+1):(cpts[j+1])])
			}
			return(tmpmean)
		}
		param.var=function(object){
			cpts=c(0,cpts(object))
			nseg=length(cpts)-1
			data=data.set(object)
			tmpvar=NULL
			for(j in 1:nseg){
				tmpvar[j]=var(data[(cpts[j]+1):(cpts[j+1])])
			}
			return(tmpvar)
		}
		param.scale=function(object,shape){
			cpts=c(0,cpts(object))
			nseg=length(cpts)-1
			data=data.set(object)
			y=c(0,cumsum(data))
			tmpscale=NULL
			for(j in 1:nseg){
				tmpscale[j]=(y[(cpts[j+1]+1)]-y[(cpts[j]+1)])/((cpts[j+1]-cpts[j]+1)*shape)
			}
			return(tmpscale)			
		}
		if(cpttype(object)=="mean"){
			param.est(object)<-list(mean=param.mean(object))
		}
		else if(cpttype(object)=="variance"){
			param.est(object)<-list(variance=param.var(object))
		}
		else if(cpttype(object)=="mean and variance"){
			if(distribution(object)=="Normal"){
				param.est(object)<-list(mean=param.mean(object),variance=param.var(object))
			}
			else if(distribution(object)=="Gamma"){
				param.est(object)<-list(scale=param.scale(object,shape=shape),shape=shape)
			}
			else if(distribution(object)=="Exponential"){
				param.est(object)<-list(rate=1/param.mean(object))
			}
			else if(distribution(object)=="Poisson"){
			  param.est(object)<-list(lambda=param.mean(object))
			}
			else{
				stop("Unknown distribution for a change in mean and variance")
			}
		}
		else{
			stop("Unknown changepoint type, must be 'mean', 'variance' or 'mean and variance'")
		}
		return(object)
	})


	setMethod("param", "cpt.reg", function(object,shape,...) {			
		param.norm=function(object){
			cpts=c(0,cpts(object))
			nseg=length(cpts)-1
			data=data.set(object)
			p=ncol(data)-1
			tmpbeta=matrix(NA,ncol=p,nrow=nseg)
			for(j in 1:nseg){
				tmpbeta[j,]=solve(t(data[(cpts[j]+1):cpts[j+1],2:(p+1)])%*%data[(cpts[j]+1):cpts[j+1],2:(p+1)],t(data[(cpts[j]+1):cpts[j+1],2:(p+1)])%*%data[(cpts[j]+1):cpts[j+1],1])
			}
			return(tmpbeta)
		}
		if(distribution(object)=="Normal"){
			param.est(object)<-list(beta=param.norm(object))
		}
		else{
			stop("Unknown distribution, must be 'Normal'")
		}
		return(object)
	})

# summary functions
	setMethod("summary","cpt",function(object){
	    cat("Changepoint type      : Change in",cpttype(object),'\n')
	    cat("Method of analysis    :",method(object),"\n")
	    cat("Assumed Distribution  :", distribution(object),"\n")
	    cat("Type of penalty       :", pen.type(object), "with value,",pen.value(object),"\n")
	    cat("Maximum no. of cpts   :", ncpts.max(object),"\n")
	    if(length(cpts(object))<=20){cat("Changepoint Locations :",cpts(object),"\n")}
	    else{cat("Number of changepoints:", length(cpts(object)),"\n")}
	})

	setMethod("summary","cpt.reg",function(object){
	    cat("Changepoint type     : Change in",cpttype(object),'\n')
	    cat("Method of analysis   :",method(object),"\n")
	    cat("Assumed Distribution :", distribution(object),"\n")
	    cat("Type of penalty      :", pen.type(object), "with value,",pen.value(object),"\n")
	    cat("Maximum no. of cpts   :", ncpts.max(object),"\n")
	    if(length(cpts(object))<=20){cat("Changepoint Locations :",cpts(object),"\n")}
	    else{cat("Number of changepoints:", length(cpts(object)),"\n")}
	})

# print functions
	setMethod("print","cpt",function(x){
	    cat("Class 'cpt' : Changepoint Object\n")
	    cat("       ~~   : S4 class containing", length(attributes(x))-1, "slots with names\n")
	    cat("             ", names(attributes(x))[1:(length(attributes(x))-1)], "\n\n")
	    cat("Created on  :", x@date, "\n\n")
	    cat("summary(.)  :\n----------\n")
	    summary(x)
	})
	setMethod("print","cpt.reg",function(x){
	    cat("Class 'cpt' : Changepoint Object\n")
	    cat("       ~~   : S4 class containing", length(attributes(x))-1, "slots with names\n")
	    cat("             ", names(attributes(x))[1:(length(attributes(x))-1)], "\n\n")
	    cat("Created on  :", x@date, "\n\n")
	    cat("summary(.)  :\n----------\n")
	    summary(x)
	})

# plot functions
	setMethod("plot","cpt",function(x,cpt.col='red',cpt.width=1,cpt.style=1,...){
		plot(data.set(x),...)
		if(cpttype(x)=="variance"){
			abline(v=cpts(x)[1:(length(cpts(x))-1)],col=cpt.col,lwd=cpt.width,lty=cpt.style)
		}
		else if(cpttype(x)=="mean"  ||  cpttype(x)=="mean and variance"){
			nseg=length(cpts(x))
			cpts=c(0,cpts(x))
			if((distribution(x)=="Normal")||(distribution(x)=="CUSUM")){
				means=param.est(x)$mean
			}
			else if(distribution(x)=="Gamma"){
				means=param.est(x)$scale*param.est(x)$shape
			}
			else if(distribution(x)=="Exponential"){
				means=1/param.est(x)$rate
			}
			else if(distribution(x)=="Poisson"){
			  means=param.est(x)$lambda
			}
			else{
				stop('Invalid Changepoint distribution type')
			}
			for(i in 1:nseg){
				segments(cpts[i]+1,means[i],cpts[i+1],means[i],col=cpt.col,lwd=cpt.width,lty=cpt.style)
			}
		}
		else{
			stop('Invalid Changepoint Type for plotting.\n Can only plot mean, variance, mean and variance')
		}
	})

	setMethod("plot","cpt.reg",function(x,cpt.col='red',cpt.width=1,cpt.style=1,...){
		if(dim(data.set(x))[2]>3){
			stop("A plot function for regression of more than one regressor is not available")
		}
		else if(dim(data.set(x))[2]==3){
			if(data.set(x)[1,2]==1){
				nseg=length(cpts(x))
				cpts=c(0,cpts(x))
				betas=param.est(x)$beta
				plot(data.set(x[,3]),data.set(x[,1]),...)
				for(i in 1:nseg){
					segments(cpts[i]+1,betas[i,1]*data.set(x)[cpts[i]+1,2]+betas[i,2]*data.set(x)[cpts[i]+1,3],cpts[i+1],betas[i,1]*data.set(x)[cpts[i+1],2]+betas[i,2]*data.set(x)[cpts[i+1],3],col=cpt.col,lwd=cpt.width,lty=cpt.style)
				}
			}
			if(data.set(x)[1,3]==1){
				nseg=length(cpts(x))
				cpts=c(0,cpts(x))
				betas=param.est(x)$beta
				plot(data.set(x[,2]),data.set(x[,1]),...)
				for(i in 1:nseg){
					segments(cpts[i]+1,betas[i,2]*data.set(x)[cpts[i]+1,3]+betas[i,1]*data.set(x)[cpts[i]+1,2],cpts[i+1],betas[i,2]*data.set(x)[cpts[i+1],3]+betas[i,1]*data.set(x)[cpts[i+1],2],col=cpt.col,lwd=cpt.width,lty=cpt.style)
				}
			}
			else{
				stop("A plot function for regression of more than one regressor is not available")
			}
		}
		else{
			nseg=length(cpts(x))
			cpts=c(0,cpts(x))
			betas=param.est(x)$beta
			plot(data.set(x[,2]),data.set(x[,1]),...)
			for(i in 1:nseg){
				segments(cpts[i]+1,betas[i,1]*data.set(x)[cpts[i]+1,2],cpts[i+1],betas[i,2]+betas[i,1]*data.set(x)[cpts[i+1],2],col=cpt.col,lwd=cpt.width,lty=cpt.style)
			}
		}
	})

# likelihood functions
	setGeneric("likelihood", function(object) standardGeneric("likelihood"))
	setMethod("likelihood", "cpt", function(object) {
		if(distribution(object)=="Normal"){
			if(cpttype(object)=="mean"){
				mll.mean=function(x2,x,n){
				  return( x2-(x^2)/n)
				}
				y2=c(0,cumsum(data.set(object)^2))
				y=c(0,cumsum(data.set(object)))
				cpts=c(0,cpts(object))
				nseg=length(cpts)-1
				tmplike=0
				for(j in 1:nseg){
			        	tmplike=tmplike+mll.mean(y2[cpts[j+1]+1]-y2[cpts[j]+1],y[cpts[j+1]+1]-y[cpts[j]+1],cpts[j+1]-cpts[j])
				}
				like=c(tmplike,tmplike+(nseg-1)*pen.value(object))
				names(like)=c("like","likepen")
			}
			else if(cpttype(object)=="variance"){
				mll.var=function(x,n){
					neg=x<=0
					x[neg==TRUE]=0.00000000001    
					return( n*(log(2*pi)+log(x/n)+1))
				}
				y2=c(0,cumsum(data.set(object)^2))
				cpts=c(0,cpts(object))
				nseg=length(cpts)-1
				tmplike=0
				for(j in 1:nseg){
					tmplike=tmplike+mll.var(y2[cpts[j+1]+1]-y2[cpts[j]+1],cpts[j+1]-cpts[j])
				}
				like=c(tmplike,tmplike+(nseg-1)*pen.value(object))
				names(like)=c("like","likepen")
			}
			else if(cpttype(object)=="mean and variance"){
				mll.meanvar=function(x2,x,n){
					sigmasq=(1/n)*(x2-(x^2)/n)
					neg=sigmasq<=0
					sigmasq[neg==TRUE]=0.00000000001
					return( n*(log(2*pi)+log(sigmasq)+1))
				}
				y2=c(0,cumsum(data.set(object)^2))
				y=c(0,cumsum(data.set(object)))
				cpts=c(0,cpts(object))
				nseg=length(cpts)-1
				tmplike=0
				for(j in 1:nseg){
					tmplike=tmplike+mll.meanvar(y2[cpts[j+1]+1]-y2[cpts[j]+1],y[cpts[j+1]+1]-y[cpts[j]+1],cpts[j+1]-cpts[j])
				}
				like=c(tmplike,tmplike+(nseg-1)*pen.value(object))
				names(like)=c("like","likepen")
			}
			else{
				stop("Unknown changepoint type, must be 'mean', 'variance' or 'mean and variance'")
			}
		}
		else if(distribution(object)=="Gamma"){
			if(cpttype(object)!="mean and variance"){
				stop("Unknown changepoint type for dist='Gamma', must be 'mean and variance'")
			}
			else{
			  mll.meanvarg=function(x,n,shape){
			    return(n*shape*log(n*shape)-n*shape*log(x))
			  }
				y=c(0,cumsum(data.set(object)))
				shape=param.est(object)$shape
				cpts=c(0,cpts(object))
				nseg=length(cpts)-1
				tmplike=0
				for(j in 1:nseg){
					tmplike=tmplike+mll.meanvarg(y[cpts[j+1]+1]-y[cpts[j]+1],cpts[j+1]-cpts[j],shape)
				}
				like=c(tmplike,tmplike+(nseg-1)*pen.value(object))
				names(like)=c("like","likepen")
			}
		}
		else if(distribution(object)=="Exponential"){
			if(cpttype(object)!="mean and variance"){
				stop("Unknown changepoint type for dist='Exponential', must be 'mean and variance'")
			}
			else{
			  mll.meanvare=function(x,n){
			    return(n*log(n)-n*log(x))
			  }
				y=c(0,cumsum(data.set(object)))
				cpts=c(0,cpts(object))
				nseg=length(cpts)-1
				tmplike=0
				for(j in 1:nseg){
					tmplike=tmplike+mll.meanvare(y[cpts[j+1]+1]-y[cpts[j]+1],cpts[j+1]-cpts[j])
				}
				like=c(tmplike,tmplike+(nseg-1)*pen.value(object))
				names(like)=c("like","likepen")
			}
		}
		else if(distribution(object)=="Poisson"){
		  if(cpttype(object)!="mean and variance"){
		    stop("Unknown changepoint type for dist='Poisson', must be 'mean and variance'")
		  }
		  else{
		    mll.meanvarp=function(x,n){
		      return(x*log(x)-x*log(n))
		    }
		    y=c(0,cumsum(data.set(object)))
		    cpts=c(0,cpts(object))
		    nseg=length(cpts)-1
		    tmplike=0
		    for(j in 1:nseg){
		      tmplike=tmplike+mll.meanvarp(y[cpts[j+1]+1]-y[cpts[j]+1],cpts[j+1]-cpts[j])
		    }
		    like=c(tmplike,tmplike+(nseg-1)*pen.value(object))
		    names(like)=c("like","likepen")
		  }
		}
		else{stop("Likelihood is only valid for distributional assumptions, not CUSUM or CSS")}
		return(like)
	})
