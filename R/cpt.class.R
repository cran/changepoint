.initcpt<-function(where){
# initialization function for cpt and cpt.reg classes
	setClass("cpt",representation(data.set="numeric", cpttype="character", method="character", 	distribution="character",pen.type="character",pen.value="numeric",cpts="numeric",ncpts.max="numeric",param.est="list",date="character"),prototype(date=date()),where=where)

	setClass("cpt.reg",representation(data.set="matrix", cpttype="character", method="character", distribution="character",pen.type="character",pen.value="numeric",cpts="numeric",ncpts.max="numeric",param.est="list",date="character"),prototype(cpttype="regression",date=date()),where=where)

# retrival functions for slots
	if(!isGeneric("data.set")) {
		if (is.function("data.set")){
			fun <- data.set
		}
		else {fun <- function(object){
				standardGeneric("data.set")
			}
		}
		setGeneric("data.set", fun,where=where)
	}
	setMethod("data.set","cpt",function(object) object@data.set,where=where)
	setMethod("data.set","cpt.reg",function(object) object@data.set,where=where)

	if(!isGeneric("cpttype")) {
		if (is.function("cpttype")){
			fun <- cpttype
		}
		else {fun <- function(object){
				standardGeneric("cpttype")
			}
		}
		setGeneric("cpttype", fun,where=where)
	}
	setMethod("cpttype","cpt",function(object) object@cpttype,where=where)
	setMethod("cpttype","cpt.reg",function(object) object@cpttype,where=where)

	if(!isGeneric("method")) {
		if (is.function("method")){
			fun <- method
		}
		else {fun <- function(object){
				standardGeneric("method")
			}
		}
		setGeneric("method", fun,where=where)
	}
	setMethod("method","cpt",function(object) object@method,where=where)
	setMethod("method","cpt.reg",function(object) object@method,where=where)
	
	if(!isGeneric("distribution")) {
		if (is.function("distribution")){
			fun <- distribution
		}
		else {fun <- function(object){
				standardGeneric("distribution")
			}
		}
		setGeneric("distribution", fun,where=where)
	}
	setMethod("distribution","cpt",function(object) object@distribution,where=where)
	setMethod("distribution","cpt.reg",function(object) object@distribution,where=where)
	
	if(!isGeneric("pen.type")) {
		if (is.function("pen.type")){
			fun <- pen.type
		}
		else {fun <- function(object){
				standardGeneric("pen.type")
			}
		}
		setGeneric("pen.type", fun,where=where)
	}
	setMethod("pen.type","cpt",function(object) object@pen.type,where=where)
	setMethod("pen.type","cpt.reg",function(object) object@pen.type,where=where)
	
	if(!isGeneric("pen.value")) {
		if (is.function("pen.value")){
			fun <- pen.value
		}
		else {fun <- function(object){
				standardGeneric("pen.value")
			}
		}
		setGeneric("pen.value", fun,where=where)
	}
	setMethod("pen.value","cpt",function(object) object@pen.value,where=where)
	setMethod("pen.value","cpt.reg",function(object) object@pen.value,where=where)
	
	if(!isGeneric("cpts")) {
		if (is.function("cpts")){
			fun <- cpts
		}
		else {fun <- function(object){
				standardGeneric("cpts")
			}
		}
		setGeneric("cpts", fun,where=where)
	}
	setMethod("cpts","cpt",function(object) object@cpts,where=where)
	setMethod("cpts","cpt.reg",function(object) object@cpts,where=where)
	
	if(!isGeneric("ncpts.max")) {
		if (is.function("ncpts.max")){
			fun <- cpts
		}
		else {fun <- function(object){
				standardGeneric("ncpts.max")
			}
		}
		setGeneric("ncpts.max", fun,where=where)
	}
	setMethod("ncpts.max","cpt",function(object) object@ncpts.max,where=where)
	setMethod("ncpts.max","cpt.reg",function(object) object@ncpts.max,where=where)

	if(!isGeneric("param.est")) {
		if (is.function("param.est")){
			fun <- param.est
		}
		else {fun <- function(object){
				standardGeneric("param.est")
			}
		}
		setGeneric("param.est", fun,where=where)
	}
	setMethod("param.est","cpt",function(object) object@param.est,where=where)
	setMethod("param.est","cpt.reg",function(object) object@param.est,where=where)
	
# ncpts function
	if(!isGeneric("ncpts")) {
		if (is.function("ncpts")){
			fun <- ncpts
		}
		else {fun <- function(object){
				standardGeneric("ncpts")
			}
		}
		setGeneric("ncpts", fun,where=where)
	}
	setMethod("ncpts","cpt",function(object) length(cpts(object))-1,where=where)
	setMethod("ncpts","cpt.reg",function(object) length(cpts(object))-1,where=where)

# replacement functions for slots
	setGeneric("data.set<-", function(object, value) standardGeneric("data.set<-"),where=where)
	setReplaceMethod("data.set", "cpt", function(object, value) {
		object@data.set <- value
		return(object)
	},where=where)
	setReplaceMethod("data.set", "cpt.reg", function(object, value) {
		object@data.set <- value
		return(object)
	},where=where)
	
	setGeneric("cpttype<-", function(object, value) standardGeneric("cpttype<-"),where=where)
	setReplaceMethod("cpttype", "cpt", function(object, value) {
		object@cpttype <- value
		return(object)
	},where=where)
	setReplaceMethod("cpttype", "cpt.reg", function(object, value) {
		object@cpttype <- value
		return(object)
	},where=where)
	
	setGeneric("method<-", function(object, value) standardGeneric("method<-"),where=where)
	setReplaceMethod("method", "cpt", function(object, value) {
		object@method <- value
		return(object)
	},where=where)
	setReplaceMethod("method", "cpt.reg", function(object, value) {
		object@method <- value
		return(object)
	},where=where)
	
	setGeneric("distribution<-", function(object, value) standardGeneric("distribution<-"),where=where)
	setReplaceMethod("distribution", "cpt", function(object, value) {
		object@distribution <- value
		return(object)
	},where=where)
	setReplaceMethod("distribution", "cpt.reg", function(object, value) {
		object@distribution <- value
		return(object)
	},where=where)
	
	setGeneric("pen.type<-", function(object, value) standardGeneric("pen.type<-"),where=where)
	setReplaceMethod("pen.type", "cpt", function(object, value) {
		object@pen.type <- value
		return(object)
	},where=where)
	setReplaceMethod("pen.type", "cpt.reg", function(object, value) {
		object@pen.type <- value
		return(object)
	},where=where)
	
	setGeneric("pen.value<-", function(object, value) standardGeneric("pen.value<-"),where=where)
	setReplaceMethod("pen.value", "cpt", function(object, value) {
		object@pen.value <- value
		return(object)
	},where=where)
	setReplaceMethod("pen.value", "cpt.reg", function(object, value) {
		object@pen.value <- value
		return(object)
	},where=where)
	
	setGeneric("cpts<-", function(object, value) standardGeneric("cpts<-"),where=where)
	setReplaceMethod("cpts", "cpt", function(object, value) {
		object@cpts <- value
		return(object)
	},where=where)
	setReplaceMethod("cpts", "cpt.reg", function(object, value) {
		object@cpts <- value
		return(object)
	},where=where)

	setGeneric("ncpts.max<-", function(object, value) standardGeneric("ncpts.max<-"),where=where)
	setReplaceMethod("ncpts.max", "cpt", function(object, value) {
		object@ncpts.max <- value
		return(object)
	},where=where)
	setReplaceMethod("ncpts.max", "cpt.reg", function(object, value) {
		object@ncpts.max <- value
		return(object)
	},where=where)
	
	setGeneric("param.est<-", function(object, value) standardGeneric("param.est<-"),where=where)
	setReplaceMethod("param.est", "cpt", function(object, value) {
		object@param.est <- value
		return(object)
	},where=where)
	setReplaceMethod("param.est", "cpt.reg", function(object, value) {
		object@param.est <- value
		return(object)
	},where=where)

# parameter functions
	setGeneric("param", function(object,...) standardGeneric("param"),where=where)
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
			else{
				stop("Unknown distribution for a change in mean and variance")
			}
		}
		else{
			stop("Unknown changepoint type, must be 'mean', 'variance' or 'mean and variance'")
		}
		return(object)
	},where=where)
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
	},where=where)

# summary functions
	setMethod("summary","cpt",function(object){
	    cat("Changepoint type      : Change in",cpttype(object),'\n')
	    cat("Method of analysis    :",method(object),"\n")
	    cat("Assumed Distribution  :", distribution(object),"\n")
	    cat("Type of penalty       :", pen.type(object), "with value,",pen.value(object),"\n")
	    cat("Maximum no. of cpts   :", ncpts.max(object),"\n")
	    if(length(cpts(object))<=20){cat("Changepoint Locations :",cpts(object),"\n")}
	    else{cat("Number of changepoints:", length(cpts(object)),"\n")}
	},where=where)
	setMethod("summary","cpt.reg",function(object){
	    cat("Changepoint type     : Change in",cpttype(object),'\n')
	    cat("Method of analysis   :",method(object),"\n")
	    cat("Assumed Distribution :", distribution(object),"\n")
	    cat("Type of penalty      :", pen.type(object), "with value,",pen.value(object),"\n")
	    cat("Maximum no. of cpts   :", ncpts.max(object),"\n")
	    if(length(cpts(object))<=20){cat("Changepoint Locations :",cpts(object),"\n")}
	    else{cat("Number of changepoints:", length(cpts(object)),"\n")}
	},where=where)

# print functions
	setMethod("print","cpt",function(x){
	    cat("Class 'cpt' : Changepoint Object\n")
	    cat("       ~~   : S4 class containing", length(attributes(x))-1, "slots with names\n")
	    cat("             ", names(attributes(x))[1:(length(attributes(x))-1)], "\n\n")
	    cat("Created on  :", x@date, "\n\n")
	    cat("summary(.)  :\n----------\n")
	    summary(x)
	},where=where)
	setMethod("print","cpt.reg",function(x){
	    cat("Class 'cpt' : Changepoint Object\n")
	    cat("       ~~   : S4 class containing", length(attributes(x))-1, "slots with names\n")
	    cat("             ", names(attributes(x))[1:(length(attributes(x))-1)], "\n\n")
	    cat("Created on  :", x@date, "\n\n")
	    cat("summary(.)  :\n----------\n")
	    summary(x)
	},where=where)

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
	},where=where)

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
	},where=where)

# likelihood functions
	setGeneric("likelihood", function(object) standardGeneric("likelihood"),where=where)
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
		else{stop("Likelihood is only valid for distributional assumptions, not CUSUM or CSS")}
		return(like)
	},where=where)

	rm(fun)
}	



