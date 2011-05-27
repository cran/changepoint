single.mean.norm.calc <-
function(data,extrainf=TRUE){
  singledim=function(data,extrainf=TRUE){
    n=length(data)
    y=c(0,cumsum(data))
    y2=c(0,cumsum(data^2))
    null=y2[n+1]-y[n+1]^2/n
    taustar=1:(n-1)
    tmp=y2[taustar+1]-y[taustar+1]^2/taustar + (y2[n+1]-y2[taustar+1]) - ((y[n+1]-y[taustar+1])^2)/(n-taustar)
    
    tau=which(tmp==min(tmp))[1]
    taulike=tmp[tau]
    if(extrainf==TRUE){
      out=c(tau,null,taulike)
      names(out)=c('cpt','null','alt')
      return(out)
    }
    else{
      return(tau)
    }
  }
    

  if(is.null(dim(data))==TRUE){
    # single data set
    cpt=singledim(data,extrainf)
    return(cpt)
  }
  else{
    rep=nrow(data)
    n=ncol(data)
    cpt=NULL
    if(extrainf==FALSE){
      for(i in 1:rep){
        cpt[i]=singledim(data[i,],extrainf)
      }
    }
    else{
      cpt=matrix(0,ncol=3,nrow=rep)
      for(i in 1:rep){
        cpt[i,]=singledim(data[i,],extrainf)
      }
      colnames(cpt)=c('cpt','null','alt')
    }
    return(cpt)
  }
}

single.mean.norm<-function(data,penalty="SIC",value=0,class=TRUE,param.estimates=TRUE){
	if(is.null(dim(data))==TRUE){ # single dataset
		n=length(data)
		if(penalty=="Asymptotic"){
			alpha=value
			alogn=(2*log(log(n)))^(-(1/2))
			blogn=(alogn^(-1))+(1/2)*alogn*log(log(log(n)))
			value=(-alogn*log(log((1-alpha+exp(-2*(pi^(1/2))*exp(blogn/alogn)))^(-1/(2*(pi^(1/2))))))+blogn)^2
		}
		tmp=single.mean.norm.calc(data,extrainf=TRUE)
		ans=decision(tmp[1],tmp[2],tmp[3],penalty,n,diffparam=1,value)
		if(class==TRUE){
			out=new("cpt")
			data.set(out)=data; cpttype(out)="mean"; method(out)="AMOC"; distribution(out)="Normal"; pen.type(out)=penalty; pen.value(out)=ans$pen; ncpts.max(out)=1
			if(ans$cpt != n){cpts(out)=c(ans$cpt,n)}
			else{cpts(out)=ans$cpt}
			if(param.estimates==TRUE){
				out=param(out)
			}
			return(out)
		}
		else{ return(ans$cpt)}
	}
	else{ 
		n=ncol(data)
		if(penalty=="Asymptotic"){
			alpha=value
			alogn=(2*log(log(n)))^(-(1/2))
			blogn=(alogn^(-1))+(1/2)*alogn*log(log(log(n)))
			value=(-alogn*log(log((1-alpha+exp(-2*(pi^(1/2))*exp(blogn/alogn)))^(-1/(2*(pi^(1/2))))))+blogn)^2
		}
		tmp=single.mean.norm.calc(data,extrainf=TRUE)
		ans=decision(tmp[,1],tmp[,2],tmp[,3],penalty,n,diffparam=1,value)
		if(class==TRUE){
			rep=nrow(data)
			out=list()
			for(i in 1:rep){
				out[[i]]=new("cpt")
				data.set(out[[i]])=data[i,];cpttype(out[[i]])="mean";method(out[[i]])="AMOC"; distribution(out[[i]])="Normal"; pen.type(out[[i]])=penalty;pen.value(out[[i]])=ans$pen;ncpts.max(out[[i]])=1
				if(ans$cpt[i] != n){cpts(out[[i]])=c(ans$cpt[i],n)}
				else{cpts(out[[i]])=ans$cpt[i]}
				if(param.estimates==TRUE){
					out[[i]]=param(out[[i]])
				}
			}
			return(out)
		}
		else{ return(ans$cpt)}
	}
}





single.var.norm.calc <-
function(data,know.mean=FALSE,mu=NA,extrainf=TRUE){
  singledim=function(data,know.mean=FALSE,mu=-1000,extrainf=TRUE){
    n=length(data)
    if((know.mean==FALSE)&(is.na(mu))){
	mu=mean(data)
    }
    y=c(0,cumsum((data-mu)^2))
    null=n*log(y[n+1]/n)
    taustar=1:(n-1)
    sigma1=y[taustar+1]/taustar
    neg=sigma1<=0
    sigma1[neg==TRUE]=1*10^(-10)
    sigman=(y[n+1]-y[taustar+1])/(n-taustar)
    neg=sigman<=0
    sigman[neg==TRUE]=1*10^(-10)
    tmp=taustar*log(sigma1) + (n-taustar)*log(sigman)
    
    tau=which(tmp==min(tmp))[1]
    taulike=tmp[tau]
    if(extrainf==TRUE){
      out=c(tau,null,taulike)
      names(out)=c('cpt','null','alt')
      return(out)
    }
    else{
      return(tau)
    }
  }
    

  if(is.null(dim(data))==TRUE){
    # single data set
    cpt=singledim(data,know.mean,mu,extrainf)
    return(cpt)
  }
  else{
    rep=nrow(data)
    n=ncol(data)
    if(length(mu)==1){
	mu=rep(mu,rep)
    }
    cpt=NULL
    if(extrainf==FALSE){
      for(i in 1:rep){
        cpt[i]=singledim(data[i,],know.mean,mu[i],extrainf)
      }
    }
    else{
      cpt=matrix(0,ncol=3,nrow=rep)
      for(i in 1:rep){
        cpt[i,]=singledim(data[i,],know.mean,mu[i],extrainf)
      }
      colnames(cpt)=c('cpt','null','alt')
    }
    return(cpt)
  }
}


single.var.norm<-function(data,penalty="SIC",value=0,know.mean=FALSE,mu=NA,class=TRUE,param.estimates=TRUE){
	if(is.null(dim(data))==TRUE){
		n=length(data)
		if(penalty=="Asymptotic"){
			alpha=value
			alogn=sqrt(2*log(log(n)))
			blogn=2*log(log(n))+ (log(log(log(n))))/2 - log(gamma(1/2))
			value=(-(log(log((1-alpha+exp(-2*exp(blogn)))^(-1/2))))/alogn + blogn/alogn)^2
		}
		tmp=single.var.norm.calc(data,know.mean,mu,extrainf=TRUE)
		ans=decision(tmp[1],tmp[2],tmp[3],penalty,n,diffparam=1,value)
		if(class==TRUE){
			out=new("cpt")
			data.set(out)=data; cpttype(out)="variance"; method(out)="AMOC"; distribution(out)="Normal"; pen.type(out)=penalty; pen.value(out)=ans$pen;ncpts.max(out)=1
			if(ans$cpt != n){cpts(out)=c(ans$cpt,n)}
			else{cpts(out)=ans$cpt}
			if(param.estimates==TRUE){
				out=param(out)
			}
			return(out)
		}
		else{ return(ans$cpt)}
	}
	else{ 
		n=ncol(data)
		if(penalty=="Asymptotic"){
			alpha=value
			alogn=sqrt(2*log(log(n)))
			blogn=2*log(log(n))+ (log(log(log(n))))/2 - log(gamma(1/2))
			value=(-(log(log((1-alpha+exp(-2*exp(blogn)))^(-1/2))))/alogn + blogn/alogn)^2
		}
		tmp=single.var.norm.calc(data,know.mean,mu,extrainf=TRUE)
		ans=decision(tmp[,1],tmp[,2],tmp[,3],penalty,n,diffparam=1,value)
		if(class==TRUE){
			rep=nrow(data)
			out=list()
			for(i in 1:rep){
				out[[i]]=new("cpt")
				data.set(out[[i]])=data[i,]; cpttype(out[[i]])="variance"; method(out[[i]])="AMOC"; distribution(out[[i]])="Normal"; pen.type(out[[i]])=penalty;pen.value(out[[i]])=ans$pen;ncpts.max(out[[i]])=1
				if(ans$cpt[i] != n){cpts(out[[i]])=c(ans$cpt[i],n)}
				else{cpts(out[[i]])=ans$cpt[i]}
				if(param.estimates==TRUE){
					out[[i]]=param(out[[i]])
				}
			}
			return(out)
		}
		else{ return(ans$cpt)}
	}
}






single.meanvar.norm.calc <-
function(data,extrainf=TRUE){
  singledim=function(data,extrainf=TRUE){
    n=length(data)
    y=c(0,cumsum(data))
    y2=c(0,cumsum((data)^2))
    null=n*log((y2[n+1]-(y[n+1]^2/n))/n)
    taustar=2:(n-1)
    sigma1=((y2[taustar+1]-(y[taustar+1]^2/taustar))/taustar)
    neg=sigma1<=0
    sigma1[neg==TRUE]=1*10^(-10)
    sigman=((y2[n+1]-y2[taustar+1])-((y[n+1]-y[taustar+1])^2/(n-taustar)))/(n-taustar)
    neg=sigman<=0
    sigman[neg==TRUE]=1*10^(-10)
    tmp=taustar*log(sigma1) + (n-taustar)*log(sigman)
    
    tau=which(tmp==min(tmp))[1]
    taulike=tmp[tau]
    if(extrainf==TRUE){
      out=c(tau,null,taulike)
      names(out)=c('cpt','null','alt')
      return(out)
    }
    else{
      return(tau)
    }
  }
    

  if(is.null(dim(data))==TRUE){
    # single data set
    cpt=singledim(data,extrainf)
    return(cpt)
  }
  else{
    rep=nrow(data)
    n=ncol(data)
    cpt=NULL
    if(extrainf==FALSE){
      for(i in 1:rep){
        cpt[i]=singledim(data[i,],extrainf)
      }
    }
    else{
      cpt=matrix(0,ncol=3,nrow=rep)
      for(i in 1:rep){
        cpt[i,]=singledim(data[i,],extrainf)
      }
      colnames(cpt)=c('cpt','null','alt')
    }
    return(cpt)
  }
}

single.meanvar.norm<-function(data,penalty="SIC",value=0,class=TRUE,param.estimates=TRUE){
	if(is.null(dim(data))==TRUE){
		n=length(data)
		if(penalty=="Asymptotic"){
			alpha=value
			alogn=sqrt(2*log(log(n)))
			blogn=2*log(log(n))+ log(log(log(n)))
			value=(-(log(log((1-alpha+exp(-2*exp(blogn)))^(-1/2))))/alogn + blogn/alogn)^2
		}
		tmp=single.meanvar.norm.calc(data,extrainf=TRUE)
		ans=decision(tmp[1],tmp[2],tmp[3],penalty,n,diffparam=2,value)
		if(class==TRUE){
			out=new("cpt")
			data.set(out)=data;cpttype(out)="mean and variance";method(out)="AMOC";distribution(out)="Normal";pen.type(out)=penalty; pen.value(out)=ans$pen;ncpts.max(out)=1
			if(ans$cpt != n){cpts(out)=c(ans$cpt,n)}
			else{cpts(out)=ans$cpt}
			if(param.estimates==TRUE){
				out=param(out)
			}
			return(out)
		}
		else{ return(ans$cpt)}
	}
	else{ 
		n=ncol(data)
		if(penalty=="Asymptotic"){
			alpha=value
			alogn=sqrt(2*log(log(n)))
			blogn=2*log(log(n))+ log(log(log(n)))
			value=(-(log(log((1-alpha+exp(-2*exp(blogn)))^(-1/2))))/alogn + blogn/alogn)^2
		}
		tmp=single.meanvar.norm.calc(data,extrainf=TRUE)
		ans=decision(tmp[,1],tmp[,2],tmp[,3],penalty,n,diffparam=2,value)
		if(class==TRUE){
			rep=nrow(data)
			out=list()
			for(i in 1:rep){
				out[[i]]=new("cpt")
				data.set(out[[i]])=data[i,];cpttype(out[[i]])="mean and variance";method(out[[i]])="AMOC";distribution(out[[i]])="Normal"; pen.type(out[[i]])=penalty;pen.value(out[[i]])=ans$pen;ncpts.max(out[[i]])=1
				if(ans$cpt[i] != n){cpts(out[[i]])=c(ans$cpt[i],n)}
				else{cpts(out[[i]])=ans$cpt[i]}
				if(param.estimates==TRUE){
					out[[i]]=param(out[[i]])
				}
			}
			return(out)
		}
		else{ return(ans$cpt)}
	}
}

