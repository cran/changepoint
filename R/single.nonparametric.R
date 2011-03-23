single.var.css.calc <-
function(data,extrainf=TRUE){
  singledim=function(data,extrainf=TRUE){
    n=length(data)
    y2=c(0,cumsum(data^2))
    taustar=1:(n-1)
		tmp=(y2[taustar+1]/y2[n+1])-taustar/n
    
		D=max(abs(tmp))
		tau=which.max(abs(tmp))
    if(extrainf==TRUE){
      out=c(tau,sqrt(n/2)*D)
      names(out)=c('cpt','test statistic')
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
      cpt=matrix(0,ncol=2,nrow=rep)
      for(i in 1:rep){
        cpt[i,]=singledim(data[i,],extrainf)
      }
      colnames(cpt)=c('cpt','test statistic')
    }
    return(cpt)
  }
}


single.var.css<-function(data,penalty="SIC",value=0,class=TRUE,param.estimates=TRUE){
	if(length(value)>1){stop('Only one dimensional penalties can be used for CSS')}
	diffparam=1
	if(is.null(dim(data))==TRUE){
		# single dataset
		n=length(data)
	}
	else{
		n=ncol(data)
	}
	if((penalty=="SIC") || (penalty=="BIC")){
		value=log(diffparam*log(n))
	}
	else if(penalty=="AIC"){
		value=log(2*diffparam)
	}
	else if(penalty=="Hannan-Quinn"){
		value=log(2*diffparam*log(log(n)))
	}
	else if(penalty=="None"){
		value=0
	}
	else if((penalty!="Manual")&&(penalty!="Asymptotic")){
		stop('Unknown Penalty')
	}
	if((penalty=="Manual")&&(is.numeric(value)==FALSE)){
		value=try(eval(parse(text=paste(value))),silent=TRUE)
		if(class(value)=='try-error'){
			stop('Your manual penalty cannot be evaluated')
		}
	}
	if(penalty=="Asymptotic"){
		if(value==0.01){value=1.628}
		else if(value==0.05){value=1.358}
		else if(value==0.1){value=1.224}
		else if(value==0.25){value=1.019}
		else if(value==0.5){value=0.828}
		else if(value==0.75){value=0.677}
		else if(value==0.9){value=0.571}
		else if(value==0.95){value=0.520}
		else{stop('Only alpha values of 0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95 are valid for CSS')}
	}
	if(is.null(dim(data))==TRUE){
		tmp=single.var.css.calc(data,extrainf=TRUE)
		ans=decision(tau=tmp[1],null=tmp[2],penalty="Manual",n=n,diffparam=1,value=value)
		if(class==TRUE){
			out=new("cpt")
			data.set(out)=data; cpttype(out)="variance"; method(out)="AMOC"; distribution(out)="CSS"; pen.type(out)=penalty; pen.value(out)=ans$pen;ncpts.max(out)=1
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
		tmp=single.var.css.calc(data,extrainf=TRUE)
		ans=decision(tau=tmp[,1],null=tmp[,2],penalty="Manual",n=n,diffparam=1,value=value)
		if(class==TRUE){
			rep=nrow(data)
			out=list()
			for(i in 1:rep){
				out[[i]]=new("cpt")
				data.set(out[[i]])=data[i,]; cpttype(out[[i]])="variance"; method(out[[i]])="AMOC"; distribution(out[[i]])="CSS"; pen.type(out[[i]])=penalty;pen.value(out[[i]])=ans$pen;ncpts.max(out[[i]])=1
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










single.mean.cusum.calc <-
function(data,extrainf=TRUE){
  singledim=function(data,extrainf=TRUE){
    n=length(data)
		ybar=mean(data)
    y=c(0,cumsum(data-ybar))
		y=y/n
    
		M=max(abs(y))
		tau=which.max(abs(y))
    if(extrainf==TRUE){
      out=c(tau,M)
      names(out)=c('cpt','test statistic')
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
      cpt=matrix(0,ncol=2,nrow=rep)
      for(i in 1:rep){
        cpt[i,]=singledim(data[i,],extrainf)
      }
      colnames(cpt)=c('cpt','test statistic')
    }
    return(cpt)
  }
}


single.mean.cusum<-function(data,penalty="Asymptotic",value=0.05,class=TRUE,param.estimates=TRUE){
	if(length(value)>1){stop('Only one dimensional penalties can be used for CUSUM')}
	if(penalty=="Asymptotic"){
		stop('No Asymptotic penalty is available for CUSUM')
	}
	if(is.null(dim(data))==TRUE){
		n=length(data)
		tmp=single.mean.cusum.calc(data,extrainf=TRUE)
		ans=decision(tau=tmp[1],null=tmp[2],penalty=penalty,n=n,diffparam=1,value=value)
		if(class==TRUE){
			out=new("cpt")
			data.set(out)=data; cpttype(out)="mean"; method(out)="AMOC"; distribution(out)="CUSUM"; pen.type(out)=penalty; pen.value(out)=ans$pen;ncpts.max(out)=1
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
		tmp=single.mean.cusum.calc(data,extrainf=TRUE)
		ans=decision(tau=tmp[,1],null=tmp[,2],penalty=penalty,n=n,diffparam=1,value=value)
		if(class==TRUE){
			rep=nrow(data)
			out=list()
			for(i in 1:rep){
				out[[i]]=new("cpt")
				data.set(out[[i]])=data[i,]; cpttype(out[[i]])="mean"; method(out[[i]])="AMOC"; distribution(out[[i]])="CUSUM"; pen.type(out[[i]])=penalty;pen.value(out[[i]])=ans$pen;ncpts.max(out[[i]])=1
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

