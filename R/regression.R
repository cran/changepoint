single.reg.norm.calc <-
function(data,extrainf=TRUE){
  singledim=function(data,extrainf=TRUE){
    n=nrow(data)
    null.beta=solve(t(data[,2:ncol(data)])%*%data[,2:ncol(data)],t(data[,2:ncol(data)])%*%data[,1])
    ymxb=data[,1]-data[,2:ncol(data)]%*%null.beta
    null=n*log(t(ymxb)%*%ymxb)
    tmp=NULL
    for(taustar in ncol(data):(n-ncol(data)+1)){
			tmp.beta1=solve(t(data[1:taustar,2:ncol(data)])%*%data[1:taustar,2:ncol(data)],t(data[1:taustar,2:ncol(data)])%*%data[1:taustar,1])
			ymxb1=data[1:taustar,1]-data[1:taustar,2:ncol(data)]%*%tmp.beta1
			tmp.betan=solve(t(data[(taustar+1):n,2:ncol(data)])%*%data[(taustar+1):n,2:ncol(data)],t(data[(taustar+1):n,2:ncol(data)])%*%data[(taustar+1):n,1])
			ymxbn=data[(taustar+1):n,1]-data[(taustar+1):n,2:ncol(data)]%*%tmp.betan
      tmp[taustar-ncol(data)]=n*log(t(ymxb1)%*%ymxb1 + t(ymxbn)%*%ymxbn)
    }
    tau=which(tmp==min(tmp))[1]+ncol(data)
    taulike=tmp[tmp==min(tmp)]
    if(extrainf==TRUE){
      return(c(tau,null,taulike))
    }
    else{
      return(tau)
    }
  }
    

  if(length(dim(data))==2){
    # single data set
    cpt=singledim(data,extrainf)
    return(cpt)
  }
  else{
    rep=dim(data)[1]
    n=dim(data)[2]
    cpt=NULL
    if(extrainf==FALSE){
      for(i in 1:rep){
        cpt[i]=singledim(data[i,,],extrainf)
      }
    }
    else{
      cpt=matrix(0,ncol=3,nrow=rep)
      for(i in 1:rep){
        cpt[i,]=singledim(data[i,,],extrainf)
      }
    }
    return(cpt)
  }
}


single.reg.norm<-function(data,penalty="SIC",value=0,class=TRUE,param.estimates=TRUE){
	if(length(dim(data))==2){
		n=nrow(data)
		if(penalty=="Asymptotic"){
			d=ncol(data)-1
			alpha=value
			top=-(log(log((1 - alpha + exp(-2*exp(2*(log(log(n)))+(d/2)*(log(log(log(n))))- log(gamma(d/2)))))^(-1/2))))  +  2*(log(log(n)))+(d/2)*(log(log(log(n))))- log(gamma(d/2))
			bottom=(2*log(log(n)))^(1/2)
			value=(top/bottom)^2
		}
		tmp=single.reg.norm.calc(data,extrainf=TRUE)
		ans=decision(tmp[1],tmp[2],tmp[3],penalty,n,diffparam=ncol(data)-1,value)
		if(class==TRUE){
			out=new("cpt.reg")
			data.set(out)=data;method(out)="AMOC";distribution(out)="Normal";pen.type(out)=penalty;pen.value(out)=ans$pen;ncpts.max(out)=1
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
		n=dim(data)[2]
		if(penalty=="Asymptotic"){
			d=dim(data)[3]-1
			alpha=value
			top=-(log(log((1 - alpha + exp(-2*exp(2*(log(log(n)))+(d/2)*(log(log(log(n))))- log(gamma(d/2)))))^(-1/2))))  +  2*(log(log(n)))+(d/2)*(log(log(log(n))))- log(gamma(d/2))
			bottom=(2*log(log(n)))^(1/2)
			value=(top/bottom)^2
		}
		tmp=single.reg.norm.calc(data,extrainf=TRUE)
		ans=decision(tmp[,1],tmp[,2],tmp[,3],penalty,n,diffparam=dim(data)[3]-1,value)
		if(class==TRUE){
			rep=nrow(data)
			out=list()
			for(i in 1:rep){
				out[[i]]=new("cpt.reg")
				data.set(out[[i]])=data[i,,];method(out[[i]])="AMOC";distribution(out[[i]])="Normal"; pen.type(out[[i]])=penalty;pen.value(out[[i]])=ans$pen;ncpts.max(out[[i]])=1
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




PELT.reg.norm=function(data,pen=0,nprune=FALSE){
  mll.reg=function(x,p){
    betaest=solve(t(x[,2:p])%*%x[,2:p],t(x[,2:p])%*%x[,1])
 	  ymxb=x[,1]-x[,2:p]%*%betaest	

		return(t(ymxb)%*%ymxb)
  }
  n=nrow(data)
	p=ncol(data)-1

  lastchangecpts=matrix(NA,nrow=n,ncol=2)
  lastchangelike=NULL
  checklist=NULL
	ncpts=NULL
	for(i in p:(2*p)){
		lastchangelike[i]=mll.reg(data[1:i,],p+1)
		lastchangecpts[i,]=c(0,i)
	}
	ncpts[1:(p-1)]=0
	ncpts[p:(2*p)]=1

  noprune=NULL
  for(tstar in (2*p+1):n){
    tmplike=NULL
    tmpt=c(checklist, tstar-p)
		for(i in 1:length(tmpt)){
	    tmplike[i]=tstar*log(lastchangelike[tmpt[i]]+mll.reg(data[(tmpt[i]+1):tstar,],p+1))+(ncpts[tmpt[i]]+1)*pen
		}
		null=tstar*log(mll.reg(data[1:tstar,],p+1))
		
    minlike=min(c(tmplike,null),na.rm=TRUE)

    if(minlike==null){
			lastchangelike[tstar]=mll.reg(data[1:tstar,],p+1)
      lastchangecpts[tstar,]=c(0,tstar)
			ncpts[tstar]=0
    }
    else{
      cpt=tmpt[tmplike==minlike][1]
			lastchangelike[tstar]=lastchangelike[cpt]+mll.reg(data[(cpt+1):tstar,],p+1)
      lastchangecpts[tstar,]=c(cpt,tstar)
			ncpts[tstar]=ncpts[cpt]+1
    }
    checklist=tmpt[tmplike<=minlike+pen]
    if(nprune==TRUE){
      noprune=c(noprune,length(checklist))
    }
  }
  if(nprune==TRUE){
    return(nprune=noprune)
  }
  else{
    fcpt=NULL
    last=n
    while(last!=0){
	fcpt=c(fcpt,lastchangecpts[last,2])
	last=lastchangecpts[last,1]
    }
    return(cpt=sort(fcpt))
  }
}





binseg.reg.norm=function(data,Q=5,pen=0){
  mll.reg=function(x,p){
    betaest=solve(t(x[,2:p])%*%x[,2:p],t(x[,2:p])%*%x[,1])
 	  ymxb=x[,1]-x[,2:p]%*%betaest	

		return(t(ymxb)%*%ymxb)
  }
  n=nrow(data)
	p=ncol(data)-1
  tau=c(0,n)
  cpt=matrix(0,nrow=2,ncol=Q)
  oldmax=1000

  for(q in 1:Q){
    lambda=rep(0,n-1)
    i=1
    st=tau[1]+1;end=tau[2]
    null=(-1/2)*(end-st+1)*log(mll.reg(data[st:end,],p+1))
    for(j in 1:(n-1)){
      if(j==end){
        st=end+1;i=i+1;end=tau[i+1]
        null=(-1/2)*(end-st+1)*log(mll.reg(data[st:end,],p+1))
      }else{
				if((j-st)<p){lambda[j]=-1*10^(100)}
				else if((end-j)<p){lambda[j]=-1*10^(100)}
				else{lambda[j]=(-1/2)*(end-st+1)*log(mll.reg(data[st:j,],p+1)+mll.reg(data[(j+1):end,],p+1))-null}
      }
    }
    k=which.max(lambda)[1]
    cpt[1,q]=k;cpt[2,q]=min(oldmax,max(lambda))
    oldmax=min(oldmax,max(lambda))
    tau=sort(c(tau,k))
  }
  op.cps=NULL
  for(i in 1:length(pen)){
    criterion=(2*cpt[2,])>=pen[i]
    if(sum(criterion)==0){
      op.cps=0
    }
    else{
      op.cps=c(op.cps,max(which((criterion)==TRUE)))
    }
  }
  return(list(cps=cpt,op.cpts=op.cps,pen=pen))
}


segneigh.reg.norm=function(data,Q=5,pen=0){
  n=nrow(data)
	p=ncol(data)-1
  all.seg=matrix(0,ncol=n,nrow=n)
  for(i in 1:(n-p)){
    for(j in (i+p):n){
	    betaest=solve(t(data[i:j,2:(p+1)])%*%data[i:j,2:(p+1)],t(data[i:j,2:(p+1)])%*%data[i:j,1])
	 	  ymxb=data[i:j,1]-data[i:j,2:(p+1)]%*%betaest	
      all.seg[i,j]=t(ymxb)%*%ymxb
    }
  }
  like.Q=matrix(0,ncol=n,nrow=Q)
  like.Q[1,]=(1:n)*log(all.seg[1,])
	like.Q[1,1:p]=1*10^{1000}
  cp=matrix(0,ncol=n,nrow=Q)
  for(q in 2:Q){
    for(j in q:n){
      like=NULL
      if((j-p-(q-1)*p)<0){v=(q-1)*p}
      else{v=((q-1)*p):(j-p)}
      like=(j)*log(exp(like.Q[q-1,v]/v)+all.seg[v+1,j])

      like.Q[q,j]= min(like,na.rm=TRUE)
      cp[q,j]=which(like==min(like,na.rm=TRUE))[1]+(q-1)*p-1
    }

  }
  cps.Q=matrix(0,ncol=Q,nrow=Q)
  for(q in 2:Q){
    cps.Q[q,1]=cp[q,n]
    for(i in 1:(q-1)){
      cps.Q[q,(i+1)]=cp[(q-i),cps.Q[q,i]]
    }
  }

  op.cps=NULL
   k=0:(Q-1)

  for(i in 1:length(pen)){
    criterion=like.Q[,n]+k*pen[i]

    op.cps=c(op.cps,which(criterion==min(criterion))-1)
  }
  return(list(cps=cps.Q,op.cpts=op.cps,pen=pen,like=criterion[op.cps+1]))
}



multiple.reg.norm=function(data,mul.method="PELT",penalty="SIC",value=0,Q=5,class=TRUE,param.estimates=TRUE){
	if(!((mul.method=="PELT")||(mul.method=="BinSeg")||(mul.method=="SegNeigh"))){
		stop("Multiple Method is not recognised")
	}
	if(is.null(dim(data))==TRUE){
		stop("Data must be a matrix with first column being the response and the remaining columns covariates")
	}
	if(length(dim(data))==2){
		# single dataset
		n=nrow(data)
		diffparam=ncol(data)-1
	}
	else{
		n=dim(data)[2]
		diffparam=dim(data)[3]-1
	}
	if((penalty=="SIC") || (penalty=="BIC")){
		value=diffparam*log(n)
	}
	else if(penalty=="AIC"){
		value=2*diffparam
	}
	else if(penalty=="Hannan-Quinn"){
		value=2*diffparam*log(log(n))
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
		alpha=value
	  top= -(log(log((1-alpha+exp(-2*exp(2*(log(log(n)))+(diffparam/2)*(log(log(log(n))))-log(gamma(diffparam/2)))))^(-1/2))))+2*(log(log(n)))+(diffparam/2)*(log(log(log(n))))-log(gamma(diffparam/2))
	  bottom=(2*log(log(n)))^(1/2)
		value=(top/bottom)^2 + 2*(log(n))
	}
	if(length(dim(data))==2){
		# single dataset
		if(mul.method=="PELT"){
			out=PELT.reg.norm(data,value,nprune=FALSE)
			cpts=out
		}
		else if(mul.method=="BinSeg"){
			out=binseg.reg.norm(data,Q,value)
			cpts=c(sort(out$cps[1,1:out$op.cpts]),n)
		}
		else if(mul.method=="SegNeigh"){
			out=segneigh.reg.norm(data,Q,value)
			if(out$op.cpts==0){cpts=n}
			else{cpts=c(sort(out$cps[out$op.cpts+1,][out$cps[out$op.cpts+1,]>0]),n)}
		}
		if(class==TRUE){
			ans=new("cpt.reg")
			data.set(ans)=data;cpttype(ans)="regression";method(ans)=mul.method; distribution(ans)="Normal";pen.type(ans)=penalty;pen.value(ans)=value;cpts(ans)=cpts
			if(mul.method=="PELT"){
				ncpts.max(ans)=Inf
			}
			else{
				ncpts.max(ans)=Q
			}
			if(param.estimates==TRUE){
				ans=param(ans)
			}
			return(ans)
		}
		else{ return(out)}
	}
	else{
		rep=dim(data)[1]
		out=list()
		if(class==TRUE){cpts=list()}
		if(mul.method=="PELT"){
			for(i in 1:rep){
				out=c(out,list(PELT.reg.norm(data[i,,],value,nprune=FALSE)))
			}
			cpts=out
		}
		else if(mul.method=="BinSeg"){
			for(i in 1:rep){
				out=c(out,list(binseg.reg.norm(data[i,,],Q,value)))
				if(class==TRUE){cpts[[i]]=c(sort(out$cps[1,1:out$op.cpts]),n)}
			}
		}
		else if(mul.method=="SegNeigh"){
			for(i in 1:rep){
				out=c(out,list(segneigh.reg.norm(data[i,,],Q,value)))
				if(class==TRUE){
					if(out$op.cpts==0){cpts[[i]]=n}
					else{cpts[[i]]=c(sort(out$cps[out$op.cpts+1,][out$cps[out$op.cpts+1,]>0]),n)}
				}
			}
		}
		if(class==TRUE){
			ans=list()
			for(i in 1:rep){
				ans[[i]]=new("cpt.reg")
				data.set(ans[[i]])=data[i,,];cpttype(ans[[i]])="regression"; method(ans[[i]])=mul.method;distribution(ans[[i]])="Normal";pen.type(ans[[i]])=penalty;pen.value(ans[[i]])=value;cpts(ans[[i]])=cpts[[i]]
				if(mul.method=="PELT"){
					ncpts.max(ans[[i]])=Inf
				}
				else{
					ncpts.max(ans[[i]])=Q
				}
				if(param.estimates==TRUE){
					ans[[i]]=param(ans[[i]])
				}
			}
			return(ans)
		}
		else{return(out)}
	}
}

