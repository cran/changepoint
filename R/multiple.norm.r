PELT.var.norm=function(data,pen=0,know.mean=FALSE,mu=-1000,nprune=FALSE){
  mll.var.EFK=function(x,n){
    neg=x<=0
    x[neg==TRUE]=0.00000000001    
    return( n*(log(2*pi)+log(x/n)+1))
  }
  if((know.mean==FALSE)&(mu==-1000)){
	mu=mean(data)
  }
  n=length(data)
  y2=c(0,cumsum((data-mu)^2))

  lastchangecpts=matrix(NA,nrow=n,ncol=2)
  lastchangelike=matrix(NA,nrow=n,ncol=2)
  checklist=NULL
  lastchangelike[1,]=c(mll.var.EFK(y2[2],1),mll.var.EFK(y2[n+1]-y2[2],n-1)+pen)
  lastchangecpts[1,]=c(0,1)
  lastchangelike[2,]=c(mll.var.EFK(y2[3],2),mll.var.EFK(y2[n+1]-y2[3],n-2)+pen)
  lastchangecpts[2,]=c(0,2)
  lastchangelike[3,]=c(mll.var.EFK(y2[4],3),mll.var.EFK(y2[n+1]-y2[4],n-3)+pen)
  lastchangecpts[3,]=c(0,3)
  noprune=NULL
  for(tstar in 4:n){
    tmplike=NULL
    tmpt=c(checklist, tstar-2)
    tmplike=lastchangelike[tmpt,1]+mll.var.EFK(y2[tstar+1]-y2[tmpt+1],tstar-tmpt)+pen
    if(tstar==n){
      lastchangelike[tstar,]=c(min(c(tmplike,mll.var.EFK(y2[tstar+1]-y2[1],tstar)),na.rm=TRUE),0)
    }
    else{
      lastchangelike[tstar,]=c(min(c(tmplike,mll.var.EFK(y2[tstar+1]-y2[1],tstar)),na.rm=TRUE),mll.var.EFK(y2[n+1]-y2[tstar+1],n-tstar)+pen)
    }
    if(lastchangelike[tstar,1]==mll.var.EFK(y2[tstar+1]-y2[1],tstar)){
      lastchangecpts[tstar,]=c(0,tstar)
    }
    else{
      cpt=tmpt[tmplike==lastchangelike[tstar,1]][1]
      lastchangecpts[tstar,]=c(cpt,tstar)
    }
    checklist=tmpt[tmplike<=lastchangelike[tstar,1]+pen]
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


PELT.mean.norm=function(data,pen=0,nprune=FALSE){
  mll.mean.EFK=function(x2,x,n){
    return( x2-(x^2)/n)
  }
  n=length(data)
  y2=c(0,cumsum(data^2))
  y=c(0,cumsum(data))

  lastchangecpts=matrix(NA,nrow=n,ncol=2)
  lastchangelike=matrix(NA,nrow=n,ncol=2)
  checklist=NULL
  lastchangelike[1,]=c(mll.mean.EFK(y2[2],y[2],1),mll.mean.EFK(y2[n+1]-y2[2],y[n+1]-y[2],n-1)+pen)
  lastchangecpts[1,]=c(0,1)
  lastchangelike[2,]=c(mll.mean.EFK(y2[3],y[3],2),mll.mean.EFK(y2[n+1]-y2[3],y[n+1]-y[3],n-2)+pen)
  lastchangecpts[2,]=c(0,2)
  lastchangelike[3,]=c(mll.mean.EFK(y2[4],y[4],3),mll.mean.EFK(y2[n+1]-y2[4],y[n+1]-y[4],n-3)+pen)
  lastchangecpts[3,]=c(0,3)
  noprune=NULL
  for(tstar in 4:n){
    tmplike=NULL
    tmpt=c(checklist, tstar-2)
    tmplike=lastchangelike[tmpt,1]+mll.mean.EFK(y2[tstar+1]-y2[tmpt+1],y[tstar+1]-y[tmpt+1],tstar-tmpt)+pen
    if(tstar==n){
      lastchangelike[tstar,]=c(min(c(tmplike,mll.mean.EFK(y2[tstar+1],y[tstar+1],tstar)),na.rm=TRUE),0)
    }
    else{
      lastchangelike[tstar,]=c(min(c(tmplike,mll.mean.EFK(y2[tstar+1],y[tstar+1],tstar)),na.rm=TRUE),mll.mean.EFK(y2[n+1]-y2[tstar+1],y[n+1]-y[tstar+1],n-tstar)+pen)
    }
    if(lastchangelike[tstar,1]==mll.mean.EFK(y2[tstar+1],y[tstar+1],tstar)){
      lastchangecpts[tstar,]=c(0,tstar)
    }
    else{
      cpt=tmpt[tmplike==lastchangelike[tstar,1]][1]
      lastchangecpts[tstar,]=c(cpt,tstar)
    }
    checklist=tmpt[tmplike<=lastchangelike[tstar,1]+pen]
    if(nprune==TRUE){
      noprune=c(noprune,length(checklist))
    }
  }
  if(nprune==TRUE){
    return(noprune)
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

PELT.meanvar.norm=function(data,pen=0,nprune=FALSE){
  mll.meanvar.EFK=function(x2,x,n){
    sigmasq=(1/n)*(x2-(x^2)/n)
    neg=sigmasq<=0
    sigmasq[neg==TRUE]=0.00000000001
    return(n*(log(2*pi)+log(sigmasq)+1))
  }
  n=length(data)
  y2=c(0,cumsum(data^2))
  y=c(0,cumsum(data))

  lastchangecpts=matrix(NA,nrow=n,ncol=2)
  lastchangelike=matrix(NA,nrow=n,ncol=2)
  checklist=NULL
  lastchangelike[1,]=c(mll.meanvar.EFK(y2[2],y[2],1),mll.meanvar.EFK(y2[n+1]-y2[2],y[n+1]-y[2],n-1)+pen)
  lastchangecpts[1,]=c(0,1)
  lastchangelike[2,]=c(mll.meanvar.EFK(y2[3],y[3],2),mll.meanvar.EFK(y2[n+1]-y2[3],y[n+1]-y[3],n-2)+pen)
  lastchangecpts[2,]=c(0,2)
  lastchangelike[3,]=c(mll.meanvar.EFK(y2[4],y[4],3),mll.meanvar.EFK(y2[n+1]-y2[4],y[n+1]-y[4],n-3)+pen)
  lastchangecpts[3,]=c(0,3)
  noprune=NULL
  for(tstar in 4:n){
    tmplike=NULL
    tmpt=c(checklist, tstar-2)
    tmplike=lastchangelike[tmpt,1]+mll.meanvar.EFK(y2[tstar+1]-y2[tmpt+1],y[tstar+1]-y[tmpt+1],tstar-tmpt)+pen
    if(tstar==n){
      lastchangelike[tstar,]=c(min(c(tmplike,mll.meanvar.EFK(y2[tstar+1],y[tstar+1],tstar)),na.rm=TRUE),0)
    }
    else{
      lastchangelike[tstar,]=c(min(c(tmplike,mll.meanvar.EFK(y2[tstar+1],y[tstar+1],tstar)),na.rm=TRUE),mll.meanvar.EFK(y2[n+1]-y2[tstar+1],y[n+1]-y[tstar+1],n-tstar)+pen)
    }
    if(lastchangelike[tstar,1]==mll.meanvar.EFK(y2[tstar+1],y[tstar+1],tstar)){
      lastchangecpts[tstar,]=c(0,tstar)
    }
    else{
      cpt=tmpt[tmplike==lastchangelike[tstar,1]][1]
      lastchangecpts[tstar,]=c(cpt,tstar)
    }
    checklist=tmpt[tmplike<=lastchangelike[tstar,1]+pen]
    if(nprune==TRUE){
      noprune=c(noprune,length(checklist))
    }
  }
  if(nprune==TRUE){
    return(noprune)
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


segneigh.var.norm=function(data,Q=5,pen=0,know.mean=FALSE,mu=-1000){
  n=length(data)
  if((know.mean==FALSE)&(mu==-1000)){
	mu=mean(data)
  }
  all.seg=matrix(0,ncol=n,nrow=n)
  for(i in 1:n){
	ssq=0
    for(j in i:n){
        m=j-i+1
        ssq=ssq+(data[j]-mu)^2
	if(ssq<=0){sigmasq=0.00000000001/m}
        else{sigmasq=ssq/m}
        all.seg[i,j]=-(m/2)*(log(2*pi)+log(sigmasq)+1)
    }
  }
  like.Q=matrix(0,ncol=n,nrow=Q)
  like.Q[1,]=all.seg[1,]
  cp=matrix(0,ncol=n,nrow=Q)
  for(q in 2:Q){
    for(j in q:n){
      like=NULL
      v=(q-1):(j-1)
      like=like.Q[q-1,v]+all.seg[v+1,j]

      like.Q[q,j]= max(like,na.rm=TRUE)
      cp[q,j]=which(like==max(like,na.rm=TRUE))[1]+(q-2)
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
    criterion=-2*like.Q[,n]+k*pen[i]

    op.cps=c(op.cps,which(criterion==min(criterion))-1)
  }
  return(list(cps=cps.Q,op.cpts=op.cps,pen=pen,like=criterion[op.cps+1]))
}


segneigh.mean.norm=function(data,Q=5,pen=0){
  n=length(data)
  all.seg=matrix(0,ncol=n,nrow=n)
  for(i in 1:n){
  	ssq=0
  	sumx=0
    for(j in i:n){
        len=j-i+1
        sumx=sumx+data[j]
        ssq=ssq+data[j]^2
        all.seg[i,j]=-0.5*(ssq-(sumx^2)/len)
    }
  }
  like.Q=matrix(0,ncol=n,nrow=Q)
  like.Q[1,]=all.seg[1,]
  cp=matrix(0,ncol=n,nrow=Q)
  for(q in 2:Q){
    for(j in q:n){
      like=NULL
      v=(q-1):(j-1)
      like=like.Q[q-1,v]+all.seg[v+1,j]

      like.Q[q,j]= max(like,na.rm=TRUE)
      cp[q,j]=which(like==max(like,na.rm=TRUE))[1]+(q-2)
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
    criterion=-2*like.Q[,n]+k*pen[i]

    op.cps=c(op.cps,which(criterion==min(criterion))-1)
  }
  return(list(cps=cps.Q,op.cpts=op.cps,pen=pen,like=criterion[op.cps+1]))
}


segneigh.meanvar.norm=function(data,Q=5,pen=0){
  n=length(data)
  all.seg=matrix(0,ncol=n,nrow=n)
  for(i in 1:n){
  	ssq=0
  	sumx=0
    for(j in i:n){
        len=j-i+1
        sumx=sumx+data[j]
        ssq=ssq+data[j]^2
        sigmasq=(1/len)*(ssq-(sumx^2)/len)
        if(sigmasq<=0){sigmasq=0.00000000001}
        all.seg[i,j]=-(len/2)*(log(2*pi)+log(sigmasq)+1)
    }
  }
  like.Q=matrix(0,ncol=n,nrow=Q)
  like.Q[1,]=all.seg[1,]
  cp=matrix(0,ncol=n,nrow=Q)
  for(q in 2:Q){
    for(j in q:n){
      like=NULL
      if((j-2-q)<0){v=q}
      else{v=(q):(j-2)}
      like=like.Q[q-1,v]+all.seg[v+1,j]

      like.Q[q,j]= max(like,na.rm=TRUE)
      cp[q,j]=which(like==max(like,na.rm=TRUE))[1]+(q-1)
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
    criterion=-2*like.Q[,n]+k*pen[i]

    op.cps=c(op.cps,which(criterion==min(criterion))-1)
  }
  return(list(cps=cps.Q,op.cpts=op.cps,pen=pen,like=criterion[op.cps+1]))
}


binseg.var.norm=function(data,Q=5,pen=0,know.mean=FALSE,mu=-1000){
  mll.var=function(x,n){
    neg=x<=0
    x[neg==TRUE]=0.00000000001    
    return( -0.5*n*(log(2*pi)+log(x/n)+1))
  }
  n=length(data)
  if((know.mean==FALSE)&(mu==-1000)){
	mu=mean(data)
  }
  y2=c(0,cumsum((data-mu)^2))
  tau=c(0,n)
  cpt=matrix(0,nrow=2,ncol=Q)
  oldmax=1000

  for(q in 1:Q){
    lambda=rep(0,n-1)
    i=1
    st=tau[1]+1;end=tau[2]
    null=mll.var(y2[end+1]-y2[st],end-st+1)
    for(j in 1:(n-1)){
      if(j==end){
        st=end+1;i=i+1;end=tau[i+1]
        null=mll.var(y2[end+1]-y2[st],end-st+1)
      }else{
        lambda[j]=mll.var(y2[j+1]-y2[st],j-st+1)+mll.var(y2[end+1]-y2[j+1],end-j)-null
      }
    }
    k=which.max(lambda)[1]
    cpt[1,q]=k;cpt[2,q]=min(oldmax,max(lambda))
    oldmax=min(oldmax,max(lambda))
    tau=sort(c(tau,k))
  }
  op.cps=NULL
  p=1:(Q-1)
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


binseg.mean.norm=function(data,Q=5,pen=0){
  mll.mean=function(x2,x,n){
    return( -0.5*(x2-(x^2)/n))
  }
  n=length(data)
  y2=c(0,cumsum(data^2))
  y=c(0,cumsum(data))
  tau=c(0,n)
  cpt=matrix(0,nrow=2,ncol=Q)
  oldmax=1000

  for(q in 1:Q){
    lambda=rep(0,n-1)
    i=1
    st=tau[1]+1;end=tau[2]
    null=mll.mean(y2[end+1]-y2[st],y[end+1]-y[st],end-st+1)
    for(j in 1:(n-1)){
      if(j==end){
        st=end+1;i=i+1;end=tau[i+1]
        null=mll.mean(y2[end+1]-y2[st],y[end+1]-y[st],end-st+1)
      }else{
        lambda[j]=mll.mean(y2[j+1]-y2[st],y[j+1]-y[st],j-st+1)+mll.mean(y2[end+1]-y2[j+1],y[end+1]-y[j+1],end-j)-null
      }
    }
    k=which.max(lambda)[1]
    cpt[1,q]=k;cpt[2,q]=min(oldmax,max(lambda)) # done so that when we do the decision later we can take the max(which(criterion==T)), rather than min(which(criterion==F))-1
    oldmax=min(oldmax,max(lambda))
    tau=sort(c(tau,k))
  }
  op.cps=NULL
  p=1:(Q-1)
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


binseg.meanvar.norm=function(data,Q=5,pen=0){
  mll.meanvar=function(x2,x,n){
    sigmasq=(1/n)*(x2-(x^2)/n)
    neg=sigmasq<=0
    sigmasq[neg==TRUE]=0.00000000001
    return(-(n/2)*(log(2*pi)+log(sigmasq)+1))
  }
  n=length(data)
  y2=c(0,cumsum(data^2))
  y=c(0,cumsum(data))
  tau=c(0,n)
  cpt=matrix(0,nrow=2,ncol=Q)
  oldmax=1000

  for(q in 1:Q){
    lambda=rep(0,n-1)
    i=1
    st=tau[1]+1;end=tau[2]
    null=mll.meanvar(y2[end+1]-y2[st],y[end+1]-y[st],end-st+1)
    for(j in 1:(n-1)){
      if(j==end){
        st=end+1;i=i+1;end=tau[i+1]
        null=mll.meanvar(y2[end+1]-y2[st],y[end+1]-y[st],end-st+1)
      }else{
	if((j-st)<2){lambda[j]=-1*10^(100)}
	else if((end-j)<2){lambda[j]=-1*10^(100)}
	else{lambda[j]=mll.meanvar(y2[j+1]-y2[st],y[j+1]-y[st],j-st+1)+mll.meanvar(y2[end+1]-y2[j+1],y[end+1]-y[j+1],end-j)-null}
      }
    }
    k=which.max(lambda)[1]
    cpt[1,q]=k;cpt[2,q]=min(oldmax,max(lambda))
    oldmax=min(oldmax,max(lambda))
    tau=sort(c(tau,k))
  }
  op.cps=NULL
  p=1:(Q-1)
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


multiple.var.norm=function(data,mul.method="PELT",penalty="SIC",value=0,Q=5,know.mean=FALSE,mu=-1000,class=TRUE,param.estimates=TRUE){
	if(!((mul.method=="PELT")||(mul.method=="BinSeg")||(mul.method=="SegNeigh"))){
		stop("Multiple Method is not recognised")
	}
	diffparam=1
	if(is.null(dim(data))==TRUE){
		# single dataset
		n=length(data)
		mu=mu[1]
	}
	else{
		n=ncol(data)
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
		alogn=sqrt(2*log(log(n)))
		blogn=2*log(log(n))+ (log(log(log(n))))/2 - log(gamma(1/2))
		value=(-(log(log((1-alpha+exp(-2*exp(blogn)))^(-1/2))))/alogn + blogn/alogn)^2
	}
	if(is.null(dim(data))==TRUE){
		# single dataset
		if(mul.method=="PELT"){
			out=PELT.var.norm(data,value,know.mean,mu,nprune=FALSE)
			cpts=out
		}
		else if(mul.method=="BinSeg"){
			out=binseg.var.norm(data,Q,value,know.mean,mu)
			cpts=c(sort(out$cps[1,1:out$op.cpts]),n)
		}
		else if(mul.method=="SegNeigh"){
			out=segneigh.var.norm(data,Q,value,know.mean,mu)
			if(out$op.cpts==0){cpts=n}
			else{cpts=c(sort(out$cps[out$op.cpts+1,][out$cps[out$op.cpts+1,]>0]),n)}
		}
		if(class==TRUE){
			ans=new("cpt")
			data.set(ans)=data;cpttype(ans)="variance";method(ans)=mul.method;distribution(ans)="Normal"; pen.type(ans)=penalty;pen.value(ans)=value;cpts(ans)=cpts
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
		rep=nrow(data)
		out=list()
		if(length(mu)!=rep){
			mu=rep(mu,rep)
		}
		if(class==TRUE){cpts=list()}
		if(mul.method=="PELT"){
			for(i in 1:rep){
				out=c(out,list(PELT.var.norm(data[i,],value,know.mean,mu[i],nprune=FALSE)))
			}
			if(class==TRUE){cpts=out}
		}
		else if(mul.method=="BinSeg"){
			for(i in 1:rep){
				out=c(out,list(binseg.var.norm(data[i,],Q,value,know.mean,mu[i])))
				if(class==TRUE){cpts[[i]]=c(sort(out$cps[1,1:out$op.cpts]),n)}
			}
		}
		else if(mul.method=="SegNeigh"){
			for(i in 1:rep){
				out=c(out,list(segneigh.var.norm(data[i,],Q,value,know.mean,mu[i])))
				if(class==TRUE){
					if(out$cps==0){cpts[[i]]=n}
					else{cpts[[i]]=c(sort(out$cps[out$op.cpts+1,][out$cps[out$op.cpts+1,]>0]),n)}
				}
			}
		}
		if(class==TRUE){
			ans=list()
			for(i in 1:rep){
				ans[[i]]=new("cpt")
				data.set(ans[[i]])=data[i,];cpttype(ans[[i]])="variance";method(ans[[i]])=mul.method;distribution(ans[[i]])="Normal"; pen.type(ans[[i]])=penalty;pen.value(ans[[i]])=value;cpts(ans[[i]])=cpts[[i]]
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


multiple.mean.norm=function(data,mul.method="PELT",penalty="SIC",value=0,Q=5,class=TRUE,param.estimates=TRUE){
	if(!((mul.method=="PELT")||(mul.method=="BinSeg")||(mul.method=="SegNeigh"))){
		stop("Multiple Method is not recognised")
	}
	diffparam=1
	if(is.null(dim(data))==TRUE){
		# single dataset
		n=length(data)
	}
	else{
		n=ncol(data)
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
			alogn=(2*log(log(n)))^(-(1/2))
			blogn=(alogn^(-1))+(1/2)*alogn*log(log(log(n)))
			value=(-alogn*log(log((1-alpha+exp(-2*(pi^(1/2))*exp(blogn/alogn)))^(-1/(2*(pi^(1/2))))))+blogn)^2
	}
	if(is.null(dim(data))==TRUE){
		# single dataset
		if(mul.method=="PELT"){
			out=PELT.mean.norm(data,value,nprune=FALSE)
			cpts=out
		}
		else if(mul.method=="BinSeg"){
			out=binseg.mean.norm(data,Q,value)
			cpts=c(sort(out$cps[1,1:out$op.cpts]),n)
		}
		else if(mul.method=="SegNeigh"){
			out=segneigh.mean.norm(data,Q,value)
			if(out$op.cpts==0){cpts=n}
			else{cpts=c(sort(out$cps[out$op.cpts+1,][out$cps[out$op.cpts+1,]>0]),n)}
		}
		if(class==TRUE){
			ans=new("cpt")
			data.set(ans)=data;cpttype(ans)="mean";method(ans)=mul.method;distribution(ans)="Normal"; pen.type(ans)=penalty;pen.value(ans)=value;cpts(ans)=cpts
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
		rep=nrow(data)
		out=list()
		if(class==TRUE){cpts=list()}
		if(mul.method=="PELT"){
			for(i in 1:rep){
				out=c(out,list(PELT.mean.norm(data[i,],value,nprune=FALSE)))
			}
			cpts=out
		}
		else if(mul.method=="BinSeg"){
			for(i in 1:rep){
				out=c(out,list(binseg.mean.norm(data[i,],Q,value)))
				if(class==TRUE){cpts[[i]]=c(sort(out$cps[1,1:out$op.cpts]),n)}
			}
		}
		else if(mul.method=="SegNeigh"){
			for(i in 1:rep){
				out=c(out,list(segneigh.mean.norm(data[i,],Q,value)))
				if(class==TRUE){
					if(out$op.cpts==0){cpts[[i]]=n}
					else{cpts[[i]]=c(sort(out$cps[out$op.cpts+1,][out$cps[out$op.cpts+1,]>0]),n)}
				}
			}
		}
		if(class==TRUE){
			ans=list()
			for(i in 1:rep){
				ans[[i]]=new("cpt")
				data.set(ans[[i]])=data[i,];cpttype(ans[[i]])="mean";method(ans[[i]])=mul.method;distribution(ans[[i]])="Normal"; pen.type(ans[[i]])=penalty;pen.value(ans[[i]])=value;cpts(ans[[i]])=cpts[[i]]
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

multiple.meanvar.norm=function(data,mul.method="PELT",penalty="SIC",value=0,Q=5,class=TRUE,param.estimates=TRUE){
	if(!((mul.method=="PELT")||(mul.method=="BinSeg")||(mul.method=="SegNeigh"))){
		stop("Multiple Method is not recognised")
	}
	diffparam=2
	if(is.null(dim(data))==TRUE){
		# single dataset
		n=length(data)
	}
	else{
		n=ncol(data)
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
			alogn=sqrt(2*log(log(n)))
			blogn=2*log(log(n))+ log(log(log(n)))
			value=(-(log(log((1-alpha+exp(-2*exp(blogn)))^(-1/2))))/alogn + blogn/alogn)^2
	}
	if(is.null(dim(data))==TRUE){
		# single dataset
		if(mul.method=="PELT"){
			out=PELT.meanvar.norm(data,value,nprune=FALSE)
			cpts=out
		}
		else if(mul.method=="BinSeg"){
			out=binseg.meanvar.norm(data,Q,value)
			cpts=c(sort(out$cps[1,1:out$op.cpts]),n)
		}
		else if(mul.method=="SegNeigh"){
			out=segneigh.meanvar.norm(data,Q,value)
			if(out$op.cpts==0){cpts=n}
			else{cpts=c(sort(out$cps[out$op.cpts+1,][out$cps[out$op.cpts+1,]>0]),n)}
		}
		if(class==TRUE){
			ans=new("cpt")
			data.set(ans)=data;cpttype(ans)="mean and variance";method(ans)=mul.method; distribution(ans)="Normal";pen.type(ans)=penalty;pen.value(ans)=value;cpts(ans)=cpts
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
		rep=nrow(data)
		out=list()
		if(class==TRUE){cpts=list()}
		if(mul.method=="PELT"){
			for(i in 1:rep){
				out=c(out,list(PELT.meanvar.norm(data[i,],value,nprune=FALSE)))
			}
			cpts=out
		}
		else if(mul.method=="BinSeg"){
			for(i in 1:rep){
				out=c(out,list(binseg.meanvar.norm(data[i,],Q,value)))
				if(class==TRUE){cpts[[i]]=c(sort(out$cps[1,1:out$op.cpts]),n)}
			}
		}
		else if(mul.method=="SegNeigh"){
			for(i in 1:rep){
				out=c(out,list(segneigh.meanvar.norm(data[i,],Q,value)))
				if(class==TRUE){
					if(out$op.cpts==0){cpts[[i]]=n}
					else{cpts[[i]]=c(sort(out$cps[out$op.cpts+1,][out$cps[out$op.cpts+1,]>0]),n)}
				}
			}
		}
		if(class==TRUE){
			ans=list()
			for(i in 1:rep){
				ans[[i]]=new("cpt")
				data.set(ans[[i]])=data[i,];cpttype(ans[[i]])="mean and variance"; method(ans[[i]])=mul.method;distribution(ans[[i]])="Normal";pen.type(ans[[i]])=penalty;pen.value(ans[[i]])=value;cpts(ans[[i]])=cpts[[i]]
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

