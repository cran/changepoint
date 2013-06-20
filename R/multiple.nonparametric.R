segneigh.var.css=function(data,Q=5,pen=0){
  n=length(data)
  if(n<4){stop('Data must have atleast 4 observations to fit a changepoint model.')}
  if(Q>((n/2)+1)){stop(paste('Q is larger than the maximum number of segments',(n/2)+1))}
  
	y2=c(0,cumsum(data^2))
	oldmax=1000

	test=NULL
  like.Q=matrix(0,ncol=n,nrow=Q)
  cp=matrix(0,ncol=n,nrow=Q)
  for(q in 2:Q){ # no of segments
    for(j in q:n){
      like=NULL
      v=(q-1):(j-1)
			if(q==2){
	      like=abs(sqrt(j/2)*(y2[v+1]/y2[j+1] -v/j))
			}
			else{
	      like=like.Q[q-1,v]+abs(sqrt((j-cp[q-1,v])/2)*((y2[v+1]-y2[cp[q-1,v]+1])/(y2[j+1]-y2[cp[q-1,v]+1]) -(v-cp[q-1,v])/(j-cp[q-1,v])))
			}
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

	op.cps=0
	flag=0
	for(q in 2:Q){
		criterion=NULL
		cpttmp=c(0,sort(cps.Q[q,1:(q-1)]),n)
		for(i in 1:(q-1)){
			criterion[i]=abs(sqrt((cpttmp[i+2]-cpttmp[i])/2)*((y2[cpttmp[i+1]+1]-y2[cpttmp[i]+1])/(y2[cpttmp[i+2]+1]-y2[cpttmp[i]+1]) -(cpttmp[i+1]-cpttmp[i])/(cpttmp[i+2]-cpttmp[i])))
			if(criterion[i]<pen){flag=1}
		}
		if(flag==1){
			break
		}
		op.cps=op.cps+1
	}
  if(op.cps==(Q-1)){warning('The number of segments identified is Q, it is advised to increase Q to make sure changepoints have not been missed.')}
  return(list(cps=cps.Q,op.cpts=op.cps,pen=pen))
}



binseg.var.css=function(data,Q=5,pen=0){
  n=length(data)
  if(n<4){stop('Data must have atleast 4 observations to fit a changepoint model.')}
  if(Q>((n/2)+1)){stop(paste('Q is larger than the maximum number of segments',(n/2)+1))}
  
  y2=c(0,cumsum(data^2))
  tau=c(0,n)
  cpt=matrix(0,nrow=2,ncol=Q)
  oldmax=Inf

  for(q in 1:Q){
    lambda=rep(0,n-1)
    i=1
    st=tau[1]+1;end=tau[2]
    for(j in 1:(n-1)){
      if(j==end){
        st=end+1;i=i+1;end=tau[i+1]
      }else{
        lambda[j]=sqrt((end-st+1)/2)*((y2[j+1]-y2[st])/(y2[end+1]-y2[st]) -(j-st+1)/(end-st+1))
      }
    }
    k=which.max(abs(lambda))
    cpt[1,q]=k;cpt[2,q]=min(oldmax,max(abs(lambda),na.rm=T))
    oldmax=min(oldmax,max(abs(lambda),na.rm=T))
    tau=sort(c(tau,k))
  }
  op.cps=NULL
  p=1:(Q-1)
  for(i in 1:length(pen)){
    criterion=(cpt[2,])>=pen[i]
    if(sum(criterion)==0){
      op.cps=0
    }
    else{
      op.cps=c(op.cps,max(which((criterion)==TRUE)))
    }
  }
  if(op.cps==Q){warning('The number of changepoints identified is Q, it is advised to increase Q to make sure changepoints have not been missed.')}
  return(list(cps=cpt,op.cpts=op.cps,pen=pen))
}


multiple.var.css=function(data,mul.method="BinSeg",penalty="SIC",pen.value=0,Q=5,class=TRUE,param.estimates=TRUE){
	if(mul.method=="PELT"){ stop("CSS does not satisfy the assumptions of PELT, use SegNeigh or BinSeg instead.") }
	else if(!((mul.method=="BinSeg")||(mul.method=="SegNeigh"))){
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
		pen.value=log(diffparam*log(n))
	}
	else if((penalty=="SIC1") || (penalty=="BIC1")){
		pen.value=log((diffparam+1)*log(n))
	}
	else if(penalty=="AIC"){
		pen.value=log(2*diffparam)
	}
	else if(penalty=="AIC1"){
		pen.value=log(2*(diffparam+1))
	}
	else if(penalty=="Hannan-Quinn"){
		pen.value=log(2*diffparam*log(log(n)))
	}
	else if(penalty=="Hannan-Quinn1"){
		pen.value=log(2*(diffparam+1)*log(log(n)))
	}
	else if(penalty=="None"){
		pen.value=0
	}
	else if((penalty!="Manual")&&(penalty!="Asymptotic")){
		stop('Unknown Penalty')
	}
	if((penalty=="Manual")&&(is.numeric(pen.value)==FALSE)){
		pen.value=try(eval(parse(text=paste(pen.value))),silent=TRUE)
		if(class(pen.value)=='try-error'){
			stop('Your manual penalty cannot be evaluated')
		}
	}
	if(penalty=="Asymptotic"){
		if(pen.value==0.01){pen.value=1.628}
		else if(pen.value==0.05){pen.value=1.358}
		else if(pen.value==0.1){pen.value=1.224}
		else if(pen.value==0.25){pen.value=1.019}
		else if(pen.value==0.5){pen.value=0.828}
		else if(pen.value==0.75){pen.value=0.677}
		else if(pen.value==0.9){pen.value=0.571}
		else if(pen.value==0.95){pen.value=0.520}
		else{stop('Only alpha values of 0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95 are valid for CSS')}
	}
	if(is.null(dim(data))==TRUE){
		# single dataset
		if(mul.method=="BinSeg"){
			out=binseg.var.css(data,Q,pen.value)
			if(out$op.cpts==0){cpts=n}
			else{cpts=c(sort(out$cps[1,1:out$op.cpts]),n)}
		}
		else if(mul.method=="SegNeigh"){
			out=segneigh.var.css(data,Q,pen.value)
			if(out$op.cpts==0){cpts=n}
			else{cpts=c(sort(out$cps[out$op.cpts+1,][out$cps[out$op.cpts+1,]>0]),n)}
		}
		if(class==TRUE){
			ans=new("cpt")
			data.set(ans)=data;cpttype(ans)="variance";method(ans)=mul.method;test.stat(ans)="CSS"; pen.type(ans)=penalty;pen.value(ans)=pen.value;cpts(ans)=cpts
			ncpts.max(ans)=Q
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
		if(mul.method=="BinSeg"){
			for(i in 1:rep){
				out=c(out,list(binseg.var.css(data[i,],Q,pen.value)))
				if(class==TRUE){
					if(out[[i]]$op.cpts==0){cpts[[i]]=n}
					else{cpts[[i]]=c(sort(out[[i]]$cps[1,1:out[[i]]$op.cpts]),n)}
				}
			}
		}
		else if(mul.method=="SegNeigh"){
			for(i in 1:rep){
				out=c(out,list(segneigh.var.css(data[i,],Q,pen.value)))
				if(class==TRUE){
					if(out[[i]]$op.cpts==0){cpts[[i]]=n}
					else{cpts[[i]]=c(sort(out[[i]]$cps[out[[i]]$op.cpts+1,][out[[i]]$cps[out[[i]]$op.cpts+1,]>0]),n)}
				}
			}
		}
		if(class==TRUE){
			ans=list()
			for(i in 1:rep){
				ans[[i]]=new("cpt")
				data.set(ans[[i]])=data[i,];cpttype(ans[[i]])="variance";method(ans[[i]])=mul.method;test.stat(ans[[i]])="CSS"; pen.type(ans[[i]])=penalty;pen.value(ans[[i]])=pen.value;cpts(ans[[i]])=cpts[[i]]
				ncpts.max(ans[[i]])=Q
				if(param.estimates==TRUE){
					ans[[i]]=param(ans[[i]])
				}
			}
			return(ans)
		}
		else{return(out)}
	}
}

















segneigh.mean.cusum=function(data,Q=5,pen=0){
  n=length(data)
  if(n<2){stop('Data must have atleast 2 observations to fit a changepoint model.')}
  if(Q>((n/2)+1)){stop(paste('Q is larger than the maximum number of segments',(n/2)+1))}
  
  y=c(0,cumsum(data))
	oldmax=1000

	test=NULL
  like.Q=matrix(0,ncol=n,nrow=Q)
  cp=matrix(0,ncol=n,nrow=Q)
  for(q in 2:Q){ # no of segments
    for(j in q:n){
      like=NULL
      v=(q-1):(j-1)
			if(q==2){
	      like=abs((y[v+1]-(v/j)*y[j+1])/j)
			}
			else{
	      like=like.Q[q-1,v]+abs(((y[v+1]-y[cp[q-1,v]+1])-((v-cp[q-1,v])/(j-cp[q-1,v]))*(y[j+1]-y[cp[q-1,v]+1]))/(j-cp[q-1,v]))
			}
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

	op.cps=0
	flag=0
	for(q in 2:Q){
		criterion=NULL
		cpttmp=c(0,sort(cps.Q[q,1:(q-1)]),n)
		for(i in 1:(q-1)){
			criterion[i]=abs(((y[cpttmp[i+1]+1]-y[cpttmp[i]+1])-((cpttmp[i+1]-cpttmp[i])/(cpttmp[i+2]-cpttmp[i]))*(y[cpttmp[i+2]+1]-y[cpttmp[i]+1]))/(cpttmp[i+2]-cpttmp[i]))
			if(criterion[i]<pen){flag=1}
		}
		if(flag==1){
			break
		}
		op.cps=op.cps+1
	}

  if(op.cps==(Q-1)){warning('The number of segments identified is Q, it is advised to increase Q to make sure changepoints have not been missed.')}
  return(list(cps=cps.Q,op.cpts=op.cps,pen=pen))
}


binseg.mean.cusum=function(data,Q=5,pen=0){
  n=length(data)
  if(n<2){stop('Data must have atleast 2 observations to fit a changepoint model.')}
  if(Q>((n/2)+1)){stop(paste('Q is larger than the maximum number of segments',(n/2)+1))}

  y=c(0,cumsum(data))
  tau=c(0,n)
  cpt=matrix(0,nrow=2,ncol=Q)
  oldmax=Inf

  for(q in 1:Q){
    lambda=rep(0,n-1)
    i=1
    st=tau[1]+1;end=tau[2]
    for(j in 1:(n-1)){
      if(j==end){
        st=end+1;i=i+1;end=tau[i+1]
      }else{
        lambda[j]=((y[j+1]-y[st])-((j-st+1)/(end-st+1))*(y[end+1]-y[st]))/(end-st+1)
      }
    }
    k=which.max(abs(lambda))
    cpt[1,q]=k;cpt[2,q]=min(oldmax,max(abs(lambda)))
    oldmax=min(oldmax,max(abs(lambda)))
    tau=sort(c(tau,k))
  }
  op.cps=NULL
  p=1:(Q-1)
  for(i in 1:length(pen)){
    criterion=(cpt[2,])>=pen[i]
    if(sum(criterion)==0){
      op.cps=0
    }
    else{
      op.cps=c(op.cps,max(which((criterion)==TRUE)))
    }
  }
  if(op.cps==Q){warning('The number of changepoints identified is Q, it is advised to increase Q to make sure changepoints have not been missed.')}
  return(list(cps=cpt,op.cpts=op.cps,pen=pen))
}


multiple.mean.cusum=function(data,mul.method="BinSeg",penalty="Asymptotic",pen.value=0.05,Q=5,class=TRUE,param.estimates=TRUE){
	if(mul.method=="PELT"){ stop("CUSUM does not satisfy the assumptions of PELT, use SegNeigh or BinSeg instead.") }
	else if(!((mul.method=="BinSeg")||(mul.method=="SegNeigh"))){
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
		pen.value=diffparam*log(n)
	}
	else if((penalty=="SIC1") || (penalty=="BIC1")){
		pen.value=(diffparam+1)*log(n)
	}
	else if(penalty=="AIC"){
		pen.value=2*diffparam
	}
	else if(penalty=="AIC1"){
		pen.value=2*(diffparam+1)
	}
	else if(penalty=="Hannan-Quinn"){
		pen.value=2*diffparam*log(log(n))
	}
	else if(penalty=="Hannan-Quinn1"){
		pen.value=2*(diffparam+1)*log(log(n))
	}
	else if(penalty=="None"){
		pen.value=0
	}
	else if((penalty!="Manual")&&(penalty!="Asymptotic")){
		stop('Unknown Penalty')
	}
	if((penalty=="Manual")&&(is.numeric(pen.value)==FALSE)){
		pen.value=try(eval(parse(text=paste(pen.value))),silent=TRUE)
		if(class(pen.value)=='try-error'){
			stop('Your manual penalty cannot be evaluated')
		}
	}
	if(penalty=="Asymptotic"){
		stop('Asymptotic penalties have not been implemented yet for CUSUM')
	}
	if(is.null(dim(data))==TRUE){
		# single dataset
		if(mul.method=="BinSeg"){
			out=binseg.mean.cusum(coredata(data),Q,pen.value)
			if(out$op.cpts==0){cpts=n}
			else{cpts=c(sort(out$cps[1,1:out$op.cpts]),n)}
		}
		else if(mul.method=="SegNeigh"){
			out=segneigh.mean.cusum(coredata(data),Q,pen.value)
			if(out$op.cpts==0){cpts=n}
			else{cpts=c(sort(out$cps[out$op.cpts+1,][out$cps[out$op.cpts+1,]>0]),n)}
		}
		if(class==TRUE){
			ans=new("cpt")
			data.set(ans)=data;cpttype(ans)="mean";method(ans)=mul.method;test.stat(ans)="CUSUM"; pen.type(ans)=penalty;pen.value(ans)=pen.value;cpts(ans)=cpts
			ncpts.max(ans)=Q
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
		if(mul.method=="BinSeg"){
			for(i in 1:rep){
				out=c(out,list(binseg.mean.cusum(data[i,],Q,pen.value)))
				if(class==TRUE){
					if(out[[i]]$op.cpts==0){cpts[[i]]=n}
					else{cpts[[i]]=c(sort(out[[i]]$cps[1,1:out[[i]]$op.cpts]),n)}
				}
			}
		}
		else if(mul.method=="SegNeigh"){
			for(i in 1:rep){
				out=c(out,list(segneigh.mean.cusum(data[i,],Q,pen.value)))
				if(class==TRUE){
					if(out[[i]]$op.cpts==0){cpts[[i]]=n}
					else{cpts[[i]]=c(sort(out[[i]]$cps[out[[i]]$op.cpts+1,][out[[i]]$cps[out[[i]]$op.cpts+1,]>0]),n)}
				}
			}
		}
		if(class==TRUE){
			ans=list()
			for(i in 1:rep){
				ans[[i]]=new("cpt")
				data.set(ans[[i]])=ts(data[i,]);cpttype(ans[[i]])="mean";method(ans[[i]])=mul.method;test.stat(ans[[i]])="CUSUM"; pen.type(ans[[i]])=penalty;pen.value(ans[[i]])=pen.value;cpts(ans[[i]])=cpts[[i]]
				ncpts.max(ans[[i]])=Q
				if(param.estimates==TRUE){
					ans[[i]]=param(ans[[i]])
				}
			}
			return(ans)
		}
		else{return(out)}
	}
}

