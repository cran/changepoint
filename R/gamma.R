single.meanvar.gamma.calc <-
function(data,shape=1,extrainf=TRUE){
  singledim=function(data,shape,extrainf=TRUE){
    n=length(data)
    y=c(0,cumsum(data))
    null=2*n*shape*log(y[n+1])-2*n*shape*log(n*shape)
    taustar=2:(n-1)
    tmp=2*taustar*shape*log(y[taustar+1]) -2*taustar*shape*log(taustar*shape) + 2*(n-taustar)*shape*log((y[n+1]-y[taustar+1]))-2*(n-taustar)*shape*log((n-taustar)*shape)
    
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
    cpt=singledim(data,shape,extrainf)
    return(cpt)
  }
  else{
    rep=nrow(data)
    n=ncol(data)
    cpt=NULL
    if(length(shape)==1){
	shape=rep(shape,rep)
    }
    if(extrainf==FALSE){
      for(i in 1:rep){
        cpt[i]=singledim(data[i,],shape[i],extrainf)
      }
    }
    else{
      cpt=matrix(0,ncol=3,nrow=rep)
      for(i in 1:rep){
        cpt[i,]=singledim(data[i,],shape[i],extrainf)
      }
      colnames(cpt)=c('cpt','null','alt')
    }
    return(cpt)
  }
}

single.meanvar.gamma<-function(data,shape=1,penalty="SIC",value=0,class=TRUE,param.estimates=TRUE){
	if(sum(data<=0)>0){stop('Gamma distribution requires positive data')}

	if(is.null(dim(data))==TRUE){
		n=length(data)
		if(penalty=="Asymptotic"){
			return('Asymptotic values for the Gamma Distribution are not defined, please choose an alternative penalty type')
		}
		tmp=single.meanvar.gamma.calc(data,shape,extrainf=TRUE)
		ans=decision(tmp[1],tmp[2],tmp[3],penalty,n,diffparam=1,value)
		if(class==TRUE){
			out=new("cpt")
			data.set(out)=data;cpttype(out)="mean and variance";method(out)="AMOC";distribution(out)="Gamma";pen.type(out)=penalty; pen.value(out)=ans$pen;ncpts.max(out)=1
			if(ans$cpt != n){cpts(out)=c(ans$cpt,n)}
			else{cpts(out)=ans$cpt}
			if(param.estimates==TRUE){
				out=param(out,shape=shape)
			}
			return(out)
		}
		else{ return(ans$cpt)}
	}
	else{ 
		n=ncol(data)
		if(penalty=="Asymptotic"){
			return('Asymptotic values for the Gamma Distribution are not defined, please choose an alternative penalty type')
		}
		tmp=single.meanvar.gamma.calc(data,shape,extrainf=TRUE)
		ans=decision(tmp[,1],tmp[,2],tmp[,3],penalty,n,diffparam=1,value)
		if(class==TRUE){
			rep=nrow(data)
			out=list()
			for(i in 1:rep){
				out[[i]]=new("cpt")
				data.set(out[[i]])=data[i,];cpttype(out[[i]])="mean and variance";method(out[[i]])="AMOC";distribution(out[[i]])="Gamma"; pen.type(out[[i]])=penalty;pen.value(out[[i]])=ans$pen;ncpts.max(out[[i]])=1
				if(ans$cpt[i] != n){cpts(out[[i]])=c(ans$cpt[i],n)}
				else{cpts(out[[i]])=ans$cpt[i]}
				if(param.estimates==TRUE){
					out[[i]]=param(out[[i]],shape=shape)
				}
			}
			return(out)
		}
		else{ return(ans$cpt)}
	}
}

# PELT.meanvar.gamma=function(data,shape=1,pen=0,nprune=FALSE){
#   mll.meanvar.EFK=function(x,n,shape){
#     return( 2*n*shape*log(x)-2*n*shape*log(n*shape))
#   }
#   n=length(data)
#   y=c(0,cumsum(data))
# 
#   lastchangecpts=matrix(NA,nrow=n,ncol=2)
#   lastchangelike=matrix(NA,nrow=n,ncol=2)
#   checklist=NULL
#   lastchangelike[1,]=c(mll.meanvar.EFK(y[2],1,shape),mll.meanvar.EFK(y[n+1]-y[2],n-1,shape)+pen)
#   lastchangecpts[1,]=c(0,1)
#   lastchangelike[2,]=c(mll.meanvar.EFK(y[3],2,shape),mll.meanvar.EFK(y[n+1]-y[3],n-2,shape)+pen)
#   lastchangecpts[2,]=c(0,2)
#   lastchangelike[3,]=c(mll.meanvar.EFK(y[4],3,shape),mll.meanvar.EFK(y[n+1]-y[4],n-3,shape)+pen)
#   lastchangecpts[3,]=c(0,3)
#   noprune=NULL
#   for(tstar in 4:n){
#     tmplike=NULL
#     tmpt=c(checklist, tstar-2)
#     tmplike=lastchangelike[tmpt,1]+mll.meanvar.EFK(y[tstar+1]-y[tmpt+1],tstar-tmpt,shape)+pen
#     if(tstar==n){
#       lastchangelike[tstar,]=c(min(c(tmplike,mll.meanvar.EFK(y[tstar+1]-y[1],tstar,shape)),na.rm=TRUE),0)
#     }
#     else{
#       lastchangelike[tstar,]=c(min(c(tmplike,mll.meanvar.EFK(y[tstar+1]-y[1],tstar,shape)),na.rm=TRUE),mll.meanvar.EFK(y[n+1]-y[tstar+1],n-tstar,shape)+pen)
#     }
#     if(lastchangelike[tstar,1]==mll.meanvar.EFK(y[tstar+1]-y[1],tstar,shape)){
#       lastchangecpts[tstar,]=c(0,tstar)
#     }
#     else{
#       cpt=tmpt[tmplike==lastchangelike[tstar,1]][1]
#       lastchangecpts[tstar,]=c(cpt,tstar)
#     }
#     checklist=tmpt[tmplike<=lastchangelike[tstar,1]+pen]
#     if(nprune==TRUE){
#       noprune=c(noprune,length(checklist))
#     }
#   }
#   if(nprune==TRUE){
#     return(nprune=noprune)
#   }
#   else{
#     fcpt=NULL
#     last=n
#     while(last!=0){
# 	fcpt=c(fcpt,lastchangecpts[last,2])
# 	last=lastchangecpts[last,1]
#     }
#     return(cpt=sort(fcpt))
#   }
# }


PELT.meanvar.gamma=function(data,shape=1,pen=0){
  # function that uses the PELT method to calculate changes in mean & variance where the segments in the data are assumed to be Exponential
	if(sum(data<=0)>0){stop('Gamma distribution requires positive data')}

  n=length(data)
  y=c(0,cumsum(data))
  error=0
  
  storage.mode(y)='double'
  cptsout=rep(0,n) # sets up null vector for changepoint answer
  storage.mode(cptsout)='integer'
  
	answer=list()
	answer[[6]]=1
	on.exit(.C("FreePELT",answer[[6]],PACKAGE='changepoint'))
	answer=.C('PELT_meanvar_gamma',y,as.integer(n),as.double(pen),cptsout,as.double(shape),as.integer(error),PACKAGE='changepoint')
	if(answer[[6]]>0){
	  print("C code error:",answer[[6]])
	  stop(call.=F)
	}
	
  return(sort(answer[[4]][answer[[4]]>0]))
}


segneigh.meanvar.gamma=function(data,shape=1,Q=5,pen=0){
	if(sum(data<=0)>0){stop('Gamma distribution requires positive data')}

  n=length(data)
  if(Q>(n/2)){stop(paste('Q is larger than the maximum number of segments',n/2))}
  all.seg=matrix(0,ncol=n,nrow=n)
  for(i in 1:n){
  	sumx=0
    for(j in i:n){
        len=j-i+1
        sumx=sumx+data[j]
        all.seg[i,j]=len*shape*log(len*shape)-len*shape*log(sumx)
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


# binseg.meanvar.gamma=function(data,shape=1,Q=5,pen=0){
#   mll.meanvar=function(x,n,shape){
#     return(n*shape*log(n*shape)-n*shape*log(x))
#   }
#   n=length(data)
#   y=c(0,cumsum(data))
#   tau=c(0,n)
#   cpt=matrix(0,nrow=2,ncol=Q)
#   oldmax=1000
# 
#   for(q in 1:Q){
#     lambda=rep(0,n-1)
#     i=1
#     st=tau[1]+1;end=tau[2]
#     null=mll.meanvar(y[end+1]-y[st],end-st+1,shape)
#     for(j in 1:(n-1)){
#       if(j==end){
#         st=end+1;i=i+1;end=tau[i+1]
#         null=mll.meanvar(y[end+1]-y[st],end-st+1,shape)
#       }else{
# 	if((j-st)<2){lambda[j]=-1*10^(100)}
# 	else if((end-j)<2){lambda[j]=-1*10^(100)}
# 	else{lambda[j]=mll.meanvar(y[j+1]-y[st],j-st+1,shape)+mll.meanvar(y[end+1]-y[j+1],end-j,shape)-null}
#       }
#     }
#     k=which.max(lambda)[1]
#     cpt[1,q]=k;cpt[2,q]=min(oldmax,max(lambda))
#     oldmax=min(oldmax,max(lambda))
#     tau=sort(c(tau,k))
#   }
#   op.cps=NULL
#   p=1:(Q-1)
#   for(i in 1:length(pen)){
#     criterion=(2*cpt[2,])>=pen[i]
#     if(sum(criterion)==0){
#       op.cps=0
#     }
#     else{
#       op.cps=c(op.cps,max(which((criterion)==TRUE)))
#     }
#   }
#   return(list(cps=cpt,op.cpts=op.cps,pen=pen))
# }

binseg.meanvar.gamma=function(data,shape=1,Q=5,pen=0){
  # function that uses the BinSeg method to calculate changes in mean & variance where the segments in the data are assumed to be Exponential
	if(sum(data<=0)>0){stop('Gamma distribution requires positive data')}

  n=length(data)
  if(Q>(n/2)){stop(paste('Q is larger than the maximum number of segments',n/2))}
  y=c(0,cumsum(data))

  storage.mode(y)='double'
  cptsout=rep(0,Q) # sets up null vector for changepoint answer
  likeout=rep(0,Q) # sets up null vector for likelihood of changepoints in cptsout
  storage.mode(cptsout)='double'
  op_cps=0
  
  answer=.C('binseg_meanvar_gamma',y,as.integer(n),as.double(pen),as.integer(Q),as.integer(cptsout),likeout,as.integer(op_cps),as.double(shape),PACKAGE='changepoint')
  return(list(cps=rbind(answer[[5]],answer[[6]]),op.cpts=answer[[7]],pen=pen))
}


multiple.meanvar.gamma=function(data,shape=1,mul.method="PELT",penalty="SIC",value=0,Q=5,class=TRUE,param.estimates=TRUE){
	if(!((mul.method=="PELT")||(mul.method=="BinSeg")||(mul.method=="SegNeigh"))){
		stop("Multiple Method is not recognised")
	}
	diffparam=1
	if(is.null(dim(data))==TRUE){
		# single dataset
		n=length(data)
		shape=shape[1]
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
  else if(penalty=='Asymptotic'){
    return('Asymptotic values for the Gamma Distribution are not defined, please choose an alternative penalty type')
  }
	if(is.null(dim(data))==TRUE){
		# single dataset
		if(mul.method=="PELT"){
			out=PELT.meanvar.gamma(data,shape,value)
			cpts=out
		}
		else if(mul.method=="BinSeg"){
			out=binseg.meanvar.gamma(data,shape,Q,value)
			if(out$op.cpts==0){cpts=n}
			else{cpts=c(sort(out$cps[1,1:out$op.cpts]),n)}
		}
		else if(mul.method=="SegNeigh"){
			out=segneigh.meanvar.gamma(data,shape,Q,value)
			if(out$op.cpts==0){cpts=n}
			else{cpts=c(sort(out$cps[out$op.cpts+1,][out$cps[out$op.cpts+1,]>0]),n)}
		}
		if(class==TRUE){
			ans=new("cpt")
			data.set(ans)=data;cpttype(ans)="mean and variance";method(ans)=mul.method; distribution(ans)="Gamma";pen.type(ans)=penalty;pen.value(ans)=value;cpts(ans)=cpts
			if(mul.method=="PELT"){
				ncpts.max(ans)=Inf
			}
			else{
				ncpts.max(ans)=Q
			}
			if(param.estimates==TRUE){
				ans=param(ans,shape)
			}
			return(ans)
		}
		else{ return(out)}
	}
	else{
		rep=nrow(data)
		out=list()
		if(length(shape)!=rep){
			shape=rep(shape,rep)
		}
		if(class==TRUE){cpts=list()}
		if(mul.method=="PELT"){
			for(i in 1:rep){
				out=c(out,list(PELT.meanvar.gamma(data[i,],shape[i],value)))
			}
			cpts=out
		}
		else if(mul.method=="BinSeg"){
			for(i in 1:rep){
				out=c(out,list(binseg.meanvar.gamma(data[i,],shape[i],Q,value)))
				if(class==TRUE){
					if(out[[i]]$op.cpts==0){cpts[[i]]=n}
					else{cpts[[i]]=c(sort(out[[i]]$cps[1,1:out[[i]]$op.cpts]),n)}
				}
			}
		}
		else if(mul.method=="SegNeigh"){
			for(i in 1:rep){
				out=c(out,list(segneigh.meanvar.gamma(data[i,],shape[i],Q,value)))
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
				data.set(ans[[i]])=data[i,];cpttype(ans[[i]])="mean and variance"; method(ans[[i]])=mul.method;distribution(ans[[i]])="Gamma";pen.type(ans[[i]])=penalty;pen.value(ans[[i]])=value;cpts(ans[[i]])=cpts[[i]]
				if(mul.method=="PELT"){
					ncpts.max(ans[[i]])=Inf
				}
				else{
					ncpts.max(ans[[i]])=Q
				}
				if(param.estimates==TRUE){
					ans[[i]]=param(ans[[i]],shape[i])
				}
			}
			return(ans)
		}
		else{return(out)}
	}
}

