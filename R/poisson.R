single.meanvar.poisson.calc <-
  function(data,extrainf=TRUE){
    singledim=function(data,extrainf=TRUE){
      n=length(data)
      y=c(0,cumsum(data))
      if(y[n+1]==0){
        null=Inf
      }
      else{
        null=2*y[n+1]*log(n) - 2*y[n+1]*log(y[n+1])
      }
      taustar=2:(n-1)
      tmp=2*log(taustar)*y[taustar+1] -2*y[taustar+1]*log(y[taustar+1]) + 2*log(n-taustar)*(y[n+1]-y[taustar+1])-2*(y[n+1]-y[taustar+1])*log((y[n+1]-y[taustar+1]))
      if(sum(is.na(tmp))!=0){
        tmp[which(is.na(tmp))]=Inf
      }
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


single.meanvar.poisson<-function(data,penalty="SIC",value=0,class=TRUE,param.estimates=TRUE){
  if((sum(data<0)>0)){stop('Poisson distribution requires positive data')}
  if(sum(as.integer(data)==data)!=length(data)){stop('Poisson distribution requires integer data')}
  
  if(is.null(dim(data))==TRUE){
    n=length(data)
    if(penalty=="Asymptotic"){
      return('Asymptotic values for the Poisson Distribution are not available yet, please choose an alternative penalty type')
    }
    tmp=single.meanvar.poisson.calc(data,extrainf=TRUE)
    ans=decision(tmp[1],tmp[2],tmp[3],penalty,n,diffparam=1,value)
    if(class==TRUE){
      out=new("cpt")
      data.set(out)=data;cpttype(out)="mean and variance";method(out)="AMOC";distribution(out)="Poisson";pen.type(out)=penalty; pen.value(out)=ans$pen;ncpts.max(out)=1
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
      return('Asymptotic values for the Poisson Distribution are not available yet, please choose an alternative penalty type')
    }
    tmp=single.meanvar.poisson.calc(data,extrainf=TRUE)
    ans=decision(tmp[,1],tmp[,2],tmp[,3],penalty,n,diffparam=1,value)
    if(class==TRUE){
      rep=nrow(data)
      out=list()
      for(i in 1:rep){
        out[[i]]=new("cpt")
        data.set(out[[i]])=data[i,];cpttype(out[[i]])="mean and variance";method(out[[i]])="AMOC";distribution(out[[i]])="Poisson"; pen.type(out[[i]])=penalty;pen.value(out[[i]])=ans$pen;ncpts.max(out[[i]])=1
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




PELT.meanvar.poisson=function(data,pen=0){
  # function that uses the PELT method to calculate changes in mean & variance where the segments in the data are assumed to be Poisson
  if((sum(data<0)>0)){stop('Poisson distribution requires positive data')}
  if(sum(as.integer(data)==data)!=length(data)){stop('Poisson distribution requires integer data')}
  n=length(data)
  y=c(0,cumsum(data))
  error=0
  
  storage.mode(y)='double'
  cptsout=rep(0,n) # sets up null vector for changepoint answer
  storage.mode(cptsout)='integer'
  
  answer=list()
  answer[[5]]=1
  on.exit(.C("FreePELT",answer[[5]],PACKAGE='changepoint'))
  answer=.C('PELT_meanvar_poisson',y,as.integer(n),as.double(pen),cptsout,as.integer(error),PACKAGE='changepoint')
  if(answer[[5]]>0){
    print("C code error:",answer[[5]])
    stop(call.=F)
  }
  
  return(sort(answer[[4]][answer[[4]]>0]))
}


segneigh.meanvar.poisson=function(data,Q=5,pen=0){
  if((sum(data<0)>0)){stop('Poisson distribution requires positive data')}
  if(sum(as.integer(data)==data)!=length(data)){stop('Poisson distribution requires integer data')}
  n=length(data)
  if(Q>(n/2)){stop(paste('Q is larger than the maximum number of segments',n/2))}
  all.seg=matrix(0,ncol=n,nrow=n)
  for(i in 1:n){
    sumx=0
    for(j in i:n){
      len=j-i+1
      sumx=sumx+data[j]
      if(sumx==0){
        all.seg[i,j]=-Inf
      }
      else{
        all.seg[i,j]=sumx*log(sumx)-sumx*log(len)
      }
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



binseg.meanvar.poisson=function(data,Q=5,pen=0){
  # function that uses the BinSeg method to calculate changes in mean & variance where the segments in the data are assumed to be Poisson
  if((sum(data<0)>0)){stop('Poisson distribution requires positive data')}
  if(sum(as.integer(data)==data)!=length(data)){stop('Poisson distribution requires integer data')}
  
  n=length(data)
  if(Q>(n/2)){stop(paste('Q is larger than the maximum number of segments',n/2))}
  y=c(0,cumsum(data))
  
  storage.mode(y)='double'
  cptsout=rep(0,Q) # sets up null vector for changepoint answer
  likeout=rep(0,Q) # sets up null vector for likelihood of changepoints in cptsout
  storage.mode(cptsout)='double'
  op_cps=0
  
  answer=.C('binseg_meanvar_poisson',y,as.integer(n),as.double(pen),as.integer(Q),as.integer(cptsout),likeout,as.integer(op_cps),PACKAGE='changepoint')
  return(list(cps=rbind(answer[[5]],answer[[6]]),op.cpts=answer[[7]],pen=pen))
}


multiple.meanvar.poisson=function(data,mul.method="PELT",penalty="SIC",value=0,Q=5,class=TRUE,param.estimates=TRUE){
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
  else if((penalty=="SIC1") || (penalty=="BIC1")){
    value=(diffparam+1)*log(n)
  }
  else if(penalty=="AIC"){
    value=2*diffparam
  }
  else if(penalty=="AIC1"){
    value=2*(diffparam+1)
  }
  else if(penalty=="Hannan-Quinn"){
    value=2*diffparam*log(log(n))
  }
  else if(penalty=="Hannan-Quinn1"){
    value=2*(diffparam+1)*log(log(n))
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
  else if(penalty=="Asymptotic"){
    return('Asymptotic values for the Poisson Distribution are not available yet, please choose an alternative penalty type')
  }
  if(is.null(dim(data))==TRUE){
    # single dataset
    if(mul.method=="PELT"){
      out=PELT.meanvar.poisson(data,value)
      cpts=out
    }
    else if(mul.method=="BinSeg"){
      out=binseg.meanvar.poisson(data,Q,value)
      if(out$op.cpts==0){cpts=n}
      else{cpts=c(sort(out$cps[1,1:out$op.cpts]),n)}
    }
    else if(mul.method=="SegNeigh"){
      out=segneigh.meanvar.poisson(data,Q,value)
      if(out$op.cpts==0){cpts=n}
      else{cpts=c(sort(out$cps[out$op.cpts+1,][out$cps[out$op.cpts+1,]>0]),n)}
    }
    if(class==TRUE){
      ans=new("cpt")
      data.set(ans)=data;cpttype(ans)="mean and variance";method(ans)=mul.method; distribution(ans)="Poisson";pen.type(ans)=penalty;pen.value(ans)=value;cpts(ans)=cpts
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
        out=c(out,list(PELT.meanvar.poisson(data[i,],value)))
      }
      cpts=out
    }
    else if(mul.method=="BinSeg"){
      for(i in 1:rep){
        out=c(out,list(binseg.meanvar.poisson(data[i,],Q,value)))
        if(class==TRUE){
          if(out[[i]]$op.cpts==0){cpts[[i]]=n}
          else{cpts[[i]]=c(sort(out[[i]]$cps[1,1:out[[i]]$op.cpts]),n)}
        }
      }
    }
    else if(mul.method=="SegNeigh"){
      for(i in 1:rep){
        out=c(out,list(segneigh.meanvar.poisson(data[i,],Q,value)))
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
        data.set(ans[[i]])=data[i,];cpttype(ans[[i]])="mean and variance"; method(ans[[i]])=mul.method;distribution(ans[[i]])="Poisson";pen.type(ans[[i]])=penalty;pen.value(ans[[i]])=value;cpts(ans[[i]])=cpts[[i]]
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
