decision<-function(tau,null,alt,penalty="SIC",n=0,diffparam=1,value=0){
	if((length(tau)!=length(null))||(length(tau)!=length(alt))){
		stop("Lengths of tau, null and alt do not match")
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
	single.decision=function(tau,null,alt,n=0,diffparam=1,value=0){
		teststat=null-alt
		if(teststat>=value){return(tau)}
		else{return(n)}
	}
	if(length(tau)==1){
		out=single.decision(tau,null,alt,n,diffparam,value)
		names(out)="cpt"
		return(list(cpt=out,pen=value))
	}
	else{
		rep=length(tau)
		out=NULL
		for(i in 1:rep){
			out[i]=single.decision(tau[i],null[i],alt[i],n,diffparam,value)
		}
		names(out)=rep("cpt",rep)
		return(list(cpt=out,pen=value))
	}
}
