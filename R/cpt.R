cpt.mean=function(data,penalty="SIC",value=0,method="AMOC",Q=5,dist="Normal",class=TRUE,param.estimates=TRUE){
	if(!((dist=="Normal")||(dist=="CUSUM"))){ stop("Invalid distribution, must be Normal or CUSUM") }
	if(dist=="Normal"){
		if(method=="AMOC"){
			return(single.mean.norm(data,penalty,value,class,param.estimates))
		}
		else if(method=="PELT" || method=="BinSeg"){
			return(multiple.mean.norm(data,mul.method=method,penalty,value,Q,class,param.estimates))
		}
		else if(method=="SegNeigh"){
			warning("SegNeigh is computationally slow, use PELT instead")
			return(multiple.mean.norm(data,mul.method=method,penalty,value,Q,class,param.estimates))
		}	
		else{
			stop("Invalid Method, must be AMOC, PELT, SegNeigh or BinSeg")
		}
	}
	else if(dist=="CUSUM"){
		if(method=="AMOC"){
			return(single.mean.cusum(data,penalty,value,class,param.estimates))
		}
		else if(method=="SegNeigh" || method=="BinSeg"){
			return(multiple.mean.cusum(data,mul.method=method,penalty,value,Q,class,param.estimates))
		}
		else{
			stop("Invalid Method, must be AMOC, SegNeigh or BinSeg")
		}
	}
}

#cpt.reg=function(data,penalty="SIC",value=0,method="AMOC",Q=5,dist="Normal",class=TRUE,param.estimates=TRUE){
#	if(dist !="Normal"){ stop("Invalid distribution, must be Normal") }
#	if(method=="AMOC"){
#		return(single.reg.norm(data,penalty,value,class,param.estimates))
#	}
#	else if(method=="PELT" || method=="BinSeg"){
#		return(multiple.reg.norm(data,mul.method=method,penalty,value,Q,class,param.estimates))
#	}
#	else if(method=="SegNeigh"){
#		warning("SegNeigh is computationally slow, use PELT instead")
#		return(multiple.reg.norm(data,mul.method=method,penalty,value,Q,class,param.estimates))
#	}
#	else{
#		stop("Invalid Method, must be AMOC, PELT, SegNeigh or BinSeg")
#	}
#}

cpt.var=function(data,penalty="SIC",value=0,know.mean=FALSE, mu=-1000,method="AMOC",Q=5,dist="Normal",class=TRUE,param.estimates=TRUE){
	if(dist =="Normal"){
		if(method=="AMOC"){
			return(single.var.norm(data,penalty,value,know.mean,mu,class,param.estimates))
		}
		else if(method=="PELT" || method=="BinSeg"){
			return(multiple.var.norm(data,mul.method=method,penalty,value,Q,know.mean,mu,class,param.estimates))
		}
		else if(method=="SegNeigh"){
			warning("SegNeigh is computationally slow, use PELT instead")
			return(multiple.var.norm(data,mul.method=method,penalty,value,Q,know.mean,mu,class,param.estimates))
		}
		else{
			stop("Invalid Method, must be AMOC, PELT, SegNeigh or BinSeg")
		}
	}
	else if(dist=="CSS"){
		if(method=="AMOC"){
			return(single.var.css(data,penalty,value,class,param.estimates))
		}
		else if(method=="PELT" || method=="SegNeigh" || method=="BinSeg"){
			return(multiple.var.css(data,mul.method=method,penalty,value,Q,class,param.estimates))
		}
		else{
			stop("Invalid Method, must be AMOC, SegNeigh or BinSeg")
		}
	}
	else{
		stop("Invalid distribution, must be Normal or CSS")
	}
}

cpt.meanvar=function(data,penalty="SIC",value=0,method="AMOC",Q=5,dist="Normal",class=TRUE,param.estimates=TRUE,shape=0){
	if(dist=="Normal"){
		if(method=="AMOC"){
			return(single.meanvar.norm(data,penalty,value,class,param.estimates))
		}
		else if(method=="PELT" || method=="BinSeg"){
			return(multiple.meanvar.norm(data,mul.method=method,penalty,value,Q,class,param.estimates))
		}
		else if(method=="SegNeigh"){
			warning("SegNeigh is computationally slow, use PELT instead")
			return(multiple.meanvar.norm(data,mul.method=method,penalty,value,Q,class,param.estimates))
		}
		else{
			stop("Invalid Method, must be AMOC, PELT, SegNeigh or BinSeg")
		}
	}
	else if(dist=="Gamma"){
		if(method=="AMOC"){
			return(single.meanvar.gamma(data,shape,penalty,value,class,param.estimates))
		}
		else if(method=="PELT" || method=="BinSeg"){
			return(multiple.meanvar.gamma(data,shape,mul.method=method,penalty,value,Q,class,param.estimates))
		}
		else if(method=="SegNeigh"){
			warning("SegNeigh is computationally slow, use PELT instead")
			return(multiple.meanvar.gamma(data,shape,mul.method=method,penalty,value,Q,class,param.estimates))
		}
		else{
			stop("Invalid Method, must be AMOC, PELT, SegNeigh or BinSeg")
		}
	}
	else if(dist=="Exponential"){
		if(method=="AMOC"){
			return(single.meanvar.exp(data,penalty,value,class,param.estimates))
		}
		else if(method=="PELT" || method=="BinSeg"){
			return(multiple.meanvar.exp(data,mul.method=method,penalty,value,Q,class,param.estimates))
		}
		else if(method=="SegNeigh"){
			warning("SegNeigh is computationally slow, use PELT instead")
			return(multiple.meanvar.exp(data,mul.method=method,penalty,value,Q,class,param.estimates))
		}
		else{
			stop("Invalid Method, must be AMOC, PELT, SegNeigh or BinSeg")
		}
	}
	else{
		stop("Invalid Distribution, must be Normal, Gamma or Exponential")
	}
}

