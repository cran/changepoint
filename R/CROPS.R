CROPS <- function(data, penalty="CROPS", pen.value, method="PELT", test.stat="Normal", class=TRUE, param.est=TRUE, minseglen, func){
  if(method != "PELT"){stop('CROPS is a valid penalty choice only if method="PELT", please change your method or your penalty.')}
  mu <- mean(data)
  sumstat=cbind(c(0,cumsum(coredata(data))),c(0,cumsum(coredata(data)^2)),cumsum(c(0,(coredata(data)-mu)^2)))
  
  switch(test.stat,
    "Normal" = {stat = "norm"},
    "Exponential" = {stat = "exp"},
    "Gamma" = {stat = "gamma"},
    "Poisson" = {stat = "poisson"},
    {stop("Only Normal, Exponential, Gamma and Poisson are valid test statistics")}
  )
  costfunc = paste0(func, ".", stat)
  
  out = range_of_penalties(sumstat, cost=costfunc, min_pen=pen.value[1], max_pen=pen.value[2], minseglen=minseglen)
  
  #calculate the matrix 
  #newmatrix = list()
  #newmatrix[[1]] = out[[2]]
  #newmatrix[[2]] = out[[1]][1,]
  #s = cbind(newmatrix[[1]], newmatrix[[2]])
  #b = t(sapply(a@cpts.full, '[', 1:max(sapply(a@cpts.full, length)))) #makes matrix with NA
  if(func=="var"){
    cpttype="variance"
  }else if(func=="meanvar"){
    cpttype="mean and variance"
  }else{
    cpttype="mean"
  }
 # browser()
  
  
  if(class==TRUE){
      ans = class_input(data=data,cpttype=cpttype, method="PELT", test.stat=test.stat, penalty=penalty, pen.value=pen.value, minseglen=minseglen, param.estimates=param.est, out=out)
#     ans=new("cpt.range")
#     data.set(ans)=data;cpttype(ans)=func;method(ans)="PELT";test.stat(ans)=test.stat; pen.type(ans)=penalty;pen.value(ans)=pen.value
#     
#     m = t(sapply(out[[2]], '[', 1:max(sapply(out[[2]], length))))
#     
#     cpts.full(ans) = m
#     pen.value.full(ans) = out[[1]][1,]
#  #   pen.value.input(ans) = pen.value
    
    return(ans)
  }else{return(out)}
  
}