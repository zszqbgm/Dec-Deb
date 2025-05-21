bias.fun=function(method.betahat,beta,repeatnum){
  betahat=rep(NA,repeatnum)
  for (t in 1:repeatnum){
    betahat[t]=method.betahat[[t]][1]
  }
  bias=mean(betahat)-beta[1]
  return('bias'=bias)
}
rmse.fun=function(method.betahat,beta,repeatnum){
  bias=rep(NA,repeatnum)
  for (t in 1:repeatnum){
    betahat=method.betahat[[t]][1]
    bias[t]=(betahat-beta[1])^2
  }
  bias.average=mean(bias)
  rmse=sqrt(bias.average)
  return(list('biasall'=bias,'bias'=bias.average,'rmse'=rmse))
}
se.fun=function(repeatnum,method.betahat){
  betahat.rep=rep(NA,repeatnum)
  for (t in 1:repeatnum){
    betahat.rep[t]=method.betahat[[t]][1]
  }
  se=sd(betahat.rep)
  return('se'=se)
}
