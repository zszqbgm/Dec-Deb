generate_sigmaE=function(struct,ratio,p){
  sigmaE=matrix(0,nrow=p,ncol=p) 
  if(struct=="indep") {
    sigmaE=diag(p) 
  } else if(struct=="ar1") {
    sigmaE=ratio^(abs(outer(1:p,1:p,"-")))
  } else if(struct=="cs") {
    sigmaE=ratio*rep(1,p)%*%t(rep(1,p)) + (1-ratio)*diag(p)
  }
  return(sigmaE)
}

generate_dataset=function(n,p,q,s,sigmaE){
  H=matrix(rnorm(n*q,mean=0,sd=1),n,q,byrow = TRUE)
  psi=matrix(rnorm(q*p,mean=1,sd=1),q,p,byrow = TRUE)
  E=rmvnorm(n,mean=rep(0,p),sigma = sigmaE)
  X=H%*%psi+E
  beta=matrix(c(rep(1,s),rep(0,p-s)),p,1,byrow = TRUE)
  xi=matrix(rnorm(q*1,mean=0,sd=1),q,1,byrow = TRUE)
  e=matrix(rnorm(n*1,mean=0,sd=1),n,1,byrow = TRUE)
  Y=X%*%beta+H%*%xi+e
  b=solve(t(psi)%*%psi+sigmaE)%*%t(psi)%*%xi
  return_list=list("X"=X,"Y"=Y,"H"=H,"psi"=psi,"xi"=xi,"beta"=beta,"E"= E,"e"=e,'b'=b)
  return(return_list)
}

generate_dataset_unconf=function(n,p,q,s,sigmaE){
  H=matrix(rnorm(n*q,mean=0,sd=1),n,q,byrow = TRUE)
  X=rmvnorm(n,mean=rep(0,p),sigma = sigmaE)
  beta=matrix(c(rep(1,s),rep(0,p-s)),p,1,byrow = TRUE)
  e=matrix(rnorm(n*1,mean=0,sd=1),n,1,byrow = TRUE)
  xi=matrix(rnorm(q*1,mean=0,sd=1),q,1,byrow = TRUE)
  Y=X%*%beta+H%*%xi+e
  return_list=list("X"=X,"Y"=Y,"beta"=beta,"E"= E,"e"=e)
  return(return_list)
}





