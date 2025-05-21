library(glmnet)
#spectral transformation
Spe_trans_Q=function(X,rho){
  n=nrow(X)
  p=ncol(X)
  UDV_list=svd(X)
  U=UDV_list$u
  D=UDV_list$d
  D_tilde =  D[floor(rho * min(n, p))] / D[1:(floor(rho * min(n, p)))]
  D_tilde = sqrt(1 - D_tilde)
  F_trim = U[, 1:(floor(rho * min(n, p)))] %*% diag(D_tilde)
  Q = diag(n) - F_trim %*% t(F_trim)
  return(Q)
}
#\hat{\beta}^{init}
Estimate_init_coeff = function(X,Y,Q,lambda1=NULL){
  Xtilde = Q %*% X
  Ytilde = Q %*% Y
  weights=apply(Xtilde,2,function(x) norm(x,type='2')/sqrt(nrow(X)))
  if (is.null(lambda1)){
    fit=cv.glmnet(Xtilde,Ytilde,penalty.factor=weights)
    lambda1=fit$lambda.min
    betahatinit=as.vector(coef(fit,lambda1)[-1])
  } else {
    fit=glmnet(Xtilde,Ytilde,lambda=lambda1,family = c("gaussian"),alpha=1,intercept = FALSE)
    betahatinit=as.vector(coef(fit))
  }
  if (all(betahatinit==0)){print('wrong')}
  return("initpara"=betahatinit)
}
#\hat{\sigma}_e
Estimate_Sigmae=function(betahatinit,Y,X,Q){
  divisor=sum(diag(Q%*%Q))
  error=(norm(Q%*%Y-Q%*%X%*%betahatinit,type='2'))^2
  sigma_e_hat=(error/divisor)^0.5
  return('Sigmaehat'=sigma_e_hat)
}
#softthreshold
softthres=function(x,mu){
  if (x>mu){return (x-mu)}
  else {if (x< (-mu)){return (x+mu)}
    else {return (0)}}
}
#\hat{\Theta}
optdual=function(X,mu,e,u,maxt1,mindif1){
  n=nrow(X)
  p=ncol(X) 
  A=t(X)%*%X/(4*n)
  if (is.null(mu)){mu=0.5*sqrt(log(p)/n)}
  eta=rep(0,p)
  eta[u]=(1-mu)/(2*A[u,u])
  t=1;dif=10
  while ((t<=maxt1)&(dif>=mindif1)) {
    eta.old=eta
    for (j in 1:p){
      a=2*(A[,j]%*%eta-A[j,j]*eta[j])-e[j]
      eta[j]=-(softthres(a,mu))/(2*A[j,j])
    }
    t=t+1
    dif=norm((eta-eta.old),type='2')
  }
  return.list=list('dual'=eta,'T'=t,'converge'=dif)
  return(return.list)
}
opttheta=function(dual,X,Q,maxt2,mindif2){
  p=ncol(X)
  n=nrow(X)
  Z=t(X)%*%Q%*%Q%*%Q%*%Q%*%X
  Sigmahat=t(X)%*%Q%*%Q%*%X/n
  theta=rep(0,p)
  dif=10;t=1
  while((t<=maxt2)&(dif>=mindif2)){
    theta.old=theta
    for (j in 1:p){
      a=Z[j,]%*%theta-Z[j,j]*theta[j]
      theta[j]=(n*Sigmahat[j,]%*%dual/2-a)/Z[j,j]
    }
    t=t+1
    dif=norm((theta-theta.old),type='2')
  }
  return.list=list('theta'=theta,'T'=t,'converge'=dif)
  return(return.list)
}
optprogram=function(X,Q,mu,maxt1,mindif1,maxt2,mindif2){
  p=ncol(X)
  n=nrow(X)
  I=diag(1,p,p)
  Theta=matrix(NA,p,p)
  for (u in 1:p){
    dual=optdual(X,mu,I[,u],u,maxt1,mindif1)$dual 
    theta=opttheta(dual,X,Q,maxt2,mindif2)$theta  
    Theta[,u]=theta
  }
  return('Theta'=t(Theta))
}
Theta_est=function(X,Q,mu,maxt1,mindif1,maxt2,mindif2){
  p=ncol(X)
  n=nrow(X)
  Xtilde=Q%*%X
  emp_cov=t(Xtilde)%*%Xtilde/n
  if (det(emp_cov)!=0){
    Theta=solve(emp_cov)
    print('inverse')
  } else{
    Theta=optprogram(X,Q,mu,maxt1,mindif1,maxt2,mindif2)
  }
  return('Theta'=Theta)
}
#deconfounded and debiasaed lasso:
decon_deb_est=function(X,Y,rho,lambda1=NULL,
                       mu,maxt1,mindif1,maxt2,mindif2,
                       alpha){
  p=ncol(X)
  n=nrow(X)
  Q=Spe_trans_Q(X,rho)
  betahat.init=Estimate_init_coeff(X,Y,Q,lambda1=NULL)
  Theta.hat=Theta_est(X,Q,mu,maxt1,mindif1,maxt2,mindif2)
  sigmaehat=Estimate_Sigmae(betahat.init,Y,X,Q)
  #coefficient est
  betahat.debiased=as.numeric(betahat.init+(Theta.hat%*%t(Q%*%X)%*%(Q%*%Y-Q%*%X%*%betahat.init))/n)
  #asy-var est
  Var=((sigmaehat^2)/(n^2))*(Theta.hat%*%t(X)%*%Q%*%Q%*%Q%*%Q%*%X%*%t(Theta.hat))
  SD=(sigmaehat/n)*(sqrt(diag(Theta.hat%*%t(X)%*%Q%*%Q%*%Q%*%Q%*%X%*%t(Theta.hat))))
  #asy-conf 
  interval=qnorm(1-(alpha/2))*SD
  CI.low=betahat.debiased-interval
  CI.up=betahat.debiased+interval
  #asy-pvalue
  pvalue=2*(1-pnorm(abs(betahat.debiased)/SD))
  
  returnList=list('Q'=Q,'betahat.init'=betahat.init,'sigmaehat'=sigmaehat,
                  'Theta.hat'=Theta.hat,
                  'betahat.debiased'=betahat.debiased,
                  'Var'=Var,'SD'=SD,
                  'interval'=interval,'CI.low'=CI.low,'CI.up'=CI.up,'pvalue'=pvalue)
  return(returnList)
}
