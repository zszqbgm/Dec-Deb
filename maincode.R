library(Matrix)
library(glmnet)
library(expm)
library(flare)
library(mvtnorm)
library(scalreg)
library(MASS)
library(RSpectra)

pro_trim_betahat=list()
guo_ddl_betahat=list()
java_dl_betahat=list()

n=
p=
struct='cs'
struct='indep'
repeatnum=1000
  
for (t in 1:repeatnum){
  sigmaE=generate_sigmaE(struct,ratio=0.5,p)
  data_list=generate_dataset_cs(n,p,q=3,s=5,sigmaE)
  X=data_list$X
  Y=data_list$Y
  
  #proposed method:rho=0.5
  mu=0.5*sqrt(log(p)/n)
  pro_trim_betahat[[t]]=decon_deb_est(X,Y,rho=0.5,lambda1=NULL,
                                      mu=mu,maxt1=200,mindif1=0.000001,maxt2=500,mindif2=0.000001,
                                      alpha=0.05)$betahat.debiased
  
  #guo_ddl:rho=0.5
  guo_ddl_betahat[[t]]=DDL(X,Y,index=c(1),rho=0.5,rhop=0.5,alpha=0.05,lambda1=NULL)$est_ddl
  
  #java:
  java_dl_betahat[[t]]=SSLasso(X,Y,alpha=0.05,lambda=NULL,mu=NULL,intercept=FALSE,resol=1.3,maxiter=200,threshold=1e-2,verbose=TRUE)
  
  print(t)
}


#######################################metrics##############################################
n=
p=
q=3
s=5
repeatnum=1000
beta=matrix(c(rep(1,s),rep(0,p-s)),p,1,byrow = TRUE)

bias.fun(pro_trim_betahat,beta,repeatnum)
bias.fun(guo_ddl_betahat,beta,repeatnum)
bias.fun(java_dl_betahat,beta,repeatnum)

rmse.fun(pro_trim_betahat,beta,repeatnum)$rmse
rmse.fun(guo_ddl_betahat,beta,repeatnum)$rmse
rmse.fun(java_dl_betahat,beta,repeatnum)$rmse

se.fun(repeatnum,pro_trim_betahat)
se.fun(repeatnum,guo_ddl_betahat)
se.fun(repeatnum,java_dl_betahat)

