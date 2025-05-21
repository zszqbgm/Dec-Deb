library(Matrix)
library(glmnet)
library(expm)
library(flare)
library(mvtnorm)
library(scalreg)
library(MASS)
library(RSpectra)

proposed_betahat_0.1rho=list()
proposed_betahat_0.3rho=list()
proposed_betahat_0.5rho=list()
proposed_betahat_0.7rho=list()
proposed_betahat_0.9rho=list()

n=250
p=300
struct='cs'
struct='indep'
repeatnum=100
for (t in 1:repeatnum){
  sigmaE=generate_sigmaE(struct,ratio=0.5,p)
  data_list=generate_dataset_cs(n,p,q=3,s=5,sigmaE)
  X=data_list$X
  Y=data_list$Y
  
  #proposed method:rho=0.1
  mu=0.5*sqrt(log(p)/n)
  proposed=decon_deb_est(X,Y,rho=0.1,lambda1=NULL,
                         mu=mu,maxt1=200,mindif1=0.000001,maxt2=500,mindif2=0.000001,
                         alpha=0.05)
  proposed_betahat_0.1rho[[t]]=proposed$betahat.debiased
  print('proposed rho=0.1')
  
  #proposed method:rho=0.3
  proposed=decon_deb_est(X,Y,rho=0.3,lambda1=NULL,
                         mu=mu,maxt1=200,mindif1=0.000001,maxt2=500,mindif2=0.000001,
                         alpha=0.05)
  proposed_betahat_0.3rho[[t]]=proposed$betahat.debiased
  print('proposed rho=0.3')
  
  #proposed method:rho=0.5
  proposed=decon_deb_est(X,Y,rho=0.5,lambda1=NULL,
                         mu=mu,maxt1=200,mindif1=0.000001,maxt2=500,mindif2=0.000001,
                         alpha=0.05)
  proposed_betahat_0.5rho[[t]]=proposed$betahat.debiased
  print('proposed rho=0.5')
  
  #proposed method:rho=0.7
  proposed=decon_deb_est(X,Y,rho=0.7,lambda1=NULL,
                         mu=mu,maxt1=200,mindif1=0.000001,maxt2=500,mindif2=0.000001,
                         alpha=0.05)
  proposed_betahat_0.7rho[[t]]=proposed$betahat.debiased
  print('proposed rho=0.7')
  
  #proposed method:rho=0.9
  proposed=decon_deb_est(X,Y,rho=0.9,lambda1=NULL,
                         mu=mu,maxt1=200,mindif1=0.000001,maxt2=500,mindif2=0.000001,
                         alpha=0.05)
  proposed_betahat_0.9rho[[t]]=proposed$betahat.debiased
  print('proposed rho=0.9')
  
  print(t)
}

#######################################metrics##############################################
n=250
p=300
q=3
s=5
repeatnum=100
beta=matrix(c(rep(1,s),rep(0,p-s)),p,1,byrow = TRUE)

bias.fun(proposed_betahat_0.1rho,beta,repeatnum)
bias.fun(proposed_betahat_0.3rho,beta,repeatnum)
bias.fun(proposed_betahat_0.5rho,beta,repeatnum)
bias.fun(proposed_betahat_0.7rho,beta,repeatnum)
bias.fun(proposed_betahat_0.9rho,beta,repeatnum)

rmse.fun(proposed_betahat_0.1rho,beta,repeatnum)$rmse
rmse.fun(proposed_betahat_0.3rho,beta,repeatnum)$rmse
rmse.fun(proposed_betahat_0.5rho,beta,repeatnum)$rmse
rmse.fun(proposed_betahat_0.7rho,beta,repeatnum)$rmse
rmse.fun(proposed_betahat_0.9rho,beta,repeatnum)$rmse

se.fun(repeatnum,proposed_betahat_0.1rho)
se.fun(repeatnum,proposed_betahat_0.3rho)
se.fun(repeatnum,proposed_betahat_0.5rho)
se.fun(repeatnum,proposed_betahat_0.7rho)
se.fun(repeatnum,proposed_betahat_0.9rho)

