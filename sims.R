######
# libraries
library(pracma)
library(dbscan)

# source(file="functions.R")

scenario=1 # 1,2,3
n=100 # 100,1000, sample size
N=1000 # number of simulation repetitions

#####
# Scenario 1
if(scenario==1)
{
  sig=0.4
  del=0.2
  alpha=beta=3
  dimen=3
  gen_cov=function(n, alpha, beta)
  {
    Z=rbeta(dimen*n,alpha,beta)
    Z=matrix(Z,ncol=dimen)
    return(Z)
  }  
  y_val=function(x)
  {
    Y=rep(0,n)
    for(i in 1:n)
      Y[i]=exp(2*cos(7*x[i,1])*sin(7*x[i,2]))*(4-8*(x[i,3]-0.5)^2)
    return(Y)
  }
  y_val_int=function(x,y,z)
  {
    exp(2*cos(7*x)*sin(7*y))*(4-8*(z-0.5)^2)
  }
  par_t= integral3(y_val_int, del, 1-del, del, 1-del, del, 1-del) 
  K0=c(1,3,6,8,10,12)
  K1=c(5,6,8,10,12,14)
}

#####
# Scenario 2
if(scenario==2)
{
  sig=0.4
  del=0
  alpha=beta=3
  dimen=3
  gen_cov=function(n, alpha, beta)
  {
    d=rbinom(dimen*n,1,1/2)
    Z=rbeta(dimen*n,alpha,beta)
    Z2=rbeta(dimen*n,1,1)
    Z[d==1]=Z2[d==1]
    Z=matrix(Z,ncol=dimen)
    return(Z)
  }  
  y_val=function(x)
  {
    Y=rep(0,n)
    for(i in 1:n)
      Y[i]=exp(2*cos(7*x[i,1])*sin(7*x[i,2]))*(4-8*(x[i,3]-0.5)^2)
    return(Y)
  }
  y_val_int=function(x,y,z)
  {
    exp(2*cos(7*x)*sin(7*y))*(4-8*(z-0.5)^2)
  }
  par_t= integral3(y_val_int, del, 1-del, del, 1-del, del, 1-del) 
  K0=c(1,2,3,5)
  K1=c(5,6,7,8)
}

#####
# Scenario 3
if(scenario==3)
{
  sig=0.2
  del=0.2
  alpha=beta=3
  dimen=2
  gen_cov=function(n, alpha, beta)
  {
    Z=rbeta(dimen*n,alpha,beta)
    Z=matrix(Z,ncol=dimen)
    return(Z)
  }  
  y_val=function(x)
  {
    Y=rep(0,n)
    for(i in 1:n)
      Y[i]=x[i,1]*x[i,2]^2-1+cos(x[i,1]/x[i,2])
    return(Y)
  }
  y_val_int=function(x,y)
  {
    x*y^2-1+cos(x/y)
  }
  par_t= integral2(y_val_int, del, 1-del, del, 1-del, sector = FALSE,
                   reltol = 1e-6)
  par_t = par_t$Q
  K0=c(1,3,4,6,8,10)
  K1=c(4,8,10,12,14,16)
}

a=rep(del,dimen)
b=rep(1-del,dimen)

#####
# Test: single run estimator, default settings
Z=gen_cov(n,alpha,beta)
Y=y_val(Z)+rnorm(n,0,sig)
expec.bias.cor(Y,Z,K=1,L=0,a=a,b=b)
expec.bias.cor(Y,Z,a=a,b=b)
par_t

#####
# Repeated Simulations

res0=matrix(rep(0,length(K0)*N),ncol=length(K0))
colnames(res0)=K0
res1=matrix(rep(0,length(K1)*N),ncol=length(K1))
colnames(res1)=K1

seed=10

for(i in 1:N)
{
  Z=gen_cov(n,alpha,beta)
  Y=y_val(Z)+rnorm(n,0,sig)
  for(j in 1:(length(K0))) res0[i,j]=expec.bias.cor(Y,Z,K=K0[j],L=0,a=a,b=b)
  for(j in 1:(length(K1))) res1[i,j]=expec.bias.cor(Y,Z,K=K1[j],a=a,b=b)
  print(i)
}

round(n^(1/2)*(apply(res0,2,mean)-par_t),4)
round(n^(1/2)*apply(res0,2,sd),4)
round(n^(1/2)*((apply(res0,2,mean)-par_t)^2+apply(res0,2,var))^(1/2),4)

round(n^(1/2)*(apply(res1,2,mean)-par_t),4)
round(n^(1/2)*apply(res1,2,sd),4)
round(n^(1/2)*((apply(res1,2,mean)-par_t)^2+apply(res1,2,var))^(1/2),4)

n
scenario


#####
# Asymptotic normality
if(scenario==1)
{
  nor_test=n^(1/2)*(res0[,1]-par_t)
  # nor_test=n^(1/2)*(res1[,3]-par_t)
  
  a=mean(nor_test)
  b=sd(nor_test)
  
  curve(dnorm(x,a,b),from=-2,to=2, ylab="", xlab="", ylim=c(0,1), main="Estimated density, L=0,K=1")
#  curve(dnorm(x,a,b),from=-2,to=2, ylab="", xlab="", ylim=c(0,1), main="Estimated density, L=1,K=6")
  lines(density(nor_test), lty=2, col="red")
  
  qqnorm(nor_test, pch=20, xlab="", ylab="", main="Estimated density, L=0,K=1")
#  qqnorm(nor_test, pch=20, xlab="", ylab="", main="Estimated density, L=1,K=6")
  qqline(nor_test, col="red")
  
  shapiro.test(nor_test)
  

}