setwd("C:\\Users\\jiaojin1\\Documents\\GitHub\\adcomp\\tmb_examples")
library(TMB)
compile("thetalog.cpp")
dyn.load(dynlib("thetalog"))

## Read data
Y <- scan("thetalog.dat", skip=3, quiet=TRUE)
data <- list(Y=Y)

## Parameter initial guess
parameters <- list(
  X = data$Y*0,
  logr0 = 0,
  logtheta = 0,
  logK = 6,
  logQ = 0,
  logR = 0
)

## Fit model
obj <- MakeADFun(data, parameters, random="X", DLL="thetalog")
newtonOption(obj, smartsearch=FALSE)

opt <- nlminb(obj$par, obj$fn, obj$gr)
rep<-sdreport(obj)
rep

###output all results including fixed and random variables
X_hat<-summary(rep,"random")
params<-summary(rep,"fixed",p.value=TRUE)
X_hat<-as.data.frame(X_hat)
x1<-X_hat$Estimate

###maximum likelihood value under estimated parameters
loglik=obj$fn()[1]
loglik

#####this number is small, it may be only calculated by part of likelihood function instead of the whole MLE.

###varify whether the above fn is the likelihood value under estimated parameters

ans=0
for(i in 1:(length(x1)-1))
{
  m<-x1[i]+exp(params[1,1])*(1.0-(exp(x1[i])/exp(params[3,1]))^exp(params[2,1]))
  ans=ans-dnorm(x1[i+1],m,sqrt(exp(params[4,1])),log=TRUE)
}
for (i in 1:length(x1))
{
  ans=ans-dnorm(data$Y[i],x1[i],sqrt(exp(params[5,1])),log=TRUE)
}
ans

##return ans almost gives the similar results as pomp MLE by vanila sequential Monte Carlo filter by pfilter. 

###conclusion: the obj$fn()[1] does not directly return the likelihood value under the mle.


#####pomp for thetalog
require(pomp)
require(ggplot2)
require(magrittr)
require(reshape2)
Y <- scan("thetalog.dat", skip=3, quiet=TRUE)
data <- list(Y=Y)

dataset=data.frame(Y=data$Y,t=seq(1,length(data$Y),1))


###version skip pointer definitation
skel <- Csnippet("Dx = r0*(1.0-pow(exp(x)/K,theta));")

partrans <- "
Tlogr0 = exp(logr0);
TlogK = exp(logK);
Tlogtheta = exp(logtheta);
TlogQ = exp(logQ);
TlogR = exp(logR);
"

paruntrans <- "
Tlogr0 = log(logr0);
TlogK = log(logK);
Tlogtheta = log(logtheta);
TlogQ = log(logQ);
TlogR = log(logR);
"

theta_rpro<-Csnippet("
                     double mean;
                     double dw;
                     dw=rnorm(0,sqrt(Q));
                     mean=(r0*(1.0-pow(exp(x)/K,theta)+dw))*dt;
                     x += mean;
                     ")
theta_rmea<-Csnippet("
                     Y=rnorm(x,sqrt(R));
                     ")

theta_dmea<-Csnippet("
                     lik = dnorm(Y, x, sqrt(R), give_log);
                     ")

pomp(dataset,times="t",t0=0,
     rprocess=euler.sim(theta_rpro,delta.t=0.5),
     rmeasure=theta_rmea,
     dmeasure=theta_dmea,
     statenames="x",
     initializer=Csnippet("x = 0;"),
    # fromEstimationScale=Csnippet(partrans),
    # toEstimationScale=Csnippet(paruntrans), 
     skeleton=vectorfield(skel),
     paramnames=c("r0","K","theta","Q","R","x.0"),params = c(r0=1,K=exp(6),theta=1,Q=1,R=1,x.0=0)
) -> theta_pomp

##use pomp function for simulation
sims1<-simulate(theta_pomp,as.data.frame=T,include.data=T)
###plot the simulated results
sims1 %>%
  melt(id=c("sim","time")) %>%
  ggplot(aes(x=time,group=sim,color=sim,y=value))+
  geom_line()+
  facet_wrap(~variable,scales="free_y")

####get the likelihood value under sequential Monte Carlo filter; MIF2 can also be used in future

pf<-pfilter(theta_pomp,params=c(r0=1,K=exp(6),theta=1,Q=1,R=1,x.0=0),Np=2000)
pf@loglik

########
########parameter estimate under MLE through mle2
##construct negative log likelihood function
library(stats4)
library(bbmle)
library(optimx)

params=c(r0 = 1,K = exp(6),theta = 1,Q = 1,R = 1)

pomp.loglik<-function(params){
  mean=trajectory(theta_pomp,params=c(params,x.0=0))
    #c(logr0=par[1],logK=par[2],logtheta=par[3],logQ=par[4],logR=par[5],x.0=0)
  sum(dnorm(Y,mean,sqrt(params["R"]),log=TRUE))
}

f1<-function(par){
  params<-c(r0=par[1],K=par[2],theta=par[3],Q=par[4],R=par[5])
  pomp.loglik(params)
}

f1(par=c(1,exp(6),1,1,1))

optimx(c(1,exp(6),1,1,1),f1,method="Nelder-Mead")


pompMLE<-optim(f1,par=c(1,exp(6),1,1,1))

nlminb(c(1,exp(6),1,1,1),f1)
  
mle2(f1,start=list(par=c(1,exp(6),1,1,1)),data=data.frame(Y))
,data=data.frame(obs(theta_pomp)))

nlminb(c(1,exp(6),1,1,1),loglik.normal)

ans=0
for(i in 1:(length(x1)-1))
{
  m<-x1[i]+exp(params[1,1])*(1.0-(exp(x1[i])/exp(params[3,1]))^exp(params[2,1]))
  ans=ans-dnorm(x1[i+1],m,sqrt(exp(params[4,1])),log=TRUE)
}
for (i in 1:length(x1))
{
  ans=ans-dnorm(data$Y[i],x1[i],sqrt(exp(params[5,1])),log=TRUE)
}
ans
#params=c(r0=1,K=exp(6),theta=1,Q=1,R=1)
loglik.normal<-function(params){
  x <- trajectory(theta_pomp1,params=params)
  sum(dnorm(x=obs(theta_pomp1),mean=x,sd=sqrt(params["R"]),log=TRUE))
}

f<-function(par){
  params<-c(r0=par[1],K=par[2],theta=par[3],Q=par[4],R=par[5],x.0=0)
  loglik.normal(params)
}

pompfit<-nlminb(c(1,exp(6),1,1,1),f)
pompfit

pompt<-system.time(nlminb(c(1,exp(6),1,1,1),f))
pompt
#> pompt
#user  system elapsed 
#0.19    0.18    0.39 




#######
########parameter estimate under MIF
########parameter estimate under MIF2
######MIF2 parameter estimate
require(doParallel)
registerDoParallel(CORES)

param.tab1<-data.frame(r0=c(0.000,2),K=c(exp(5.5),exp(6.5)),theta=c(0.0001,2),Q=c(0.0001,2),R=c(0.0001,2),x.0=c(0.0001,2))
rownames(param.tab1)=c("lower.bound","upper.bound")

tic <- Sys.time()
mpar <- foreach(
  i=1:2,
  .packages=c('pomp'),
  .inorder=FALSE) %dopar% 
  {
    Sys.sleep(i*.1)
    NMIF<-20#0 ##update for flux
    NP<-1000#0 ##update for flux
    METHOD="mif2"
    source("Daphnia_1.R")
    param.tab <- param.tab1
    LV.pars <- c("r0","K","theta","Q","R","x.0")
    LV.ivps <- c("x.0")
    LV.rw.sd<- rw.sd(r0=0.02, K=0.02,theta=0.02,Q=0.02,R=0.02,x.0=ivp(0.1))
    
    LV.hyperparams <-
      list(min=unlist(param.tab["lower.bound",]),max=unlist(param.tab["upper.bound",]))
    
    LV.rprior <- function(hyperparams, ...)
    {
      r<-runif(length(hyperparams$min),min=hyperparams$min,max=hyperparams$max)
      names(r) <- names(hyperparams$min)
      return(r)
    }
    set.seed(8100+i)
    Sys.sleep(i*0.1)
    th.draw <-LV.rprior(LV.hyperparams)
    m<-try(mif2(theta_pomp1,
                Nmif=NMIF,
                start=th.draw, # we will initialize
                rw.sd=LV.rw.sd,
                Np=NP,
                cooling.type='geometric',
                cooling.fraction= 0.3,
                max.fail=200,
                transform=TRUE
    ))
    list(pomp=m,start=th.draw)
  }

require(magrittr)
require(plyr)
require(ggplot2)

??eta_prof<-readRDS("eta-profile1.rds")
eta_prof %<>%
  subset(nfail.max==0) %>%
  mutate(eta=exp(signif(log(eta),5))) %>%
  ddply(~eta,subset,rank(-loglik)<=2)

eta_prof<-subset(eta_prof,loglik>(-260))

eta_prof %>%
  ggplot(aes(x=log(eta),y=loglik))+
  geom_point() +
  geom_smooth(method="loess",span=0.25) + labs(y="Log-likelihood",x="Log(eta)") + 
  geom_vline(xintercept = c(log(0.0424),log(0.45))) + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),panel.background = element_blank())




#####another version pomp using pointer(clapse may due to pointer usage)
theta_rprocess<-Csnippet("
                         double *x=&X1, mean;
                         int i, n=(int)N;
                         for (i=0;i< n; i++) {
                         mean=x[i]+r0*(1.0-pow(exp(x[i])/K,theta));
                         x[i+1]=rnorm(mean,sqrt(Q));}
                         ")

rInit <- Csnippet("
                  double *x = &X1;
                  int i, n = (int) N;
                  for (i=0; i < n; i++) x[i] = 0.0;
                  ")

#theta_dmeas <- Csnippet("lik = dnorm(Y,X,sqrt(R),give_log);")
#theta_rmeas <- Csnippet("Y = rnorm(X,sqrt(R));")
parameter<-c(N=200,r0=1,K=1,theta=1,Q=1,R=1)
statX = c(sprintf("X%d",1:length(data$Y)))

pomp(dataset,times="t",t0=0,
     rprocess=euler.sim(theta_rprocess,delta.t=1),
     measurement.model=Y~norm(-X1,sqrt(R)),
     statenames=statX,
     paramnames=c("N","r0","K","theta","Q","R"),params = parameter,initializer = rInit
) -> theta_pomp2

sims2 <-simulate(theta_pomp2,as.data.frame=TRUE,include.data=TRUE)
pompt2<-system.time(simulate(theta_pomp2,as.data.frame=TRUE,include.data=TRUE))
pompt2

# rdf
sims2 %>%
  melt(id=c("sim","time")) %>%
  ggplot(aes(x=time,group=sim,color=sim,y=value))+
  geom_line()+
  facet_wrap(~variable,scales="free_y")


