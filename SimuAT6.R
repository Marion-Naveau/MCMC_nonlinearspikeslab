rm(list=ls())

library(ggplot2)
library(nlme)
library(cowplot)
library(glmnet)
library(mvnfast)
library(doParallel)
library(microbenchmark)
library(tidyverse)
library(nimble)

S=30                         #nb de simulations
n=200                          #nb d'individus
J=10                            #nb de repetitions
p=500                          #nb de covariables
sigma2=30                      #variance residuelle des y
Gamma2=200                     #variance inter-individuelle des phi_i
gamma2=1                       #variance pour simulation composantes de V_i selon N(0,gamma2) ind?pendantes
mu=1200                         #intercept
beta=c(100,50,20,rep(0,p-3))  #vecteur d'effets fixes
betatild=c(mu,beta)
c0=200                         #parametre de la fct g
c2=300                         #parametre de la fct g, doit etre different de 0
psi=c(c0,c2)

g <- function(phi_i,psi,t_ij){
  return(psi[1]/(1+exp(-(t_ij-phi_i)/psi[2])))
}

nu0=0.01

ncore = 10
cl = makeCluster(ncore)
registerDoParallel(cl)

nb_iter_finish=0
save(nb_iter_finish,file="Nb_iter_finish_AT1.Rdata")

Temps_Tot<-foreach(s = 1:S, .packages = c("nlme","mvnfast","doParallel","ggplot2","cowplot","glmnet","microbenchmark","tidyverse","nimble")) %dopar% {

  print(s)
  # Simulation du modèle
  set.seed(s)
  phi=rep(0,n)
  Vtild=matrix(NA,nrow=n,ncol=p+1)
  Vtild[,1]=rep(1,n)
  for (i in 1:n) {
    # simulation de V_i
    Vtild[i,2:(p+1)]=rmvn(1,mu=rep(0,p),sigma=diag(1,p))
  }
  Vtild[,c(2:(p+1))]=scale(Vtild[,c(2:(p+1))])
  # simulation de phi_i
  for (i in 1:n){
    phi[i]=sum(betatild*Vtild[i,]) + rnorm(1,mean=0,sd=sqrt(Gamma2))
  }
  V=Vtild[,-1]
  # Simulation des données observées y

  Y=rep(0,n*J)                    #vecteur des y_ij avec y_ij=Y[(i-1)*J+j]
  t=rep(seq(150,3000,length=J),n) #vecteur des x_ij avec x_ij=x[(i-1)*J+j]
  for (i in 1:n){
    for (j in 1:J) {
      Y[(i-1)*J+j]=g(phi_i=phi[i],psi=psi,t_ij=t[(i-1)*J+j])+rnorm(1,mean=0,sd=sqrt(sigma2))
    }
  }

  Id=rep(c(1:n),each=J)
  data=data.frame(Id,Y,t)
  id=as.matrix(Id)
  data_input = data %>%
    as_tibble() %>%
    left_join(V %>%
                as_tibble() %>%
                rowid_to_column("Id")) %>%
    rename(times = t, i = Id)

  # Valeurs initiales des paramètres à estimer avec le MCMC-SAEM
  betatild_init=c(1400,rep(100,10),rep(1,p-10))
  sigma2_init=100
  Gamma2_init=5000

  eta_init=rep(400,2)
  Omega_init=rep(20,2) #Omega=(w_1^2,w_2^2,w_3^2)

  #Hyperparamètres
  nu0=0.01             #parametre spike
  nu1=12000         #parametre slab
  nu_Gamma=1        #parametre prior Gamma2
  lb_Gamma=1        #parametre prior Gamma2
  nu_sigma=1        #parametre prior sigma2
  lb_sigma=1        #parametre prior sigma2
  a=1               #parametre prior theta
  b=p               #parametre prior theta
  sigma2_mu=3000^2   #parametre prior mu
  rho2=rep(1200,2)               #parametre prior eta
  nu_omega=1                  #parametre prior omega_m
  lb_omega=1               #parametre prior omega_m

  #Initialisation
  tau2=c(rep(c(rep(1,39),0.9),12),rep(1,20))
  tau=c(0.98,tau2)
  param_init=list(beta=betatild_init,alpha=0.5,Gamma2=Gamma2_init,sigma2=sigma2_init,eta=eta_init,Omega=Omega_init)
  hyperparam=list(nu0=nu0,nu1=nu1,nu_Gamma=nu_Gamma,lb_Gamma=lb_Gamma,nu_sigma=nu_sigma,lb_sigma=lb_sigma,a=a,b=b,sigma2_mu=sigma2_mu,rho2=rho2,nu_omega=nu_omega,lb_omega=lb_omega,tau=tau)

  temps=microbenchmark(SAEM_MAP(niter, nburnin, niterMH_phi,niterMH_psi,Y,t,id, covariables = Vtild, param_init = param_init,hyperparam = hyperparam),times=5,unit="s")

  load("Nb_iter_finish_AT1.Rdata")
  nb_iter_finish=nb_iter_finish+1
  save(nb_iter_finish,file="Nb_iter_finish_AT1.Rdata")

  temps
}
stopCluster(cl)

save(Temps_Tot,file="Res_simuAT1.Rdata")

load("Res_simuAT1.Rdata")

temps_min=rep(0,S)
for (s in 1:S){
  temps_min[s]=summary(Temps_Tot[[s]])$min
}
min=min(temps_min)
