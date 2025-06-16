######### Section 5.2 Poisson Random Effects ##########

### SETUP 

# The Poisson random effects model is of the following form: 
# y_{ij} | \eta_i \sim \text{Poi}(e^{\eta_i}), j = 1, \ldots, J, \\  
# \eta_i | \mu \sim N(\mu,1), i = 1, \ldots, I, \\
# \mu \sim N(0, \sigma_\mu^2).

# Considered Schemes: Euler-Maruyama, semi-implicit Euler, Skew-Symmetric
# Case 1: Initialise at truth
# Case 2: Warm start


## Load functions
source('poisson_functions.R')

## Setup up the model and data

### scenarios 1, 2 and 3 refer to the ones in the paper Livingstone&Zanella(2020)
scenario<-1 # set to 1, 2, or 3 (1=easiest, 3=hardest)
n1<-50  # number of groups 
n<-n1+1  #number of parameters
N<-n1*5  # number of observations
vmu=1/(10^2) # prior precision for global means
blk.ind<-rep(c(1:n1),each=5) # group membership of the observations
if (scenario==1){
  va=1/(1^2) # prior precision for random effects
  true_mu<-5 # define data-generating mu
}
if (scenario==2){
  va=1/(3^2) # prior precision for random effects
  true_mu<-5 # define data-generating mu
}
if (scenario==3){
  va=1/(3^2) # prior precision for random effects
  true_mu<-10 # define data-generating mu
}
set_nested_pois_target(n1=n1,vmu=vmu,va=va,true_mu=true_mu,N=N,blk.ind=blk.ind)
stopifnot(!any(is.na(y.vec)))


## Case 1: Initialise at truth
mse_mu_ULA <- c()
mse_mu_UB <- c()
mse_mu_impULA <- c()

iter <- 100
for (j in 1:iter){
  
  mu <- 5 # right from the truth 
  print(j)
  print(paste("starting mu (from prior) is ",mu))
  start_x<-c(mu,mu+rnorm(n1,0,1/sqrt(va)))
  set_nested_pois_target(n1=n1,vmu=vmu,va=va,true_mu=true_mu,N=N,blk.ind=blk.ind)
  
  T<-5*10^4
  sigma_vec<-c(0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05)
  
  for(sig in sigma_vec){
    print(sig)
    output_ULA<-grad_MCMC(T=T,sigma=sig,n=n,method="ULA",start_x=start_x)
    mse_mu_ULA <- c(mse_mu_ULA,(mean(output_ULA$x_samples[,1]) - true_mu)^2)
    
    output_UB<-grad_MCMC(T=T,sigma=sig,n=n,method="UB",start_x=start_x)
    mse_mu_UB <- c(mse_mu_UB,(mean(output_UB$x_samples[,1]) - true_mu)^2)
    
    output_impULA<-try(grad_MCMC(T=T,sigma=sig,n=n,method="ImpULA",start_x=start_x))
    if (inherits(output_impULA,"try-error")){
      mse_mu_impULA <- c(mse_mu_impULA,NA)
    }else{
      mse_mu_impULA <- c(mse_mu_impULA,(mean(output_impULA$x_samples[,1]) - true_mu)^2)
    }
  }

}
mse_final_mu_ULA <- t(matrix(mse_mu_ULA,ncol=iter))
mse_final_mu_UB <- t(matrix(mse_mu_UB,ncol=iter))
mse_final_mu_impULA <- t(matrix(mse_mu_impULA,ncol=iter))

sim1_true_start_ULA <- mse_final_mu_ULA
sim1_true_start_UB <- mse_final_mu_UB
sim1_true_start_impULA <- mse_final_mu_impULA



## Case 2: Warm start
mse_mu_ULA <- c()
mse_mu_UB <- c()
mse_mu_impULA <- c()

iter <- 100
for (j in 1:iter){
  mu<-rnorm(1,5,1/sqrt(vmu)) # warm start
  print(j)
  print(paste("starting mu (from prior) is ",mu))
  start_x<-c(mu,mu+rnorm(n1,0,1/sqrt(va)))
  set_nested_pois_target(n1=n1,vmu=vmu,va=va,true_mu=true_mu,N=N,blk.ind=blk.ind)
  
  T<-6*10^4
  sigma_vec<-c(0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05)
  
  for(sig in sigma_vec){
    print(sig)
    output_ULA<-grad_MCMC(T=T,sigma=sig,n=n,method="ULA",start_x=start_x)
    mse_mu_ULA <- c(mse_mu_ULA,(mean(output_ULA$x_samples[10001:60000,1]) - true_mu)^2)
    
    output_UB<-grad_MCMC(T=T,sigma=sig,n=n,method="UB",start_x=start_x)
    mse_mu_UB <- c(mse_mu_UB,(mean(output_UB$x_samples[10001:60000,1]) - true_mu)^2)
    
    output_impULA<-try(grad_MCMC(T=T,sigma=sig,n=n,method="ImpULA",start_x=start_x))
    if (inherits(output_impULA,"try-error")){
      mse_mu_impULA <- c(mse_mu_impULA,NA)
    }else{
      mse_mu_impULA <- c(mse_mu_impULA,(mean(output_impULA$x_samples[,1]) - true_mu)^2)
    }
  }
  
}

mse_final_mu_ULA <- t(matrix(mse_mu_ULA,ncol=iter))
mse_final_mu_UB <- t(matrix(mse_mu_UB,ncol=iter))
mse_final_mu_impULA <- t(matrix(mse_mu_impULA,ncol=iter))

sim1_warm_start_ULA <- mse_final_mu_ULA
sim1_warm_start_UB <- mse_final_mu_UB
sim1_warm_start_impULA <- mse_final_mu_impULA


save(sim1_true_start_ULA,
     sim1_true_start_UB,
     sim1_true_start_impULA,
     sim1_warm_start_ULA,
     sim1_warm_start_UB,
     sim1_warm_start_impULA,
     file = "poisson.RData")