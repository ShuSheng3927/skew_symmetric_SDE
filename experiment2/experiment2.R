## LOAD FUNCTIONS FROM GITHUB REPOSITORY
source('functions.R')

# scenarios 1, 2 and 3 refer to the ones in the paper Livingstone&Zanella(2020)
scenario<-1 # set to 1, 2, or 3 (1=easiest, 3=hardest)

#define target (Hierarchical Poisson regression example, with hyperparameters depending on scenario)
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
print(paste("scenario=",scenario))
# generate synth data
set_nested_pois_target(n1=n1,vmu=vmu,va=va,true_mu=true_mu,N=N,blk.ind=blk.ind)
stopifnot(!any(is.na(y.vec)))





## ABOVE: general set up
## BELOW: experiment 2




# MSE Comparison Experiment  SIM1

## start from truth, SIM1a (no burn-in required)

# saved data 
sim1_from_truth_ULA <- mse_final_mu_ULA
sim1_from_truth_UB <- mse_final_mu_UB


mse_mu_ULA <- c()
mse_mu_UB <- c()

iter <- 100
for (j in 1:iter){
  
  # Sample starting point (common to all algorithms)
  # mu<-rnorm(1,5,1/sqrt(vmu)) # warm start
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
    #plot(output_ULA$x_samples[,1],main=sig,ylab="ULA")
    mse_mu_ULA <- c(mse_mu_ULA,(mean(output_ULA$x_samples[,1]) - true_mu)^2)
    
    output_UB<-grad_MCMC(T=T,sigma=sig,n=n,method="UB",start_x=start_x)
    #plot(output_UB$x_samples[,1],main=sig,ylab="UB")
    mse_mu_UB <- c(mse_mu_UB,(mean(output_UB$x_samples[,1]) - true_mu)^2)
  }

}

mse_final_mu_ULA <- t(matrix(mse_mu_ULA,ncol=iter))
mse_final_mu_UB <- t(matrix(mse_mu_UB,ncol=iter))


# sim1a-UB
plot(colMeans(sim1_from_truth_UB),xaxt="n",ylab="MSE",xlab="Step Size",main="UB",type="b")
axis(1, at=1:10,labels=sigma_vec)
# sim1a-ULA
plot(colMeans(sim1_from_truth_ULA),xaxt="n",ylab="MSE",xlab="Step Size",main="ULA",type="b")
axis(1, at=1:10,labels=sigma_vec)
# sim1a-compare
plot(colMeans(sim1_from_truth_UB),xaxt="n",ylab="MSE",xlab="Step Size",main="UB vs. ULA",ylim=c(0,max(colMeans(sim1_from_truth_UB),colMeans(sim1_from_truth_ULA)[1:5])),type="b")
axis(1, at=1:10,labels=sigma_vec)
legend("topright",c("UB","ULA"),col=c("black","red"),pch=1)
points(colMeans(sim1_from_truth_ULA),col="red",type="b")


## warm start,  SIM1b (need to remove burn-in)

# saved data
sim1_warm_start_ULA <- mse_final_mu_ULA
sim1_warm_start_UB <- mse_final_mu_UB


mse_mu_ULA <- c()
mse_mu_UB <- c()

iter <- 100
for (j in 1:iter){
  
  # Sample starting point (common to all algorithms)
  mu<-rnorm(1,5,1/sqrt(vmu)) # warm start
  # mu <- 5 # right from the truth 
  print(j)
  print(paste("starting mu (from prior) is ",mu))
  start_x<-c(mu,mu+rnorm(n1,0,1/sqrt(va)))
  set_nested_pois_target(n1=n1,vmu=vmu,va=va,true_mu=true_mu,N=N,blk.ind=blk.ind)
  
  T<-6*10^4
  sigma_vec<-c(0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05)
  
  for(sig in sigma_vec){
    print(sig)
    output_ULA<-grad_MCMC(T=T,sigma=sig,n=n,method="ULA",start_x=start_x)
    #plot(output_ULA$x_samples[,1],main=sig,ylab="ULA")
    mse_mu_ULA <- c(mse_mu_ULA,(mean(output_ULA$x_samples[10001:60000,1]) - true_mu)^2)
    
    output_UB<-grad_MCMC(T=T,sigma=sig,n=n,method="UB",start_x=start_x)
    #plot(output_UB$x_samples[,1],main=sig,ylab="UB")
    mse_mu_UB <- c(mse_mu_UB,(mean(output_UB$x_samples[10001:60000,1]) - true_mu)^2)
  }
  
}

mse_final_mu_ULA <- t(matrix(mse_mu_ULA,ncol=iter))
mse_final_mu_UB <- t(matrix(mse_mu_UB,ncol=iter))


# sim1b-UB
plot(colMeans(sim1_warm_start_UB),xaxt="n",ylab="MSE",xlab="Step Size",main="UB",type="b")
axis(1, at=1:10,labels=sigma_vec)
# sim1b-ULA
plot(colMeans(sim1_warm_start_ULA),xaxt="n",ylab="MSE",xlab="Step Size",main="ULA",type="b")
axis(1, at=1:10,labels=sigma_vec)
# sim1b-compare
plot(colMeans(sim1_warm_start_UB),xaxt="n",ylab="MSE",xlab="Step Size",main="UB vs. ULA",ylim=c(0,max(colMeans(sim1_warm_start_UB),colMeans(sim1_warm_start_ULA)[1:5])),type="b")
axis(1, at=1:10,labels=sigma_vec)
legend("topright",c("UB","ULA"),col=c("black","red"),pch=1)
points(colMeans(sim1_warm_start_ULA),col="red",type="b")








# variance and quantile of UB over stepsize SIM2
sim2_var_UB <- final_var_mu_UB
sim2_quantile_UB <- final_quantile_mu_UB

# variance of mu
var_mu_UB <- c()
iter <- 100
for (j in 1:iter){
  
  # Sample starting point (common to all algorithms)
  vmu=1/(10^2) # prior precision for global means
  mu<-rnorm(1,5,1/sqrt(vmu)) # warm start
  true_var <- 1/sqrt(vmu)
  # mu <- 5 # right from the truth 
  print(j)
  print(paste("starting mu (from prior) is ",mu))
  start_x<-c(mu,mu+rnorm(n1,0,1/sqrt(va)))
  set_nested_pois_target(n1=n1,vmu=vmu,va=va,true_mu=true_mu,N=N,blk.ind=blk.ind)
  
  T<-5.5*10^4
  sigma_vec <- c(1,2,3,4,5)
  
  for(sig in sigma_vec){
    output_UB<-grad_MCMC(T=T,sigma=sig,n=n,method="UB",start_x=start_x)
    sample_var <- sum((output_UB$x_samples[5001:T,1] - mean(output_UB$x_samples[5001:T,1]))^2)/(T-5001)
    var_mu_UB <- c(var_mu_UB, sample_var)
  }
  
}

final_var_mu_UB <- t(matrix(var_mu_UB,ncol=iter))

plot(colMeans(final_var_mu_UB),xaxt="n",ylab="Sample Variance",xlab="Step Size",main="",type="b")
axis(1, at=1:5,labels=sigma_vec)


# quantile
quantile_mu_UB <- c()

iter <- 100
for (j in 1:iter){
  
  # Sample starting point (common to all algorithms)
  vmu=1/(10^2) # prior precision for global means
  mu<-rnorm(1,5,1/sqrt(vmu)) # warm start
  true_mu_90 <- qnorm(0.9,5,1/sqrt(vmu)) # truth of 90 percentile of mu
  # mu <- 5 # right from the truth 
  print(paste("starting mu (from prior) is ",mu))
  start_x<-c(mu,mu+rnorm(n1,0,1/sqrt(va)))
  set_nested_pois_target(n1=n1,vmu=vmu,va=va,true_mu=true_mu,N=N,blk.ind=blk.ind)
  
  T<-5.5*10^4
  sigma_vec <- c(1,2,3,4,5)
  
  for(sig in sigma_vec){
    output_UB<-grad_MCMC(T=T,sigma=sig,n=n,method="UB",start_x=start_x)
    sample_quant <- quantile(output_UB$x_samples[5001:T,1],0.9)
    quantile_mu_UB <- c(quantile_mu_UB,sample_quant)
  }
  
}

#mse_final_mu_ULA <- t(matrix(mse_mu_ULA,ncol=iter))
final_quantile_mu_UB <- t(matrix(quantile_mu_UB,ncol=iter))


plot(colMeans(final_quantile_mu_UB),xaxt="n",ylab="Sample 90 Percentile",xlab="Step Size",main="",type="b")
axis(1, at=1:5,labels=sigma_vec)


