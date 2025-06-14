# this file contains the R functions used by the other main files to 
# reproduce the simulations reported in the paper

# required packages

require(mvtnorm)

###################################################
####### ISOTROPIC PROPOSALS - NO ADAPTATION #######
###################################################
# this function implements RWM/MALA/Barker with isotropic proposal and no adaptation
grad_MCMC<-function(T=T,sigma=sigma,n=n,method=NULL,start_x){
  stopifnot(method=="MALA"|method=="RWM"|method=="Barker"|method=="ULA"|method=="UB"|method=="ImpULA")
  # specify proposal and acceptance probability
  if(method=="MALA"){
    log_q_ratio<-log_q_ratio_mala
    rprop<-rmala
    ma <- 1
  }
  if(method=="RWM"){
    log_q_ratio<-log_q_ratio_rw
    rprop<-rrw
    ma <- 1
  }
  if(method=="Barker"){
    log_q_ratio<-log_q_ratio_barker
    rprop<-rbarker
    ma <- 1
  }
  if(method=="ULA"){
    log_q_ratio<-log_q_ratio_mala
    rprop<-rmala
    ma <- 0
  }
  if(method=="UB"){
    log_q_ratio<-log_q_ratio_barker
    rprop<-rbarker
    ma <- 0
  }
  if(method=="ImpULA"){
    
    x<-start_x
    t<-1
    x_samples<-matrix(NA,nrow = T,ncol = n)
    x_samples[1,]<-x
    
    # mcmc iterations
    for (t in 2:T){
      # propose new state
      y<-rimp_mala(g_prime,x,sigma)
      x <- y
      # store samples
      x_samples[t,]<-x
    }
    return(list(x_samples=x_samples))
  }
  # initialize sampler and output objects
  x<-start_x
  t<-1
  x_samples<-matrix(NA,nrow = T,ncol = n)
  x_samples[1,]<-x
  # mcmc iterations
  for (t in 2:T){
    
    # propose new state
    y<-x+rprop(g_prime(x),sigma)
    
    # metropolis adjustment
    if(ma == 1){
      # compute acceptance rate
      ap<-min(1,exp(log_f_ratio(x,y)+log_q_ratio(x,y,sigma)))
      if (runif(1)<ap){# accept/reject
        x<-y
      }
    }
    if(ma == 0){
      x <- y
    }
    
    # store samples
    x_samples[t,]<-x
  }
  return(list(x_samples=x_samples))
}
## MALA PROPOSAL AND ACCEPTANCE PROBABILITY (ISOTROPIC)
log_q_ratio_mala<-function(x,y,sigma){return(sum(
  -(x-y-sigma^2*g_prime(y)/2)^2/(2*sigma^2)+
    (y-x-sigma^2*g_prime(x)/2)^2/(2*sigma^2)
))}
rmala <-function(c.,sigma.){
  return(
    sigma.*rnorm(n=length(c.))+sigma.^2*c./2
  )}
rimp_mala <- function(g_prime, start_x,sigma){
  noise = rnorm(n=length(start_x))
  tol = 1e-3
  max_iter = 500
  theta = 0.2
  
  fixed_point_x = start_x
  new_x_option = start_x + (1-theta) * g_prime(start_x) * sigma**2 / 2 + theta * g_prime(fixed_point_x) * sigma**2 / 2 + sigma * noise
  counter = 0
  while(norm(new_x_option - fixed_point_x,"2") > tol){
      if (counter < max_iter){
        fixed_point_x = new_x_option
        new_x_option =  start_x + (1-theta) * g_prime(start_x) * sigma**2 / 2 + theta * g_prime(fixed_point_x) * sigma**2 / 2 + sigma * noise
        counter = counter + 1
      }else{
        break
    }
  }
  return(new_x_option)
}

## RANDOM WALK PROPOSAL AND ACCEPTANCE PROBABILITY (ISOTROPIC)
log_q_ratio_rw<-function(x,y,sigma){return(0)}
rrw<-function(c,sigma){return(sigma*rnorm(n=length(c)))}
## BARKER PROPOSAL AND ACCEPTANCE PROBABILITY (ISOTROPIC)
log_q_ratio_barker<-function(x,y,sigma){
  beta1<-  c(-g_prime(y)*(x-y))
  beta2<-  c(-g_prime(x)*(y-x))
  return(sum(# compute acceptance with log_sum_exp trick for numerical stability
    -(pmax(beta1,0)+log1p(exp(-abs(beta1))))+
      (pmax(beta2,0)+log1p(exp(-abs(beta2))))
  ))
}
rbarker <-function(c,sigma){
  z<-sigma*rnorm(n=length(c)) # rnorm could be replaced with any symm kernel
  b<-2*(runif(n=length(c))< 1/(1+exp(-c*z)))-1
  return(z*b)
}
###################################################################
### DEFINE THE 4 SCENARIOS USED IN HETEROGENEOUS TARGETS SECTION  #####
###################################################################
set_target<-function(scenario){
  print(paste("Set target as in scenario ",scenario))
  if(scenario==1){#Gaussian with one small scale
    n<<-100 # number of dimensions
    Sigma_targ<<-diag(rep(1,n))
    Sigma_targ[1,1]<<-0.01^2
    Prec_target<<-solve(Sigma_targ)
    log_f_ratio<<-function(x,y){return(-0.5*(
      matrix(y,nrow = 1,ncol = n)%*%Prec_target%*%matrix(y,nrow = n,ncol = 1)-
        matrix(x,nrow = 1,ncol = n)%*%Prec_target%*%matrix(x,nrow = n,ncol = 1)
    ))}
    # g_prime computes the gradient of the log_target
    g_prime<<-function(x){return(-c(matrix(x,nrow = 1,ncol = n)%*%Prec_target))}
  }
  if(scenario==2){#Gaussian with random scales
    n<<-100
    scales<<-exp(rnorm(n,mean = 0,sd = 1))
    variances<<-scales^2
    Sigma_targ<<-diag(variances)
    Prec_target<<-solve(Sigma_targ)
    log_f_ratio<<-function(x,y){return(-0.5*(
      matrix(y,nrow = 1,ncol = n)%*%Prec_target%*%matrix(y,nrow = n,ncol = 1)-
        matrix(x,nrow = 1,ncol = n)%*%Prec_target%*%matrix(x,nrow = n,ncol = 1)
    ))}
    g_prime<<-function(x){return(-c(matrix(x,nrow = 1,ncol = n)%*%Prec_target))}
  }
  if(scenario==3){#Hyperbolic with random scales
    n<<-100
    delta2<<-0.1
    scales<<-exp(rnorm(n,mean = 0,sd = 1))
    log_f_ratio<<-function(x,y){return(-sum(
      (sqrt(delta2+(y/scales)^2)-sqrt(delta2+(x/scales)^2))
    ))}
    g_prime<<-function(x){
      x_resc<- x/scales
      grad<- -x_resc/sqrt(delta2+x_resc^2)
      return(grad/scales)
    }
  }
  if(scenario==4){#Skew-normal with random scales
    n<<-100
    scales<<-exp(rnorm(n,mean = 0,sd = 1))
    alpha<<-4
    log_f_ratio<<-function(x,y){
      return(sum(
        dnorm(y/scales,0,1, log = TRUE)+pnorm(alpha*y/scales,0,1, log.p = TRUE)-
          (dnorm(x/scales,0,1, log = TRUE)+pnorm(alpha*x/scales,0,1, log.p = TRUE))
      ))
    }
    g_prime<<-function(x){
      x_resc<- x/scales
      grad<- -x_resc+
        alpha*exp(dnorm(alpha*x_resc,0,1, log = TRUE)-pnorm(alpha*x_resc,0,1, log.p = TRUE))
      return(grad/scales)
    }
  }
}

###################################################################
### DEFINE THE HIERARCHICAL POISSON REGRESSION TARGET  #####
###################################################################
set_nested_pois_target<-function(n1=10,
                                 vmu=1/(10^2), # prior precision for global means
                                 va=1/(2^2), # prior precision for random effects
                                 true_mu=5,
                                 N=NULL, # num of obs
                                 blk.ind=NULL
){
  stopifnot(length(blk.ind)==N)
  n<<-1+n1# number of parameters to sample
  par_names<<-rep(NA,n)
  par_names[1]<<-c("mu")
  ind1<<-1+1:n1
  par_names[ind1]<<-c("eta")
  ## generate data
  true_eta<<-true_mu+rnorm(n1,mean = 0,sd = 1/sqrt(va))
  blk.list<<-lapply(X = c(1:n1),FUN = function(i){which(blk.ind==i)})
  y.vec<<-c(rpois(n = N,lambda = exp(true_eta[blk.ind])))
  # target and gradient
  log_f_ratio<<-function(x,y){
    return(
      log_f(y)-log_f(x)
    )}
  log_f<<-function(x){
    mu<-x[1]
    eta<-x[ind1]
    log_lambda<-eta[blk.ind]
    return(
      -mu^2*vmu/2-sum((eta-mu)^2)*va/2+  #log prior
        sum(-exp(log_lambda)+y.vec*(log_lambda)) #log likelihoods
    )
  }
  g_prime<<-function(x,sigma_0=sigma0){
    mu<-x[1]
    eta<-x[ind1]
    lambda.vec<-exp(eta[blk.ind])
    ly_diff<-lambda.vec-y.vec
    return(c(-mu*vmu-sum(mu-eta)*va,#grad_mu
             -(eta-mu)*va-vapply(blk.list,FUN = function(ii){sum(ly_diff[ii])},FUN.VALUE = 1)#grad_eta
    ))
  }
}