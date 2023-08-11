l2 <- function(x){
  return(norm(x,type="2"))
}

euler_scheme <- function(mu,sigma,delta,time, start){
  # the Euler-Maruyama scheme for a (one-dim) autonomous SDE
  # mu       drift term of the SDE
  # sigma    volatility term of the SDE
  # delta    step size 
  # time     time of the simulation T = N * delta
  # start    starting point 
  
  step <- floor(time / delta)
  curr <- start
  result <- c(curr)
  for(i in 1:step){
    new <- curr + delta * mu(curr) + sigma(curr) * rnorm(1,mean=0,sd=sqrt(delta))
    result <- c(result, new)
    curr <- new 
  }
  return(result)
}

tamed_euler_scheme <- function(mu,sigma,delta,time, start){
  # the Tamed Euler-Maruyama scheme for a (one-dim) autonomous SDE
  # mu       drift term of the SDE
  # sigma    volatility term of the SDE
  # delta    step size 
  # time     time of the simulation T = N * delta
  # start    starting point 
  
  step <- floor(time / delta)
  curr <- start
  result <- c(curr)
  for(i in 1:step){
    new <- curr + delta * mu(curr) / (1 + delta * l2(mu(curr)) )+ sigma(curr) * rnorm(1,mean=0,sd=sqrt(delta))
    result <- c(result, new)
    curr <- new 
  }
  return(result)
}

barker_scheme <- function(mu, sigma, delta, prob, time, start){
  # the Barker scheme for a (one-dim) autonomous SDE
  # mu       drift term of the SDE
  # sigma    volatility term of the SDE
  # delta    step size 
  # prob     probability function for injection of skewness
  # time     time of the simulation T = N * delta
  # start    starting point 
  
  step <- floor(time / delta)
  curr <- start
  result <- c(curr)
  for(i in 1:step){
    xi <- rnorm(1) * sqrt(delta) * sigma(curr)
    b_prob <- prob(mu=mu,sigma=sigma,delta=delta,curr=curr,xi=xi)
    if (runif(1) < b_prob){
      b <- 1
    }
    else{
      b <- -1 
    }
    new <- curr + b*xi 
    result <- c(result, new)
    curr <- new 
  }
  return(result)
}

prob_cauchy <- function(mu, sigma, delta, curr, xi){
  # the Cauchy CDF probability function for injecting skewness
  # mu       drift term of the SDE
  # sigma    volatility term of the SDE
  # delta    step size 
  # curr     current position 
  if(sigma(curr) == 0){
    exponential <- "inf"
    prob <- 1
  }
  else{
    exp_power <- (2 * xi * mu(curr) / sigma(curr) )/ sigma(curr)
    prob <- 1/(1+exp(-exp_power))
  }
  return(prob)
}



# Ornstein-Uhlenbeck 
# dX_t = theta(mu - X_t) dt + sigma dW_t
ornstein_uhlenbeck_drift <- function(curr){
  mu <- 0
  theta <- 1 
  return(theta*(mu-curr))
}
ornstein_uhlenbeck_volatility <- function(curr){
  sigma <- sqrt(2)
  return(sigma)
}

mu <- 0
theta <- 1 
sigma <- sqrt(2)
time <- 5



# Scheme setup 
delta_vec <- c(exp(-3),exp(-2.5),exp(-2),exp(-1.5),exp(-1))
start <- 1
iter <- 10000


# running the scemes 
yt_barker <- c()
yt_euler <- c()
yt_tamed_euler <- c()

for (i in 1:iter){
  for (delta in delta_vec){
    barker_sample <- barker_scheme(mu = ornstein_uhlenbeck_drift,sigma = ornstein_uhlenbeck_volatility, delta = delta, prob = prob_cauchy, time = time, start = start)
    yt_barker <- c(yt_barker, barker_sample[length(barker_sample)])
    
    euler_sample <- euler_scheme(mu = ornstein_uhlenbeck_drift,sigma = ornstein_uhlenbeck_volatility, delta = delta, time = time, start = start)
    yt_euler <- c(yt_euler, euler_sample[length(euler_sample)])
    
    tamed_euler_sample <- tamed_euler_scheme(mu = ornstein_uhlenbeck_drift,sigma = ornstein_uhlenbeck_volatility, delta = delta, time = time, start = start)
    yt_tamed_euler <- c(yt_tamed_euler, tamed_euler_sample[length(tamed_euler_sample)])
  }
}

final_yt_barker <- t(matrix(yt_barker,ncol=iter))
final_yt_euler <- t(matrix(yt_euler,ncol=iter))
final_yt_tamed_euler <- t(matrix(yt_tamed_euler,ncol=iter))



# compute the errors 
expected_solution <- start*exp(-theta * time) + mu * (1- exp(-theta*time))

f <- function(x){
  return(x^2)
}

log_error_barker <- log(colMeans(f(final_yt_barker)) - f(expected_solution))
log_error_euler <- log(colMeans(f(final_yt_euler)) - f(expected_solution))
log_error_tamed_euler <- log(colMeans(f(final_yt_tamed_euler)) - f(expected_solution))


# generate the plot
plot(delta_vec,log_error_barker,col="red",type="b",ylim=c(min(log_error_barker,log_error_euler,log_error_tamed_euler),max(log_error_barker,log_error_euler,log_error_tamed_euler)),xlab="Step Size",ylab="log f(error)",xaxt="n",main="f(x) = x^2",log='x')
points(delta_vec,log_error_euler,col="black",type="b")
points(delta_vec,log_error_tamed_euler,col="blue",type="b")
axis(1, at=delta_vec,labels=c('exp(-3)','exp(-2.5)','exp(-2)','exp(-1.5)','exp(-1)'))
legend("topleft", legend=c("Euler","Tamed Euler","Barker"), col=c("black","blue","red"), lty=c(1,1), cex=1)






