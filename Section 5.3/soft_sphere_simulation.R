######### Section 5.3 Soft Spheres ##########

### SETUP 

# The particle dynamics are governed by the stochastic differential equation
# dY_t^{(i)} = 4B(\beta - Y_t^{(i)}) \|Y_t^{(i)} - \beta \|^2 \, dt + \frac{A}{Nr^2} \sum_{j=1}^N (Y_t^{(i)} - Y_t^{(j)})e^{-\| Y_t^{(i)} - Y_t^{(j)} \|^2/2r^2} \, dt + \sqrt{2D}\,dW_t,
# for i = 1, 2, ..., N where \beta is the position of the trap, A > 0 is the strength of the repulsion between spheres with radius r, and B > 0 is the strength of the trap. We set \beta = (0,0)^T, D = 0.25, A = 30, r = 0.15 and N = 50.

# B = 0.1, 0.2, ..., 1
# stpesize = 0.1, 0.2, ..., 1
# step = 10
# iteration = 100

# Considered Schemes: Euler-Maruyama, semi-implicit Euler, Skew-Symmetric

# Count the frequency of numerical unstability


l2 <- function(a,b){
  return(norm(a-b,type="2"))
}

soft_sphere_drift <- function(x_t,B){
  # drift for soft spheres in a 2d anharmonic trap model
  # A > 0   strength of the repulsion between spheres
  # B > 0   strength of the trap
  # r       size of the sphere
  # beta_t  a moving trap
  # x_t     position of particles, N*dim matrix
  
  A <- 10
  r <- 0.5
  N <- nrow(x_t)
  beta_t <- c(0,0)
  
  output <- c()
  for (i in 1:N){
    soft_sphere <- 0
    for(j in 1:N){
      soft_sphere <- soft_sphere + (x_t[i,]-x_t[j,])*exp(-l2(x_t[i,],x_t[j,])^2/2/r^2)
    }
    drift_i <- 4*B*(beta_t - x_t[i,])*l2(beta_t,x_t[i,])^2 + A/N/r^2*soft_sphere
    output <- c(output, drift_i)
  }
  return(t(matrix(output,nrow=2)))
}

soft_sphere_vol <- function(x_t){
  D <- 0.25
  return(sqrt(2*D))
}



EM <- function(curr, drift, volatility, stepsize, B){
  curr + drift(curr, B) * stepsize + matrix(rnorm(length(curr)), ncol=2) * sqrt(stepsize) * volatility(curr)
}

impEM <- function(curr, drift, volatility, stepsize, B,theta = 0.2){
  tol = 1e-3
  max_iter = 500
  
  noise = matrix(rnorm(length(curr)), ncol=2)
  fixed_point <- curr
  fixed_point_new <- curr + (1-theta) * drift(curr, B)*stepsize + theta*drift(fixed_point, B) * stepsize + noise * sqrt(stepsize) * volatility(curr)
  counter = 0
  
  while (norm(fixed_point_new - fixed_point, "2") > tol){
    if (counter < max_iter){
      fixed_point <- fixed_point_new
      fixed_point_new <- curr + (1-theta) * drift(curr, B)*stepsize + theta*drift(fixed_point, B) * stepsize +  noise * sqrt(stepsize) * volatility(curr)
      counter = counter + 1
    }else{
      break
    }
  }
  return (fixed_point)
}

cauchy_prob <- function(x, mu, sigma, xi,B){
  return(1/(1 + exp(-2*xi*mu(x,B) / sigma(x)^2)))
}
Barker <- function(curr, drift, volatility, stepsize,B){
  xi <- matrix(rnorm(length(curr), 0, sd= sqrt(stepsize) * volatility(curr)),ncol=2)
  b <- matrix(c(rbinom(length(curr), 1, cauchy_prob(curr,drift,volatility,xi,B))),ncol=2)
  now <- curr + b*xi
  return(now)
}



B_vals <- seq(0.1, 1, by = 0.1)
stepsize_vals <- seq(0.1, 1, by = 0.1)

results <- data.frame()

steps = 10
iter = 100

N = 50
run_simulation <- function(method_fn, curr, stepsize, B, steps, method_name, theta = NULL) {
  for (i in 1:steps) {
    if (method_name == "impEM") {
      new_pos <- try(impEM(curr, soft_sphere_drift, soft_sphere_vol, stepsize, B, theta), silent = TRUE)
      if (inherits(new_pos, "try-error") || any(is.infinite(new_pos))) return(TRUE)
    } else {
      new_pos <- method_fn(curr, soft_sphere_drift, soft_sphere_vol, stepsize, B)
      if (any(is.infinite(new_pos))) return(TRUE)
    }
    curr <- new_pos
  }
  return(FALSE)
}

for (B in B_vals) {
  for (stepsize in stepsize_vals) {
    cat("Running B =", B, "stepsize =", stepsize, "\n")
    
    count_EM <- 0
    count_Barker <- 0
    count_impEM_theta02 <- 0
    count_impEM_theta1 <- 0
    
    for (j in 1:iter) {
      set.seed(j)
      curr <- cbind(runif(N, -1, 1), runif(N, -1, 1))
      if (run_simulation(EM, curr, stepsize, B, steps, "EM")) {
        count_EM <- count_EM + 1
      }
      
      curr <- cbind(runif(N, -1, 1), runif(N, -1, 1))
      if (run_simulation(Barker, curr, stepsize, B, steps, "Barker")) {
        count_Barker <- count_Barker + 1
      }
      
      curr <- cbind(runif(N, -1, 1), runif(N, -1, 1))
      if (run_simulation(NULL, curr, stepsize, B, steps, "impEM", theta = 0.2)) {
        count_impEM_theta02 <- count_impEM_theta02 + 1
      }
      
      curr <- cbind(runif(N, -1, 1), runif(N, -1, 1))
      if (run_simulation(NULL, curr, stepsize, B, steps, "impEM", theta = 1)) {
        count_impEM_theta1 <- count_impEM_theta1 + 1
      }
    }
    
    # Store results
    results <- rbind(results, data.frame(
      B = B,
      stepsize = stepsize,
      EM_unstable = count_EM,
      impEM_theta02_unstable = count_impEM_theta02,
      impEM_theta1_unstable = count_impEM_theta1,
      Barker_unstable = count_Barker
    ))
  }
}

save(results, "soft_sphere.RData")