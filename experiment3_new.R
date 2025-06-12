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

cauchy_prob <- function(x, mu, sigma, xi,B){
  return(1/(1 + exp(-2*xi*mu(x,B) / sigma(x)^2)))
}



EM <- function(curr, drift, volatility, stepsize, B){
  curr + drift(curr, B) * stepsize + matrix(rnorm(length(curr)), ncol=2) * sqrt(stepsize) * volatility(curr)
}
impEM <- function(curr, drift, volatility, stepsize, B){
  tol = 1e-3
  max_iter = 20
  
  noise = matrix(rnorm(length(curr)), ncol=2)
  fixed_point <- curr
  fixed_point_new <- curr + drift(fixed_point, B) * stepsize + noise * sqrt(stepsize) * volatility(curr)
  counter = 0
  
  while (norm(fixed_point_new - fixed_point, "2") > tol){
    #cat("Position: ", norm(fixed_point_new,"2"), "\n")
    if (counter < max_iter){
      fixed_point <- fixed_point_new
      fixed_point_new <- curr + drift(fixed_point, B) * stepsize +  noise * sqrt(stepsize) * volatility(curr)
      counter = counter + 1
    }else{
      break
    }
  }
  return (fixed_point)
}

Barker <- function(curr, drift, volatility, stepsize,B){
  xi <- matrix(rnorm(length(curr), 0, sd= sqrt(stepsize) * volatility(curr)),ncol=2)
  b <- matrix(c(rbinom(length(curr), 1, cauchy_prob(curr,drift,volatility,xi,B))),ncol=2)
  now <- curr + b*xi
  return(now)
}



stepsize = 0.2
B = 0.5
steps = 20
iter = 10

N = 50

unstability_EM = 0
for (j in 1:iter){
  set.seed(j)
  x_pos <- runif(N,-1,1)
  y_pos <- runif(N,-1,1)
  curr <- cbind(x_pos, y_pos)
  for (i in 1:steps){
    new_pos = EM(curr, soft_sphere_drift, soft_sphere_vol, stepsize, B)
    if (sum(is.infinite(new_pos)) > 0){
      unstability_EM = unstability_EM + 1
      break
    }else{
      curr = new_pos
    }
  }
}
unstability_EM


unstability_Barker = 0
for (j in 1:iter){
  set.seed(j)
  x_pos <- runif(N,-1,1)
  y_pos <- runif(N,-1,1)
  curr <- cbind(x_pos, y_pos)
  for (i in 1:steps){
    new_pos = Barker(curr, soft_sphere_drift, soft_sphere_vol, stepsize, B)
    if (sum(is.infinite(new_pos)) > 0){
      unstability_Barker = unstability_Barker + 1
      break
    }else{
      curr = new_pos
    }
  }
}
unstability_Barker



unstability_impEM = 0
for (j in 1:iter){
  set.seed(j)
  x_pos <- runif(N,-1,1)
  y_pos <- runif(N,-1,1)
  curr <- cbind(x_pos, y_pos)
  for (i in 1:steps){
    new_pos = try(impEM(curr, soft_sphere_drift, soft_sphere_vol, stepsize, B))
    if (inherits(new_pos, "try-error") | sum(is.infinite(new_pos)) > 0){
      unstability_impEM = unstability_impEM + 1
      break
    }else{
      curr = new_pos
    }
  }
}
unstability_impEM


