######### Section 5.1-(b) Stochastic Ginzburg-Landau ##########

### SETUP 

# SDE: dX_t = [(eta-alpha^2/2)X_t - lambda X_t^3] dt + alpha X_t dW_t 
# alpha = {0.5, 2}
# X_0 = {0.1, 1, 10}
# stepsize = {0.001, 0.005, 0.010, 0.020, 0.050}
# T = 5
# Monte Carlo Iteration = 100,000

# Considered Schemes: Euler-Maruyama, tamed Euler, Skew-Symmetric

# The considered SDE admits closed form solution, but is hard to be used to calculate its mean. Instead, we simulate Brownian motion paths which are then used to estimate X_T, so one could use Monte Carlo to estimate the mean. Afterwards, we compare the weak error of the mean across the Monte Carlo iterations for each of the three schemes. 

library(parallel)

termT = 5
dt = 0.0001
tt = seq(0,termT, dt)
eta = 0
lambda = 1

# Parameter grid
alphas <- c(0.5, 2)
starts <- c(0.1, 1, 10)

error_total = c() 

# Loop over alpha and start
for (alpha in alphas) {
  for (start in starts) {
    
    # Define model functions for current alpha
    drift <- function(x) (eta+alpha^2/2)*x - lambda*x^3
    volatility <- function(x) alpha*x
    
    EM <- function(x, drift, volatility, stepsize) {
      x + drift(x) * stepsize + volatility(x) * rnorm(1) * sqrt(stepsize)
    }
    
    tamed_EM <- function(x, drift, volatility, stepsize) {
      x + drift(x) / (1 + abs(drift(x))) * stepsize + volatility(x) * rnorm(1) * sqrt(stepsize)
    }
    
    Barker <- function(x, drift, volatility, stepsize) {
      noise <- rnorm(1)
      flip_prob <- plogis(sqrt(stepsize) * noise * 2 / volatility(x) * drift(x))
      flip <- 2 * rbinom(1, 1, flip_prob) - 1
      x + sqrt(stepsize) * noise * flip * volatility(x)
    }
    
    # Simulation function
    simulate_once <- function(args) {
      seed <- args[1]
      time <- 5
      stepsize <- args[2]
      iter <- time / stepsize
      set.seed(seed)
      
      curr_EM <- curr_Barker <- curr_tamed_EM <- start
      
      for (i in 1:iter) {
        curr_EM <- EM(curr_EM, drift, volatility, stepsize)
        curr_Barker <- Barker(curr_Barker, drift, volatility, stepsize)
        curr_tamed_EM <- tamed_EM(curr_tamed_EM, drift, volatility, stepsize)
      }
      
      list(
        seed = seed,
        stepsize = stepsize,
        EM = curr_EM,
        Barker = curr_Barker,
        TamedEM = curr_tamed_EM
      )
    }
    
    # Build argument list
    argument_values <- list()
    for (seed in 1:100000) { 
      for (Delta in c(0.001, 0.005, 0.010, 0.020, 0.050)) {
        argument_values <- append(argument_values, list(c(seed, Delta)))
      }
    }
    
    # Run in parallel
    results_total <- mclapply(argument_values, simulate_once, mc.cores = detectCores(logical = TRUE)-1)
    
    # Convert results
    results_total <- do.call(rbind, lapply(results_total, function(x) {
      data.frame(
        seed = x$seed,
        stepsize = x$stepsize,
        EM = x$EM,
        Barker = x$Barker,
        TamedEM = x$TamedEM
      )
    }))
    
    # Compute theoretical moments
    iter = 100000
    one_xt <- function(i) {
      set.seed(i)  
      
      ww_gap <- c(0, rnorm(termT / dt, 0, sqrt(dt)))
      ww     <- cumsum(ww_gap)
      
      integrand <- exp(2 * eta * tt + 2 * alpha * ww)
      
      xt <- start * exp(eta * termT + alpha * ww[length(ww)]) /
        sqrt(
          1 + 2 * start^2 * lambda * (
            sum(integrand) * dt -
              (integrand[1] + integrand[length(integrand)]) * dt / 2
          )
        )
      
      return(xt)
    }
    n_cores <- detectCores(logical = TRUE)-1 
    xt_list <- unlist(
      mclapply(
        X         = 1:iter,
        FUN       = one_xt,
        mc.cores  = n_cores-1
      )
    )
    
    true_mean <- mean(xt_list)
    true_std <- sd(xt_list)
    
    stepsizes <- sort(unique(results_total$stepsize))
    
    # Initialize storage
    error_summary <- data.frame(
      stepsize = stepsizes,
      start = start,
      alpha = alpha,
      EM_mean_error = NA_real_,
      Barker_mean_error = NA_real_,
      TamedEM_mean_error = NA_real_,
      EM_sd_error = NA_real_,
      Barker_sd_error = NA_real_,
      TamedEM_sd_error = NA_real_
    )
    
    # Compute error statistics
    for (i in seq_along(stepsizes)) {
      s <- stepsizes[i]
      subset_data <- results_total[results_total$stepsize == s, ]
      
      error_summary$EM_mean_error[i]      <- abs(mean(subset_data$EM) - true_mean)
      error_summary$Barker_mean_error[i]  <- abs(mean(subset_data$Barker) - true_mean)
      error_summary$TamedEM_mean_error[i] <- abs(mean(subset_data$TamedEM) - true_mean)
      
      error_summary$EM_sd_error[i]        <- abs(sd(subset_data$EM) - true_std)
      error_summary$Barker_sd_error[i]    <- abs(sd(subset_data$Barker) - true_std)
      error_summary$TamedEM_sd_error[i]   <- abs(sd(subset_data$TamedEM) - true_std)
    }
    
    # Save result for this parameter combination
    error_total = rbind(error_total, error_summary)
  }
}

save(error_total, file = "ginzburg_landau_data.RData")
