library(parallel)

# Parameter grid
alphas <- c(0.1, 1, 10)
starts <- c(0.1, 1, 10)

error_total = c() 

# Loop over alpha and start
for (alpha in alphas) {
  for (start in starts) {
    
    # Define model functions for current alpha
    drift <- function(x) -x
    volatility <- function(x) alpha * x
    
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
      for (Delta in c(0.001, 0.005, 0.01, 0.1, 0.2)) {
        argument_values <- append(argument_values, list(c(seed, Delta)))
      }
    }
    
    # Run in parallel
    results_total <- mclapply(argument_values, simulate_once, mc.cores = 10)
    
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
    true_mean <- start * exp(-5)
    true_std <- start^2 * exp(-2 * 5) * (exp(5 * alpha^2) - 1)
    
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

save(error_total, file = "montecarlo_results_aggregated.RData")
