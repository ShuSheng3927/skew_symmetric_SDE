# Experiment 4 

### BACKGROUND

# The goal of this experiment is to investigate the relative performance of EM, Tamed EM, and SS schemes for SDEs with different | mu / sigma | ratios. 
# It is hypothesised that when the ratio is large, SS will underperform due to its skewing update mechanism. When the ratio is moderate or small, SS should perform reasonably well despite having small volatility magnitude. 

### SETUP 

# SDE: dX_t = -X_t dt + alpha X_t dW_t # alpha = {0.1, 1, 10}
# X_0 = {0.1, 1, 10}
# stepsize = {0.001, 0.005, 0.01, 0.1, 0.2}
# T = 5
# Monte Carlo Iteration = 100,000

# Considered Schemes: Euler-Maruyama 

# The considered linear SDE admits closed form solution, thus we compare our simulated results with the analytically obtained ground truth. 

# We compare the weak error of the mean and standard deviation across the Monte Carlo iterations for each of the three schemes. 


load("./montecarlo_results_aggregated_full.RData")


error_summary = error_total[error_total$start == 1 & error_total$alpha == 1, -2]
error_summary = error_summary[-2]

stepsize = error_summary$stepsize


par(mfrow=c(1,2))
ymin = min(error_summary[,2:4])
ymax = max(error_summary[,2:4])

plot(stepsize, error_summary$EM_mean_error,log='xy',type="b",pch=1,xlab="Stepsize", ylab="Weak Error",ylim=c(ymin, ymax),main="Weak Error of Mean")
points(stepsize, error_summary$Barker_mean_error, type="b",pch=2, col="red")
points(stepsize, error_summary$TamedEM_mean_error, type="b", pch=3, col="blue")
legend("center", legend=c("EM", "SS", "TamedEM"), pch=c(1,2,3), col=c("black", "red", "blue"))


ymin = min(error_summary[,5:7])
ymax = max(error_summary[,5:7])

plot(stepsize, error_summary$EM_sd_error,log='xy',type="b",pch=1,xlab="Stepsize", ylab="Weak Error",ylim=c(ymin, ymax),main="Weak Error of Standard Deviation")
points(stepsize, error_summary$Barker_sd_error, type="b",pch=2, col="red")
points(stepsize, error_summary$TamedEM_sd_error, type="b", pch=3, col="blue")
legend("center", legend=c("EM", "SS", "TamedEM"), pch=c(1,2,3), col=c("black", "red", "blue"))