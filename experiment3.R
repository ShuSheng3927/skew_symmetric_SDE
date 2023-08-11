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

# Library for plots
library(ggplot2)
library(hrbrthemes)



# general simulation setup
B_vec <- rep(seq(0.1,1,0.1),10)
stepsize_vec <- c(rep(0.1,10),rep(0.2,10),rep(0.3,10),rep(0.4,10),rep(0.5,10),rep(0.6,10),rep(0.7,10),rep(0.8,10),rep(0.9,10),rep(1,10))
iter <- 100
steps <- 10
N <- 50

# simulation for euler 
blowup_euler <- data.frame(stepsize_vec, B_vec, rep(0,100))
colnames(blowup_euler) <- c("stepsize", "B", "blowup")

for(i in seq(0.1,1,0.1)){
  for(j in seq(0.1,1,0.1)){
    stepsize <- i
    B <- j
    for (l in 1:iter){
      x_pos <- runif(N,-1,1)
      y_pos <- runif(N,-1,1)
      results_euler <- data.frame(1:N,x_pos,y_pos,rep(0,N))
      colnames(results_euler) <- c("id","x","y","time") 
      t <- 0
      for (k in 1:steps){
        if(sum(c(is.infinite(results_euler$x),is.infinite(results_euler$y))) == 0){
          curr <- as.matrix(results_euler[1:N+N*(k-1),2:3])
          now <- curr + stepsize*soft_sphere_drift(curr,B) + soft_sphere_vol(curr)*matrix(rnorm(2*N,0,sqrt(stepsize)),ncol=2)
          t <- t + stepsize 
          new <- cbind(1:N,now,rep(t,N))
          colnames(new) <- c("id","x","y","time") 
          rownames(new) <- 1:N+N*k
          results_euler <- rbind(results_euler,new)
        }
        else{
          index <- i / 0.1 * 10 - 10 + j / 0.1
          if(index > 50 & index < 60){
            index <- index+1 # to bypass an issue of R
          }
          blowup_euler$blowup[index] <- blowup_euler$blowup[index] + 1
          break
        }
      }
    }
  }
}
blowup_euler$blowup <- blowup_euler$blowup/iter

# generate the heatmap for euler
ggplot(blowup_euler, aes(B, stepsize, fill= blowup)) + 
  geom_tile(color = "white",lwd = 0.2,linetype = 1) +
  theme_ipsum() + 
  geom_text(aes(label = blowup), color = "white", size = 4) +
  xlab("B") + 
  ylab("Stepsize") + 
  labs(fill="Blowup Prob.") + 
  ggtitle("Euler Maruyama") + 
  scale_y_continuous(breaks = seq(0.1,1,0.1)) +
  scale_x_continuous(breaks = seq(0.1,1,0.1)) +
  coord_fixed() + 
  guides(fill = guide_colourbar(barwidth = 0.5,barheight = 15))



# simualtion for barker
blowup_barker <- data.frame(stepsize_vec, B_vec, rep(0,100))
colnames(blowup_barker) <- c("stepsize", "B", "blowup")

for(i in seq(0.1,1,0.1)){
  for(j in seq(0.1,1,0.1)){
    stepsize <- i
    B <- j
    for (l in 1:iter){
      x_pos <- runif(N,-1,1)
      y_pos <- runif(N,-1,1)
      results_barker <- data.frame(1:N,x_pos,y_pos,rep(0,N))
      colnames(results_barker) <- c("id","x","y","time") 
      time <- 0
      for (k in 1:steps){
        if(sum(c(is.infinite(results_barker$x),is.infinite(results_barker$y))) == 0){
          curr <- as.matrix(results_barker[1:N+N*(k-1),2:3])
          xi <- matrix(rnorm(2*N, 0, sd= sqrt(stepsize) * soft_sphere_vol(curr)),ncol=2)
          b <- matrix(c(rbinom(2*N, 1, cauchy_prob(curr,soft_sphere_drift,soft_sphere_vol,xi,B))),ncol=2)
          now <- curr + b*xi
          time <- time + stepsize 
          new <- cbind(1:N,now,rep(time,N))
          colnames(new) <- c("id","x","y","time") 
          rownames(new) <- 1:N+N*k
          results_barker <- rbind(results_barker,new)
        }
        else{
          index <- i / 0.1 * 10 - 10 + j / 0.1
          if(index > 50 & index < 60){
            index <- index+1 # to bypass an issue of R
          }
          blowup_barker$blowup[index] <- blowup_barker$blowup[index] + 1
          break
        }
      }
    }
  }
}
blowup_barker$blowup <- blowup_barker$blowup/iter

# generate the heatmap for barker
ggplot(blowup_barker, aes(B, stepsize, fill= blowup)) + 
  geom_tile(color = "white",lwd = 0.2,linetype = 1) +
  theme_ipsum() + 
  geom_text(aes(label = blowup), color = "white", size = 4) +
  xlab("B") + 
  ylab("Stepsize") + 
  labs(fill="Blowup Prob.") + 
  ggtitle("Unadjusted Barker") + 
  scale_y_continuous(breaks = seq(0.1,1,0.1)) +
  scale_x_continuous(breaks = seq(0.1,1,0.1)) +
  coord_fixed() + 
  guides(fill = guide_colourbar(barwidth = 0.5,barheight = 15))






