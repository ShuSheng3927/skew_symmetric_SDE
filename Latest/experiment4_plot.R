######### Experiment 4 ##########

### BACKGROUND

# The goal of this experiment is to investigate the relative performance of EM, Tamed EM, and SS schemes for SDEs with different | mu / sigma | ratios. 
# It is hypothesised that when the ratio is large, SS will underperform due to its skewing update mechanism. When the ratio is moderate or small, SS should perform reasonably well despite having small volatility magnitude. 

### SETUP 

# SDE: dX_t = - X_t dt + alpha X_t dW_t 
# alpha = {0.5, 1, 2}
# X_0 = {0.1, 1, 10}
# stepsize = {0.001, 0.005, 0.01, 0.1, 0.2, 0.3}
# T = 5
# Monte Carlo Iteration = 100,000

# Considered Schemes: Euler-Maruyama, tamed Euler, Skew-Symmetric

# The considered linear SDE admits closed form solution, thus we compare our simulated results with the analytically obtained ground truth. 

# We compare the weak error of the mean and standard deviation across the Monte Carlo iterations for each of the three schemes. 

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

load("./experiment4_data.RData")

p <- ggplot(experiment4_data, aes(x = stepsize, y = MeanError,
                           color = Method, shape = Method, linetype = Method)) +
  geom_line() +
  geom_point(size = 2) +
  scale_x_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  scale_y_log10() +
  facet_grid(rows = vars(alpha), cols = vars(start),
             labeller = label_bquote(rows = .(paste("a =", alpha)),
                                     cols = .(paste("x =", start)))) +
  scale_color_manual(values = c(
    "Euler-Maruyama" = "black",
    "Tamed Euler" = "blue",
    "Skew-Symmetric" = "red"
  )) +
  scale_shape_manual(values = c(
    "Euler-Maruyama" = 16,
    "Tamed Euler" = 17,
    "Skew-Symmetric" = 15
  )) +
  scale_linetype_manual(values = c(
    "Euler-Maruyama" = "solid",
    "Tamed Euler" = "dashed",
    "Skew-Symmetric" = "dotted"
  )) +
  labs(x = "Stepsize", y = "Weak Error", color = NULL,
       shape = NULL, linetype = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    legend.margin = margin(t = -15),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

print(p)

