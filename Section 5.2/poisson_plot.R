library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(patchwork)

load("./poisson.RData")

stepsizes <- c(0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05)

# Define appearance
method_colors <- c(
  "Euler-Maruyama" = "black",
  "Skew-Symmetric" = "red",
  "Implicit Euler" = "purple"
)
method_shapes <- c(
  "Euler-Maruyama" = 16,
  "Skew-Symmetric" = 15,
  "Implicit Euler" = 17
)
method_linetypes <- c(
  "Euler-Maruyama" = "solid",
  "Skew-Symmetric" = "dotted",
  "Implicit Euler" = "dashed"
)

# Log10 breaks
log_breaks <- 10^round(seq(log10(0.005), log10(0.05), by = 0.2), 2)

# ===== True Start plot =====
mean_ULA <- colMeans(sim1_true_start_ULA, na.rm = FALSE)
mean_UB <- colMeans(sim1_true_start_UB, na.rm = FALSE)
mean_impULA <- colMeans(sim1_true_start_impULA, na.rm = FALSE)

df_equilibrium <- data.frame(
  Stepsize = rep(stepsizes, times = 3),
  MSE = c(mean_ULA, mean_UB, mean_impULA),
  Method = rep(c("Euler-Maruyama", "Skew-Symmetric", "Implicit Euler"), each = length(stepsizes)),
  Start = "True Start"
)

# ===== Warm Start plot =====
mean_ULA <- colMeans(sim1_warm_start_ULA, na.rm = FALSE)
mean_UB <- colMeans(sim1_warm_start_UB, na.rm = FALSE)
mean_impULA <- colMeans(sim1_warm_start_impULA, na.rm = FALSE)

df_warm <- data.frame(
  Stepsize = rep(stepsizes, times = 3),
  MSE = c(mean_ULA, mean_UB, mean_impULA),
  Method = rep(c("Euler-Maruyama", "Skew-Symmetric", "Implicit Euler"), each = length(stepsizes)),
  Start = "Warm Start"
)

# Combine data
df_all <- bind_rows(df_equilibrium, df_warm)

# ===== Plot =====
p <- ggplot(df_all, aes(x = Stepsize, y = MSE, color = Method, shape = Method, linetype = Method)) +
  geom_line() +
  geom_point(size = 2) +
  scale_x_log10(
    breaks = log_breaks,
    labels = trans_format("log10", math_format(10^.x))
  ) +
  scale_y_log10() +
  scale_color_manual(values = method_colors) +
  scale_shape_manual(values = method_shapes) +
  scale_linetype_manual(values = method_linetypes) +
  labs(
    x = "Stepsize",
    y = "Mean Squared Error",
    color = NULL,
    shape = NULL,
    linetype = NULL
  ) +
  facet_wrap(~Start, nrow = 1, labeller = label_value) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    legend.margin = margin(t = -15),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

print(p)
