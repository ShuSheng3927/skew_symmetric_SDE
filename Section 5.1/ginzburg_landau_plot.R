library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

load("./ginzburg_landau_data.RData")

ginzburg_landau_data <- error_total %>%
  pivot_longer(
    cols = c(EM_mean_error, Barker_mean_error, TamedEM_mean_error),
    names_to  = "Method",
    values_to = "MeanError"
  ) %>%
  mutate(
    Method = recode(Method,
                    EM_mean_error      = "Euler-Maruyama",
                    TamedEM_mean_error = "Tamed Euler",
                    Barker_mean_error  = "Skew-Symmetric"
    ),
    Method = factor(Method,
                    levels = c("Euler-Maruyama", "Tamed Euler", "Skew-Symmetric")),
    start = factor(start,
                   levels = c(0.5, 1, 5),
                   labels = c("0.5", "1", "5")
    ),
    alpha = factor(alpha,
                   levels = c(0.5, 2),
                   labels = c("0.5", "2")
    )
  ) %>% 
  select(-EM_sd_error, -Barker_sd_error, -TamedEM_sd_error)

p <- ggplot(ginzburg_landau_data, aes(x = stepsize, y = MeanError,
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
                                     cols = .(paste("x =", start))), scales = 'free_y') +
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
  labs(x = "Step-Size", y = "Weak Error", color = NULL,
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

