library(ggplot2)
library(tidyr)
library(dplyr)

# Load and process
load("./soft_sphere.RData")

df_long <- results %>%
  pivot_longer(
    cols = -c(B, stepsize),
    names_to = "method",
    values_to = "unstable_count"
  ) %>%
  mutate(
    unstable_count = unstable_count / 100,
    stepsize = factor(stepsize),
    B = factor(B),
    method = recode(method,
                    "Barker_unstable" = "Skew-Symmetric",
                    "EM_unstable" = "Euler-Maruyama",
                    "impEM_theta02_unstable" = "Implicit Euler",
                    "impEM_theta1_unstable" = "Implicit Euler (theta = 1)")
  ) %>%
  filter(method != "Implicit Euler (theta = 1)")

# Shared fill range
fill_limits <- range(df_long$unstable_count)

# Plot
ggplot(df_long, aes(x = stepsize, y = B, fill = unstable_count)) +
  geom_tile(color = "grey90", linewidth = 0.3) +
  geom_text(aes(label = round(unstable_count, 2)), size = 3, color = "black") +
  facet_wrap(~ method, nrow = 1) +
  scale_fill_gradient(
    low = "#ffffff",
    high = "#f66b6b",
    limits = fill_limits,
    name = "Frequency",
    guide = guide_colorbar(
      barheight = unit(5, "cm")  # <-- increase this as needed
    )
  ) +
  labs(
    x = "Stepsize",
    y = "B"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold"),
    legend.position = "right",
    legend.title = element_text(angle = 0, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
