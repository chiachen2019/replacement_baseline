library(tidyverse)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(patchwork)

initial_blocking_results <- read.csv("../output/output_initial_blocking_results_2026-10-02.csv")

partial_vector_control_scenario_results <- read.csv("../output/output_partial_2026-10-02.csv")

partial_vector_control_scenario_results |>
  filter(
    blocking == 0.5,
#    time > 365 * 8 -5,
    time < 5,
    K0_red == 0.5, 
    K0_red_post_rel == 0, 
    blocking == 0.5
  ) |> select(time, A_w, K0_red, K0_red_post_rel)

# ============================================================
# 1. Baseline source-reduction scenarios
# ============================================================

sr_data <- partial_vector_control_scenario_results |>
  filter(
    blocking == 0.5,
    time > 365 * 8,
    time < 365 * 8.5
  ) |>
  mutate(
    time = time - 365 * 8,
    prop = A_w / (A_w + A_u),
    K0_red = factor(K0_red),
    
    retention = factor(
      K0_red_post_rel,
      levels = c(0, 0.25, 0.5, 0.75),
      labels = c(
        "0% of baseline",
        "25% of baseline",
        "50% of baseline",
        "75% of baseline"
      )
    )
  )


# Au
p1 <- ggplot(
  sr_data,
  aes(time, A_u, colour = K0_red)
) +
  geom_line(linewidth = 0.6) +
  facet_wrap(
    ~ retention,
    ncol = 2
  ) +
  labs(
    title = "(A) Uninfected adults",
    x = NULL,
    y = expression(A[u]),
    colour = "Baseline source\nreduction"
  ) +
  theme_classic(base_size = 12)


# Aw
p2 <- ggplot(
  sr_data,
  aes(time, A_w, colour = K0_red)
) +
  facet_wrap(
    ~ retention,
    ncol = 2
  ) +
  geom_line(linewidth = 0.6) +
  labs(
    title = expression("(B) " * italic(Wolbachia) * "-infected adults"), 
    x = NULL,
    y = expression(A[w]),
    colour = "Baseline source\nreduction"
  ) +
  theme_classic(base_size = 12)


# Wolbachia proportion
p3 <- ggplot(
  sr_data,
  aes(time, prop, colour = K0_red)
) +
  geom_line(linewidth = 0.6) +
  geom_hline(
    yintercept = 0.99,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  facet_wrap(
    ~ retention,
    ncol = 2
  ) +
  scale_y_continuous(
    labels = scales::percent_format(),
    limits = c(0, 1)
  ) +
  labs(
    title = expression("(C) " * italic(Wolbachia) * "-infected frequency"), 
    x = NULL,
    y = expression(
      italic(Wolbachia)~"frequency"
    ),
    colour = "Baseline source\nreduction"
  ) +
  theme_classic(base_size = 12)


# ============================================================
# 2. Initial-population scenarios
# ============================================================

initial_data <- initial_blocking_results |>
  filter(
    blocking == 0.5,
    time < 365 * 0.5,
    initial_proportion < 1
  ) |>
  mutate(
    prop = A_w / (A_w + A_u),
    initial_proportion = factor(initial_proportion)
  )


# Au
p4 <- ggplot(
  initial_data,
  aes(
    time,
    A_u,
    colour = initial_proportion
  )
) +
  geom_line(linewidth = 0.6) +
  labs(
    title = "(D) Uninfected adults",
    x = "Time (days)",
    y = expression(A[u]),
    colour = "Initial\nproportion"
  ) +
  theme_classic(base_size = 12)


# Aw
p5 <- ggplot(
  initial_data,
  aes(
    time,
    A_w,
    colour = initial_proportion
  )
) +
  geom_line(linewidth = 0.6) +
  labs(
    title = expression("(E) " * italic(Wolbachia) * "-infected adults"), 
    x = "Time (days)",
    y = expression(A[w]),
    colour = "Initial\nproportion"
  ) +
  theme_classic(base_size = 12)


# Wolbachia proportion
p6 <- ggplot(
  initial_data,
  aes(
    time,
    prop,
    colour = initial_proportion
  )
) +
  geom_line(linewidth = 0.6) +
  geom_hline(
    yintercept = 0.99,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  scale_y_continuous(
    labels = scales::percent_format(),
    limits = c(0, 1)
  ) +
  labs(
    title = expression("(F) " * italic(Wolbachia) * "-infected frequency"), 
    x = "Time (days)",
    y = expression(
      italic(Wolbachia)~"frequency"
    ),
    colour = "Initial\nproportion"
  ) +
  theme_classic(base_size = 12)


# ============================================================
# 3. Combine into 2 x 3 figure
# ============================================================

(p1 + p2 + p3) /
  (p4 + p5 + p6) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    plot.title = element_text(
      size = 12
    ),
    axis.title = element_text(
      face = "bold"
    )
  )

ggsave("est_time_series.pdf", width = 15, height = 8)
