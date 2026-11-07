library(tidyverse)
library(ggplot2)
library(ggpubr)
library(cowplot)


initial_blocking_results <- read.csv("../output/output_initial_blocking_results_2026-10-02.csv")
vector_control_only_results <- read.csv("../output/vector_control_only_results_2026-09-28.csv")
no_intervention <- read.csv("../output/no_release_results_2026-09-28.csv")

source("plots_func.R")

## time to reach equilbrium
no_intervention_equilibrium <- no_intervention |>
  filter(
    time < 365,
    A_u > 0.99 * equilibrium_A_u
  ) |>
  group_by(initial_proportion) |>
  slice_min(time, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(initial_proportion, time) |>
  rename(time_no_intervention = time)


wolbachia_equilibrium <- initial_blocking_results |>
  filter(
    time < 365,
    A_u + A_w > 0.99 * equilibrium_A_u
  ) |>
  group_by(initial_proportion) |>
  slice_min(time, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(initial_proportion, time) |>
  rename(time_wolbachia = time)


equilibrium_time <- no_intervention_equilibrium |>
  full_join(
    wolbachia_equilibrium,
    by = "initial_proportion"
  ) |> mutate(perc_shortened = 
                (time_wolbachia-time_no_intervention)/time_no_intervention)

write.csv(equilibrium_time, 
          file = paste0("../output/equilbrium_time_", 
                                          Sys.Date(), ".csv"))

##time series 

initial_prop_plot <- c(0.1, 0.3, 0.7,  0.9)
blocking_plot     <- c(0.7, 0.9)
source_red_plot   <- c(0.5, 0.9)

VC_total_all <- bind_rows(
  
  # No intervention
  no_intervention |>
    filter(initial_proportion %in% initial_prop_plot) |>
    transmute(
      initial_proportion,
      time,
      VC_total,
      scenario = "No intervention"
    ),
  
  # Wolbachia
  initial_blocking_results |>
    filter(
      initial_proportion %in% initial_prop_plot,
      blocking %in% blocking_plot
    ) |>
    transmute(
      initial_proportion,
      time,
      VC_total,
      scenario = paste0("Blocking effect ", blocking)
    ),
  
  # Source reduction
  vector_control_only_results |>
    filter(
      initial_proportion %in% initial_prop_plot,
      K0_red %in% source_red_plot
    ) |>
    transmute(
      initial_proportion,
      time,
      VC_total,
      scenario = paste0("Source reduction ", K0_red)
    )
)


##
VC_total_time_series <- ggplot(
  VC_total_all  |> filter(time < 365*2.5)  |>
    mutate(
      scenario = factor(
        scenario,
        levels = c(
          "No intervention",
          "Blocking effect 0.7",
          "Blocking effect 0.9",
          "Source reduction 0.5",
          "Source reduction 0.9"
        )
      )
    ),
  aes(
    x = time,
    y = VC_total,
    colour = scenario,
    group = scenario
  )
) +
  facet_wrap(
    ~ initial_proportion,
    nrow = 1,
    labeller = labeller(
      initial_proportion = function(x)
        paste0("Initial population = ", as.numeric(x))
    )
  ) +
  scale_colour_manual(
    values = c(
      "No intervention" = "black",
      
      "Blocking effect 0.7" = "#FCA5A5",
      "Blocking effect 0.9" = "#991B1B",
      
      "Source reduction 0.5" = "#3B82F6",
      "Source reduction 0.9" = "#1E3A8A"
    )
  ) +
  scale_y_continuous() +
  geom_line(
    linewidth =1,
    alpha = 0.7
  ) +
  labs(
    x = "Time (days)",
    y = "Total vectorial capacity",
    colour = NULL
  ) +
  theme_classic() +
  theme(
    legend.position = "left"
  ) + geom_vline(
    xintercept = 365,
    linetype = "dashed",
    linewidth = 0.6,
    colour = "grey40"
  )

VC_total_time_series

#heat map 
# ============================================================
# Wolbachia
# ============================================================

#first time established 
time_to_est_plot <- initial_blocking_results |>  
  filter(
    blocking == 0.5,
    A_w / (A_w + A_u) > 0.99, 
    initial_proportion > 0.01 & initial_proportion <1 
  ) |>
  group_by(
    initial_proportion,
  ) |>
  summarise(
    first_time = min(time),
    .groups = "drop"
  ) |> ggplot (aes(x = initial_proportion, y= first_time)) +
  geom_line() +
  geom_point() +
  labs(
    x = "Initial population",
    y = expression("Days to >99% " * italic("Wolbachia-Aedes")),
    title = "Time to establishment" 
  ) + 
  theme_classic()



wolb_period <- initial_blocking_results |>
  filter(
    blocking %in% c(0.5, 0.6, 0.7, 0.8, 0.9, 1),
    initial_proportion < 1
  ) |>
  mutate(
    period = case_when(
      time > 1 & time <= 365 ~ "Days 1–365",
      time > 365 ~ "Days 365 +",
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(period)) |>
  group_by(
    initial_proportion,
    blocking,
    period
  ) |>
  summarise(
    VC_relative_mean = mean(
      VC_relative,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  mutate(
    intervention = paste0(
      "Wolbachia blocking = ",
      blocking
    )
  )



# ============================================================
# Vector control
# ===========================================================

source_reduction_period <- vector_control_only_results |>
  filter(
    K0_red %in% c(0.4, 0.5, 
                  0.6, 0.7, 0.8, 0.9) &
      initial_proportion > 0.01 
    &
      initial_proportion < 1
  ) |> mutate (
    period = case_when(
      time > 1 & time <= 365 ~ "Days 1–365",
      time > 365  ~ "Days 365 +",
      TRUE ~ NA_character_
    ) ) |>
  filter(!is.na(period)) |>
  group_by(
    initial_proportion,
    K0_red,
    period
  ) |>
  summarise(
    VC_relative_mean = mean(
      VC_relative,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  mutate(
    intervention = paste0(
      "Vector control K0 red = ",
      K0_red
    )
  )


# ============================================================
# Combine
# ============================================================
intervention_period_data <- bind_rows(
  
  wolb_period |>
    mutate(
      intervention_type = "Wolbachia",
      intervention_level = blocking
    ) |>
    select(
      initial_proportion,
      period,
      VC_relative_mean,
      intervention_type,
      intervention_level
    ),
  
  source_reduction_period |>
    mutate(
      intervention_type = "Source reduction",
      intervention_level = K0_red
    ) |>
    select(
      initial_proportion,
      period,
      VC_relative_mean,
      intervention_type,
      intervention_level
    )
  
) |>
  mutate(
    period_label = factor(
      period,
      levels = c(
        "Days 1–365",
        "Days 365 +"
      )
    ),
    
    intervention_type = factor(
      intervention_type,
      levels = c(
        "Wolbachia",
        "Source reduction"
      )
    )
  )


# ============================================================
# no intervention
# ============================================================
no_intervention <- no_intervention |>
  filter(
    initial_proportion < 1
  ) |>
  mutate(
    period = case_when(
      time >= 1 & time <= 365 ~ "Days 1–365",
      time > 365 ~ "Days 365 +",
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(period)) |>
  group_by(
    initial_proportion,
    period
  ) |>
  summarise(
    VC_relative_mean = mean(
      VC_relative,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  mutate(
    period = factor(
      period,
      levels = c(
        "Days 1–365",
        "Days 365 +"
      )
    )
  )

# ============================================================
# plot
# ============================================================


wolb_1yr <- intervention_period_data |>
  filter(
    period_label == "Days 1–365",
    intervention_type == "Wolbachia"
  )

wolb_later <- intervention_period_data |>
  filter(
    period_label == "Days 365 +",
    intervention_type == "Wolbachia"
  )

source_1yr <- intervention_period_data |>
  filter(
    period_label == "Days 1–365",
    intervention_type == "Source reduction"
  )

source_later <- intervention_period_data |>
  filter(
    period_label == "Days 365 +",
    intervention_type == "Source reduction"
  )

p1 <- make_heatmap(
  wolb_1yr,
  title_text = expression(italic(Wolbachia) ~ "— Days 1–365"),
  x_lab_text = "Blocking effect")

p2 <- make_heatmap(
  source_1yr,
  title_text =  "Source reduction — Days 1–365",
  x_lab_text ="Source reduction intensity"
)

p3 <- make_heatmap(
  wolb_later,
  title_text = expression(italic(Wolbachia) ~ "— Days 365+"),
  x_lab_text = "Blocking effect"
)

p4 <- make_heatmap(
  source_later,
  title_text = "Source reduction — Days 365+",
  x_lab_text ="Source reduction intensity"
)

p5 <- ggplot(
  no_intervention,
  aes(
    x = period,
    y = factor(initial_proportion),
    fill = VC_relative_mean
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.5
  ) +
  geom_text(
    aes(label = sprintf("%.2f", VC_relative_mean)),
    size = 3
  ) +
  scale_fill_gradient2(
    low = "darkgreen",
    mid = "white",
    high = "darkorange",
    midpoint = 1,
    limits = c(0, 4),
    oob = scales::squish,
    name = "Fold change\nvectorial capacity"
  ) +
  labs(
    x = "Period",
    y = "Initial population",
    title = "No intervention"
  ) +
  theme_classic(
    base_size = 10
  )


# Extract one legend
legend <- get_legend(
  p1 + theme(legend.position = "right")
)

# Remove legends from individual plots
p5_noleg <- p5 + theme(legend.position = "none")
p1_noleg <- p1 + theme(legend.position = "none")
p2_noleg <- p2 + theme(legend.position = "none")
p3_noleg <- p3 + theme(legend.position = "none")
p4_noleg <- p4 + theme(legend.position = "none")


# Top row
top_row <- ggarrange(
  time_to_est_plot,
  p1_noleg,
  p3_noleg,
  ncol = 3,
  labels = c("(A)", "(C)", "(E)")
)


# Bottom-left cell:
# legend on top, p5 underneath and smaller
bottom_left <- ggarrange(
  legend,
  p5_noleg,
  ncol = 2,
  widths = c(0.3, 0.7),
  labels = c("", "(B)")
)

# Bottom row
bottom_row <- ggarrange(
  bottom_left,
  p2_noleg,
  p4_noleg,
  ncol = 3,
  labels = c("", "(D)", "(F)")
)


time_series <- ggarrange(
  VC_total_time_series,
  ncol = 1,
  labels = "(G)"
)

# Final figure
combined_heatmap <- ggarrange(
  top_row,
  bottom_row,
  time_series, 
  ncol = 1,
  nrow = 3,
  heights = c(1, 1, 1)
)

combined_heatmap

ggsave(paste0("initial_com_approach_", 
              Sys.Date(),".pdf"),
       width = 13, 
       height = 10)
