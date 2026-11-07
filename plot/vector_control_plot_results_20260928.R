library(ggplot2)
library(ggpubr)
library(dplyr)
library(scales)


source("../R/vector_control_para_20260928.r")
source("plots_func.r")

partial_vector_control_scenario_results <- read.csv("../output/output_partial_2026-10-02.csv")

#1. heatmap, #one year since release started ---------
hp1 <- plot_vc_heatmap(
  partial_vector_control_scenario_results,
  parameters,
  post_control = 0
)

hp2 <- plot_vc_heatmap(
  partial_vector_control_scenario_results,
  parameters,
  post_control = 0.25
)

hp3 <- plot_vc_heatmap(
  partial_vector_control_scenario_results,
  parameters,
  post_control = 0.5
)

hp4 <- plot_vc_heatmap(
  partial_vector_control_scenario_results,
  parameters,
  post_control = 0.75
)

combined_heat_plot <- ggarrange(
  hp1, hp2,
  hp3, hp4,
  ncol = 2,
  nrow = 2,
  common.legend = TRUE,
  legend = "right",
  labels = c("(A)", "(B)", "(C)", "(D)")
)
combined_heat_plot

ggsave(file = paste0("heatmap_vector_control_x_blocking_", 
                     Sys.Date(),".pdf"),
       width = 10, 
       height = 10)


#2. mosquito results 
mosquito_data <- partial_vector_control_scenario_results |>  
  mutate(K0_red_post_rel = K0_red_post_rel *100) |>
  filter(
    blocking == 0.5,
    time > 365 * 8,
    time < 365 * 8 + 365,
    A_w / (A_w + A_u) > 0.99
  ) |>
  group_by(
    K0_red,
    K0_red_post_rel
  ) |>
  summarise(
    first_time = min(time) - 365 * 8,
    .groups = "drop"
  ) 

write.csv(mosquito_data, file = paste0("time_to_est_", 
Sys.Date(),".csv"))




mosquito_plot <- mosquito_data|>
  ggplot(
    aes(
      x = K0_red,
      y = first_time,
      colour = factor(K0_red_post_rel),
      group = K0_red_post_rel
    )
  ) +
  geom_line() +
  geom_point() +
  scale_colour_discrete(
    labels = function(x) paste0(x, "%")
  ) + scale_y_continuous(
    limits = c(70, 80)
  ) +  
  labs(
    x = "Baseline source reduction intensity",
    y = expression("Days to >99% " * italic("Wolbachia-Aedes")),
    colour = "Relative post-release\nsource reduction",
    title = expression(bold("(A)") * " Time to establishment") 
  ) + 
  theme_classic()

mosquito_plot <- mosquito_plot +
  theme(
    plot.margin = margin(5, 5, 5, 5),
    legend.position = "bottom"
  ) 
mosquito_column <- ggarrange(
  NULL,
  mosquito_plot,
  NULL,
  ncol = 1,
  nrow = 3,
  heights = c(0.5, 1, 0.5)
)

final_plot <- ggarrange(
  mosquito_column,
  combined_heat_plot,
  ncol = 2,
  widths = c(1, 2)
)

final_plot

ggsave(file = paste0("heatmap_vector_control_x_blocking_time_to_establishment_", 
                     Sys.Date(),".pdf"),
       width = 17, 
       height = 9)
