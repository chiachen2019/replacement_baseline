library(tidyverse)
library(ggplot2)
library(ggpubr)
library(colorspace)  
library(viridis)  

##-----1. plot time series ------
multiple_release <- read.csv("temporal_output_all2026-01-30.csv") 
target_release <- read.csv("temporal_output_80goal_all2026-01-30.csv") 


data <- multiple_release |> 
  mutate(release = "Continuous") |> 
  bind_rows(target_release |> mutate(release = "Threshold")) |>
  mutate(baseline = case_when(
    Condition %in% c("Condition_1", "Condition_4", "Condition_7") ~ 10,
    Condition %in% c("Condition_2", "Condition_5", "Condition_8") ~ 20,
    Condition %in% c("Condition_3", "Condition_6", "Condition_9") ~ 90)) |>
  mutate(start_day = case_when(
    Condition %in% c("Condition_1", "Condition_2", "Condition_3") ~ 28-20,
    Condition %in% c("Condition_4", "Condition_5", "Condition_6") ~ 21-20,
    Condition %in% c("Condition_7", "Condition_8", "Condition_9") ~ 50-20)) |> 
  filter(start_day == 1)


colors <- c(
  "10" = "darkorange",
  "20" = "purple",
  "90" = "blue"
) 

theme_black <- theme(
  text = element_text(color = "black"),
  axis.text = element_text(color = "black"),
  axis.title = element_text(color = "black"),
  strip.text = element_text(color = "black"),
  legend.text = element_text(color = "black"),
  legend.title = element_text(color = "black")
)

timeseries_a <- data |> filter(Day > 19) |> 
  ggplot() +
  geom_line(aes(x = Day-20, y = ProportionInfected, 
                col = as.factor(baseline))) +
  facet_wrap(~ release, nrow = 2) +
  geom_hline(yintercept = 1) +
  scale_color_manual(
    name = "Baseline population ratio %",
    values = colors
  ) + 
  theme_classic() +
  labs(y = "Wolbachia establishment", 
       x= "Time (day)") +theme_black

timeseries_b <- data |> 
  filter(Day > 19) |> 
  ggplot() +
  geom_line(aes(x = Day-20, y = VC, col = as.factor(baseline))) +
  facet_wrap(~ release, nrow = 2) +
  scale_color_manual(
    name = "Baseline population ratio %",
    values = colors
  ) + 
  theme_classic() +
  labs(y = "Vectorial capacity", 
       x= "Time (day)")+theme_black

timeseries_c<-data |> filter(Day > 19) |> 
  ggplot() + 
  geom_line(aes(x= Day-20, y = perc_diff_vc, col = as.factor(baseline))) +
  facet_wrap(~release, nrow = 2) +
  scale_color_manual(
    name = "Baseline population ratio %",
    values = colors
  ) + 
  geom_hline(yintercept = 1, col = "red") + 
  theme_classic()+
  labs(y = "Fold change in vectorial capacity", 
       x= "Time (day)")+theme_black


timeseries_plot <- ggarrange(timeseries_a, timeseries_b, timeseries_c,
                             ncol = 3, nrow = 1,
                             labels = c("(A)", "(B)", "(C)"), common.legend = TRUE) 
timeseries_plot
ggsave(paste0("timeseries_plot_", Sys.Date(),".pdf"), 
       timeseries_plot, width = 8, height = 6)


###----2. varying blocking effect and baseline population ratio --------
vary_para <- bind_rows(read.csv("final_results_goal802026-01-30.csv") |> 
                         mutate(release = "Threshold"), 
                       read.csv("final_results2026-01-30.csv") |> 
                         mutate(release = "Continuous")) 

#exact key results 
vary_para |> filter(final_perc_diff_vc<0.5 & bi == 0.2) |> group_by(release) |> 
  summarise(min(uninfected_t_infected)) 

vary_para |> filter(final_perc_diff_vc<0.5 & bi == 0.5) |> group_by(release) |> 
  summarise(min(uninfected_t_infected)) 


# generate plots continuous 

pdf(paste0("varying_para_", Sys.Date(), ".pdf"), width = 8, height = 6)
ggplot(read.csv("final_results2026-01-30.csv"), aes(x = final_results[, 3], y = final_results[, 4])) +
  geom_point(aes(#color = factor(I_a_initial), 
    fill = factor(reduction_transmission)), 
    shape = 21, size = 3, stroke = 0.5) +
  scale_fill_viridis_d() +  # Viridis for bi
  labs(x = "Baseline population ratio (%)", 
       y = "Fold change in vectorial capacity",
       fill = "Reduction in transmission") +
  geom_hline(yintercept = 1, color = "red", size = 1) +  # Horizontal line at y = 0
  theme_minimal() +
  theme(legend.position = "top",
        legend.box = "vertical") +  # Stack legends in vertical layout
  guides(fill = guide_legend(nrow = 1),  # Arrange fill legend in 1 row
         color = guide_legend(nrow = 1)) +
  ylim(c(0, 8)) +
  ggtitle ("(A) Continuous strategy")
dev.off()

#### generate plots threshold 
pdf(paste0("varying_para_80goal", Sys.Date(),".pdf"), width = 8, height = 6)
ggplot(read.csv("final_results_goal802026-01-30.csv"), aes(x = final_results[, 3], y = final_results[, 4])) +
  geom_point(aes( 
    fill = factor(reduction_transmission)), 
    shape = 21, size = 3, stroke = 0.5) +
  scale_fill_viridis_d() +  
  labs(x = "Baseline population ratio (%)", 
       y = "Fold change in vectorial capacity",
       fill = "Reduction in transmission") +
  geom_hline(yintercept = 1, color = "red", size = 1) + 
  theme_minimal() +
  theme(legend.position = "top",
        legend.box = "vertical") +  
  guides(fill = guide_legend(nrow = 1),  
         color = guide_legend(nrow = 1)) +
  ylim(c(0, 7)) + 
  ggtitle ("(B) Threshold strategy")
dev.off()


### --- 3. spatial results -----
spatial <- read.csv("spatial_result_2026-01-30.csv")

#extract key results
spatial|> filter(Condition == "Condition_4") |> filter(Day == 20) |> 
  mutate(baseline = Uninfected/46184)

spatial|> filter(Condition == "Condition_4") |> 
  filter(Day > 20+365) |> 
  group_by(Habitat) |> summarise(mean(perc_diff_vc),
                                 q025_vc   = quantile(perc_diff_vc, 0.025, na.rm = TRUE),
                                 q975_vc   = quantile(perc_diff_vc, 0.975, na.rm = TRUE))
1-0.371

spatial |>
  filter(Condition == "Condition_4", Habitat %in% c(1, 2, 4, 5), 
         Day > 21) |>
  summarise(
    mean_vc   = mean(perc_diff_vc, na.rm = TRUE),
    q025_vc   = quantile(perc_diff_vc, 0.025, na.rm = TRUE),
    q975_vc   = quantile(perc_diff_vc, 0.975, na.rm = TRUE)
  )

###generate plots 
plot_wolbachia_tile <- function(data,
                                condition = "Condition_1",
                                day_points = c(1, 27, 100, 150, 200, 250, 300, 350, 400, 450)) {
  
  data |>
    filter(Condition == condition,
           day %in% day_points) |>
    ggplot(aes(x = as.factor(Habitat),
               y = as.factor(day),
               fill = ProportionInfected,
               color = Target)) +
    geom_tile(width = 0.9, height = 0.9) +
    scale_fill_gradient(low = "white", high = "black",
                        name = "Prop Wolbachia") +
    scale_color_manual(values = c("Target" = "#D5006D", "Non-target" = "black")) +
    labs(x = "Location", y = " ", fill = "Wolbachia Proportion") +
    scale_y_discrete(labels = function(x) paste0("Day = ", x)) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      legend.position = "right"
    ) +
    guides(color = "none")
}

plot_FCVC_tile <- function(data,
                           condition = "Condition_1",
                           day_points = c(1, 27, 100, 150, 200, 250, 300, 350, 400, 450)) {
  
  ggplot(combined_df |> 
           filter(Condition == condition) |> 
           filter(day %in% day_points),
         aes(x = as.factor(Habitat), y = as.factor(day),
             fill = perc_diff_vc, color = Target)) +
    geom_tile(width = 0.9, height = 0.9) +
    scale_fill_gradient2(low = "green", high = "red", midpoint = 1, limits = c(0, 4), ) +  # Grayscale color based on VC
    scale_color_manual(values = c("Target" = "#D5006D", "Non-target" = "black")) +
    labs(x = "Location", y = " ", fill = "Fold change \nin vectorial capacity") +
    scale_y_discrete(labels = function(x) paste0("Day = ", x))  +  # Format y-axis labels
    theme_minimal() +
    theme(panel.grid = element_blank(),
          legend.position = "right") +
    guides(color = "none")
}


p1 <- plot_wolbachia_tile(combined_df,
                          condition = "Condition_4",
                          day_points = c(1, 50, 100, 150, 200, 250, 300, 350, 400, 450))
p2 <- plot_FCVC_tile(combined_df,
                     condition = "Condition_4",
                     day_points = c(1, 50, 100, 150, 200, 250, 300, 350, 400, 450))   

final_plot <- ggarrange(
  p1,
  p2, 
  labels = c("(A)", "(B)"),
  ncol = 2, legend = "bottom",
  font.label  = list(size = 12, face = "bold"))

final_plot
ggsave(paste0("spatial_plot", Sys.Date(), ".pdf"), 
       final_plot, width = 8, height = 4)



