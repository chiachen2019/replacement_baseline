source("turelli2020_space_func_20251208.R")

### ----- 1. Load parameters -----
param_file <- "./spatial_parameters.csv"
param_df <- read.csv(param_file, stringsAsFactors = FALSE)

required_cols <- c(
  "vi", "Fi", "tau", "vu", "Fu", "s_h", "days", "K_0",
  "I_a_initial", "U_a_initial_neighbor", "U_a_initial_target", 
  "ai", "bi", "ci",
  "eti", "au", "bu", "cu", "etu", "m_jki", "m_jku","N", "release_site", "release_start"
)

param_df <- param_df[, required_cols]
param_matrix <- as.matrix(param_df)

### ---- 2. Run simulation ---
simulation_results <- simulate_wolbachia_space_matrix(param_matrix)

## ---- 3. Combine and export CSV ----

combined_df <- do.call(rbind, lapply(seq_along(simulation_results), function(i) {
  mat <- simulation_results[[i]]
  df <- as.data.frame(mat)
  df$Condition <- paste0("Condition_", i)
  df
}))

combined_df <- combined_df[, c(ncol(combined_df), 1:(ncol(combined_df)-1))]


start_lookup <- param_df |> 
  select(release_start) |> 
  mutate(Condition = c("Condition_1", "Condition_2", "Condition_3", "Condition_4"))

combined_df <- combined_df |>
  left_join(start_lookup, by = "Condition") |>
  group_by(Habitat, Condition) |>
  mutate(
    initial_avg_vc = mean(VC[Day < release_start], na.rm = TRUE),
    perc_diff_vc   = VC / initial_avg_vc
  ) |>
  ungroup()


combined_df$day <- combined_df$Day - 20
combined_df$day <- factor(combined_df$day, levels = sort(unique(combined_df$day), decreasing = TRUE))
combined_df$Target <- ifelse(combined_df$Habitat == 3, "Target", "Non-target")

combined_df|> filter(Condition == "Condition_4") |> 
  filter(Day == 20) |> 
  mutate(baseline = Uninfected/46184) |> select(baseline)

combined_df|> filter(Condition == "Condition_4") |> filter(Day > 20+365) |> 
  group_by(Habitat) |> summarise(mean(perc_diff_vc))

write.csv(combined_df, paste0("spatial_result_", Sys.Date(),".csv"))


