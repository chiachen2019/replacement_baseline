library(dplyr)
library(tidyr)
library(ggplot2)
source("turelli2020_vc_80percentgoal_func_20251103.R")

## ----- 1. Load parameters ----
param_file <- "./parameters.csv"
param_df <- read.csv(param_file, stringsAsFactors = FALSE)
required_cols <- c(
  "vi", "Fi", "tau", "vu", "Fu", "s_h", "days", "K_0",
  "I_a_initial", "U_a_initial", 
  "ai", "bi", "ci",
  "eti", "au", "bu", "cu", "etu", "release_start"
)

param_df <- param_df[, required_cols]
param_matrix <- as.matrix(param_df)

## ---- 2. Run simulations ----

simulation_results <- simulate_multiple_wolbachia_80percent_matrix(param_matrix)

## ---- 3. Add fold-change in VC (perc_diff_vc) ----

# Use the release_start from the *first* row 
release_start_day <- param_df$release_start[1]

simulation_results <- lapply(simulation_results, function(df) {
  initial_avg_vc <- mean(df$VC[df$Day < release_start_day], na.rm = TRUE)
  df$perc_diff_vc <- (df$VC) / initial_avg_vc
  df
})


## ---- 4. Combine and export CSV ----

combined_df <- do.call(rbind, lapply(seq_along(simulation_results), function(i) {
  mat <- simulation_results[[i]]
  df <- as.data.frame(mat)
  df$Condition <- paste0("Condition_", i)
  df
}))

combined_df <- combined_df[, c(ncol(combined_df), 1:(ncol(combined_df)-1))]

write.csv(
  combined_df,
  paste0("temporal_output_80goal_all", Sys.Date(), ".csv"),
  row.names = FALSE
)
