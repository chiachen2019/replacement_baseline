source("turelli2020_vc_func_20251103.R")

#----1. func to run simple simulations -------
simulate_with_varying_params <- function(param_matrix, param_names, values_list) {
  baseline_params <- param_matrix[1, ]
  
  param_grid <- expand.grid(values_list)
  colnames(param_grid) <- param_names
  
  varied_param_matrix_list <- list()
  
  for (i in 1:nrow(param_grid)) {
    vals <- param_grid[i, ]
    params <- baseline_params
    params[param_names] <- as.numeric(vals)  
    varied_param_matrix_list[[i]] <- params  
  }
  
  varied_param_matrix <- do.call(rbind, varied_param_matrix_list)
  
  results <- simulate_multiple_wolbachia_matrix(varied_param_matrix)
  
  return(list(results = results, param_grid = param_grid))
}


param_matrix <- matrix(c(
  0.8431209, 0.269*0.85, 19, 0.85, 0.269, 1, 365*1.5, 1740*60,
  3500, 4618.347,
  0.2147, 1*0.3, 0.9702003, 9.50689,
  0.2147, 1, 0.9702003, 9.50689, 21 
), ncol = 19, byrow = TRUE)

colnames(param_matrix) <- c("vi", "Fi", "tau", "vu", "Fu", "s_h", "days", "K_0",
                            "I_a_initial", "U_a_initial", 
                            "ai", "bi", "ci", "eti", 
                            "au", "bu", "cu", "etu", "release_start")

#### Define multiple parameters to vary 
param_names <- c("U_a_initial", "bi")
values_list <- list(seq(46183.47*0.05,(46183.47), by = 46183.47 * 0.025), 
                    seq(0, 0.5, by = 0.1))  
####-----2. Run the simulation -----
sim_results <- simulate_with_varying_params(param_matrix, param_names, values_list)

####-----3. Extract results and parameter combinations -------
results <- sim_results$results
param_grid <- sim_results$param_grid

# Extract the necessary values from results
release_start_value <- param_matrix[1, "release_start"]

uninfected_t_infected <- sapply(results, function(res) {
  mean(res$Uninfected[res$Day < release_start_value], na.rm = TRUE)
})  
final_infected <- sapply(results, function(res) res$Infected[nrow(res)])
final_ProportionInfected <- sapply(results, function(res) res$ProportionInfected[nrow(res)])
final_perc_diff_vc <- sapply(results, function(res) {
  res$perc_diff_vc <- res$VC / mean(res$VC[res$Day < release_start_value])
  res$perc_diff_vc[length(res$perc_diff_vc)]
})

final_results <- cbind(param_grid, 
                       uninfected_t_infected = uninfected_t_infected / 46183.47 * 100, 
                       final_perc_diff_vc = final_perc_diff_vc)

final_results <- as.data.frame(final_results)
bi <- unique(final_results$bi)
final_results$reduction_transmission = 1- final_results$bi

write.csv(final_results, file = paste0("final_results", Sys.Date(),".csv"))




