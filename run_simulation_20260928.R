library(ggplot2)
library(tidyverse)

# ============================================================
# EXISTING VECTOR CONTROL
# ============================================================
#
# Objective:
# Evaluate how relaxing or maintaining existing source
# reduction during Wolbachia replacement affects vectorial
# capacity.
#
# Simulation setup:
# - Source reduction begins in year 3.
# - Wolbachia releases begin in year 8.
# - Weekly releases continue for one year.
# - Each release is proportional to the mean wild-type
#   adult abundance during the source-reduction period.
#
# Parameters varied:
# - K0_red_values:
#   Baseline reduction in aquatic carrying capacity (10-90%).
#
# - post_control_values:
#   Fraction of baseline source reduction retained
#   following Wolbachia introduction (0-75%).
#
# - blocking_values:
#   Reduction in mosquito-to-human transmission due
#   to Wolbachia infection (50-100%).
#
# Outcome:
# Relative vectorial capacity compared with mean vectorial
# capacity during the pre-release source-reduction period.
#
# Total scenarios: 9 x 4 x 6 = 216
#
# ============================================================


source("R/parameters/vector_control_para_20260928.r")
source("R/vector_control_func_20260811.r")
source("R/calculate_vectorial_capacity_func_20260928.r")
source("R/run_k0_blocking_scenarios_func_20260928.r")

partial_vector_control_scenario_results <- run_k0_blocking_scenarios(
  
  parameters = parameters,
  state = state,
  times = times,
  K0_red_values = c(0.1, 0.2, 0.3, 0.4,
                    0.5, 0.6, 0.7, 0.8, 0.9),
  post_control_values = c(0, 0.25, 0.5, 0.75),
  blocking_values = c(0.5, 0.6, 0.7, 0.8, 0.9, 1)
)

partial_vector_control_scenario_results|> 
  write.csv(file = paste0("output/output_partial_", Sys.Date(),".csv"))




# ============================================================
# WOLBACHIA REPLACEMENT IN INITIALLY LOW
# MOSQUITO POPULATIONS
# ============================================================
#
# Objective:
# Evaluate the effectiveness of Wolbachia replacement in
# recently invaded areas where mosquito populations are
# initially below their ecological equilibrium.
#
# Simulation setup:
# - No conventional source reduction is applied.
# - Initial wild-type mosquito abundance varies from
#   10% to 100% of its equilibrium population.
# - Wolbachia releases begin on day 1 and continue
#   weekly for one year.
# - Release numbers are proportional to the initial
#   wild-type adult mosquito population.
#
# Parameters varied:
# - Initial population: 10-100% of equilibrium
# - Viral blocking efficacy: 50-100%
#
# Outcome:
# Fold change in vectorial capacity relative to the
# initial population before Wolbachia introduction.
#
# ============================================================

source("R/parameters/no_vector_control_para_20260928.r")
source("R/vector_control_func_20260811.r")
source("R/calculate_vectorial_capacity_func_20260928.r")
source("R/run_initial_blocking_scenarios_func_20260928.R")

initial_blocking_results <- run_initial_blocking_scenarios(
  parameters = parameters,
  times = times,
  initial_proportions = c(0.1, 0.20, 0.30, 0.40, 0.5, 0.6, 0.7, 0.8, 0.9, 1.00),
  blocking_values = c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
)


initial_blocking_results |> write.csv(file = paste0("output/output_initial_blocking_results_", Sys.Date(),".csv"))



# ============================================================
# NO INTERVENTION (BASELINE REFERENCE)
# ============================================================
#
# Objective:
# Evaluate mosquito population growth and changes in
# vectorial capacity in the absence of any intervention.
#
# Simulation setup:
# - No Wolbachia-infected mosquitoes are released.
# - No conventional source reduction is applied.
# - Initial wild-type mosquito abundance varies from
#   10% to 100% of its equilibrium population.
# - Mosquito populations evolve naturally towards
#   ecological equilibrium under seasonal conditions.
#
# Parameters:
# - Initial population: 10-100% of equilibrium
# - Release ratio: 0 (no Wolbachia releases)
# - Viral blocking efficacy: 0 (not applicable)
#
# Outcome:
# Fold change in vectorial capacity relative to the
# initial population, without intervention.
#
# ============================================================

no_release <- run_initial_blocking_scenarios(
  parameters = parameters,
  times = times,
  initial_proportions = c(
    0.1, 0.2, 0.3, 0.4, 0.5,
    0.6, 0.7, 0.8, 0.9, 1.0
  ),
  blocking_values = 0,
  release_ratio = 0
  
)

no_release |> write.csv(file = paste0("output/no_release_results_", Sys.Date(),".csv"))


# ============================================================
# CONVENTIONAL VECTOR CONTROL ONLY IN INITIALLY LOW
# MOSQUITO POPULATIONS
# ============================================================
#
# Objective:
# Evaluate the effectiveness of conventional source
# reduction as an alternative to Wolbachia replacement
# in initially low mosquito populations.
#
# Simulation setup:
# - No Wolbachia-infected mosquitoes are released.
# - Source reduction begins on day 1 and is maintained
#   throughout the simulation period.
# - Initial wild-type mosquito abundance varies from
#   1% to 100% of its equilibrium population.
# - Source reduction decreases the aquatic carrying
#   capacity without directly removing adult mosquitoes.
#
# Parameters varied:
# - Initial population: 1-100% of equilibrium
# - Source reduction intensity: 10-90%
# - Release ratio: 0 (no Wolbachia releases)
# - Viral blocking efficacy: 0 (not applicable)
#
# Outcome:
# Fold change in vectorial capacity relative to the
# initial population before source reduction.
#
# ============================================================

source("R/parameters/vector_control_para_20260928.r")
source("R/vector_control_func_20260811.r")
source("R/calculate_vectorial_capacity_func_20260928.r")
source("R/run_initial_blocking_scenarios_func_20260928.R")

# Start vector control at day 1
parameters[["t_control_start"]] <- 1

# Vector control does NOT change at the nominal release time
# because there is actually no Wolbachia release
parameters[["K0_red_post_rel"]] <- 1


K0_red_values <- c(
  0.1, 0.2, 0.3, 0.4, 0.5,
  0.6, 0.7, 0.8, 0.9
)

vector_control_only_list <- list()

for (k in K0_red_values) {
  
  parameters_k <- parameters
  parameters_k[["K0_red"]] <- k
  
  temp_results <- run_initial_blocking_scenarios(
    parameters = parameters_k,
    times = times,
    initial_proportions = c(
      0.1, 0.2, 0.3, 0.4, 0.5,
      0.6, 0.7, 0.8, 0.9, 1.0
    ),
    # No Wolbachia, so blocking is irrelevant
    blocking_values = 0,
    release_ratio = 0
  )
  
  # Record vector-control intensity
  temp_results$K0_red <- k
  vector_control_only_list[[as.character(k)]] <- temp_results
}


# Combine all K0 scenarios
vector_control_only_results <- bind_rows(
  vector_control_only_list
)

vector_control_only_results |>
  write.csv(
    file = paste0("output/vector_control_only_results_",Sys.Date(),".csv"
    ))

