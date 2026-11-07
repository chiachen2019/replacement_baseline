run_k0_blocking_scenarios <- function(
    parameters,
    state,
    times,
    K0_red_values = c(0, 0.25, 0.50, 0.75, 0.90),
    post_control_values = c(0, 0.5),
    blocking_values = c(0.5, 0.7, 0.8, 0.9, 1.0),
    release_ratio = 1,
    release_interval = 7,
    baseline_type = "control_period"
) {
  
  library(deSolve)
  
  scenario_grid <- expand.grid(
    K0_red = K0_red_values,
    K0_red_post_rel = post_control_values,
    blocking = blocking_values
  )
  
  all_results <- list()
  
  for (i in seq_len(nrow(scenario_grid))) {
    
    this_K0_red <- scenario_grid$K0_red[i]
    this_post_control <- scenario_grid$K0_red_post_rel[i]
    this_blocking <- scenario_grid$blocking[i]
    
    this_parameters <- parameters
    this_parameters["K0_red"] <- this_K0_red
    this_parameters["K0_red_post_rel"] <- this_post_control
    
    
    # ------------------------------------------------------
    # 1. Simulate up to the start of Wolbachia release
    # ------------------------------------------------------
    
    pre_times <- times[
      times <= this_parameters["t_release_start"]
    ]
    
    out_pre <- ode(
      y = state,
      times = pre_times,
      func = vector_control_wolbachia_model,
      parms = this_parameters
    )
    
    out_pre_df <- as.data.frame(out_pre)
    
    
    # ------------------------------------------------------
    # 2. Mean Au during vector-control period
    # ------------------------------------------------------
    
    control_idx <- (
      out_pre_df$time >= this_parameters["t_control_start"] &
        out_pre_df$time <  this_parameters["t_release_start"]
    )
    
    mean_Au_control <- mean(
      out_pre_df$A_u[control_idx],
      na.rm = TRUE
    )
    
    
    # ------------------------------------------------------
    # 3. Calculate release number
    # ------------------------------------------------------
    
    release_number <- release_ratio * mean_Au_control
    
    
    # ------------------------------------------------------
    # 4. Create scenario-specific release events
    # ------------------------------------------------------
    
    release_times <- seq(
      from = this_parameters["t_release_start"],
      to   = this_parameters["t_release_stop"],
      by   = release_interval
    )
    
    release_events <- data.frame(
      var = "A_w",
      time = release_times,
      value = release_number,
      method = "add"
    )
    
    
    # ------------------------------------------------------
    # 5. Run full simulation
    # ------------------------------------------------------
    
    this_times <- sort(unique(c(times, release_times)))
    
    out <- ode(
      y = state,
      times = this_times,
      func = vector_control_wolbachia_model,
      parms = this_parameters,
      events = list(data = release_events)
    )
    
    out_df <- as.data.frame(out)
    
    
    # ------------------------------------------------------
    # 6. Calculate vectorial capacity
    # ------------------------------------------------------
    
    vc_df <- calculate_vectorial_capacity(
      out = out_df,
      parameters = this_parameters,
      blocking = this_blocking,
      baseline_type = baseline_type
    )
    
    
    # ------------------------------------------------------
    # 7. Store scenario information
    # ------------------------------------------------------
    
    vc_df$K0_red <- this_K0_red
    vc_df$K0_red_post_rel <- this_post_control
    vc_df$blocking <- this_blocking
    
    vc_df$mean_Au_control <- mean_Au_control
    vc_df$release_number <- release_number
    vc_df$release_ratio <- release_ratio
    
    all_results[[i]] <- vc_df
  }
  
  results <- do.call(rbind, all_results)
  
  rownames(results) <- NULL
  
  return(results)
}
