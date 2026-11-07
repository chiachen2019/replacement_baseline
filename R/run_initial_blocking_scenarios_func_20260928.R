run_initial_blocking_scenarios <- function(
    parameters,
    times,
    initial_proportions = c(0.10, 0.25, 0.50, 0.75, 1.00),
    blocking_values = c(0.5, 0.7, 0.8, 0.9, 1.0),
    release_ratio = 1,
    release_interval = 7,
    baseline_type = c("initial")
) {
  
  library(deSolve)
  
  # ------------------------------------------------------------
  # 1. Calculate Wolbachia-free, non-seasonal equilibrium
  # ------------------------------------------------------------
  
  b_u     <- parameters[["b_u"]]
  K0      <- parameters[["K0"]]
  alpha_u <- parameters[["alpha_u"]]
  d_L_u   <- parameters[["d_L_u"]]
  d_A_u   <- parameters[["d_A_u"]]
  
  equilibrium_term <- 1 -
    d_A_u * (alpha_u + d_L_u) /
    (b_u * alpha_u)
  
  if (!is.finite(equilibrium_term) || equilibrium_term <= 0) {
    stop(
      "A positive uninfected equilibrium does not exist. ",
      "Check b_u, alpha_u, d_L_u and d_A_u."
    )
  }
  
  L_u_eq <- K0 * equilibrium_term
  A_u_eq <- alpha_u / d_A_u * L_u_eq
  
  
  # ------------------------------------------------------------
  # 2. Create scenario combinations
  # ------------------------------------------------------------
  
  scenario_grid <- expand.grid(
    initial_proportion = initial_proportions,
    blocking = blocking_values,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  all_results <- vector(
    mode = "list",
    length = nrow(scenario_grid)
  )
  
  
  # ------------------------------------------------------------
  # 3. Run scenarios
  # ------------------------------------------------------------
  
  for (i in seq_len(nrow(scenario_grid))) {
    
    this_initial_proportion <-
      scenario_grid$initial_proportion[i]
    
    this_blocking <-
      scenario_grid$blocking[i]
    
    
    # ----------------------------------------------------------
    # Initial population
    # ----------------------------------------------------------
    
    this_state <- c(
      L_u = this_initial_proportion * L_u_eq,
      A_u = this_initial_proportion * A_u_eq,
      L_w = 0,
      A_w = 0
    )
    
    # Safety check
    if (!identical(
      names(this_state),
      c("L_u", "A_u", "L_w", "A_w")
    )) {
      stop(
        "Incorrect state names: ",
        paste(names(this_state), collapse = ", ")
      )
    }
    
    
    # ----------------------------------------------------------
    # Wolbachia release
    # Each release = 1.5 × initial Au
    # ----------------------------------------------------------
    
    initial_A_u <- unname(this_state["A_u"])
    
    release_number <- release_ratio * initial_A_u
    
    release_times <- seq(
      from = parameters["t_release_start"],
      to   = parameters["t_release_stop"],
      by   = release_interval
    )
    
    release_events <- data.frame(
      var = "A_w",
      time = release_times,
      value = release_number,
      method = "add"
    )
    
    
    # ----------------------------------------------------------
    # Run model
    # ----------------------------------------------------------
    
    this_times <- sort(unique(c(times, release_times)))
    
    out <- ode(
      y = this_state,
      times = this_times,
      func = vector_control_wolbachia_model,
      parms = parameters,
      events = list(data = release_events)
    )
    
    out_df <- as.data.frame(out)
    
    
    # ----------------------------------------------------------
    # Calculate vectorial capacity
    # ----------------------------------------------------------
    
    vc_df <- calculate_vectorial_capacity(
      out = out_df,
      parameters = parameters,
      blocking = this_blocking,
      baseline_type = baseline_type
    )
    
    
    # ----------------------------------------------------------
    # Scenario identifiers
    # ----------------------------------------------------------
    
    vc_df$initial_proportion <- this_initial_proportion
    vc_df$initial_percentage <- 100 * this_initial_proportion
    vc_df$blocking <- this_blocking
    
    
    # ----------------------------------------------------------
    # Initial and equilibrium values
    # ----------------------------------------------------------
    
    vc_df$initial_L_u <- unname(this_state["L_u"])
    vc_df$initial_A_u <- initial_A_u
    
    vc_df$equilibrium_L_u <- L_u_eq
    vc_df$equilibrium_A_u <- A_u_eq
    
    
    # ----------------------------------------------------------
    # Release information
    # ----------------------------------------------------------
    
    vc_df$release_ratio <- release_ratio
    vc_df$release_number <- release_number
    
    
    all_results[[i]] <- vc_df
  }
  
  
  # ------------------------------------------------------------
  # 4. Combine results
  # ------------------------------------------------------------
  
  results <- do.call(rbind, all_results)
  
  rownames(results) <- NULL
  
  return(results)
}
