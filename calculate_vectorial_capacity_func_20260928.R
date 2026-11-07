calculate_vectorial_capacity <- function(
    out,
    parameters,
    human_pop = 3000,
    biting_rate = 0.2146,
    b_u = 1,
    c_u = 0.97,
    eip = 9.5,
    blocking = 0.8,
    relative_to_baseline = TRUE,
    baseline_type = c("initial", "control_period")
) {
  
  # Match the selected baseline option
  baseline_type <- match.arg(baseline_type)
  
  out_df <- as.data.frame(out)
  
  # Check required columns
  required_columns <- c("time", "A_u", "A_w")
  
  missing_columns <- setdiff(required_columns, names(out_df))
  
  if (length(missing_columns) > 0) {
    stop(
      "The output is missing the following required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  # Extract parameter values as unnamed scalars
  t_control_start   <- parameters[["t_control_start"]]
  t_release_start   <- parameters[["t_release_start"]]
  d_A_u           <- parameters[["d_A_u"]]
  c_d_A_w         <- parameters[["c_d_A_w"]]
  
  
  # Effective adult mortality rates
  d_A_u <- parameters[["d_A_u"]]
  c_d_A_w <- parameters[["c_d_A_w"]]
  
  d_A_u_eff <- d_A_u
  d_A_w_eff <- d_A_u * (1 + c_d_A_w)
  
  # Convert continuous mortality rates to daily survival probabilities
  p_u <- exp(-d_A_u_eff)
  p_w <- exp(-d_A_w_eff)
  
  # Adult mosquitoes per human
  m_u <- 0.5*out_df$A_u / human_pop
  m_w <- 0.5*out_df$A_w / human_pop
  
  # Wolbachia-associated reduction in vector competence
  b_w <- b_u * (1 - blocking)
  c_w <- c_u
  
  # Vectorial capacity of uninfected mosquitoes
  VC_u <- (
    m_u *
      biting_rate^2 *
      b_u *
      c_u *
      p_u^eip /
      (-log(p_u))
  )
  
  # Vectorial capacity of Wolbachia-infected mosquitoes
  VC_w <- (
    m_w *
      biting_rate^2 *
      b_w *
      c_w *
      p_w^eip /
      (-log(p_w))
  )
  
  # Total vectorial capacity
  VC_total <- VC_u + VC_w
  
  result <- data.frame(
    time = out_df$time,
    A_u = out_df$A_u,
    A_w = out_df$A_w,
    p_u = p_u,
    p_w = p_w,
    VC_u = VC_u,
    VC_w = VC_w,
    VC_total = VC_total
  )
  
  # ----------------------------------------------------------
  # Calculate relative vectorial capacity
  # ----------------------------------------------------------
  
  if (relative_to_baseline) {
    
    if (baseline_type == "initial") {
      
      # Use the first simulation time as the baseline
      initial_time <- min(result$time, na.rm = TRUE)
      
      baseline_idx <- result$time == initial_time
      
    } else if (baseline_type == "control_period") {
      
      # Use the mean VC during the vector-control period
      baseline_idx <- (
        result$time >= t_control_start &
          result$time < t_release_start
      )
      
      if (!any(baseline_idx, na.rm = TRUE)) {
        stop(
          "No simulation times fall within the vector-control period. ",
          "Check t_control_start, t_release_start, and the simulation times, ",
          "or use baseline_type = 'initial'."
        )
      }
    }
    
    baseline_VC <- mean(
      result$VC_total[baseline_idx],
      na.rm = TRUE
    )
    
    if (!is.finite(baseline_VC) || baseline_VC <= 0) {
      stop(
        "The selected baseline vectorial capacity is zero, ",
        "negative, NA, or non-finite."
      )
    }
    
    result$baseline_type <- baseline_type
    result$baseline_VC <- baseline_VC
    result$VC_relative <- result$VC_total / baseline_VC
  }
  
  return(result)
}