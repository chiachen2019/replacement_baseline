simulate_multiple_wolbachia_80percent_matrix <- function(param_matrix) {
  all_results <- list()
  
  for (i in 1:nrow(param_matrix)) {
    params <- as.list(param_matrix[i, ])
    names(params) <- colnames(param_matrix)
    
    vi <- params$vi
    Fi <- params$Fi
    tau <- params$tau
    vu <- params$vu
    Fu <- params$Fu
    s_h <- params$s_h
    days <- params$days
    K_0 <- params$K_0
    I_a_initial <- params$I_a_initial
    U_a_initial <- params$U_a_initial
    ai <- params$ai
    bi <- params$bi
    ci <- params$ci
    eti <- params$eti
    au <- params$au
    bu <- params$bu
    cu <- params$cu
    etu <- params$etu
    
    # New parameter: release_start (default 28 if not provided)
    release_start <- if (!is.null(params$release_start)) as.integer(params$release_start) else 28
    release_start <- max(1L, min(as.integer(release_start), as.integer(days)))  # clamp to [1, days]
    
    
    I_a <- numeric(days)
    U_a <- numeric(days)
    P <- numeric(days)
    lambda_u <- numeric(days)
    lambda_i <- numeric(days)
    
    if (release_start > 1) I_a[1:(release_start - 1)] <- 0
    I_a[release_start] <- I_a_initial
    U_a[1:21] <- U_a_initial
    P[1] <- I_a[1] / (I_a[1] + U_a[1])
    lambda_u[1] <- 0
    lambda_i[1] <- 0
    
    release_allowed <- TRUE  # Track whether releases are still allowed
    
    t_start <- max(as.integer(tau) + 1L, 2L) 
    
    for (t in t_start:(days - 1)) {
      if (t >= 21 & t < release_start) {
        U_a[t + 1] <- vu * U_a[t] + Fu * (1 - (s_h * (P[t - tau]))) * U_a[t - tau] * (1 - (I_a[t] + U_a[t]) / K_0)
      }
      
      # Add release only if still allowed
      if (t > release_start && (t - release_start) %% 7 == 0 && release_allowed) {
        I_a[t] <- I_a[t] + I_a_initial
      }
      
      if (t >= release_start) {
        I_a[t + 1] <- vi * I_a[t] + Fi * I_a[t - tau] * (1 - (I_a[t] + U_a[t]) / K_0)
        U_a[t + 1] <- vu * U_a[t] + Fu * (1 - (s_h * (P[t - tau]))) * U_a[t - tau] * (1 - (I_a[t] + U_a[t]) / K_0)
      }
      
      # Update proportion and growth
      P[t + 1] <- I_a[t + 1] / (I_a[t + 1] + U_a[t + 1])
      lambda_u[t + 1] <- U_a[t + 1] / U_a[t]
      lambda_i[t + 1] <- I_a[t + 1] / I_a[t]
      
      # Check if threshold exceeded, and stop future releases permanently
      if (release_allowed && P[t + 1] > 0.8) {
        release_allowed <- FALSE
      }
    }
    
    VCi <- (I_a * ai^2 * bi * ci * exp(-(-log(vi)) * eti)) / -log(vi)
    VCu <- (U_a * au^2 * bu * cu * exp(-(-log(vu)) * etu)) / -log(vu)
    VC <- (VCi + VCu) / 3000
    
    results <- data.frame(
      Day = 1:days,
      Infected = I_a,
      Uninfected = U_a,
      ProportionInfected = P,
      lambda_u = lambda_u,
      lambda_i = lambda_i,
      VCi = VCi,
      VCu = VCu,
      VC = VC,
      Simulation = i
    )
    
    all_results[[i]] <- results
  }
  
  return(all_results)
}
