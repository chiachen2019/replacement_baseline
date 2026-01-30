simulate_wolbachia_space_matrix <- function(param_matrix) {
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
    ai <- params$ai
    bi <- params$bi
    ci <- params$ci
    eti <- params$eti
    au <- params$au
    bu <- params$bu
    cu <- params$cu
    etu <- params$etu
    
    m_jku <- params$m_jku
    m_jki <- params$m_jki
    N <- params$N
    release_site<- params$release_site
    
    I_a_initial<- params$I_a_initial
    U_a_initial_neighbor<- params$U_a_initial_neighbor
    U_a_initial_target<- params$U_a_initial_target
    
    # New parameter: release_start (default 28 if not provided)
    release_start <- if (!is.null(params$release_start)) as.integer(params$release_start) else 28
    release_start <- max(1L, min(as.integer(release_start), as.integer(days)))  # clamp to [1, days]
    
    
    # Initial state
    I_a <- matrix(0, nrow = days, ncol = N)
    U_a <- matrix(0, nrow = days, ncol = N)
    P <- matrix(0, nrow = days, ncol = N)
    
    if (release_start > 1) I_a[1:(release_start - 1)] <- 0
    I_a[release_start, release_site] <- I_a_initial
    
    U_a[1:21, ] <- U_a_initial_neighbor
    U_a[1:21, release_site] <- U_a_initial_target
    
    
    P <- matrix(0, nrow = days, ncol = N)
    VCi <- matrix(0, nrow = days, ncol = N)
    VCu <- matrix(0, nrow = days, ncol = N)
    VC <- matrix(0, nrow = days, ncol = N)
    
    t_start <- max(as.integer(tau) + 1L, 2L) 
    
    for (t in t_start:(days - 1)) {
      for (j in 1:N) {
        # Dispersal terms (to left and right neighbors)
        dispersal_I <- 0
        dispersal_U <- 0
        
        if (j > 1) {
          dispersal_I <- dispersal_I + m_jki * (I_a[t, j - 1] - I_a[t, j])
          dispersal_U <- dispersal_U + m_jku * (U_a[t, j - 1] - U_a[t, j])
        }
        if (j < N) {
          dispersal_I <- dispersal_I + m_jki * (I_a[t, j + 1] - I_a[t, j])
          dispersal_U <- dispersal_U + m_jku * (U_a[t, j + 1] - U_a[t, j])
        }
        
        # Population dynamics with dispersal
        if (t >= 21 & t < release_start) {
          U_a[t + 1, j] <- vu * U_a[t, j] + 
            Fu * (1 - (s_h * P[t - tau, j])) * U_a[t - tau, j] * (1 - (I_a[t, j] + U_a[t, j]) / K_0) + 
            dispersal_U
        }
        
        if (t > release_start & t %% 7 == 0 & j == release_site) {
          I_a[t, j] <- I_a[t, j] + I_a_initial
        }
        
        if (t >= release_start) {
          I_a[t + 1, j] <- vi * I_a[t, j] + Fi * I_a[t - tau, j] * (1 - (I_a[t, j] + U_a[t, j]) / K_0) + 
            dispersal_I
          U_a[t + 1, j] <- vu * U_a[t, j] + Fu * (1 - (s_h * P[t - tau, j])) * U_a[t - tau, j] * (1 - (I_a[t, j] + U_a[t, j]) / K_0) + 
            dispersal_U
        }
        
        total_pop <- I_a[t + 1, j] + U_a[t + 1, j]
        P[t + 1, j] <- I_a[t + 1, j] / total_pop
        
      }
      
    }
    
    VCi <- (I_a * ai^2 * bi * ci * exp(-(-log(vi)) * eti)) / -log(vi) 
    VCu <- (U_a * au^2 * bu * cu * exp(-(-log(vu)) * etu)) / -log(vu)
    VC <- (VCi + VCu)/3000
    
    results_long <- data.frame(
      expand.grid(
        Day = 1:days,
        Habitat = 1:N
      ),
      Infected = as.vector(I_a),
      Uninfected = as.vector(U_a),
      ProportionInfected = as.vector(P),
      VCi = as.vector(VCi),
      VCu = as.vector(VCu),
      VC = as.vector(VC),
      Simulation = i
    )
    
    all_results[[i]] <- results_long
  }
  
  return(all_results)
}





