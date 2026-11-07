vector_control_wolbachia_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # 1. Vector control
    if (t < t_control_start) {
      
      K0_red_eff <- 0
      
    } else if (t < t_release_start) {
      
      K0_red_eff <- K0_red
      
    } else {
      
      K0_red_eff <- K0_red_post_rel * K0_red
      
    }
    
    K0_eff <- K0 * (1 - K0_red_eff)
    
    
    # Wolbachia life-history parameters
    b_w_base   <- b_u * (1 - c_b_w)
    alpha_w    <- alpha_u * (1 - c_alpha_w)
    d_L_w <- d_L_u * (1 + c_d_L_w)
    d_A_w <- d_A_u * (1 + c_d_A_w)
    
    # Seasonal carrying capacity
    K_t <- K0_eff * (1 + eta * cos(omega * t))
    
    total_L <- L_u + L_w
    total_A <- A_u + A_w
    
    mating_frac_u <- ifelse(total_A > 0, A_u / total_A, 0)
    
    deriv_L_u <- b_u * A_u * mating_frac_u *
      (1 - total_L / K_t) -
      alpha_u * L_u -
      d_L_u * L_u
    
    deriv_A_u <- alpha_u * L_u -
      d_A_u * A_u
    
    deriv_L_w <- b_w_base * A_w *
      (1 - total_L / K_t) -
      alpha_w * L_w -
      d_L_w * L_w
    
    deriv_A_w <- alpha_w * L_w -
      d_A_w * A_w
    
    list(c(
      deriv_L_u,
      deriv_A_u,
      deriv_L_w,
      deriv_A_w
    ))
  })
}
