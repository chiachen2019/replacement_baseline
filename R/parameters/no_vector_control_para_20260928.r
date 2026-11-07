state <- c(
  L_u = 1000, A_u = 1000,
  L_w = 0, A_w = 0
)

parameters <- c(
  b_u       = 10, #birth
  K0        = 30000, #carrying capacity 
  eta       = 0.2,
  omega     = 2 * pi / 365,
  alpha_u   = 1/10, #development rate
  d_L_u     = 1/20, #death rate
  d_A_u     = 1/14, #death rate
  
  #wolbachia cost
  c_b_w     = 0.1, #birth cost 
  c_alpha_w = 0, 
  c_d_L_w = 0.15, #mortality cost
  c_d_A_w   = 0.15, #mortality cost 
  
  # vector control
  t_control_start = Inf,
  t_control_stop  = Inf,
  K0_red          = 0.80,   # 80% reduction in K0
  dL_add          = 0.00,
  dA_add          = 0.00,
  
  # Wolbachia release
  t_release_start = 1,
  t_release_stop  = 365*1,
  release_ratio  = 1
)

release_interval <- 7

release_times <- seq(
  from = parameters["t_release_start"],
  to   = parameters["t_release_stop"],
  by   = release_interval
)

release_events <- data.frame(
  var = "A_w",
  time = release_times,
  value = unname(parameters["release_ratio"]),
  method = "add"
)

times <- seq(0, 365 * 15, by = 1)
times <- sort(unique(c(times, release_times)))
