plot_vc_scenario <- function(
    scenario_results,
    parameters,
    post_control = 0
) {
  
  library(dplyr)
  library(ggplot2)
  
  # Select post-release vector-control scenario
  plot_df <- scenario_results |>
    filter(
      abs(K0_red_post_rel - post_control) < 1e-8
    ) |>
    mutate(
      vector_control_intensity = paste0(
        "Source reduction = ",
        K0_red * 100,
        "%"
      )
    )
  
  
  ggplot(
    plot_df,
    aes(
      x = time,
      y = VC_relative,
      group = interaction(blocking, K0_red)
    )
  ) +
    
    # Before Wolbachia release: black line
    geom_line(
      data = plot_df |>
        filter(time <= parameters[["t_release_start"]]),
      color = "black",
      linewidth = 0.8,
      alpha = 1
    ) +
    
    # After Wolbachia release: coloured by blocking
    geom_line(
      data = plot_df |>
        filter(time >= parameters[["t_release_start"]]),
      aes(color = factor(blocking)),
      linewidth = 0.8,
      alpha = 1
    ) +
    
    # Vector-control start
    geom_vline(
      xintercept = parameters[["t_control_start"]],
      linetype = "dotted",
      color = "black"
    ) +
    
    # Wolbachia release start/stop
    geom_vline(
      xintercept = parameters[["t_release_start"]],
      linetype = "dashed",
      color = "red"
    ) +
    
    geom_vline(
      xintercept = parameters[["t_release_stop"]],
      linetype = "dashed",
      color = "red"
    ) +
    
    geom_hline(
      yintercept = 1,
      linetype = 2
    ) +
    
    facet_wrap(
      ~ vector_control_intensity,
      nrow = 3,
      scales = "free_y"
    ) +
    
    labs(
      x = "Time",
      y = "Fold change vectorial capacity",
      color = "Blocking",
      title = paste0(
        "Relative post-release source reduction = ",
        post_control * 100,
        "%"
      )
    ) +
    
    theme_classic(base_size = 8)
}


plot_vc_heatmap <- function(
    scenario_results,
    parameters,
    post_control = 0
) {
  
  library(dplyr)
  library(ggplot2)
  
  heatmap_df <- scenario_results |>
    filter(
      abs(K0_red_post_rel - post_control) < 1e-8,
      time >= parameters[["t_release_stop"]]
    ) |>
    group_by(K0_red, blocking) |>
    summarise(
      VC_relative_mean = mean(
        VC_relative,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
  
  
  ggplot(
    heatmap_df,
    aes(
      y = factor(K0_red),
      x = factor(blocking),
      fill = VC_relative_mean
    )
  ) +
    
    geom_tile(
      color = "white"
    ) +
    
    geom_text(
      aes(
        label = round(
          VC_relative_mean,
          2
        )
      ),
      size = 3
    ) +
    
    scale_fill_gradient2(
      low = "darkgreen",
      mid = "white",
      high = "darkorange",
      midpoint = 1,
      limits = c(0, 4),
      oob = scales::squish,
      breaks = c(0, 1, 2, 3, 4),
      name = "Fold change in \nvectorial capacity"
    ) +
    
    labs(
      x = expression(italic(Wolbachia) ~ "blocking effect"),
      y = "Baseline source reduction intensity",
      title = paste0(
        "Post release source reduction reduced to ",
        post_control * 100,
        "% of baseline"
      )
    ) +
    
    theme_classic(
      base_size = 10
    )
}

make_vc_spatial_heatmap <- function(data, panel_label, show_legend = TRUE) {
  
  p <- ggplot(
    data |> filter(time > 8 * 365, blocking < 1),
    aes(
      x = factor(patch),
      y = time - 8 * 365,
      fill = VC_relative
    )
  ) +
    geom_tile() +
    scale_fill_gradient2(
      low = "#003300",
      mid = "white",
      high = "darkorange",
      midpoint = 1,
      limits = c(0, 5),
      oob = scales::squish,
      name = "Fold change\nVC"
    ) +
    facet_grid(
      blocking ~ panel_lab,
      labeller = labeller(
        blocking = function(x) x,
        panel_lab = function(x) x
      )
    ) +
    labs(
      x = "Patch",
      y = "Days since release starts"
    ) +
    theme_classic(base_size = 9)
  
  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  }
  
  ggdraw() +
    draw_plot(
      p,
      x = 0.04,
      y = 0,
      width = 0.90,
      height = 0.90
    ) +
    draw_label(
      panel_label,
      x = 0.35,
      y = 0.985,
      size = 10,
      fontface = "bold"
    ) +
    draw_label(
      "Source reduction",
      x = 0.50,
      y = 0.94,
      size = 10
    ) +
    draw_label(
      "Blocking effect",
      x = 0.785,
      y = 0.48,
      angle = 270,
      size = 10
    )
}

plot_initial_population_scenario <- function(
    initial_blocking_results,
    parameters,
    exclude_initial = 1.0,
    max_time = 365 * 10
) {
  
  library(dplyr)
  library(ggplot2)
  
  plot_df <- initial_blocking_results |>
    filter(time < max_time) |>
    mutate(
      initial_population = factor(
        initial_proportion,
        levels = c(
          0.10, 0.20, 0.30, 0.40, 0.50,
          0.60, 0.70, 0.80, 0.90, 1.00
        ),
        ordered = TRUE
      ),
      blocking = factor(
        blocking,
        levels = sort(unique(blocking))
      ),
      initial_proportion_label = paste0(
        "Initial population = ",
        initial_proportion
      )
    ) |>
    filter(
      !initial_proportion %in% exclude_initial
    )
  
  
  ggplot(
    plot_df,
    aes(
      x = time,
      y = VC_relative,
      group = interaction(
        blocking,
        initial_proportion_label
      )
    )
  ) +
    
    # Before Wolbachia release
    geom_line(
      data = plot_df |>
        filter(
          time <= parameters[["t_release_start"]]
        ),
      color = "black",
      linewidth = 0.8,
      alpha = 1
    ) +
    
    # After Wolbachia release
    geom_line(
      data = plot_df |>
        filter(
          time >= parameters[["t_release_start"]]
        ),
      aes(
        color = blocking
      ),
      linewidth = 0.8,
      alpha = 1
    ) +
    
    # Wolbachia release start
    geom_vline(
      xintercept = parameters[["t_release_start"]],
      linetype = "dashed",
      color = "red"
    ) +
    
    # Wolbachia release stop
    geom_vline(
      xintercept = parameters[["t_release_stop"]],
      linetype = "dashed",
      color = "red"
    ) +
    
    # Baseline VC
    geom_hline(
      yintercept = 1,
      linetype = "dashed",
      linewidth = 0.5
    ) +
    
    facet_wrap(
      ~ initial_proportion_label,
      nrow = 3,
      scales = "free_y"
    ) +
    
    labs(
      x = "Time",
      y = "Fold change vectorial capacity",
      color = "Blocking"
    ) +
    
    theme_classic(
      base_size = 8
    )
}


plot_initial_population_heatmap <- function(
    initial_blocking_results,
    parameters
) {
  
  library(dplyr)
  library(ggplot2)
  
  heatmap_df <- initial_blocking_results |>
    mutate(
      period = case_when(
        time < parameters[["t_release_start"]] ~ "pre-release",
        
        time >= parameters[["t_release_start"]] &
          time < parameters[["t_release_stop"]] ~ "release_period",
        
        time >= parameters[["t_release_stop"]] ~ "post-release"
      )
    ) |>
    
    filter(period == "post-release") |>
    
    group_by(
      initial_proportion,
      blocking
    ) |>
    
    summarise(
      VC_relative_mean = mean(
        VC_relative,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) |>
    
    mutate(
      initial_population = factor(
        initial_proportion,
        levels = c(
          0.10, 0.20, 0.30, 0.40, 0.50,
          0.60, 0.70, 0.80, 0.90, 1.00
        ),
        ordered = TRUE
      ),
      
      blocking = factor(
        blocking,
        levels = sort(unique(blocking))
      )
    )
  
  
  ggplot(
    heatmap_df,
    aes(
      x = blocking,
      y = initial_population,
      fill = VC_relative_mean
    )
  ) +
    
    geom_tile(
      color = "white",
      linewidth = 0.5
    ) +
    
    geom_text(
      aes(
        label = sprintf(
          "%.2f",
          VC_relative_mean
        )
      ),
      size = 3
    ) +
    
    scale_fill_gradient2(
      low = "darkgreen",
      mid = "white",
      high = "darkorange",
      midpoint = 1,
      limits = c(0, 4),
      oob = scales::squish,
      breaks = c(0, 1, 2, 3, 4),
      name = "Fold change in\nvectorial capacity"
    ) +
    
    labs(
      x = expression(
        italic(Wolbachia) ~ "blocking effect"
      ),
      y = "Initial population"
    ) +
    
    theme_classic(
      base_size = 12
    )
}

# helper function
make_heatmap <- function(data, title_text, x_lab_text) {
  
  ggplot(
    data,
    aes(
      y = factor(initial_proportion),
      x = factor(intervention_level),
      fill = VC_relative_mean
    )
  ) +
    geom_tile(
      color = "white",
      linewidth = 0.5
    ) +
    geom_text(
      aes(label = sprintf("%.2f", VC_relative_mean)),
      size = 3
    ) +
    scale_fill_gradient2(
      low = "darkgreen",
      mid = "white",
      high = "darkorange",
      midpoint = 1,
      limits = c(0, 4),
      oob = scales::squish,
      name = "Fold change\nvectorial capacity"
    ) +
    labs(
      y = "Initial population",
      x = x_lab_text,
      title = title_text
    ) +
    theme_classic(
      base_size = 10
    )
}
