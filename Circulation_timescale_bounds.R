# Function to calculate circulation time bounds for your parameter space
calculate_circulation_bounds <- function() {
  
  # Define parameter ranges
  k_range <- c(1e-14, 1e-12, 1e-10, 1e-8)  # Extended range including extremes
  layer_thickness_range <- c(10, 100, 1000, 10000, 194200)  # Including full core
  
  # Store results
  results <- data.frame(
    k = numeric(),
    layer_thickness = numeric(),
    circulation_time = numeric()
  )
  
  config <- EnceladusConfig$new()
  
  # Test all combinations
  for(k in k_range) {
    for(thickness in layer_thickness_range) {
      # Create and run simulation
      sim <- Simulator$new(k = k, layer_thickness = thickness, config = config)
      sim$run_simulation(timesteps = 100)  # Reduced for speed
      
      # Get circulation time
      metrics <- sim$get_summary_metrics()
      circ_time <- metrics$circulation_time_years
      
      results <- rbind(results, data.frame(
        k = k,
        layer_thickness = thickness,
        circulation_time = circ_time
      ))
      
      cat(sprintf("k=%.0e, thickness=%g m: %.1f years\n", k, thickness, circ_time))
    }
  }
  
  # Find bounds
  min_circulation <- min(results$circulation_time, na.rm = TRUE)
  max_circulation <- max(results$circulation_time, na.rm = TRUE)
  
  cat(sprintf("\nCirculation time bounds:\n"))
  cat(sprintf("Minimum: %.1f years\n", min_circulation))
  cat(sprintf("Maximum: %.1f years\n", max_circulation))
  
  return(list(
    min_time = min_circulation,
    max_time = max_circulation,
    results = results
  ))
}

# Run the analysis
bounds <- calculate_circulation_bounds()

# Create the adaptive grouping function
create_adaptive_grouping <- function(min_circ_time, max_circ_time, 
                                     min_layer_frac = 0.001, max_layer_frac = 1.0) {
  
  function(layer_fraction, circulation_time_years) {
    # Normalize both parameters to [0,1], with clamping
    layer_norm <- pmax(0, pmin(1, (layer_fraction - min_layer_frac) / (max_layer_frac - min_layer_frac)))
    time_norm <- pmax(0, pmin(1, (circulation_time_years - min_circ_time) / (max_circ_time - min_circ_time)))
    
    # Bilinear scaling: 1 year (min) to 1000 years (max) time steps
    target_time_step <- 1 + (layer_norm * time_norm * 999)
    
    # Calculate required grouping
    grouping <- max(1, floor(circulation_time_years / target_time_step))
    
    return(grouping)
  }
}

# Use the empirically determined bounds
adaptive_grouping_fn <- create_adaptive_grouping(bounds$min_time, bounds$max_time)
