sweep_results <- list()

for (i in 1:nrow(parameter_sweep)) {
  print(sprintf("\nRunning scenario %d of %d:", i, nrow(parameter_sweep)))
  print(parameter_sweep[i, ])

  tryCatch({
    sim_i <- Simulator$new(
      k = parameter_sweep$k[i],
      layer_thickness = parameter_sweep$layer_thickness[i],
      config = enceladus
    )

    sim_i$run_simulation()
    metrics_i <- sim_i$get_summary_metrics()

    phreeqc_i    <- PhreeqcIntegrator$new(sim_i)
    phreeqc_file <- NA
    if (generate_phreeqc) {
      phreeqc_file <- phreeqc_i$generate_phreeqc_input_tidal_kinetic(
        mode             = cycle_mode,
        n_circulations   = n_circulations,
        n_orbital_cycles = n_orbital_cycles
      )
    }

    sweep_results[[i]] <- list(
      params = data.frame(
        k = parameter_sweep$k[i],
        layer_thickness = parameter_sweep$layer_thickness[i]
      ),
      metrics = data.frame(
        mean_flow_rate_Ls         = metrics_i$mean_flow_rate_Ls,
        circulation_time_years    = metrics_i$circulation_time_years,
        complete_circulation_time = metrics_i$complete_circulation_time
      ),
      phreeqc_params = phreeqc_i$phreeqc_params,
      phreeqc_file   = phreeqc_file
    )
  }, error = function(e) {
    warning(sprintf("Error in scenario %d: %s", i, e$message))
    sweep_results[[i]] <<- list(
      params = data.frame(
        k = parameter_sweep$k[i],
        layer_thickness = parameter_sweep$layer_thickness[i]
      ),
      metrics = data.frame(
        mean_flow_rate_Ls         = NA,
        circulation_time_years    = NA,
        complete_circulation_time = NA
      ),
      phreeqc_params = list(),
      phreeqc_file   = NA
    )
  })
}
