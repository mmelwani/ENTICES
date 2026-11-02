# k ~ 1e-14 m2 at phi0 = 0.32: the Kozeny-Carman grain diameter for it. Run to the 1e8 yr backstop (PI 2026-11-02):
# the ocean circulation time (~5e8 yr) is longer than t_end, so this case is not at steady state by design.
case <- list(inherits = "reference", name = "k_1e-14",
             grain_diameter_m = 5.040e-6,   # m; kozeny_carman_d(1e-14, 0.32)
             ss_tol = NULL,                 # no steady-state stop
             label = "not at steady state by design")
