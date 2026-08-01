------------------------------------------------------------
-- Simulation: DEGRACE high-res IC, L=1024 Mpc/h, 2048^3
-- Same seed (2026) as the 1024^3 2LPTic IC -> identical
-- large-scale phases. Copy to input.lua before running
-- (the code always reads "input.lua" from the cwd).
------------------------------------------------------------
box_1d      = 1024.0
Nmesh       = 2048
Npart_1D    = 2048
random_seed = 2026

------------------------------------------------------------
-- Initial conditions
------------------------------------------------------------
lpt_order = 2          -- 1 = Zel'dovich (1LPT), 2 = 2LPT
z_ini     = 49.0

------------------------------------------------------------
-- Input / Output
------------------------------------------------------------
pofk_path_at_zini = "pofk_bli_z49.txt"

output_prefix = "snap/IC_Np1d_2048_L_1024"

------------------------------------------------------------
-- Cosmology
------------------------------------------------------------
OmegaM      = 0.3
OmegaLambda = 1.0 - OmegaM
h           = 0.7

------------------------------------------------------------
-- Growth rates (enter velocities only, not positions)
-- f(z=49) = OmegaM(z)^0.55 = 0.99999 for this cosmology,
-- so the EdS value 1.0 is accurate to 1e-5 here.
------------------------------------------------------------
growth_rate1 = 1.0
growth_rate2 = 2.0 * growth_rate1

------------------------------------------------------------
-- Numerical parameters
------------------------------------------------------------
buffer_factor             = 1.25
interpolation_method      = "CIC"
density_assignment_method = "CIC"
interlacing               = true

------------------------------------------------------------
-- Random field
------------------------------------------------------------
fix_amplitude     = false
ic_reverse_phases = false

------------------------------------------------------------
-- Outputs
------------------------------------------------------------
write_pofk = true

pofk_output_path = "snap/pofks_L1024_Np2048.txt"

write_gadget = true
