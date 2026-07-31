------------------------------------------------------------
-- Validation vs DEGRACE low-res 2LPTic IC
-- /cosma7/data/dp004/bl267/Runs/DEGRACE/ICs/IC_data/L1024/Node_002/ics.*
-- Same config as its 2LPTic.param: L=1024, Nmesh=Nsample=1024,
-- seed 2026, z=49; input spectrum = the exact effective
-- P(k,z=49) of that run (PK_2LPTic_z00_002.dat shape renormalized
-- to sigma8=0.8000650, /38.9492^2; checked against inputspec_ics.txt
-- to 2e-5).
------------------------------------------------------------
box_1d      = 1024.0
Nmesh       = 1024
Npart_1D    = 1024
random_seed = 2026

lpt_order = 2
z_ini     = 49.0

pofk_path_at_zini = "pofk_lowres_exact_z49.txt"

output_prefix = "snap/IC_ours_Np1024_L1024"

OmegaM      = 0.3
OmegaLambda = 1.0 - OmegaM
h           = 0.7

-- f(z=49) = 1 to 1e-5 for this cosmology
growth_rate1 = 1.0
growth_rate2 = 2.0 * growth_rate1

buffer_factor             = 1.25
interpolation_method      = "CIC"
density_assignment_method = "CIC"
interlacing               = true

fix_amplitude     = false
ic_reverse_phases = false

write_pofk = true
pofk_output_path = "snap/pofks_ours_1024.txt"

write_gadget = true
