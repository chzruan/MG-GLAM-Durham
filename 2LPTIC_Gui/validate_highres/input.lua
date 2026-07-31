------------------------------------------------------------
-- Validation vs Gui's 2048^3 FML IC (HEFT set on cosma8):
-- /cosma8/data/dp203/bl267/Projects/Ongoing/HEFT/ICs/IC_highres/
--   IC_Np1d_2048_L_1024_2LPT.{0..255}
-- Same config: L=1024, Nmesh=Npart=2048, seed 2026, z=49,
-- input spectrum pofk_bli_z49.txt (HEFT normalization,
-- x1.0181 in P above the DEGRACE low-res normalization).
------------------------------------------------------------
box_1d      = 1024.0
Nmesh       = 2048
Npart_1D    = 2048
random_seed = 2026

lpt_order = 2
z_ini     = 49.0

pofk_path_at_zini = "pofk_bli_z49.txt"

output_prefix = "snap/IC_ours_Np2048_L1024"

OmegaM      = 0.3
OmegaLambda = 1.0 - OmegaM
h           = 0.7

growth_rate1 = 1.0
growth_rate2 = 2.0 * growth_rate1

buffer_factor             = 1.25
interpolation_method      = "CIC"
density_assignment_method = "CIC"
interlacing               = true

fix_amplitude     = false
ic_reverse_phases = false

write_pofk = true
pofk_output_path = "snap/pofks_ours_2048.txt"

write_gadget = true
