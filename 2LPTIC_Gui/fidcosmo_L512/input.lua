------------------------------------------------------------
-- Fiducial-cosmology 2LPTIC IC to compare with the GLAM
-- fid_LCDM_L512Np2048Ng4096 Run1-5 suite (same box, Np, and
-- input P(k); differs only in IC method: 2LPT @ z=49 vs
-- GLAM ZA @ z=100, and in realization seed).
-- Input: fid PkTable.dat (CAMB z=0, sigma8=0.8159) rescaled
-- by D(z=49)/D(0) = 1/39.2136 for Om=0.3089/OL=0.6911,
-- with GLAM's own k^-3 extension to k=100.
------------------------------------------------------------
box_1d      = 512.0
Nmesh       = 2048
Npart_1D    = 2048
random_seed = 2026

lpt_order = 2
z_ini     = 49.0

pofk_path_at_zini = "pofk_fid_z49.txt"

output_prefix = "snap/IC_fid_Np2048_L512"

OmegaM      = 0.3089
OmegaLambda = 1.0 - OmegaM
h           = 0.677

-- f(z=49) = 0.99999 for this cosmology
growth_rate1 = 1.0
growth_rate2 = 2.0 * growth_rate1

buffer_factor             = 1.25
interpolation_method      = "CIC"
density_assignment_method = "CIC"
interlacing               = true

fix_amplitude     = false
ic_reverse_phases = false

write_pofk = true
pofk_output_path = "snap/pofks_fid_2048.txt"

write_gadget = true
