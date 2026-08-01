#include <FML/ComputePowerSpectra/ComputePowerSpectrum.h>
#include <FML/FFTWGrid/FFTWGrid.h>
#include <FML/LPT/DisplacementFields.h>
#include <FML/MPIParticles/MPIParticles.h>
#include <FML/MemoryLogging/MemoryLogging.h>
#include <FML/ParticleTypes/SimpleParticle.h>
#include <FML/RandomFields/GaussianRandomField.h>
#include <FML/GadgetUtils/GadgetUtils.h>
#include <FML/Global/Global.h>
#include <FML/LuaFileParser/LuaFileParser.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

//=====================================================
// Simple 1LPT/2LPT IC generator controlled by a Lua file.
//
// Required Lua parameters:
//   box_1d       = 1024.0
//   random_seed  = 2026
//   Nmesh        = 2048
//   Npart_1D     = 2048
//   pofk_path    = "pofk_bli_z49.txt"
//   output_prefix = "/path/to/IC_Np1d_2048_L_1024_LPT"
//   lpt_order     = 2   -- 1 for Zel'dovich/1LPT, 2 for 2LPT
//
// Useful optional Lua parameters:
//   z_ini                     = 49.0
//   OmegaM                    = 0.3
//   h                         = 0.7
//   growth_rate1              = 0.99059529
//   growth_rate2              = 2.0 * growth_rate1   -- or set explicitly
//   buffer_factor             = 1.5
//   interpolation_method      = "CIC"
//   density_assignment_method = "PQS"
//   interlacing               = true
//   fix_amplitude             = false
//   ic_reverse_phases         = false
//   write_pofk                = true
//   pofk_output_path          = "/path/to/pofks.txt"
//   write_gadget              = true
//   lpt_order                 = 2
//=====================================================

constexpr int Ndim = 3;

template <class T>
using MPIParticles = FML::PARTICLE::MPIParticles<T>;
template <int N>
using FFTWGrid = FML::GRID::FFTWGrid<N>;
using Particle = HEFTParticle<Ndim>;
using RandomGenerator = FML::RANDOM::RandomGenerator;
using GSLRandomGenerator = FML::RANDOM::GSLRandomGenerator;
using GadgetWriter = FML::FILEUTILS::GADGET::GadgetWriter;

int main() {
#ifdef MEMORY_LOGGING
    auto * mem = FML::MemoryLog::get();
#endif

#ifndef USE_GSL
#error "This IC generator requires USE_GSL for GSLRandomGenerator. Compile with USE_GSL."
#endif

#ifdef GADGET_LONG_INT_IDS
    if (FML::ThisTask == 0)
        std::cout << "Using 64-bit Gadget IDs\n";
#else
    if (FML::ThisTask == 0)
        std::cout << "Using 32-bit Gadget IDs\n";
#endif

    const std::string parameter_file = std::string("input.lua");
    FML::FILEUTILS::LuaFileParser lfp(parameter_file);

    //============================================================
    // Read parameters from Lua
    //============================================================
    const double box = lfp.read_double("box_1d", 1024.0, lfp.required);
    const int random_seed = lfp.read_int("random_seed", 2026, lfp.required);
    const int Nmesh = lfp.read_int("Nmesh", 0, lfp.required);
    const int Npart_1D = lfp.read_int("Npart_1D", 0, lfp.required);
    const std::string pofk_path_at_zini = lfp.read_string("pofk_path_at_zini", "", lfp.required);
    const std::string output_prefix = lfp.read_string("output_prefix", "", lfp.required);
    const int lpt_order = lfp.read_int("lpt_order", 2, lfp.optional);
    const bool use_2LPT = (lpt_order == 2);

    const double z_ini = lfp.read_double("z_ini", 49.0, lfp.optional);
    const double a_ini = 1.0 / (1.0 + z_ini);
    const double OmegaM = lfp.read_double("OmegaM", 0.3, lfp.optional);
    const double OmegaLambda = lfp.read_double("OmegaLambda", 1.0 - OmegaM, lfp.optional);
    const double h = lfp.read_double("h", 0.7, lfp.optional);

    const double buffer_factor = lfp.read_double("buffer_factor", 1.5, lfp.optional);
    const std::string interpolation_method = lfp.read_string("interpolation_method", "CIC", lfp.optional);
    const std::string density_assignment_method = lfp.read_string("density_assignment_method", "PQS", lfp.optional);
    const bool interlacing = lfp.read_bool("interlacing", true, lfp.optional);

    const bool fix_amplitude = lfp.read_bool("fix_amplitude", false, lfp.optional);
    const bool ic_reverse_phases = lfp.read_bool("ic_reverse_phases", false, lfp.optional);

    const bool write_pofk = lfp.read_bool("write_pofk", true, lfp.optional);
    const std::string pofk_output_path = lfp.read_string("pofk_output_path", output_prefix + "_pofks.txt", lfp.optional);
    const bool write_gadget = lfp.read_bool("write_gadget", true, lfp.optional);

    const double HoverH0 = std::sqrt(OmegaM / (a_ini * a_ini * a_ini) + OmegaLambda);
    const double growth_rate1 = lfp.read_double("growth_rate1", 1.0, lfp.optional);
    const double growth_rate2 = lfp.read_double("growth_rate2", 2.0 * growth_rate1, lfp.optional);
    const double vfac_1LPT = a_ini * a_ini * HoverH0 * growth_rate1;
    const double vfac_2LPT = a_ini * a_ini * HoverH0 * growth_rate2;

    if (Nmesh <= 0 or Npart_1D <= 0)
        throw std::runtime_error("Nmesh and Npart_1D must be positive");
    if (box <= 0.0)
        throw std::runtime_error("box_1d must be positive");
    if (z_ini <= -1.0)
        throw std::runtime_error("z_ini must be larger than -1");
    if (lpt_order != 1 and lpt_order != 2)
        throw std::runtime_error("lpt_order must be either 1 or 2");

#ifndef GADGET_LONG_INT_IDS
    const long double npart_total_check = std::pow(static_cast<long double>(Npart_1D), static_cast<int>(Ndim));
    if (npart_total_check > static_cast<long double>(std::numeric_limits<unsigned int>::max()) and FML::ThisTask == 0) {
        std::cerr << "WARNING: Npart_1D^3 exceeds 32-bit Gadget IDs. "
                  << "Compile with GADGET_LONG_INT_IDS to avoid wrapped IDs.\n";
    }
#endif

    if (FML::ThisTask == 0) {
        std::cout << "#=====================================================================\n";
        std::cout << " Start of Lua parameters\n";
        std::cout << "#=====================================================================\n";
        std::cout << "parameter_file                      = " << parameter_file << "\n";
        std::cout << "box_1d                             = " << box << "\n";
        std::cout << "Nmesh                              = " << Nmesh << "\n";
        std::cout << "Npart_1D                           = " << Npart_1D << "\n";
        std::cout << "pofk_path_at_zini                  = " << pofk_path_at_zini << "\n";
        std::cout << "output_prefix                       = " << output_prefix << "\n";
        std::cout << "lpt_order                          = " << lpt_order << "\n";
        std::cout << "z_ini                              = " << z_ini << "\n";
        std::cout << "OmegaM                             = " << OmegaM << "\n";
        std::cout << "OmegaLambda                        = " << OmegaLambda << "\n";
        std::cout << "h                                  = " << h << "\n";
        std::cout << "interpolation_method               = " << interpolation_method << "\n";
        std::cout << "density_assignment_method          = " << density_assignment_method << "\n";
        std::cout << "interlacing                        = " << interlacing << "\n";
        std::cout << "fix_amplitude                       = " << fix_amplitude << "\n";
        std::cout << "ic_reverse_phases                  = " << ic_reverse_phases << "\n";
        std::cout << "#=====================================================================\n";
        std::cout << " End of Lua parameters\n";
        std::cout << "#=====================================================================\n";
    }

    //============================================================
    // Read P(k)
    //============================================================
    if (FML::ThisTask == 0) {
        std::cout << "#=====================================================\n";
        std::cout << " Reading and splining P(k)\n";
        std::cout << "#=====================================================\n";
    }

    std::ifstream fp(pofk_path_at_zini.c_str());
    if (!fp.is_open())
        throw std::runtime_error("Cannot open P(k) file: " + pofk_path_at_zini);

    std::vector<double> logk, logpofk;
    double previous_k = -1.0;
    for (;;) {
        double kin, pofkin;
        fp >> kin >> pofkin;
        if (fp.eof())
            break;
        if (!fp)
            throw std::runtime_error("Malformed line while reading P(k) file: " + pofk_path_at_zini);
        if (kin <= 0.0 or pofkin <= 0.0)
            throw std::runtime_error("P(k) file must contain positive k and P(k)");
        if (previous_k > 0.0 and kin <= previous_k)
            throw std::runtime_error("P(k) k-array must be strictly increasing for the spline");
        previous_k = kin;
        logk.push_back(std::log(kin));
        logpofk.push_back(std::log(pofkin));
    }
    if (logk.size() < 2)
        throw std::runtime_error("Need at least two P(k) samples to build a spline");

    FML::INTERPOLATION::SPLINE::Spline logpofk_spline(logk, logpofk);

    auto Pofk_of_kBox_over_volume = [&](double kBox) {
        const double k_phys = kBox / box;
        const double volume = std::pow(box, Ndim);
        return std::exp(logpofk_spline(std::log(k_phys))) / volume;
    };

    //============================================================
    // Set up random field
    //============================================================
    std::shared_ptr<RandomGenerator> rng = std::make_shared<GSLRandomGenerator>();
    rng->set_seed(random_seed);

    if (FML::ThisTask == 0) {
        std::cout << "#=====================================================\n";
        std::cout << " Creating delta_L_fourier\n";
        std::cout << "#=====================================================\n";
    }

    const auto nextra = FML::INTERPOLATION::get_extra_slices_needed_for_density_assignment(interpolation_method);
    FFTWGrid<Ndim> delta(Nmesh, nextra.first, nextra.second);

    FML::RANDOM::GAUSSIAN::generate_gaussian_random_field_fourier(
        delta,
        rng.get(),
        Pofk_of_kBox_over_volume,
        fix_amplitude);
    delta.set_grid_status_real(false);

    // Reverse phases? For pair-fixed ICs this maps delta(k) -> -delta(k).
    // We do it explicitly here instead of relying on the random-field helper,
    // so the Lua flag ic_reverse_phases has a well-defined effect.
    if (ic_reverse_phases) {
        if (FML::ThisTask == 0) {
            std::cout << "#=====================================================\n";
            std::cout << " Reversing Fourier phases: delta_L(k) -> -delta_L(k)\n";
            std::cout << "#=====================================================\n";
        }

        const int Local_nx = delta.get_local_nx();

#ifdef USE_OMP
#pragma omp parallel for
#endif
        for (int islice = 0; islice < Local_nx; islice++) {
            for (auto && fourier_index : delta.get_fourier_range(islice, islice + 1)) {
                delta.set_fourier_from_index(
                    fourier_index,
                    FML::GRID::FloatType(-1.0) * delta.get_fourier_from_index(fourier_index));
            }
        }
    }

    //=====================================================
    // Generate LPT potentials and displacement fields
    //=====================================================
    if (FML::ThisTask == 0) {
        std::cout << "#=====================================================\n";
        std::cout << " Computing " << lpt_order << "LPT displacement fields\n";
        std::cout << "#=====================================================\n";
    }

    FFTWGrid<Ndim> phi_1LPT;
    FML::COSMOLOGY::LPT::compute_1LPT_potential_fourier(delta, phi_1LPT);

    std::array<FFTWGrid<Ndim>, Ndim> Psi_1LPT_vector;
    FML::COSMOLOGY::LPT::from_LPT_potential_to_displacement_vector<Ndim>(phi_1LPT, Psi_1LPT_vector);
    phi_1LPT.free();

    std::array<FFTWGrid<Ndim>, Ndim> Psi_2LPT_vector;
    if (use_2LPT) {
        FFTWGrid<Ndim> phi_2LPT;
        FML::COSMOLOGY::LPT::compute_2LPT_potential_fourier(delta, phi_2LPT);
        FML::COSMOLOGY::LPT::from_LPT_potential_to_displacement_vector<Ndim>(phi_2LPT, Psi_2LPT_vector);
        phi_2LPT.free();
    }

    //=====================================================
    // Make particle grid and interpolate displacements
    //=====================================================
    MPIParticles<Particle> part;
    part.create_particle_grid(Npart_1D, buffer_factor, FML::xmin_domain, FML::xmax_domain);
    part.info();

    // create_particle_grid only sets pos: record the unperturbed lattice
    // position in q so the particle IDs computed below are correct
    {
        auto * pp = part.get_particles_ptr();
#ifdef USE_OMP
#pragma omp parallel for
#endif
        for (size_t ind = 0; ind < part.get_npart(); ind++) {
            auto * pos = pp[ind].get_pos();
            auto * q = pp[ind].get_q();
            for (int idim = 0; idim < Ndim; idim++)
                q[idim] = pos[idim];
        }
    }

    std::array<std::vector<FML::GRID::FloatType>, Ndim> displacements_1LPT;
    FML::INTERPOLATION::interpolate_grid_vector_to_particle_positions<Ndim, Particle>(
        Psi_1LPT_vector, part.get_particles_ptr(), part.get_npart(), displacements_1LPT, interpolation_method);
    for (int idim = 0; idim < Ndim; idim++)
        Psi_1LPT_vector[idim].free();

    std::array<std::vector<FML::GRID::FloatType>, Ndim> displacements_2LPT;
    if (use_2LPT) {
        FML::INTERPOLATION::interpolate_grid_vector_to_particle_positions<Ndim, Particle>(
            Psi_2LPT_vector, part.get_particles_ptr(), part.get_npart(), displacements_2LPT, interpolation_method);
        for (int idim = 0; idim < Ndim; idim++)
            Psi_2LPT_vector[idim].free();
    }

    //=====================================================
    // Add displacements and velocities
    //=====================================================
    double max_disp_1LPT = 0.0;
    double max_disp_2LPT = 0.0;
    double max_vel_1LPT = 0.0;
    double max_vel_2LPT = 0.0;

    auto * part_ptr = part.get_particles_ptr();
    const long long Np = static_cast<long long>(Npart_1D);

#ifdef USE_OMP
#pragma omp parallel for reduction(max : max_disp_1LPT, max_disp_2LPT, max_vel_1LPT, max_vel_2LPT)
#endif
    for (size_t ind = 0; ind < part.get_npart(); ind++) {
        auto * pos = part_ptr[ind].get_pos();
        auto * vel = part_ptr[ind].get_vel();
        auto * q = part_ptr[ind].get_q();

        long long ix = static_cast<long long>(std::floor(q[0] * Np));
        long long iy = static_cast<long long>(std::floor(q[1] * Np));
        long long iz = static_cast<long long>(std::floor(q[2] * Np));
        ix = std::min(std::max(ix, 0LL), Np - 1);
        iy = std::min(std::max(iy, 0LL), Np - 1);
        iz = std::min(std::max(iz, 0LL), Np - 1);

        const long long id = 1 + ix + Np * (iy + Np * iz);
        part_ptr[ind].set_id(id);

        for (int idim = 0; idim < Ndim; idim++) {
            const double dpos_1LPT = displacements_1LPT[idim][ind];
            const double dpos_2LPT = use_2LPT ? static_cast<double>(displacements_2LPT[idim][ind]) : 0.0;

            pos[idim] += dpos_1LPT + dpos_2LPT;
            vel[idim] = vfac_1LPT * dpos_1LPT + vfac_2LPT * dpos_2LPT;

            max_disp_1LPT = std::max(max_disp_1LPT, std::fabs(dpos_1LPT));
            max_disp_2LPT = std::max(max_disp_2LPT, std::fabs(dpos_2LPT));
            max_vel_1LPT = std::max(max_vel_1LPT, std::fabs(vfac_1LPT * dpos_1LPT));
            max_vel_2LPT = std::max(max_vel_2LPT, std::fabs(vfac_2LPT * dpos_2LPT));

            while (pos[idim] >= 1.0)
                pos[idim] -= 1.0;
            while (pos[idim] < 0.0)
                pos[idim] += 1.0;
        }
    }

    for (int idim = 0; idim < Ndim; idim++) {
        displacements_1LPT[idim].clear();
        displacements_1LPT[idim].shrink_to_fit();
        if (use_2LPT) {
            displacements_2LPT[idim].clear();
            displacements_2LPT[idim].shrink_to_fit();
        }
    }

    FML::MaxOverTasks(&max_disp_1LPT);
    FML::MaxOverTasks(&max_disp_2LPT);
    FML::MaxOverTasks(&max_vel_1LPT);
    FML::MaxOverTasks(&max_vel_2LPT);

    if (FML::ThisTask == 0) {
        std::cout << "Maximum 1LPT displacement: " << max_disp_1LPT * Nmesh << " grid cells\n";
        if (use_2LPT)
            std::cout << "Maximum 2LPT displacement: " << max_disp_2LPT * Nmesh << " grid cells\n";
        std::cout << "Maximum 1LPT velocity: " << max_vel_1LPT * 100.0 * box / a_ini << " km/s peculiar\n";
        if (use_2LPT)
            std::cout << "Maximum 2LPT velocity: " << max_vel_2LPT * 100.0 * box / a_ini << " km/s peculiar\n";
    }

    part.communicate_particles();

    //=====================================================
    // Power-spectrum check
    //=====================================================
    if (write_pofk) {
        if (FML::ThisTask == 0) {
            std::cout << "#=====================================================\n";
            std::cout << " Writing pofks.txt\n";
            std::cout << "#=====================================================\n";
        }
        const auto nlr = FML::INTERPOLATION::get_extra_slices_needed_for_density_assignment(density_assignment_method);
        FFTWGrid<Ndim> delta_from_part(Nmesh, nlr.first, nlr.second);
        FML::CORRELATIONFUNCTIONS::PowerSpectrumBinning<Ndim> pofk_cross(Nmesh / 2);
        FML::CORRELATIONFUNCTIONS::PowerSpectrumBinning<Ndim> pofk_part(Nmesh / 2);
        FML::CORRELATIONFUNCTIONS::PowerSpectrumBinning<Ndim> pofk_ini(Nmesh / 2);

        FML::INTERPOLATION::particles_to_fourier_grid(part.get_particles_ptr(),
                                                      part.get_npart(),
                                                      part.get_npart_total(),
                                                      delta_from_part,
                                                      density_assignment_method,
                                                      interlacing);
        FML::INTERPOLATION::deconvolve_window_function_fourier<Ndim>(delta_from_part, density_assignment_method);
        FML::CORRELATIONFUNCTIONS::bin_up_power_spectrum(delta_from_part, pofk_part);
        FML::CORRELATIONFUNCTIONS::bin_up_power_spectrum(delta, pofk_ini);
        FML::CORRELATIONFUNCTIONS::bin_up_cross_power_spectrum(delta, delta_from_part, pofk_cross);

        pofk_part.scale(box);
        pofk_ini.scale(box);
        pofk_cross.scale(box);

        if (FML::ThisTask == 0) {
            std::ofstream fout(pofk_output_path.c_str());
            if (!fout.is_open())
                throw std::runtime_error("Cannot open pofk output file: " + pofk_output_path);
            fout << "# k[h/Mpc] P_linear(k,zini)[(Mpc/h)^3] P_particle_" << lpt_order << "LPT(k,zini)[(Mpc/h)^3] r_cross P_ini_over_input\n";
            for (int i = 0; i < pofk_ini.n; i++) {
                const double k = pofk_ini.kbin[i];
                const double r_cross = pofk_cross.pofk[i] / std::sqrt(pofk_part.pofk[i] * pofk_ini.pofk[i]);
                const double p_ini_over_input = pofk_ini.pofk[i] / std::exp(logpofk_spline(std::log(k)));
                fout << k << " " << pofk_ini.pofk[i] << " " << pofk_part.pofk[i] << " "
                     << r_cross << " " << p_ini_over_input << "\n";
            }
        }
    }

    //=========================================================
    // Write Gadget files
    //=========================================================
    if (write_gadget) {
        if (FML::ThisTask == 0) {
            std::cout << "#=====================================================\n";
            std::cout << " Writing Gadget2 snapshot\n";
            std::cout << "#=====================================================\n";
        }
        GadgetWriter gw;
        const int nfiles = FML::NTasks;
        const double pos_norm = box;
        // Gadget-1 stores u = v_peculiar/sqrt(a) in km/s; internal velocities are
        // a^2 (H/H0) f Psi in box units, so v_pec[km/s] = vel * 100 * box / a
        // (same convention as COLASolver's output_gadget)
        const double vel_norm = 100.0 * box / (a_ini * std::sqrt(a_ini));
        gw.write_gadget_single(output_prefix + "." + std::to_string(FML::ThisTask),
                               part.get_particles_ptr(),
                               part.get_npart(),
                               part.get_npart_total(),
                               nfiles,
                               a_ini,
                               box,
                               OmegaM,
                               OmegaLambda,
                               h,
                               pos_norm,
                               vel_norm);
    }

    part.free();

#ifdef MEMORY_LOGGING
    mem->print();
#endif
}
