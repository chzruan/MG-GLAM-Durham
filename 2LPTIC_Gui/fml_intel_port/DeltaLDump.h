#pragma once
//=====================================================================
// delta_L dump for MG-GLAM (item I7 of docs/degrace_companion_implementation_plan.md).
//
// Writes the low-k modes of the Gaussian field that seeds the ICs,
// normalised to linear theory at a = 1, in the exact layout and
// convention of MG-GLAM's deltam/mom_* products so that the two can
// be cross-correlated mode by mode:
//   delta_k = N^-3 sum_x delta(x) exp(-i k.x), k = (2 pi/Box)(ix,iy,iz),
//   grid value j at x = j Box/N (no cell-centre shift), cube
//   |ix|,|iy| <= n, 0 <= iz <= n, n = floor(kc Box/2pi), stored as
//   complex<float> C[ix+n][iy+n][iz] (iz fastest), k = 0 set to zero,
//   raw little-endian, no record markers, plus a key = value .meta file.
// FML's FFTWGrid already uses this normalisation (fftw_r2c divides by
// N^3) and its half-complex axis is z, so the copy is index-exact.
//=====================================================================
#include <FML/FFTWGrid/FFTWGrid.h>
#include <FML/Global/Global.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef USE_MPI
#include <mpi.h>
#endif

namespace DELTAL {

    // Linear growth factor of the flat matter+Lambda background, D(1) = 1
    // (Simpson in ln a; the same integral as MG-GLAM's GrowthFactor).
    inline double growth_raw(double a, double Om, double OL) {
        const int nint = 4000;
        const double x0 = std::log(1e-8);
        const double h = (std::log(a) - x0) / nint;
        double s = 0.0;
        for (int i = 0; i <= nint; i++) {
            const double x = x0 + h * i;
            const double f = std::exp(2.5 * x) / std::pow(Om + OL * std::exp(3.0 * x), 1.5);
            s += (i == 0 || i == nint) ? f : ((i % 2) ? 4.0 * f : 2.0 * f);
        }
        return std::sqrt(Om / (a * a * a) + OL) * s * h / 3.0;
    }
    inline double growth(double a, double Om, double OL) { return growth_raw(a, Om, OL) / growth_raw(1.0, Om, OL); }

    struct Meta {
        std::string pofk_file, generator, lattice;
        int seed = 0, lpt_order = 0, Npart_1D = 0;
        bool fix_amplitude = false, reverse_phases = false;
        double h = 0.0, OmegaM = 0.0, OmegaLambda = 0.0, z_ini = 0.0;
    };

    // delta: the Gaussian field in Fourier space at a_ini (FML convention).
    template <int N>
    void write_cube(FML::GRID::FFTWGrid<N> & delta, double box, double a_ini, double kc,
                    const std::string & prefix, const Meta & m) {
        static_assert(N == 3, "delta_L dump is written for 3 dimensions");
        const int Nmesh = delta.get_nmesh();
        int n = int(kc * box / (2.0 * M_PI));
        if (n > Nmesh / 2 - 1) {
            if (FML::ThisTask == 0)
                std::cout << "[deltaL] n = " << n << " limited to Nmesh/2-1 = " << Nmesh / 2 - 1 << "\n";
            n = Nmesh / 2 - 1;
        }
        const long long n1 = 2LL * n + 1, nz = n + 1;
        const long long ncube = n1 * n1 * nz;
        // D is computed identically on every task, so each task can normalise its own slab and the
        // gather below is a sum with exactly one non-zero contribution per element: exact in float.
        const double D = growth(a_ini, m.OmegaM, m.OmegaLambda);
        std::vector<std::complex<float>> cube(ncube, std::complex<float>(0.0f, 0.0f));   // 8 B per mode

        const int local_nx = delta.get_local_nx();
        const long long x_start = delta.get_local_x_start();
        const long long nover2p1 = Nmesh / 2 + 1;
        for (int ixl = 0; ixl < local_nx; ixl++) {
            const long long ix = x_start + ixl;
            const long long kx = (ix <= Nmesh / 2) ? ix : ix - Nmesh;
            if (kx < -n || kx > n) continue;
            for (long long iy = 0; iy < Nmesh; iy++) {
                const long long ky = (iy <= Nmesh / 2) ? iy : iy - Nmesh;
                if (ky < -n || ky > n) continue;
                for (long long iz = 0; iz <= n; iz++) {
                    const long long index = (static_cast<long long>(ixl) * Nmesh + iy) * nover2p1 + iz;
                    const auto v = delta.get_fourier_from_index(index);
                    cube[iz + nz * ((ky + n) + n1 * (kx + n))] =
                        std::complex<float>(float(double(v.real()) / D), float(double(v.imag()) / D));
                }
            }
        }
#ifdef USE_MPI
        // chunked so that the int element count cannot overflow for large k_c
        float * buf = reinterpret_cast<float *>(cube.data());
        const long long ntot = 2 * ncube, chunk = 1LL << 28;
        for (long long off = 0; off < ntot; off += chunk) {
            const int cnt = int(std::min(chunk, ntot - off));
            if (FML::ThisTask == 0)
                MPI_Reduce(MPI_IN_PLACE, buf + off, cnt, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
            else
                MPI_Reduce(buf + off, nullptr, cnt, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
        }
#endif
        if (FML::ThisTask != 0) return;
        cube[0 + nz * ((0 + n) + n1 * (0 + n))] = std::complex<float>(0.0f, 0.0f);   // k = 0
        double sum2 = 0.0;
        for (long long i = 0; i < ncube; i++) sum2 += std::norm(std::complex<double>(cube[i]));
        const auto & out = cube;

        const std::string fbin = prefix + "_deltaL.bin", fmeta = prefix + "_deltaL.meta";
        std::ofstream fb(fbin, std::ios::binary);
        if (!fb) throw std::runtime_error("cannot open " + fbin);
        fb.write(reinterpret_cast<const char *>(out.data()), out.size() * sizeof(std::complex<float>));
        fb.close();

        std::ofstream fm(fmeta);
        if (!fm) throw std::runtime_error("cannot open " + fmeta);
        fm << std::setprecision(10);
        fm << "product = deltaL\n";
        fm << "format = raw little-endian complex64 (re, im float32), no record markers\n";
        fm << "layout = C(0:n,-n:n,-n:n) with the first index iz; numpy: fromfile(dtype=complex64).reshape(2n+1,2n+1,n+1) -> [ix+n, iy+n, iz]\n";
        fm << "fourier_norm = delta_k = N^-3 sum_x delta(x) exp(-i k.x), k = (2 pi/Box)(ix,iy,iz), P(k) = Box^3 <|delta_k|^2>; grid value j at x = j Box/N\n";
        fm << "modes = cube |ix|,|iy| <= n, 0 <= iz <= n (Hermitian-redundant iz = 0 modes kept); k = 0 set to zero\n";
        fm << "normalisation = Gaussian field drawn at z_ini from P_lin(k, z_ini) divided by D(a_ini): linear theory at a = 1 (flat matter+Lambda growth, D(1) = 1)\n";
        fm << "n = " << n << "\n";
        fm << "k_c = " << kc << "\n";
        fm << "N_mesh = " << Nmesh << "\n";
        fm << "N_dep = " << Nmesh << "\n";
        fm << "N_part_1D = " << m.Npart_1D << "\n";
        fm << "Box = " << box << "\n";
        fm << "units = Mpc/h\n";
        fm << "z_ini = " << m.z_ini << "\n";
        fm << "a_ini = " << a_ini << "\n";
        fm << "D_gen(a_ini) = " << D << "\n";
        fm << "Omega_m = " << m.OmegaM << "\n";
        fm << "Omega_L = " << m.OmegaLambda << "\n";
        fm << "hubble = " << m.h << "\n";
        fm << "seed = " << m.seed << "\n";
        fm << "lpt_order = " << m.lpt_order << "\n";
        fm << "fix_amplitude = " << (m.fix_amplitude ? "true" : "false") << "\n";
        fm << "reverse_phases = " << (m.reverse_phases ? "true" : "false") << "\n";
        fm << "pofk_file = " << m.pofk_file << "\n";
        fm << "lattice = " << m.lattice << "\n";
        fm << "axis_map = ix, iy, iz follow the Gadget x, y, z of the IC particles\n";
        fm << "generator = " << m.generator << "\n";
        fm << "rms_deltaL_cube = " << std::sqrt(sum2 / double(ncube)) << "\n";
        fm.close();
        std::cout << "[deltaL] wrote " << fbin << ": n = " << n << ", " << ncube << " modes, D(a_ini) = " << D
                  << ", rms |delta_L| in the cube = " << std::sqrt(sum2 / double(ncube)) << "\n";
    }
}   // namespace DELTAL
