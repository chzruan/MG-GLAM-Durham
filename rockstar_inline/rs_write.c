/* rs_write.c: publication of CATALOGS/Rockstar.<step>.<box>.bin + .meta.
 * Part of the MG-GLAM in-line Rockstar adapter; GNU GPL v3 or later.
 *
 * Both files are written as hidden temporaries in the catalogue directory and
 * renamed into place: first the old .meta is removed, then the .bin is renamed,
 * then the .meta, so a .meta always describes the .bin next to it and marks a
 * complete catalogue (an interrupted run leaves at most a .bin without .meta,
 * which check_run.py and glamio.pack refuse, and the rerun replaces). Skipped or
 * failed outputs publish an empty .bin and a .meta stating why. */
#define _POSIX_C_SOURCE 200809L
#include <errno.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include "rs_inline_internal.h"

struct column { const char *name, *type; size_t offset; const char *unit; };

#define H(f) offsetof(struct rs_record, h.f)
#define R(f) offsetof(struct rs_record, f)
static const struct column columns[] = {
    {"id", "<i8", H(id), "row number in this catalogue (= Rockstar halo id)"},
    {"x", "<f4", H(pos[0]), "Mpc/h comoving, in [0, Box]"},
    {"y", "<f4", H(pos[1]), "Mpc/h comoving, in [0, Box]"},
    {"z", "<f4", H(pos[2]), "Mpc/h comoving, in [0, Box]"},
    {"vx", "<f4", H(pos[3]), "km/s physical peculiar (Rockstar halo velocity, out_*.list VX)"},
    {"vy", "<f4", H(pos[4]), "km/s physical peculiar"},
    {"vz", "<f4", H(pos[5]), "km/s physical peculiar"},
    {"corevel_x", "<f4", H(corevel[0]), "km/s physical peculiar (core velocity)"},
    {"corevel_y", "<f4", H(corevel[1]), "km/s"},
    {"corevel_z", "<f4", H(corevel[2]), "km/s"},
    {"bulkvel_x", "<f4", H(bulkvel[0]), "km/s physical peculiar (bulk velocity)"},
    {"bulkvel_y", "<f4", H(bulkvel[1]), "km/s"},
    {"bulkvel_z", "<f4", H(bulkvel[2]), "km/s"},
    {"m", "<f4", H(m), "Msun/h, MASS_DEFINITION, all assigned particles (out_*.list Mvir_all)"},
    {"r", "<f4", H(r), "kpc/h comoving, radius of MASS_DEFINITION (from m)"},
    {"child_r", "<f4", H(child_r), "Rockstar internal"},
    {"vmax_r", "<f4", H(vmax_r), "Rockstar internal"},
    {"mgrav", "<f4", H(mgrav), "Msun/h, MASS_DEFINITION, bound particles (out_*.list Mvir when BOUND_PROPS = 1)"},
    {"vmax", "<f4", H(vmax), "km/s physical"},
    {"rvmax", "<f4", H(rvmax), "kpc/h comoving, radius of vmax"},
    {"rs", "<f4", H(rs), "kpc/h comoving, NFW scale radius"},
    {"klypin_rs", "<f4", H(klypin_rs), "kpc/h comoving, scale radius from vmax and mass"},
    {"vrms", "<f4", H(vrms), "km/s physical"},
    {"Jx", "<f4", H(J[0]), "(Msun/h) (Mpc/h) km/s physical"},
    {"Jy", "<f4", H(J[1]), "(Msun/h) (Mpc/h) km/s physical"},
    {"Jz", "<f4", H(J[2]), "(Msun/h) (Mpc/h) km/s physical"},
    {"energy", "<f4", H(energy), "Rockstar internal (total energy used for the spin)"},
    {"spin", "<f4", H(spin), "Peebles spin"},
    {"alt_m_1", "<f4", H(alt_m[0]), "Msun/h, MASS_DEFINITION2, all assigned particles"},
    {"alt_m_2", "<f4", H(alt_m[1]), "Msun/h, MASS_DEFINITION3"},
    {"alt_m_3", "<f4", H(alt_m[2]), "Msun/h, MASS_DEFINITION4"},
    {"alt_m_4", "<f4", H(alt_m[3]), "Msun/h, MASS_DEFINITION5"},
    {"Xoff", "<f4", H(Xoff), "kpc/h comoving"},
    {"Voff", "<f4", H(Voff), "km/s physical"},
    {"b_to_a", "<f4", H(b_to_a), "axis ratio (MASS_DEFINITION)"},
    {"c_to_a", "<f4", H(c_to_a), "axis ratio (MASS_DEFINITION)"},
    {"Ax", "<f4", H(A[0]), "kpc/h comoving, largest shape axis"},
    {"Ay", "<f4", H(A[1]), "kpc/h comoving"},
    {"Az", "<f4", H(A[2]), "kpc/h comoving"},
    {"b_to_a2", "<f4", H(b_to_a2), "axis ratio (500c)"},
    {"c_to_a2", "<f4", H(c_to_a2), "axis ratio (500c)"},
    {"A2x", "<f4", H(A2[0]), "kpc/h comoving, largest shape axis (500c)"},
    {"A2y", "<f4", H(A2[1]), "kpc/h comoving"},
    {"A2z", "<f4", H(A2[2]), "kpc/h comoving"},
    {"bullock_spin", "<f4", H(bullock_spin), "Bullock spin"},
    {"kin_to_pot", "<f4", H(kin_to_pot), "T/|U| (Newtonian potential)"},
    {"m_pe_b", "<f4", H(m_pe_b), "Msun/h, pseudo-evolution mass (Behroozi)"},
    {"m_pe_d", "<f4", H(m_pe_d), "Msun/h, pseudo-evolution mass (Diemer)"},
    {"halfmass_radius", "<f4", H(halfmass_radius), "kpc/h comoving"},
    {"num_p", "<i8", H(num_p), "particles assigned to the halo (out_*.list Np)"},
    {"num_child_particles", "<i8", H(num_child_particles), "Rockstar internal"},
    {"p_start", "<i8", H(p_start), "-1 (no member list is written)"},
    {"desc", "<i8", H(desc), "-1 (no merger tree)"},
    {"flags", "<i8", H(flags), "Rockstar internal"},
    {"n_core", "<i8", H(n_core), "Rockstar internal"},
    {"min_pos_err", "<f4", H(min_pos_err), "Rockstar internal"},
    {"min_vel_err", "<f4", H(min_vel_err), "Rockstar internal"},
    {"min_bulkvel_err", "<f4", H(min_bulkvel_err), "Rockstar internal"},
    {"inertia_xx", "<f4", H(inertia_tensor[0]), "inertia tensor (MASS_DEFINITION), as MPI-Rockstar OUTPUT_INERTIA_TENSOR"},
    {"inertia_yy", "<f4", H(inertia_tensor[1]), ""},
    {"inertia_zz", "<f4", H(inertia_tensor[2]), ""},
    {"inertia_xy", "<f4", H(inertia_tensor[3]), ""},
    {"inertia_yz", "<f4", H(inertia_tensor[4]), ""},
    {"inertia_zx", "<f4", H(inertia_tensor[5]), ""},
    {"inertia500c_xx", "<f4", H(inertia_tensor2[0]), "inertia tensor (500c)"},
    {"inertia500c_yy", "<f4", H(inertia_tensor2[1]), ""},
    {"inertia500c_zz", "<f4", H(inertia_tensor2[2]), ""},
    {"inertia500c_xy", "<f4", H(inertia_tensor2[3]), ""},
    {"inertia500c_yz", "<f4", H(inertia_tensor2[4]), ""},
    {"inertia500c_zx", "<f4", H(inertia_tensor2[5]), ""},
    {"chi2", "<f4", H(chi2), "NFW fit chi^2 (OUTPUT_NFW_CHI2)"},
    {"pid", "<i8", R(pid), "parent row as Rockstar's find_parents (radius r_strict with the strict pass, else r); -1: host"},
    {"upid", "<i8", R(upid), "root of the pid chain; -1: host"},
    {"ps_parent", "<i8", R(ps_parent), "Rockstar's phase-space parent row (sub_of); -1: none or not written"},
    {"fof_group", "<i8", R(fof_group), "canonical FoF group index of the halo"},
    {"m_strict", "<f4", R(m_strict), "Msun/h, MASS_DEFINITION, every particle inside (STRICT_SO_MASSES)"},
    {"r_strict", "<f4", R(r_strict), "kpc/h comoving, radius of m_strict"},
    {"alt_m_strict_1", "<f4", R(alt_m_strict[0]), "Msun/h, MASS_DEFINITION2, every particle inside"},
    {"alt_m_strict_2", "<f4", R(alt_m_strict[1]), "Msun/h, MASS_DEFINITION3"},
    {"alt_m_strict_3", "<f4", R(alt_m_strict[2]), "Msun/h, MASS_DEFINITION4"},
    {"alt_m_strict_4", "<f4", R(alt_m_strict[3]), "Msun/h, MASS_DEFINITION5"},
    {"strict_crossed", "<i4", R(strict_crossed), "bit 0: MASS_DEFINITION threshold reached by the strict pass; bits 1-4: the others"},
    {"pad", "<i4", R(pad), "0"},
};
#undef H
#undef R

int rs_columns_text(char *buf, size_t cap)
{
    size_t used = 0;
    int n = (int)(sizeof columns / sizeof columns[0]);
    int w = snprintf(buf, cap, "columns = %d\n", n);
    if (w < 0 || (size_t)w >= cap) return -1;
    used = (size_t)w;
    for (int i = 0; i < n; i++) {
        w = snprintf(buf + used, cap - used, "column_%02d = %s %s %zu %s\n", i + 1, columns[i].name,
                     columns[i].type, columns[i].offset, columns[i].unit);
        if (w < 0 || (size_t)w >= cap - used) return -1;
        used += (size_t)w;
    }
    return (int)used;
}

static int write_file(const char *path, const void *data, size_t bytes, char *err, int errlen)
{
    FILE *f = fopen(path, "wb");
    if (!f) {
        snprintf(err, errlen, "cannot create %s: %s", path, strerror(errno));
        return -1;
    }
    const char *c = data;
    size_t left = bytes;
    while (left > 0) {
        size_t chunk = left > (64u << 20) ? (64u << 20) : left;
        if (fwrite(c, 1, chunk, f) != chunk) {
            snprintf(err, errlen, "short write to %s: %s", path, strerror(errno));
            fclose(f);
            unlink(path);
            return -1;
        }
        c += chunk;
        left -= chunk;
    }
    if (fflush(f) != 0 || fsync(fileno(f)) != 0 || fclose(f) != 0) {
        snprintf(err, errlen, "cannot finish %s: %s", path, strerror(errno));
        unlink(path);
        return -1;
    }
    return 0;
}

int rs_write_catalogue(const char *dir, int step, int box, const struct rs_record *rec, int64_t nrows,
                       const char *meta, char *err, int errlen)
{
    char fin_bin[4096], fin_meta[4096], tmp_bin[4096], tmp_meta[4096];
    if (mkdir(dir, 0755) != 0 && errno != EEXIST) {
        snprintf(err, errlen, "cannot create %s: %s", dir, strerror(errno));
        return -1;
    }
    snprintf(fin_bin, sizeof fin_bin, "%s/Rockstar.%04d.%04d.bin", dir, step, box);
    snprintf(fin_meta, sizeof fin_meta, "%s/Rockstar.%04d.%04d.meta", dir, step, box);
    snprintf(tmp_bin, sizeof tmp_bin, "%s/.Rockstar.%04d.%04d.bin.tmp.%ld", dir, step, box, (long)getpid());
    snprintf(tmp_meta, sizeof tmp_meta, "%s/.Rockstar.%04d.%04d.meta.tmp.%ld", dir, step, box, (long)getpid());
    if (write_file(tmp_bin, rec, (size_t)nrows * sizeof(struct rs_record), err, errlen) != 0) return -1;
    if (write_file(tmp_meta, meta, strlen(meta), err, errlen) != 0) {
        unlink(tmp_bin);
        return -1;
    }
    if (unlink(fin_meta) != 0 && errno != ENOENT) {
        snprintf(err, errlen, "cannot remove the old %s: %s", fin_meta, strerror(errno));
        unlink(tmp_bin); unlink(tmp_meta);
        return -1;
    }
    if (rename(tmp_bin, fin_bin) != 0) {
        snprintf(err, errlen, "cannot rename to %s: %s", fin_bin, strerror(errno));
        unlink(tmp_bin); unlink(tmp_meta);
        return -1;
    }
    if (rename(tmp_meta, fin_meta) != 0) {
        snprintf(err, errlen, "cannot rename to %s: %s", fin_meta, strerror(errno));
        unlink(tmp_meta);
        return -1;
    }
    return 0;
}
