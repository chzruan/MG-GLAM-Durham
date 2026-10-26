/* rs_write.c: publication of CATALOGS/Rockstar.<step>.<box>.bin + .meta.
 * Part of the MG-GLAM in-line Rockstar adapter; GNU GPL v3 or later.
 *
 * Both files are written as hidden temporaries in the catalogue directory and
 * renamed into place: first the old .meta is removed, then the .bin is renamed,
 * then the .meta, so a .meta always describes the .bin next to it and marks a
 * complete catalogue. An interrupted run leaves a .bin without .meta and possibly
 * the temporaries .Rockstar.<step>.<box>.{bin,meta}.tmp; check_run.py and
 * glamio.pack refuse both, and a rerun of the step replaces both (the temporary
 * names are fixed). Skipped or failed outputs publish an empty .bin and a .meta
 * stating why (rs_inline.c keeps an existing complete catalogue instead). */
#define _POSIX_C_SOURCE 200809L
#include <errno.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include "rs_inline_internal.h"

#if !defined(__BYTE_ORDER__) || __BYTE_ORDER__ != __ORDER_LITTLE_ENDIAN__
#error "the catalogue is documented as little-endian records"
#endif

/* unit_unbound: the description when BOUND_PROPS = 0, for the columns Rockstar computes in its bound
 * pass when BOUND_PROPS = 1 (properties.c:457-490); NULL: the same either way */
struct column { const char *name, *type; size_t offset; const char *unit, *unit_unbound; };

#define H(f) offsetof(struct rs_record, h.f)
#define R(f) offsetof(struct rs_record, f)
static const struct column columns[] = {
    {"id", "<i8", H(id), "row number in this catalogue (= Rockstar halo id)", NULL},
    {"x", "<f4", H(pos[0]), "Mpc/h comoving, in [0, Box)", NULL},
    {"y", "<f4", H(pos[1]), "Mpc/h comoving, in [0, Box)", NULL},
    {"z", "<f4", H(pos[2]), "Mpc/h comoving, in [0, Box)", NULL},
    {"vx", "<f4", H(pos[3]), "km/s physical peculiar (Rockstar halo velocity, out_*.list VX)", NULL},
    {"vy", "<f4", H(pos[4]), "km/s physical peculiar", NULL},
    {"vz", "<f4", H(pos[5]), "km/s physical peculiar", NULL},
    {"corevel_x", "<f4", H(corevel[0]), "km/s physical peculiar (core velocity)", NULL},
    {"corevel_y", "<f4", H(corevel[1]), "km/s", NULL},
    {"corevel_z", "<f4", H(corevel[2]), "km/s", NULL},
    {"bulkvel_x", "<f4", H(bulkvel[0]), "km/s physical peculiar (bulk velocity)", NULL},
    {"bulkvel_y", "<f4", H(bulkvel[1]), "km/s", NULL},
    {"bulkvel_z", "<f4", H(bulkvel[2]), "km/s", NULL},
    {"m", "<f4", H(m), "Msun/h, MASS_DEFINITION, all assigned particles (out_*.list Mvir_all)", NULL},
    {"r", "<f4", H(r), "kpc/h comoving, radius of MASS_DEFINITION from the bound mass mgrav (BOUND_PROPS = 1)",
     "kpc/h comoving, radius of MASS_DEFINITION from m"},
    {"child_r", "<f4", H(child_r), "Rockstar internal", NULL},
    {"vmax_r", "<f4", H(vmax_r), "Rockstar internal", NULL},
    {"mgrav", "<f4", H(mgrav), "Msun/h, MASS_DEFINITION, bound particles (out_*.list Mvir when BOUND_PROPS = 1)", NULL},
    {"vmax", "<f4", H(vmax), "km/s physical", NULL},
    {"rvmax", "<f4", H(rvmax), "kpc/h comoving, radius of vmax", NULL},
    {"rs", "<f4", H(rs), "kpc/h comoving, NFW scale radius", NULL},
    {"klypin_rs", "<f4", H(klypin_rs), "kpc/h comoving, scale radius from vmax and mass", NULL},
    {"vrms", "<f4", H(vrms), "km/s physical", NULL},
    {"Jx", "<f4", H(J[0]), "(Msun/h) (Mpc/h) km/s physical", NULL},
    {"Jy", "<f4", H(J[1]), "(Msun/h) (Mpc/h) km/s physical", NULL},
    {"Jz", "<f4", H(J[2]), "(Msun/h) (Mpc/h) km/s physical", NULL},
    {"energy", "<f4", H(energy), "Rockstar internal (total energy used for the spin)", NULL},
    {"spin", "<f4", H(spin), "Peebles spin", NULL},
    {"alt_m_1", "<f4", H(alt_m[0]), "Msun/h, MASS_DEFINITION2, bound particles (BOUND_PROPS = 1); 0: not reached",
     "Msun/h, MASS_DEFINITION2, all assigned particles; 0: not reached"},
    {"alt_m_2", "<f4", H(alt_m[1]), "Msun/h, MASS_DEFINITION3, as alt_m_1", NULL},
    {"alt_m_3", "<f4", H(alt_m[2]), "Msun/h, MASS_DEFINITION4, as alt_m_1", NULL},
    {"alt_m_4", "<f4", H(alt_m[3]), "Msun/h, MASS_DEFINITION5, as alt_m_1", NULL},
    {"Xoff", "<f4", H(Xoff), "kpc/h comoving", NULL},
    {"Voff", "<f4", H(Voff), "km/s physical", NULL},
    {"b_to_a", "<f4", H(b_to_a), "axis ratio (MASS_DEFINITION)", NULL},
    {"c_to_a", "<f4", H(c_to_a), "axis ratio (MASS_DEFINITION)", NULL},
    {"Ax", "<f4", H(A[0]), "kpc/h comoving, largest shape axis", NULL},
    {"Ay", "<f4", H(A[1]), "kpc/h comoving", NULL},
    {"Az", "<f4", H(A[2]), "kpc/h comoving", NULL},
    {"b_to_a2", "<f4", H(b_to_a2), "axis ratio (MASS_DEFINITION4; 0 when alt_m_3 = 0)", NULL},
    {"c_to_a2", "<f4", H(c_to_a2), "axis ratio (MASS_DEFINITION4)", NULL},
    {"A2x", "<f4", H(A2[0]), "kpc/h comoving, largest shape axis (MASS_DEFINITION4)", NULL},
    {"A2y", "<f4", H(A2[1]), "kpc/h comoving", NULL},
    {"A2z", "<f4", H(A2[2]), "kpc/h comoving", NULL},
    {"bullock_spin", "<f4", H(bullock_spin), "Bullock spin", NULL},
    {"kin_to_pot", "<f4", H(kin_to_pot), "T/|U| (Newtonian potential)", NULL},
    {"m_pe_b", "<f4", H(m_pe_b), "Msun/h, pseudo-evolution mass (Behroozi)", NULL},
    {"m_pe_d", "<f4", H(m_pe_d), "Msun/h, pseudo-evolution mass (Diemer)", NULL},
    {"halfmass_radius", "<f4", H(halfmass_radius), "kpc/h comoving", NULL},
    {"num_p", "<i8", H(num_p), "particles assigned to the halo (out_*.list Np)", NULL},
    {"num_child_particles", "<i8", H(num_child_particles), "Rockstar internal", NULL},
    {"p_start", "<i8", H(p_start), "-1 (no member list is written)", NULL},
    {"desc", "<i8", H(desc), "-1 (no merger tree)", NULL},
    {"flags", "<i8", H(flags), "Rockstar internal", NULL},
    {"n_core", "<i8", H(n_core), "Rockstar internal", NULL},
    {"min_pos_err", "<f4", H(min_pos_err), "Rockstar internal", NULL},
    {"min_vel_err", "<f4", H(min_vel_err), "Rockstar internal", NULL},
    {"min_bulkvel_err", "<f4", H(min_bulkvel_err), "Rockstar internal", NULL},
    {"inertia_xx", "<f4", H(inertia_tensor[0]), "inertia tensor (MASS_DEFINITION), as MPI-Rockstar OUTPUT_INERTIA_TENSOR", NULL},
    {"inertia_yy", "<f4", H(inertia_tensor[1]), "", NULL},
    {"inertia_zz", "<f4", H(inertia_tensor[2]), "", NULL},
    {"inertia_xy", "<f4", H(inertia_tensor[3]), "", NULL},
    {"inertia_yz", "<f4", H(inertia_tensor[4]), "", NULL},
    {"inertia_zx", "<f4", H(inertia_tensor[5]), "", NULL},
    {"inertia500c_xx", "<f4", H(inertia_tensor2[0]), "inertia tensor (MASS_DEFINITION4, 500c by default)", NULL},
    {"inertia500c_yy", "<f4", H(inertia_tensor2[1]), "", NULL},
    {"inertia500c_zz", "<f4", H(inertia_tensor2[2]), "", NULL},
    {"inertia500c_xy", "<f4", H(inertia_tensor2[3]), "", NULL},
    {"inertia500c_yz", "<f4", H(inertia_tensor2[4]), "", NULL},
    {"inertia500c_zx", "<f4", H(inertia_tensor2[5]), "", NULL},
    {"chi2", "<f4", H(chi2), "NFW fit chi^2 (OUTPUT_NFW_CHI2)", NULL},
    {"pid", "<i8", R(pid), "parent row as Rockstar's find_parents with the radius r (ties in vmax by row); -1: host", NULL},
    {"upid", "<i8", R(upid), "root of the pid chain; -1: host", NULL},
    {"ps_parent", "<i8", R(ps_parent), "Rockstar's phase-space parent row (sub_of); -1: none or not written", NULL},
    {"fof_group", "<i8", R(fof_group), "canonical FoF group index of the halo", NULL},
    {"m_strict", "<f4", R(m_strict), "Msun/h, MASS_DEFINITION, every particle inside (STRICT_SO_MASSES)", NULL},
    {"r_strict", "<f4", R(r_strict), "kpc/h comoving, radius of m_strict", NULL},
    {"alt_m_strict_1", "<f4", R(alt_m_strict[0]), "Msun/h, MASS_DEFINITION2, every particle inside; 0: not reached", NULL},
    {"alt_m_strict_2", "<f4", R(alt_m_strict[1]), "Msun/h, MASS_DEFINITION3", NULL},
    {"alt_m_strict_3", "<f4", R(alt_m_strict[2]), "Msun/h, MASS_DEFINITION4", NULL},
    {"alt_m_strict_4", "<f4", R(alt_m_strict[3]), "Msun/h, MASS_DEFINITION5", NULL},
    {"strict_crossed", "<i4", R(strict_crossed), "bit 0: MASS_DEFINITION threshold reached by the strict pass; bits 1-4: the others", NULL},
    {"pad", "<i4", R(pad), "0", NULL},
};
#undef H
#undef R

int rs_columns_text(char *buf, size_t cap, int bound_props)
{
    size_t used = 0;
    int n = (int)(sizeof columns / sizeof columns[0]);
    int w = snprintf(buf, cap, "columns = %d\n", n);
    if (w < 0 || (size_t)w >= cap) return -1;
    used = (size_t)w;
    for (int i = 0; i < n; i++) {
        const char *unit = (!bound_props && columns[i].unit_unbound) ? columns[i].unit_unbound : columns[i].unit;
        w = snprintf(buf + used, cap - used, "column_%02d = %s %s %zu %s\n", i + 1, columns[i].name,
                     columns[i].type, columns[i].offset, unit);
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
    snprintf(tmp_bin, sizeof tmp_bin, "%s/.Rockstar.%04d.%04d.bin.tmp", dir, step, box);
    snprintf(tmp_meta, sizeof tmp_meta, "%s/.Rockstar.%04d.%04d.meta.tmp", dir, step, box);
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

int rs_catalogue_ok(const char *dir, int step, int box)
{
    char path[4096], line[512];
    snprintf(path, sizeof path, "%s/Rockstar.%04d.%04d.meta", dir, step, box);
    FILE *f = fopen(path, "r");
    int ok = 0;
    if (!f) return 0;
    while (fgets(line, sizeof line, f))
        if (!strncmp(line, "status = ", 9)) {
            ok = !strncmp(line + 9, "ok", 2) && (line[11] == '\n' || line[11] == 0);
            break;
        }
    fclose(f);
    return ok;
}
