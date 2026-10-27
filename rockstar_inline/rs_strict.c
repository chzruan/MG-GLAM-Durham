/* rs_strict.c: strict spherical-overdensity masses (Rockstar's STRICT_SO_MASSES)
 * of the in-line catalogue, from every particle of the box.
 * Part of the MG-GLAM in-line Rockstar adapter; GNU GPL v3 or later.
 *
 * As MPI-Rockstar's io_bgc2.c:279-365: for each written halo, all particles
 * within 1.1 max_halo_radius of the centre (periodic distances as
 * square_dist_from_center), sorted by distance; walking outwards with radii
 * floored at FORCE_RES, the enclosed mass is the strict mass of each definition
 * whenever the mean density exceeds its threshold, so the outermost crossing
 * wins; r is recomputed from the strict mass. A definition whose threshold is
 * never reached keeps the kernel's mass, as in Rockstar; strict_crossed says
 * which were reached. Particles come from GLAM's arrays through the FoF cell
 * index (perm + offsets) instead of a tree. */
#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include <string.h>

#include "config_vars.h"
#include "groupies.h"
#include "rs_inline_internal.h"

#define BGC2_R (1.1 / 1e3)      /* io_bgc2.h: search radius in units of max_halo_radius (kpc/h) */

static int cmp_float(const void *a, const void *b)
{
    float c = *(const float *)a, d = *(const float *)b;
    return c < d ? -1 : (c > d);
}

/* io_bgc2.c square_dist_from_center */
static inline float square_dist(const float a[3], const float c[3])
{
    float ds = 0, dx = 0;
    for (int i = 0; i < 3; i++) {
        dx = fabs(a[i] - c[i]);
        if (PERIODIC && (dx > (BOX_SIZE / 2.0)))
            dx = BOX_SIZE - dx;
        ds += dx * dx;
    }
    return ds;
}

int rs_strict_pass(int64_t n, struct rs_record *rec, const float *x, const float *y, const float *z,
                   const struct rs_units *u, const struct rs_cellindex *cells)
{
    int failed = 0;
    const int64_t nc = cells->nc;
    const double cs = cells->cs;
#pragma omp parallel copyin(particle_thresh_dens, particle_rvir_dens, particle_rvir_dens_z0, min_dens_index, \
                            dynamical_time, scale_dx)
    {
        double dens_thresh[5];
        for (int i = 0; i < 5; i++) dens_thresh[i] = particle_thresh_dens[i] * (4.0 * M_PI / 3.0) * PARTICLE_MASS;
        float *d2 = NULL;
        int64_t cap = 0;
#pragma omp for schedule(dynamic, 8)
        for (int64_t i = 0; i < n; i++) {
            if (failed) continue;
            struct halo *h = &rec[i].h;
            float R = BGC2_R * max_halo_radius(h);
            float R2 = R * R;
            int64_t lo[3], hi[3];
            for (int a = 0; a < 3; a++) {
                lo[a] = (int64_t)floor(((double)h->pos[a] - R) / cs);
                hi[a] = (int64_t)floor(((double)h->pos[a] + R) / cs);
                if (hi[a] - lo[a] + 1 >= nc) { lo[a] = 0; hi[a] = nc - 1; }
            }
            int64_t np = 0;
            for (int64_t ix = lo[0]; ix <= hi[0]; ix++)
            for (int64_t iy = lo[1]; iy <= hi[1]; iy++)
            for (int64_t iz = lo[2]; iz <= hi[2]; iz++) {
                int64_t c = ((((ix % nc) + nc) % nc) * nc + (((iy % nc) + nc) % nc)) * nc + (((iz % nc) + nc) % nc);
                for (int64_t k = cells->offsets[c]; k < cells->offsets[c + 1]; k++) {
                    int64_t pi = cells->perm[k];
                    float pp[3] = {rs_pos(u, x[pi]), rs_pos(u, y[pi]), rs_pos(u, z[pi])};
                    float ds = square_dist(pp, h->pos);
                    if (ds < R2) {
                        if (np == cap) {
                            int64_t nc2 = cap ? 2 * cap : 4096;
                            float *t = realloc(d2, (size_t)nc2 * sizeof(float));
                            if (!t) {
#pragma omp atomic write
                                failed = 1;
                                goto next_halo;
                            }
                            d2 = t;
                            cap = nc2;
                        }
                        d2[np++] = ds;
                    }
                }
            }
            qsort(d2, (size_t)np, sizeof(float), cmp_float);
            double total_mass = 0;
            float m = h->m, alt[4];
            int32_t crossed = 0;
            memcpy(alt, h->alt_m, sizeof alt);
            for (int64_t j = 0; j < np; j++) {
                float r = sqrt(d2[j]);
                if (r < FORCE_RES) r = FORCE_RES;
                total_mass += PARTICLE_MASS;
                float cur_dens = total_mass / (r * r * r);
                for (int k = 1; k < 5; k++)
                    if (cur_dens > dens_thresh[k]) {
                        alt[k - 1] = total_mass;
                        crossed |= 1 << k;
                    }
                if (cur_dens > dens_thresh[0]) {
                    m = total_mass;
                    crossed |= 1;
                }
            }
            rec[i].m_strict = m;
            rec[i].r_strict = cbrt(m / dens_thresh[0]) * 1e3;
            memcpy(rec[i].alt_m_strict, alt, sizeof alt);
            rec[i].strict_crossed = crossed;
        next_halo:;
        }
        free(d2);
    }
    return failed ? -1 : 0;
}
