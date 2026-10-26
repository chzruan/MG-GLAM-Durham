/* rs_glue.c: the symbols Rockstar's per-group kernel needs from the parts of
 * MPI-Rockstar the in-line adapter does not link (rockstar.c, fun_times.c,
 * io/io_generic.c), and the assertion hook of shadow/assert.h.
 *
 * Part of the MG-GLAM in-line Rockstar adapter; GNU GPL v3 or later.
 * convert_and_sort_core_particles and sort_core_particles are copied verbatim
 * from MPI-Rockstar src/fun_times.c (commit 0ae386f; Rockstar (c) 2011 Peter
 * Behroozi, MPI-Rockstar (c) 2024 Tokuue et al., GPLv3). */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "particle.h"
#include "halo.h"
#include "groupies.h"
#include "error.h"
#include "io/io_generic.h"

/* rockstar.c: the particle array the kernel indexes; the adapter points it at its
 * array-of-structs buffer of FoF-group members before the kernel runs */
struct particle *p = NULL;

/* io/io_generic.c: the optional per-halo analysis hook (unused) */
AHG analyze_halo_generic = NULL;

/* fun_times.c, verbatim */
int sort_core_particles(const void *a, const void *b) {
    const struct particle *c = a;
    const struct particle *d = b;
    if (c->pos[0] < d->pos[0])
        return -1;
    if (c->pos[0] > d->pos[0])
        return 1;
    return 0;
}

void convert_and_sort_core_particles(struct halo *h, struct particle *hp,
                                     float max_r, int64_t *n_core) {
    int64_t i, j, np = h->num_p;
    float   ds, dx;
    for (i = 0; i < np; i++) {
        for (j = 0, dx = 0; j < 3; j++) {
            ds = h->pos[j] - hp[i].pos[j];
            dx += ds * ds;
        }
        hp[i].pos[0] = dx;
        if (max_r && (dx > max_r * max_r)) {
            struct particle tmp = hp[i];
            hp[i]               = hp[np - 1];
            hp[np - 1]          = tmp;
            np--;
            i--;
        }
    }
    qsort(hp, np, sizeof(struct particle), sort_core_particles);
    if (n_core)
        *n_core = np;
    for (i = 0; i < h->num_p; i++)
        hp[i].pos[0] = p[hp[i].id].pos[0];
}

/* fun_times.c: with temporal halo finding off no previous halos are loaded, and
 * Rockstar's function returns exactly this (its first test: !num_prev_halos) */
float find_previous_mass(struct halo *h, struct particle *hp, int64_t *best_num_p,
                         float max_r, const struct HaloInfo *haloinfo) {
    (void)h; (void)hp; (void)max_r; (void)haloinfo;
    *best_num_p = 0;
    return 0;
}

/* shadow/assert.h: a failed assertion takes Rockstar's error path (a throw) */
void rs_assert_fail(const char *expr, const char *file, int line) {
    fprintf(stderr, "[rockstar_inline] assertion failed: %s (%s:%d)\n", expr, file, line);
    rockstar_abort_impl(134, file, line);
}
