/* rs_parents.c: which halos are distinct (the catalogue keeps only those), as
 * Rockstar's find_parents (MPI-Rockstar src/parents.c) decides it: halos in order
 * of decreasing vmax; each takes the halos whose centres lie inside its radius r
 * (periodic, capped at Box/2.01) and that have a smaller radius as children, a
 * later (lower-vmax) host overwriting an earlier one; a halo without a parent is
 * distinct. Unlike parents.c, halos of equal vmax are taken in row order, so the
 * result does not depend on the C library's qsort. The radius is the kernel's r
 * (that of mgrav when BOUND_PROPS = 1), as in Rockstar's out_*.list without
 * STRICT_SO_MASSES; the strict radius is not used, because a subhalo's strict
 * sphere also counts its host's particles.
 * Part of the MG-GLAM in-line Rockstar adapter; GNU GPL v3 or later. */
#include <stdint.h>
#include <stdlib.h>

#include "check_syscalls.h"
#include "config_vars.h"

struct rs_ph {
    float   pos[3];
    float   radius;      /* kpc/h comoving */
    float   vmax;
    int64_t id, parent;
};

#define FAST3TREE_TYPE   struct rs_ph
#define FAST3TREE_PREFIX RS_PARENTS
#include "fast3tree.c"

#include "rs_inline_internal.h"

/* the call's allocations, released by rs_parents_cleanup also after a Rockstar error (rs_guard.cpp) */
static struct rs_ph *ph = NULL;
static int64_t *order = NULL;
static struct fast3tree *tree = NULL;
static struct fast3tree_results *nearest = NULL;

void rs_parents_cleanup(void)
{
    if (nearest) fast3tree_results_free(nearest);
    if (tree) fast3tree_free(&tree);
    free(order);
    free(ph);
    nearest = NULL;
    tree = NULL;
    order = NULL;
    ph = NULL;
}

static int cmp_order(const void *a, const void *b)
{
    const struct rs_ph *c = ph + *(const int64_t *)a, *d = ph + *(const int64_t *)b;
    if (c->vmax != d->vmax) return c->vmax > d->vmax ? -1 : 1;
    return c->id < d->id ? -1 : (c->id > d->id);
}

int rs_find_parents(int64_t n, const struct halo *h, int64_t *pid)
{
    for (int64_t i = 0; i < n; i++) pid[i] = -1;
    if (n < 2) return 0;
    ph = malloc((size_t)n * sizeof(struct rs_ph));
    order = malloc((size_t)n * sizeof(int64_t));
    if (!ph || !order) {
        rs_parents_cleanup();
        return -1;
    }
    for (int64_t i = 0; i < n; i++) {
        for (int d = 0; d < 3; d++) ph[i].pos[d] = h[i].pos[d];
        ph[i].radius = h[i].r;
        ph[i].vmax = h[i].vmax;
        ph[i].id = i;
        ph[i].parent = -1;
    }
    tree = fast3tree_init(n, ph);                  /* reorders ph: rows are mapped back through id */
    _fast3tree_set_minmax(tree, 0, BOX_SIZE);
    for (int64_t i = 0; i < n; i++) order[i] = i;
    qsort(order, (size_t)n, sizeof(int64_t), cmp_order);
    const float max_dist = BOX_SIZE / 2.01;        /* parents.c:36 */
    nearest = fast3tree_results_init();
    for (int64_t i = 0; i < n; i++) {
        struct rs_ph *h1 = ph + order[i];
        float range = h1->radius * 1.0e-3;         /* kpc/h -> Mpc/h in double, as parents.c:63 */
        if (max_dist < range) range = max_dist;
        fast3tree_find_sphere_periodic(tree, nearest, h1->pos, range);
        for (int64_t j = 0; j < nearest->num_points; j++) {
            struct rs_ph *h2 = nearest->points[j];
            if (h2->radius < h1->radius) h2->parent = h1->id;
        }
    }
    for (int64_t i = 0; i < n; i++) pid[ph[i].id] = ph[i].parent;
    rs_parents_cleanup();
    return 0;
}
