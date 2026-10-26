/* rs_fof.h: parallel exact periodic friends-of-friends of the in-line Rockstar
 * adapter (the parent groups handed to Rockstar's per-group kernel).
 *
 * Part of the MG-GLAM in-line Rockstar adapter; GNU GPL v3 or later, see
 * rs_convert.h.
 *
 * Linking (as Rockstar's own test, rockstar.c:59, fast3tree.c:405): particles
 * are friends when the periodic minimum-image distance of their float positions
 * (rs_pos of rs_convert.h) satisfies d^2 < r^2 in float, r = (float)(spacing*b).
 * The result is exact (every such pair is linked; Rockstar's FoF in contrast
 * merges groups within 2r around particles with > 16 neighbours) and does not
 * depend on the thread count or scheduling: union-find roots are the smallest
 * sorted index of each group, and the groups come in a canonical order (members
 * by particle index, groups by size descending, then smallest particle index).
 *
 * Particles are sorted by coarse cell (about Nrow/2 per side) and, inside each
 * coarse cell, by fine cell (side >= r (1 + 1e-3), so linked pairs always sit in
 * the same or adjacent fine cells even after float rounding of wrapped
 * differences). The coarse-cell index (perm + cell offsets) is kept after the FoF
 * for the strict spherical-overdensity pass. */
#ifndef RS_FOF_H
#define RS_FOF_H

#include <stdint.h>

#include "rs_convert.h"

#ifdef __cplusplus
extern "C" {
#endif

/* the shared distance test (rs_fof.c and the serial reference fof_ref.c) */
static inline float rs_dwrap(float d, float box, float half)
{
    if (d > half) d -= box;
    else if (d < -half) d += box;
    return d;
}

/* (dx*dx + dy*dy) + dz*dz, each operation rounded to float: the files that use it are compiled with
 * -ffp-contract=off, so no fused multiply-add changes the rounding */
static inline int rs_linked(const float *a, const float *b, float box, float half, float r2)
{
    float dx = rs_dwrap(a[0] - b[0], box, half);
    float dy = rs_dwrap(a[1] - b[1], box, half);
    float dz = rs_dwrap(a[2] - b[2], box, half);
    return dx * dx + dy * dy + dz * dz < r2;
}

#define RS_FOF_CELL_MARGIN 1e-3     /* cell sides >= r (1 + margin), and >= r + 4 ulp(Box) (checked) */

/* return codes of rs_fof_run */
enum { RS_FOF_OK = 0, RS_FOF_ENOMEM = -1, RS_FOF_EGRID = -2, RS_FOF_ECOORD = -3, RS_FOF_EGROUPS = -4 };
const char *rs_fof_strerror(int rc);

struct rs_cellindex {               /* particles sorted by coarse cell, kept for the strict SO pass */
    int64_t  nc;                    /* coarse cells per side */
    double   cs;                    /* coarse cell side, Mpc/h */
    int64_t *perm;                  /* sorted position k -> particle index (n entries) */
    int64_t *offsets;               /* nc^3 + 1 entries: particles of cell c are perm[offsets[c]..offsets[c+1]) */
};

struct rs_fof_groups {
    int64_t  ngroups;               /* groups with >= nmin members, canonical order */
    int64_t  nmembers;              /* particles in them */
    int64_t  largest;
    int64_t *offset;                /* ngroups + 1: canonical member layout */
    int64_t *member;                /* nmembers particle indices, group by group, each group by index
                                       (NULL unless requested: the adapter scatters straight into its buffer) */
    int32_t *group_of_sorted;       /* n entries (sorted position k -> group or -1; fewer than 2^31 groups,
                                       checked); NULL after rs_fof_release_labels */
    int64_t  nbad;                  /* particles whose coordinates are not finite or more than one box outside
                                       [1, Ngrid+1) (rs_fof_run then returns RS_FOF_ECOORD) */
    double   rss_link_gib;          /* VmRSS sampled when every array of the linking stage is in use */
};

/* Runs the FoF on n particles (GLAM grid units, converted by u) with the OpenMP
 * default thread count. b is the linking length in units of the mean spacing
 * `spacing` (Mpc/h); groups of at least nmin members are kept. want_members: also
 * fill groups->member. Returns RS_FOF_OK, or an RS_FOF_E* code (rs_fof_strerror)
 * with everything released. */
int rs_fof_run(int64_t n, const float *x, const float *y, const float *z,
               const struct rs_units *u, int64_t nrow, double spacing, double b,
               int64_t nmin, int want_members,
               struct rs_cellindex *cells, struct rs_fof_groups *groups);

/* the float linking length r = (float)(spacing * b) and its square */
float rs_fof_linking_length(double spacing, double b);

void rs_fof_release_labels(struct rs_fof_groups *g);    /* frees group_of_sorted */
void rs_fof_free_groups(struct rs_fof_groups *g);
void rs_fof_free_cells(struct rs_cellindex *c);

/* bytes the FoF needs at its peak for n particles (the memory pre-check), and the
 * bytes it still holds when it returns (cell index and group labels) */
double rs_fof_peak_bytes(int64_t n, int64_t nrow, int nthreads);
double rs_fof_retained_bytes(int64_t n, int64_t nrow);

/* VmRSS of this process in GiB (0 if unknown) */
double rs_vmrss_gib(void);

#ifdef __cplusplus
}
#endif

#endif
