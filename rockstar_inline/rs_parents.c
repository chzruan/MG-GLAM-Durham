/* rs_parents.c: pid / upid of the in-line catalogue with Rockstar's own
 * find_parents (MPI-Rockstar src/parents.c, included as the template it is):
 * halos in order of decreasing vmax; each takes the halos whose centres lie
 * inside its radius (periodic) and that have a smaller radius as children, a
 * later (lower-vmax) host overwriting an earlier one. upid is the root of the
 * resulting chain (radii strictly increase along it, so it ends).
 * Part of the MG-GLAM in-line Rockstar adapter; GNU GPL v3 or later. */
#include <stdint.h>
#include <stdlib.h>

#include "check_syscalls.h"
#include "config_vars.h"

struct rs_ph {
    float   pos[3];
    float   radius;      /* kpc/h, as the Rvir column find_parents reads */
    float   vmax;
    int64_t id, parent;
};

#define FAST3TREE_TYPE   struct rs_ph
#define FAST3TREE_PREFIX RS_PARENTS
#include "fast3tree.c"

static struct rs_ph *rs_ph_list = NULL;
#define GROUP_LIST rs_ph_list
#define RADIUS     radius
#include "parents.c"

#include "rs_inline_internal.h"

int rs_find_parents(int64_t n, struct rs_record *rec, int use_strict_radius)
{
    for (int64_t i = 0; i < n; i++) rec[i].pid = rec[i].upid = -1;
    if (n < 2) return 0;
    rs_ph_list = malloc((size_t)n * sizeof(struct rs_ph));
    if (!rs_ph_list) return -1;
    for (int64_t i = 0; i < n; i++) {
        for (int d = 0; d < 3; d++) rs_ph_list[i].pos[d] = rec[i].h.pos[d];
        rs_ph_list[i].radius = use_strict_radius ? rec[i].r_strict : rec[i].h.r;
        rs_ph_list[i].vmax = rec[i].h.vmax;
        rs_ph_list[i].id = i;
        rs_ph_list[i].parent = -1;
    }
    find_parents(n);                         /* reorders rs_ph_list (tree build): map back by id */
    for (int64_t i = 0; i < n; i++) rec[rs_ph_list[i].id].pid = rs_ph_list[i].parent;
    free(rs_ph_list);
    rs_ph_list = NULL;
    for (int64_t i = 0; i < n; i++) {
        int64_t u = rec[i].pid, steps = 0;
        if (u < 0) continue;
        while (rec[u].pid >= 0 && steps++ < n) u = rec[u].pid;
        rec[i].upid = u;
    }
    return 0;
}
