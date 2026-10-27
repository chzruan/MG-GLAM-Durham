/* fof_ref.c: a deliberately simple serial periodic friends-of-friends, the
 * independent reference (test T2) of the adapter's parallel FoF, and the tool
 * that measures the grouped fraction f for the memory estimates.
 *
 * Part of the MG-GLAM in-line Rockstar adapter; GNU GPL v3 or later, see
 * rs_convert.h.
 *
 * Definition (as Rockstar's linking test, rockstar.c:59 / fast3tree.c:405):
 * particles i and j are friends if the periodic minimum-image distance of their
 * float32 positions satisfies d^2 < r^2, with r a float. Groups of at least
 * nmin particles are kept (fof.c:240-252 drops groups below MIN_HALO_PARTICLES).
 * Algorithm: particles sorted by cell (cells of side >= r (1 + RS_FOF_CELL_MARGIN),
 * as the adapter's FoF, so float rounding of wrapped differences cannot hide a
 * pair two cells apart), each cell linked with itself and its 13 forward
 * neighbours, union-find with path halving. The pair test is rs_linked of
 * rs_fof.h, the adapter's own. */
#include "fof_ref.h"
#include "../rs_fof.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int64_t find(int64_t *parent, int64_t i)
{
    while (parent[i] != i) {
        parent[i] = parent[parent[i]];
        i = parent[i];
    }
    return i;
}

static void unite(int64_t *parent, int64_t a, int64_t b)
{
    a = find(parent, a);
    b = find(parent, b);
    if (a == b) return;
    if (a < b) parent[b] = a; else parent[a] = b;      /* root = smallest index */
}

static int cmp_key(const void *pa, const void *pb)
{
    const int64_t *a = pa, *b = pb;           /* pairs (key, index): order by key, then index */
    if (a[0] != b[0]) return a[0] < b[0] ? -1 : 1;
    return a[1] < b[1] ? -1 : (a[1] > b[1]);
}

/* first position in keys[0..n) (stride 2) with key >= k */
static int64_t lower(const int64_t *kv, int64_t n, int64_t k)
{
    int64_t lo = 0, hi = n;
    while (lo < hi) {
        int64_t mid = lo + (hi - lo) / 2;
        if (kv[2 * mid] < k) lo = mid + 1; else hi = mid;
    }
    return lo;
}

int64_t fof_reference(int64_t n, const float *x, const float *y, const float *z,
                      float box, float r, int64_t nmin, int64_t *label,
                      struct fof_stats *st)
{
    const int64_t nc = (int64_t)floor((double)box / ((double)r * (1.0 + RS_FOF_CELL_MARGIN)));
    if (nc < 3) { fprintf(stderr, "fof_reference: box too small for r\n"); return -1; }
    const double inv = (double)nc / (double)box;
    const float r2 = r * r, half = 0.5f * box;
    int64_t *kv = malloc((size_t)n * 2 * sizeof(int64_t));
    int64_t *parent = malloc((size_t)n * sizeof(int64_t));
    if (!kv || !parent) { fprintf(stderr, "fof_reference: out of memory\n"); free(kv); free(parent); return -1; }
    for (int64_t i = 0; i < n; i++) {
        int64_t cx = (int64_t)(x[i] * inv), cy = (int64_t)(y[i] * inv), cz = (int64_t)(z[i] * inv);
        if (cx >= nc) cx = nc - 1;
        if (cy >= nc) cy = nc - 1;
        if (cz >= nc) cz = nc - 1;
        kv[2 * i] = (cx * nc + cy) * nc + cz;
        kv[2 * i + 1] = i;
        parent[i] = i;
    }
    qsort(kv, (size_t)n, 2 * sizeof(int64_t), cmp_key);
    static const int off[13][3] = {{0,0,1},{0,1,-1},{0,1,0},{0,1,1},{1,-1,-1},{1,-1,0},{1,-1,1},
                                   {1,0,-1},{1,0,0},{1,0,1},{1,1,-1},{1,1,0},{1,1,1}};
    int64_t a0 = 0;
    while (a0 < n) {
        const int64_t key = kv[2 * a0];
        int64_t a1 = a0;
        while (a1 < n && kv[2 * a1] == key) a1++;
        for (int64_t p = a0; p < a1; p++) {                      /* within the cell */
            int64_t i = kv[2 * p + 1];
            for (int64_t q = p + 1; q < a1; q++) {
                int64_t j = kv[2 * q + 1];
                const float pi[3] = {x[i], y[i], z[i]}, pj[3] = {x[j], y[j], z[j]};
                if (rs_linked(pi, pj, box, half, r2)) unite(parent, i, j);
            }
        }
        const int64_t cx = key / (nc * nc), cy = (key / nc) % nc, cz = key % nc;
        for (int o = 0; o < 13; o++) {                           /* with the 13 forward neighbours */
            int64_t nx = (cx + off[o][0] + nc) % nc, ny = (cy + off[o][1] + nc) % nc, nz = (cz + off[o][2] + nc) % nc;
            int64_t nkey = (nx * nc + ny) * nc + nz;
            int64_t b0 = lower(kv, n, nkey);
            if (b0 >= n || kv[2 * b0] != nkey) continue;
            int64_t b1 = b0;
            while (b1 < n && kv[2 * b1] == nkey) b1++;
            for (int64_t p = a0; p < a1; p++) {
                int64_t i = kv[2 * p + 1];
                for (int64_t q = b0; q < b1; q++) {
                    int64_t j = kv[2 * q + 1];
                    const float pi[3] = {x[i], y[i], z[i]}, pj[3] = {x[j], y[j], z[j]};
                    if (rs_linked(pi, pj, box, half, r2)) unite(parent, i, j);
                }
            }
        }
        a0 = a1;
    }
    free(kv);
    /* group sizes (stored at the roots), then labels: root index for kept groups, -1 otherwise */
    int64_t *size = calloc((size_t)n, sizeof(int64_t));
    if (!size) { fprintf(stderr, "fof_reference: out of memory\n"); free(parent); return -1; }
    for (int64_t i = 0; i < n; i++) size[find(parent, i)]++;
    struct fof_stats s = {0};
    for (int64_t i = 0; i < n; i++) {
        if (size[i] >= nmin) {
            s.ngroups++;
            s.nmembers += size[i];
            if (size[i] > s.largest) s.largest = size[i];
        }
    }
    if (label)
        for (int64_t i = 0; i < n; i++) {
            int64_t root = find(parent, i);       /* the root is the smallest member index */
            label[i] = size[root] >= nmin ? root : -1;
        }
    s.fraction = (double)s.nmembers / (double)n;
    if (st) *st = s;
    free(size);
    free(parent);
    return s.ngroups;
}
