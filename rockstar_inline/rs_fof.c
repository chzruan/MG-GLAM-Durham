/* rs_fof.c: parallel exact periodic friends-of-friends; see rs_fof.h.
 * Part of the MG-GLAM in-line Rockstar adapter; GNU GPL v3 or later.
 *
 * Steps (all deterministic):
 *  1. counting sort of the particle indices by x-slab of the coarse grid
 *     (per-chunk histograms, chunks in index order);
 *  2. per slab: stable counting sort by (y, z) coarse cell, then a stable sort
 *     of each coarse cell by fine cell; this fills perm and the cell offsets;
 *  3. positions (float, Mpc/h) and fine-cell ids in sorted order;
 *  4. linking: every unordered pair of same or adjacent fine cells is visited
 *     once (by the cell with the smaller key), lock-free union-find on the sorted
 *     index with the larger root linked under the smaller one, so each root ends
 *     as the smallest sorted index of its group whatever the thread schedule;
 *  5. groups: sizes, the kept groups (>= nmin) and their canonical order. The
 *     positions' memory is reused for the int64 group sizes and the int32 labels,
 *     so step 5 never needs more than step 4.
 * Coordinates that rs_pos cannot map (not finite, or more than one box outside
 * GLAM's range) stop the FoF in step 1 (RS_FOF_ECOORD). */
#define _POSIX_C_SOURCE 200809L
#include "rs_fof.h"

#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

float rs_fof_linking_length(double spacing, double b) { return (float)(spacing * b); }

const char *rs_fof_strerror(int rc)
{
    switch (rc) {
    case RS_FOF_OK:      return "ok";
    case RS_FOF_ENOMEM:  return "out of memory";
    case RS_FOF_EGRID:   return "no cell grid: the box is too small for the linking length, or float positions "
                                "too coarse for it";
    case RS_FOF_ECOORD:  return "particle coordinates not finite or outside GLAM's range";
    case RS_FOF_EGROUPS: return "2^31 groups or more";
    default:             return "unknown error";
    }
}

double rs_vmrss_gib(void)
{
    FILE *f = fopen("/proc/self/status", "r");
    char line[256];
    double kb = 0;
    if (!f) return 0;
    while (fgets(line, sizeof line, f))
        if (!strncmp(line, "VmRSS:", 6)) {
            kb = atof(line + 6);
            break;
        }
    fclose(f);
    return kb / 1048576.0;
}

struct grid {
    int64_t nc, m;          /* coarse cells per side, fine cells per coarse side */
    double  cs, inv_cs;     /* coarse side and its inverse */
    double  fs, inv_fs;     /* fine side and its inverse */
};

static inline int64_t clampi(int64_t v, int64_t lo, int64_t hi) { return v < lo ? lo : (v > hi ? hi : v); }

/* coarse cell (cx, cy, cz) and the local fine cell of a position (Mpc/h, float) */
static inline void cell_of(const struct grid *g, const float p[3], int64_t c[3], int64_t f[3])
{
    for (int d = 0; d < 3; d++) {
        c[d] = clampi((int64_t)((double)p[d] * g->inv_cs), 0, g->nc - 1);
        f[d] = clampi((int64_t)(((double)p[d] - (double)c[d] * g->cs) * g->inv_fs), 0, g->m - 1);
    }
}

static inline void pos_of(const struct rs_units *u, const float *x, const float *y, const float *z,
                          int64_t i, float p[3])
{
    p[0] = rs_pos(u, x[i]);
    p[1] = rs_pos(u, y[i]);
    p[2] = rs_pos(u, z[i]);
}

/* ---- lock-free union-find on int64 parents (roots only ever point to smaller indices) ---- */
static inline int64_t uf_find(int64_t *parent, int64_t i)
{
    for (;;) {
        int64_t p = __atomic_load_n(&parent[i], __ATOMIC_RELAXED);
        if (p == i) return i;
        int64_t gp = __atomic_load_n(&parent[p], __ATOMIC_RELAXED);
        if (gp != p) __atomic_compare_exchange_n(&parent[i], &p, gp, 0, __ATOMIC_RELAXED, __ATOMIC_RELAXED);
        i = gp;
    }
}

static inline void uf_union(int64_t *parent, int64_t a, int64_t b)
{
    for (;;) {
        a = uf_find(parent, a);
        b = uf_find(parent, b);
        if (a == b) return;
        if (a > b) { int64_t t = a; a = b; b = t; }
        int64_t expected = b;
        if (__atomic_compare_exchange_n(&parent[b], &expected, a, 0, __ATOMIC_ACQ_REL, __ATOMIC_RELAXED))
            return;
    }
}

static int setup_grid(struct grid *g, double box, float r, int64_t nrow)
{
    double rmin = (double)r * (1.0 + RS_FOF_CELL_MARGIN);
    /* a wrapped float difference is rounded by up to ulp(Box)/2: the margin must absorb it */
    const float boxf = (float)box;
    const double ulp = (double)nextafterf(boxf, INFINITY) - (double)boxf;
    if (rmin < (double)r + 4.0 * ulp) rmin = (double)r + 4.0 * ulp;
    int64_t ncmax = (int64_t)floor(box / rmin);
    if (ncmax < 3) return -1;
    int64_t nc = nrow / 2;
    if (nc < 3) nc = 3;
    if (nc > ncmax) nc = ncmax;
    if (nc > 2048) nc = 2048;                      /* offsets of nc^3 int64: <= 64 GiB */
    g->nc = nc;
    g->cs = box / (double)nc;
    g->inv_cs = (double)nc / box;
    int64_t m = (int64_t)floor(g->cs / rmin);
    if (m < 1) m = 1;
    if (m > 40) m = 40;                            /* fine ids fit uint16 (m^3 <= 64000) */
    g->m = m;
    g->fs = g->cs / (double)m;
    g->inv_fs = (double)m / g->cs;
    return 0;
}

double rs_fof_peak_bytes(int64_t n, int64_t nrow, int nthreads)
{
    int64_t nc = nrow / 2 < 3 ? 3 : nrow / 2;
    if (nc > 2048) nc = 2048;
    const double dn = (double)n, offsets = 8.0 * (double)nc * nc * nc, T = nthreads > 0 ? nthreads : 1;
    /* step 2: perm 8 + positions 12 + fine 2 per particle and the cell offsets, plus per thread a slab
     * buffer (tperm 8 + tkey 4 + tfine 2 + ffine 2 = 16 B per slab particle, grown by 1.25, slabs up to 3x
     * the mean) and the (y, z) and fine-cell counters */
    double step2 = 22.0 * dn + offsets + T * (16.0 * 1.25 * 3.0 * dn / (double)nc + 8.0 * (double)nc * nc + 8.0 * 64001.0);
    /* step 4 (the scratch of step 2 is gone): + parent 8 per particle; step 5 holds less */
    double step4 = 30.0 * dn + offsets;
    return 1.15 * (step2 > step4 ? step2 : step4);     /* allocator and page granularity */
}

void rs_fof_free_cells(struct rs_cellindex *c)
{
    free(c->perm);
    free(c->offsets);
    c->perm = c->offsets = NULL;
}

void rs_fof_release_labels(struct rs_fof_groups *g)
{
    free(g->group_of_sorted);
    g->group_of_sorted = NULL;
}

void rs_fof_free_groups(struct rs_fof_groups *g)
{
    free(g->offset);
    free(g->member);
    free(g->group_of_sorted);
    memset(g, 0, sizeof *g);
}

/* first position in fine[lo, hi) with fine id >= v (fine ids are sorted inside a coarse cell) */
static inline int64_t lower_u16(const uint16_t *fine, int64_t lo, int64_t hi, int64_t v)
{
    while (lo < hi) {
        int64_t mid = lo + (hi - lo) / 2;
        if ((int64_t)fine[mid] < v) lo = mid + 1; else hi = mid;
    }
    return lo;
}

/* stable insertion sort of (perm, fine) pairs by fine id (small cells) */
static void isort_fine(int64_t *perm, uint16_t *fine, int64_t n)
{
    for (int64_t i = 1; i < n; i++) {
        int64_t pv = perm[i];
        uint16_t fv = fine[i];
        int64_t j = i - 1;
        while (j >= 0 && fine[j] > fv) {
            perm[j + 1] = perm[j];
            fine[j + 1] = fine[j];
            j--;
        }
        perm[j + 1] = pv;
        fine[j + 1] = fv;
    }
}

struct gsort { int64_t size, minidx, gp; };

static int cmp_groups(const void *pa, const void *pb)
{
    const struct gsort *a = pa, *b = pb;
    if (a->size != b->size) return a->size > b->size ? -1 : 1;
    return a->minidx < b->minidx ? -1 : (a->minidx > b->minidx);
}

static int cmp_i64(const void *pa, const void *pb)
{
    int64_t a = *(const int64_t *)pa, b = *(const int64_t *)pb;
    return a < b ? -1 : (a > b);
}

/* the canonical group order (a total order, so the result does not depend on the algorithm): chunks
 * sorted in parallel, then merged pairwise; plain qsort when small or without scratch memory */
static void sort_groups(struct gsort *gs, int64_t ng)
{
    int nchunk = 1;
    while (nchunk < omp_get_max_threads() && nchunk < 64) nchunk *= 2;
    struct gsort *tmp = (ng >= 65536 && nchunk > 1) ? malloc((size_t)ng * sizeof(struct gsort)) : NULL;
    if (!tmp) {
        qsort(gs, (size_t)ng, sizeof(struct gsort), cmp_groups);
        return;
    }
#pragma omp parallel for schedule(static, 1)
    for (int c = 0; c < nchunk; c++) {
        int64_t lo = ng * c / nchunk, hi = ng * (c + 1) / nchunk;
        qsort(gs + lo, (size_t)(hi - lo), sizeof(struct gsort), cmp_groups);
    }
    struct gsort *src = gs, *dst = tmp;
    for (int width = 1; width < nchunk; width *= 2) {
#pragma omp parallel for schedule(dynamic, 1)
        for (int c = 0; c < nchunk; c += 2 * width) {
            int cm = c + width < nchunk ? c + width : nchunk, ch = c + 2 * width < nchunk ? c + 2 * width : nchunk;
            int64_t i = ng * c / nchunk, mid = ng * cm / nchunk, j = mid, hi = ng * ch / nchunk, o = i;
            while (i < mid && j < hi) dst[o++] = cmp_groups(src + j, src + i) < 0 ? src[j++] : src[i++];
            while (i < mid) dst[o++] = src[i++];
            while (j < hi) dst[o++] = src[j++];
        }
        struct gsort *t = src; src = dst; dst = t;
    }
    if (src != gs) memcpy(gs, src, (size_t)ng * sizeof(struct gsort));
    free(tmp);
}

int rs_fof_run(int64_t n, const float *x, const float *y, const float *z,
               const struct rs_units *u, int64_t nrow, double spacing, double b,
               int64_t nmin, int want_members,
               struct rs_cellindex *cells, struct rs_fof_groups *groups)
{
    memset(cells, 0, sizeof *cells);
    memset(groups, 0, sizeof *groups);
    const float r = rs_fof_linking_length(spacing, b), r2 = r * r;
    const float boxf = u->box_f, half = 0.5f * u->box_f;
    struct grid g;
    if (setup_grid(&g, u->box, r, nrow) != 0) {
        fprintf(stderr, "rs_fof: no cell grid for box %.6g and linking length %.6g\n", u->box, r);
        return RS_FOF_EGRID;
    }
    const int64_t nc = g.nc, nc2 = nc * nc, nc3 = nc2 * nc, m = g.m, m3 = m * m * m;
    const int T = omp_get_max_threads();

    int64_t *perm = malloc((size_t)n * sizeof(int64_t));
    int64_t *offsets = malloc((size_t)(nc3 + 1) * sizeof(int64_t));
    int64_t *hist = calloc((size_t)T * (size_t)nc, sizeof(int64_t));
    int64_t *slab = malloc((size_t)(nc + 1) * sizeof(int64_t));
    if (!perm || !offsets || !hist || !slab) goto fail_early;

    /* ---- 1. by x-slab: chunk t (indices [n t/T, n (t+1)/T)) is always handled as unit t ---- */
    int64_t nbad = 0;
#pragma omp parallel for schedule(static, 1) reduction(+ : nbad)
    for (int t = 0; t < T; t++) {
        int64_t lo = n * t / T, hi = n * (t + 1) / T;
        int64_t *h = hist + (int64_t)t * nc;
        for (int64_t i = lo; i < hi; i++) {
            float p[3]; int64_t c[3], f[3];
            if (!rs_pos_ok(u, x[i]) || !rs_pos_ok(u, y[i]) || !rs_pos_ok(u, z[i])) nbad++;
            pos_of(u, x, y, z, i, p);
            cell_of(&g, p, c, f);
            h[c[0]]++;
        }
    }
    if (nbad > 0) {
        free(perm); free(offsets); free(hist); free(slab);
        groups->nbad = nbad;
        fprintf(stderr, "rs_fof: %lld particles with coordinates outside GLAM's range\n", (long long)nbad);
        return RS_FOF_ECOORD;
    }
    {
        int64_t acc = 0;
        for (int64_t s = 0; s < nc; s++) {
            slab[s] = acc;
            for (int t = 0; t < T; t++) {
                int64_t v = hist[(int64_t)t * nc + s];
                hist[(int64_t)t * nc + s] = acc;
                acc += v;
            }
        }
        slab[nc] = acc;
    }
#pragma omp parallel for schedule(static, 1)
    for (int t = 0; t < T; t++) {
        int64_t lo = n * t / T, hi = n * (t + 1) / T;
        int64_t *h = hist + (int64_t)t * nc;
        for (int64_t i = lo; i < hi; i++) {
            float p[3]; int64_t c[3], f[3];
            pos_of(u, x, y, z, i, p);
            cell_of(&g, p, c, f);
            perm[h[c[0]]++] = i;
        }
    }
    free(hist);
    hist = NULL;

    /* ---- 2. per slab: by (y, z) cell, then by fine cell; positions and fine ids in sorted order ---- */
    float *spos = malloc((size_t)n * 3 * sizeof(float) + 8);        /* + 8: the int64 sizes of step 5 */
    uint16_t *fine = malloc((size_t)n * sizeof(uint16_t));
    int failed = !spos || !fine;
    if (!failed) {
#pragma omp parallel
        {
            int64_t cap = 0;
            int64_t *tperm = NULL, *cnt = malloc((size_t)(nc2 + 1) * sizeof(int64_t));
            int32_t *tkey = NULL;
            uint16_t *tfine = NULL, *ffine = NULL;
            int64_t *fcnt = malloc((size_t)(m3 + 1) * sizeof(int64_t));
            if (!cnt || !fcnt) __atomic_store_n(&failed, 1, __ATOMIC_RELAXED);
#pragma omp for schedule(dynamic, 1)
            for (int64_t s = 0; s < nc; s++) {
                if (__atomic_load_n(&failed, __ATOMIC_RELAXED)) continue;
                const int64_t s0 = slab[s], ns = slab[s + 1] - s0;
                if (ns > cap) {
                    cap = ns + ns / 4 + 1024;
                    int64_t *a = realloc(tperm, (size_t)cap * sizeof(int64_t));
                    if (a) tperm = a;
                    int32_t *k = realloc(tkey, (size_t)cap * sizeof(int32_t));
                    if (k) tkey = k;
                    uint16_t *f1 = realloc(tfine, (size_t)cap * sizeof(uint16_t));
                    if (f1) tfine = f1;
                    uint16_t *f2 = realloc(ffine, (size_t)cap * sizeof(uint16_t));
                    if (f2) ffine = f2;
                    if (!a || !k || !f1 || !f2) {
                        __atomic_store_n(&failed, 1, __ATOMIC_RELAXED);
                        cap = 0;
                        continue;
                    }
                }
                memset(cnt, 0, (size_t)(nc2 + 1) * sizeof(int64_t));
                for (int64_t q = 0; q < ns; q++) {
                    float p[3]; int64_t c[3], f[3];
                    pos_of(u, x, y, z, perm[s0 + q], p);
                    cell_of(&g, p, c, f);
                    tkey[q] = (int32_t)(c[1] * nc + c[2]);
                    tfine[q] = (uint16_t)((f[0] * m + f[1]) * m + f[2]);
                    cnt[tkey[q] + 1]++;
                }
                for (int64_t k = 0; k < nc2; k++) cnt[k + 1] += cnt[k];
                for (int64_t k = 0; k < nc2; k++) offsets[s * nc2 + k] = s0 + cnt[k];
                for (int64_t q = 0; q < ns; q++) {           /* stable scatter by (y, z) cell */
                    int64_t d = cnt[tkey[q]]++;
                    tperm[d] = perm[s0 + q];
                    ffine[d] = tfine[q];
                }
                /* each (y, z) cell by fine id, stably */
                for (int64_t k = 0; k < nc2; k++) {
                    int64_t a0 = offsets[s * nc2 + k] - s0;
                    int64_t a1 = (k + 1 < nc2 ? offsets[s * nc2 + k + 1] : slab[s + 1]) - s0;
                    int64_t na = a1 - a0;
                    if (na <= 1) continue;
                    if (na < 64) {
                        isort_fine(tperm + a0, ffine + a0, na);
                    } else {
                        memset(fcnt, 0, (size_t)(m3 + 1) * sizeof(int64_t));
                        for (int64_t q = a0; q < a1; q++) fcnt[ffine[q] + 1]++;
                        for (int64_t f = 0; f < m3; f++) fcnt[f + 1] += fcnt[f];
                        /* scatter into the tail of tfine/tkey as scratch (they are free now) */
                        int64_t *dst = perm + s0 + a0;             /* perm[s0+a0 ..] is free to reuse */
                        for (int64_t q = a0; q < a1; q++) {
                            int64_t d = fcnt[ffine[q]]++;
                            dst[d] = tperm[q];
                            tfine[a0 + d] = ffine[q];
                        }
                        memcpy(tperm + a0, dst, (size_t)na * sizeof(int64_t));
                        memcpy(ffine + a0, tfine + a0, (size_t)na * sizeof(uint16_t));
                    }
                }
                if (ns > 0) {
                    memcpy(perm + s0, tperm, (size_t)ns * sizeof(int64_t));
                    memcpy(fine + s0, ffine, (size_t)ns * sizeof(uint16_t));
                }
                for (int64_t q = 0; q < ns; q++) pos_of(u, x, y, z, perm[s0 + q], spos + 3 * (s0 + q));
            }
            free(tperm); free(tkey); free(tfine); free(ffine); free(cnt); free(fcnt);
        }
    }
    offsets[nc3] = n;
    free(slab);
    slab = NULL;
    if (failed) { free(spos); free(fine); goto fail_early; }

    /* ---- 4. linking ---- */
    int64_t *parent = malloc((size_t)n * sizeof(int64_t));
    if (!parent) { free(spos); free(fine); goto fail_early; }
#pragma omp parallel for schedule(static)
    for (int64_t k = 0; k < n; k++) parent[k] = k;
    groups->rss_link_gib = rs_vmrss_gib();             /* every array of the linking stage is touched now */

    const int64_t nf = nc * m;          /* fine cells per side over the whole box */
#pragma omp parallel for schedule(dynamic, 64)
    for (int64_t c = 0; c < nc3; c++) {
        const int64_t c0 = offsets[c], c1 = offsets[c + 1];
        if (c0 == c1) continue;
        const int64_t cx = c / nc2, cy = (c / nc) % nc, cz = c % nc;
        int64_t a0 = c0;
        while (a0 < c1) {                                   /* runs of one fine cell (sorted within the cell) */
            const int64_t fa = fine[a0];
            int64_t a1 = a0 + 1;
            while (a1 < c1 && fine[a1] == fa) a1++;
            const int64_t gx = cx * m + fa / (m * m), gy = cy * m + (fa / m) % m, gz = cz * m + fa % m;
            const int64_t keya = c * m3 + fa;
            for (int64_t dx = -1; dx <= 1; dx++)
            for (int64_t dy = -1; dy <= 1; dy++)
            for (int64_t dz = -1; dz <= 1; dz++) {
                const int64_t hx = (gx + dx + nf) % nf, hy = (gy + dy + nf) % nf, hz = (gz + dz + nf) % nf;
                const int64_t cb = ((hx / m) * nc + hy / m) * nc + hz / m;
                const int64_t fb = ((hx % m) * m + hy % m) * m + hz % m;
                const int64_t keyb = cb * m3 + fb;
                if (keyb < keya) continue;                  /* that pair belongs to the other fine cell */
                int64_t b0, b1;
                if (keyb == keya) {
                    b0 = a0;
                    b1 = a1;
                } else {
                    b0 = lower_u16(fine, offsets[cb], offsets[cb + 1], fb);
                    b1 = lower_u16(fine, b0, offsets[cb + 1], fb + 1);
                    if (b0 == b1) continue;
                }
                for (int64_t i = a0; i < a1; i++) {
                    const float *pi = spos + 3 * i;
                    int64_t ri = uf_find(parent, i);        /* stale after other unions, but still i's set */
                    for (int64_t j = (keyb == keya ? i + 1 : b0); j < b1; j++)
                        if (rs_linked(pi, spos + 3 * j, boxf, half, r2)) {
                            int64_t rj = uf_find(parent, j);
                            if (rj != ri) {                 /* rj == ri: ri is a root and j is in its set */
                                uf_union(parent, ri, rj);
                                ri = uf_find(parent, ri);
                            }
                        }
                }
            }
            a0 = a1;
        }
    }

    free(fine);

    /* ---- 5. groups: the positions' memory (12 n + 8 bytes) now holds the int32 labels in its first 4 n
     *      bytes and the int64 group sizes (then provisional ids) after them ---- */
    int32_t *label = (int32_t *)spos;
    int64_t *cnt = (int64_t *)((char *)spos + (((size_t)n * 4 + 7) & ~(size_t)7));
#pragma omp parallel for schedule(static)
    for (int64_t k = 0; k < n; k++) {                                     /* full compression: no unions now */
        __atomic_store_n(&parent[k], uf_find(parent, k), __ATOMIC_RELAXED);
        cnt[k] = 0;
    }
#pragma omp parallel for schedule(static)
    for (int64_t k = 0; k < n; k++) {
#pragma omp atomic update
        cnt[parent[k]]++;
    }
    /* kept roots in sorted-index order get provisional ids (chunked prefix, deterministic) */
    int64_t *chunkcount = calloc((size_t)T + 1, sizeof(int64_t));
    if (!chunkcount) { free(parent); free(spos); goto fail_early; }
#pragma omp parallel for schedule(static, 1)
    for (int t = 0; t < T; t++) {
        int64_t lo = n * t / T, hi = n * (t + 1) / T, c = 0;
        for (int64_t k = lo; k < hi; k++) if (parent[k] == k && cnt[k] >= nmin) c++;
        chunkcount[t + 1] = c;
    }
    for (int t = 0; t < T; t++) chunkcount[t + 1] += chunkcount[t];
    const int64_t ng = chunkcount[T];
    if (ng > INT32_MAX) {
        free(chunkcount); free(parent); free(spos);
        free(perm); free(offsets);
        memset(cells, 0, sizeof *cells);
        return RS_FOF_EGROUPS;
    }
    struct gsort *gs = malloc((size_t)(ng > 0 ? ng : 1) * sizeof(struct gsort));
    int64_t *canon = malloc((size_t)(ng > 0 ? ng : 1) * sizeof(int64_t));
    groups->offset = malloc((size_t)(ng + 1) * sizeof(int64_t));
    if (!gs || !canon || !groups->offset) {
        free(gs); free(canon); free(chunkcount); free(parent); free(spos);
        free(groups->offset); groups->offset = NULL; goto fail_early;
    }
#pragma omp parallel for schedule(static, 1)
    for (int t = 0; t < T; t++) {
        int64_t lo = n * t / T, hi = n * (t + 1) / T, id = chunkcount[t];
        for (int64_t k = lo; k < hi; k++) {
            if (parent[k] != k) continue;
            if (cnt[k] >= nmin) {
                gs[id].size = cnt[k];
                gs[id].minidx = INT64_MAX;
                gs[id].gp = id;
                cnt[k] = id++;
            } else {
                cnt[k] = -1;
            }
        }
    }
    free(chunkcount);
    /* label every sorted position with its provisional group, and the smallest particle index per group */
#pragma omp parallel for schedule(static)
    for (int64_t k = 0; k < n; k++) {
        int64_t gp = cnt[parent[k]];
        label[k] = (int32_t)gp;
        if (gp >= 0) {
            int64_t pi = perm[k];
            int64_t cur = __atomic_load_n(&gs[gp].minidx, __ATOMIC_RELAXED);
            while (pi < cur && !__atomic_compare_exchange_n(&gs[gp].minidx, &cur, pi, 0,
                                                             __ATOMIC_RELAXED, __ATOMIC_RELAXED)) {}
        }
    }
    free(parent);
    {
        int32_t *shrunk = realloc(spos, (size_t)(n > 0 ? n : 1) * sizeof(int32_t));    /* keeps the labels */
        if (shrunk) label = shrunk;
    }
    sort_groups(gs, ng);
    groups->ngroups = ng;
    groups->offset[0] = 0;
    int64_t largest = 0;
    for (int64_t gi = 0; gi < ng; gi++) {
        canon[gs[gi].gp] = gi;
        groups->offset[gi + 1] = groups->offset[gi] + gs[gi].size;
        if (gs[gi].size > largest) largest = gs[gi].size;
    }
    groups->nmembers = groups->offset[ng];
    groups->largest = largest;
    free(gs);
#pragma omp parallel for schedule(static)
    for (int64_t k = 0; k < n; k++)
        if (label[k] >= 0) label[k] = (int32_t)canon[label[k]];
    free(canon);
    groups->group_of_sorted = label;

    if (want_members) {
        groups->member = malloc((size_t)(groups->nmembers > 0 ? groups->nmembers : 1) * sizeof(int64_t));
        int64_t *cursor = malloc((size_t)(ng > 0 ? ng : 1) * sizeof(int64_t));
        if (!groups->member || !cursor) { free(cursor); rs_fof_free_groups(groups); goto fail_early; }
        memcpy(cursor, groups->offset, (size_t)ng * sizeof(int64_t));
#pragma omp parallel for schedule(static)
        for (int64_t k = 0; k < n; k++) {
            int64_t gi = label[k];
            if (gi < 0) continue;
            int64_t d;
#pragma omp atomic capture
            d = cursor[gi]++;
            groups->member[d] = perm[k];
        }
        free(cursor);
#pragma omp parallel for schedule(dynamic, 16)
        for (int64_t gi = 0; gi < ng; gi++)
            qsort(groups->member + groups->offset[gi], (size_t)(groups->offset[gi + 1] - groups->offset[gi]),
                  sizeof(int64_t), cmp_i64);
    }
    cells->nc = nc;
    cells->cs = g.cs;
    cells->perm = perm;
    cells->offsets = offsets;
    return RS_FOF_OK;

fail_early:
    free(perm); free(offsets); free(hist); free(slab);
    memset(cells, 0, sizeof *cells);
    fprintf(stderr, "rs_fof: out of memory\n");
    return RS_FOF_ENOMEM;
}
