/* test_fof.c: the adapter's FoF (rs_fof.c) against an O(n^2) brute force on
 * adversarial particle sets at the production geometry (Box 1024 Mpc/h, Ngrid 4096,
 * spacing 0.5 Mpc/h, b = 0.28: fine cells of 1/7 Mpc/h; coarse cells of 2 Mpc/h
 * instead of 1, so the cell offsets take 1 GB instead of 8). The brute
 * force has its own conversion and distance test (the definition: float
 * differences, minimum image, (dx*dx + dy*dy) + dz*dz < r*r in float), so it does
 * not share code with rs_fof.c. Also checked: the canonical group order, bitwise
 * identical output with 1 to 4 threads, and the error codes for bad coordinates
 * and an impossible grid.
 *   test_fof            prints the checks and exits 0 when all pass
 * Part of the MG-GLAM in-line Rockstar adapter; GNU GPL v3 or later. */
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rs_fof.h"

#define BOX 1024.0
#define NGRID 4096
#define NROW 1024                    /* sets the coarse grid only (Nrow/2 cells per side) */
#define SPACING 0.5                  /* Mpc/h, production's 1024/2048 */
#define NMAX 40000

static float gx[NMAX], gy[NMAX], gz[NMAX];
static int64_t np = 0;
static uint64_t rng = 0x9E3779B97F4A7C15ULL;

static double urand(void)
{
    rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17;
    return (double)(rng >> 11) * (1.0 / 9007199254740992.0);
}

static double gauss(void)
{
    double a = urand(), b = urand();
    return sqrt(-2.0 * log(a > 1e-300 ? a : 1e-300)) * cos(2.0 * M_PI * b);
}

/* a particle at a position in Mpc/h (any value), stored as GLAM grid coordinates in [1, Ngrid+1) */
static void add(double x, double y, double z)
{
    double q[3] = {x, y, z};
    float g[3];
    for (int d = 0; d < 3; d++) {
        double v = fmod(q[d], BOX);
        if (v < 0) v += BOX;
        g[d] = (float)(v / BOX * NGRID + 1.0);
        if (g[d] >= (float)(NGRID + 1)) g[d] -= (float)NGRID;
    }
    gx[np] = g[0]; gy[np] = g[1]; gz[np] = g[2]; np++;
}

/* a pair at distance s apart in a random direction around (x, y, z) */
static void add_pair(double x, double y, double z, double s)
{
    double u[3] = {gauss(), gauss(), gauss()}, l = sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
    add(x, y, z);
    add(x + s * u[0] / l, y + s * u[1] / l, z + s * u[2] / l);
}

/* ---- the brute force: own conversion, own predicate, serial union-find ---- */
static float to_mpc(float xg)
{
    const float boxf = (float)BOX;
    double x = ((double)xg - 1.0) * (BOX / NGRID);
    if (x >= BOX) x -= BOX;
    if (x < 0) x += BOX;
    float f = (float)x;
    return f < boxf ? f : 0.0f;
}

static int friends(const float *a, const float *b, float r2)
{
    const float boxf = (float)BOX, half = 0.5f * boxf;
    float d[3];
    for (int k = 0; k < 3; k++) {
        float t = a[k] - b[k];
        if (t > half) t = t - boxf;
        else if (t < -half) t = t + boxf;
        d[k] = t;
    }
    float s = d[0] * d[0];
    s = s + d[1] * d[1];
    s = s + d[2] * d[2];
    return s < r2;
}

static int64_t root(int64_t *par, int64_t i)
{
    while (par[i] != i) i = par[i] = par[par[i]];
    return i;
}

static void brute(float r, int64_t *lab)
{
    float (*q)[3] = malloc((size_t)np * sizeof *q);
    for (int64_t i = 0; i < np; i++) { q[i][0] = to_mpc(gx[i]); q[i][1] = to_mpc(gy[i]); q[i][2] = to_mpc(gz[i]); }
    const float r2 = r * r;
    for (int64_t i = 0; i < np; i++) lab[i] = i;
    for (int64_t i = 0; i < np; i++)
        for (int64_t j = i + 1; j < np; j++)
            if (friends(q[i], q[j], r2)) {
                int64_t a = root(lab, i), b = root(lab, j);
                if (a != b) lab[a > b ? a : b] = a < b ? a : b;
            }
    for (int64_t i = 0; i < np; i++) lab[i] = root(lab, i);          /* the smallest index of the group */
    free(q);
}

struct out { struct rs_cellindex c; struct rs_fof_groups g; };

static int run(int threads, double b, int64_t nmin, struct out *o)
{
    struct rs_units u;
    rs_units_init(&u, BOX, NGRID, 1.0, 0.01);
    omp_set_num_threads(threads);
    return rs_fof_run(np, gx, gy, gz, &u, NROW, SPACING, b, nmin, 1, &o->c, &o->g);
}

static int same(const struct out *a, const struct out *b)
{
    if (a->g.ngroups != b->g.ngroups || a->g.nmembers != b->g.nmembers || a->c.nc != b->c.nc) return 0;
    if (memcmp(a->g.offset, b->g.offset, (size_t)(a->g.ngroups + 1) * 8)) return 0;
    if (memcmp(a->g.member, b->g.member, (size_t)a->g.nmembers * 8)) return 0;
    if (memcmp(a->g.group_of_sorted, b->g.group_of_sorted, (size_t)np * 4)) return 0;
    if (memcmp(a->c.perm, b->c.perm, (size_t)np * 8)) return 0;
    const int64_t nc3 = a->c.nc * a->c.nc * a->c.nc;
    return !memcmp(a->c.offsets, b->c.offsets, (size_t)(nc3 + 1) * 8);
}

static void release(struct out *o)
{
    rs_fof_free_groups(&o->g);
    rs_fof_free_cells(&o->c);
}

int main(void)
{
    const double b = 0.28, r = (double)rs_fof_linking_length(SPACING, b);
    /* background */
    for (int i = 0; i < 4000; i++) add(urand() * BOX, urand() * BOX, urand() * BOX);
    /* a cloud on the corner, wrapping in all three dimensions */
    for (int i = 0; i < 3000; i++) add(0.3 * gauss(), 0.3 * gauss(), 0.3 * gauss());
    /* pairs at r (1 +- 2e-4) across the x = 0 / Box face and across fine-cell boundaries near x = 1000 */
    for (int i = 0; i < 2000; i++) {
        double s = r * (1.0 + (i % 2 ? 2e-4 : -2e-4));
        add_pair(BOX - 0.5 * s * urand(), urand() * BOX, urand() * BOX, s);
        add_pair(1000.0 + (i % 7) / 7.0 + 1e-4 * gauss(), urand() * BOX, urand() * BOX, s);
    }
    /* the edges of GLAM's coordinate range: x_grid = 1, Ngrid + 1 and their neighbouring floats */
    const float edges[] = {1.0f, nextafterf(1.0f, 2.0f), (float)(NGRID + 1), nextafterf((float)(NGRID + 1), 0.0f),
                           (float)NGRID, nextafterf((float)NGRID, 0.0f)};
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++) {
            gx[np] = edges[i]; gy[np] = edges[j]; gz[np] = edges[(i + j) % 6]; np++;
        }
    /* a dense blob, and a chain at 0.9 r across the periodic boundary */
    double cx = urand() * BOX, cy = urand() * BOX, cz = urand() * BOX;
    for (int i = 0; i < 5000; i++) add(cx + 0.5 * gauss(), cy + 0.5 * gauss(), cz + 0.5 * gauss());
    for (int i = 0; i < 300; i++) add(BOX - 10.0 + 0.9 * r * i, 500.0, 500.0);

    int fails = 0;
    int64_t *lab = malloc((size_t)np * sizeof(int64_t));
    brute((float)r, lab);
    int64_t *size = calloc((size_t)np, sizeof(int64_t));
    for (int64_t i = 0; i < np; i++) size[lab[i]]++;
    printf("test_fof: %lld particles, r = %.9g Mpc/h\n", (long long)np, r);

    const int64_t nmins[] = {2, 10};
    for (int t = 0; t < 2; t++) {
        const int64_t nmin = nmins[t];
        struct out ref;
        if (run(1, b, nmin, &ref) != RS_FOF_OK) { printf("  rs_fof_run failed\n"); return 1; }
        int64_t want = 0, wantm = 0;
        for (int64_t i = 0; i < np; i++) if (lab[i] == i && size[i] >= nmin) { want++; wantm += size[i]; }
        int64_t bad = 0, order = 0;
        for (int64_t gi = 0; gi < ref.g.ngroups; gi++) {
            const int64_t a = ref.g.offset[gi], e = ref.g.offset[gi + 1], first = ref.g.member[a];
            if (lab[first] != first || size[first] != e - a) bad++;            /* the group's smallest index */
            for (int64_t k = a; k < e; k++) {
                if (lab[ref.g.member[k]] != first) bad++;
                if (k > a && ref.g.member[k] <= ref.g.member[k - 1]) order++;
            }
            if (gi > 0) {
                const int64_t sp = a - ref.g.offset[gi - 1];
                if (e - a > sp || (e - a == sp && first <= ref.g.member[ref.g.offset[gi - 1]])) order++;
            }
        }
        int ok = ref.g.ngroups == want && ref.g.nmembers == wantm && bad == 0 && order == 0;
        printf("  nmin %lld: groups %lld (brute force %lld), members %lld (%lld), mismatches %lld, order %lld: %s\n",
               (long long)nmin, (long long)ref.g.ngroups, (long long)want, (long long)ref.g.nmembers,
               (long long)wantm, (long long)bad, (long long)order, ok ? "ok" : "FAIL");
        fails += !ok;
        for (int th = 2; th <= 4; th++) {
            struct out o;
            int same_ok = run(th, b, nmin, &o) == RS_FOF_OK && same(&ref, &o);
            printf("  nmin %lld, %d threads: %s\n", (long long)nmin, th, same_ok ? "identical" : "DIFFERENT");
            fails += !same_ok;
            release(&o);
        }
        release(&ref);
    }

    /* error codes: a non-finite coordinate, and a box too small for the linking length */
    {
        struct out o;
        float keep = gx[0];
        gx[0] = NAN;
        int rc = run(2, b, 10, &o);
        gx[0] = keep;
        int ok = rc == RS_FOF_ECOORD && o.g.nbad == 1;
        printf("  NaN coordinate: %s (%s)\n", ok ? "ok" : "FAIL", rs_fof_strerror(rc));
        fails += !ok;
        struct rs_units u;
        rs_units_init(&u, 0.3, NGRID, 1.0, 0.01);
        rc = rs_fof_run(np, gx, gy, gz, &u, NROW, SPACING, b, 10, 1, &o.c, &o.g);
        ok = rc == RS_FOF_EGRID;
        printf("  box of 0.3 Mpc/h: %s (%s)\n", ok ? "ok" : "FAIL", rs_fof_strerror(rc));
        fails += !ok;
    }
    free(lab);
    free(size);
    printf("test_fof: %s\n", fails ? "FAIL" : "PASS");
    return fails ? 1 : 0;
}
