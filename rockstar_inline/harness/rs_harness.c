/* rs_harness.c: test harness of the MG-GLAM in-line Rockstar adapter.
 *
 * Part of the MG-GLAM in-line Rockstar adapter; GNU GPL v3 or later, see
 * rs_convert.h. The INTERNAL particle header below mirrors
 * struct binary_output_header of MPI-Rockstar (src/io/io_internal.h, GPLv3).
 *
 * Reads a GLAM PM snapshot, converts it with rs_convert.h (the conversion the
 * in-line finder uses) and
 *   info     DIR STEP                        prints the header and unit checks
 *   internal DIR STEP OUTDIR NFILES [opts]   writes NFILES Rockstar INTERNAL
 *            particle files OUTDIR/glam.<block>.rbin plus OUTDIR/rockstar_ref.cfg
 *            for standalone MPI-Rockstar (the reference of the in-line finder)
 *   inline   DIR STEP CONFIG OUTDIR [--shift F]   the in-line finder on the snapshot, exactly as
 *            GLAM calls it (rs_inline_run): OUTDIR/Rockstar.<step>.<box>.bin + .meta;
 *            --shift F moves every particle by F*Box first (periodic test T3);
 *            --fail kernel|assert|alloc injects that failure into a first call, then
 *            runs a clean one in the same process (error-path test); --repeat N calls it N
 *            times in one process and prints VmRSS before each (residue test)
 *   fofcmp   DIR STEP [--no-reference]       test T2: the adapter's parallel FoF (rs_fof.c)
 *            against the serial reference, exact equality of the partitions
 *   fof      DIR STEP [B]                    serial reference FoF (fof_ref.c) with linking
 *            length B (default 0.28) times the mean spacing: groups of >= 10, the
 *            grouped fraction f and the largest group
 *   options: --shuffle SEED   random particle order across and within the files
 *            --shift F        add F*Box to every coordinate (periodic test)
 *            --force-res R    FORCE_RES of the config (comoving Mpc/h, default 0.25)
 *            --strict         STRICT_SO_MASSES = 1 in the config */
#define _FILE_OFFSET_BITS 64
#define _POSIX_C_SOURCE 200809L
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include <omp.h>

#include "../rs_convert.h"
#include "../rs_fof.h"
#include "../rs_inline.h"
#include "fof_ref.h"
#include "pm_read.h"

#define ROCKSTAR_MAGIC        (uint64_t)0xfadedacec0c0d0d0
#define PARTICLE_TYPE_FULL    2
#define VERSION_MAX_SIZE      12
#define CRITICAL_DENSITY      2.77519737e11   /* (Msun/h) / (Mpc/h)^3, Rockstar universal_constants.h */

struct rs_particle { int64_t id; float pos[6]; };            /* Rockstar struct particle */

struct binary_output_header {                                 /* Rockstar io_internal.h */
    uint64_t magic;
    int64_t  snap, chunk;
    float    scale, Om, Ol, h0;
    float    bounds[6];
    int64_t  num_halos, num_particles;
    float    box_size, particle_mass;
    int64_t  particle_type;
    int32_t  format_revision;
    char     rockstar_version[VERSION_MAX_SIZE];
    int32_t  add_flag;
    char     unused[256 - VERSION_MAX_SIZE - 4 * 12 - 4 * 2 - 8 * 6];
};
_Static_assert(sizeof(struct binary_output_header) == 256, "INTERNAL header must be 256 bytes");
_Static_assert(sizeof(struct rs_particle) == 32, "Rockstar particle must be 32 bytes");

struct snap {
    struct pm_header h;
    float *x, *y, *z, *vx, *vy, *vz;
};

static void die(const char *msg) { fprintf(stderr, "rs_harness: %s\n", msg); exit(1); }

static void load(const char *dir, int step, struct snap *s)
{
    if (pm_read_header(dir, step, &s->h) != 0) exit(1);
    size_t n = (size_t)s->h.nparticles;
    float **a[6] = {&s->x, &s->y, &s->z, &s->vx, &s->vy, &s->vz};
    for (int c = 0; c < 6; c++)
        if (!(*a[c] = malloc(n * sizeof(float)))) die("out of memory");
    int64_t bad = pm_read_particles(dir, step, &s->h, s->x, s->y, s->z, s->vx, s->vy, s->vz);
    if (bad < 0) exit(1);
    if (bad > 0) fprintf(stderr, "rs_harness: warning: %lld coordinates outside [1, Ngrid+1)\n", (long long)bad);
}

static double particle_mass(const struct pm_header *h)
{
    return CRITICAL_DENSITY * h->om * h->box * h->box * h->box / (double)h->nparticles;
}

static void print_info(const struct snap *s)
{
    const struct pm_header *h = &s->h;
    struct rs_units u;
    rs_units_init(&u, h->box, h->ngrid, h->aexpn, h->astep);
    printf("title         %.45s\n", h->title);
    printf("step          %d   a = %.6f   da = %.6f   z = %.5f\n", h->istep, h->aexpn, h->astep, 1.0 / h->aexpn - 1.0);
    printf("a_v           %.6f  (velocity epoch AEXPN - ASTEP/2)\n", u.a_vel);
    printf("Nrow %d  Ngrid %d  N %lld  Box %.3f Mpc/h  realisation %d\n", h->nrow, h->ngrid,
           (long long)h->nparticles, h->box, h->nrealization);
    printf("Om %.5f  OmL %.5f  h %.5f  m_p %.6e Msun/h  spacing %.6f Mpc/h\n", h->om, h->oml, h->hubble,
           particle_mass(h), cbrt(particle_mass(h) / (h->om * CRITICAL_DENSITY)));
    double pmin[3] = {1e30, 1e30, 1e30}, pmax[3] = {-1e30, -1e30, -1e30}, v2[3] = {0, 0, 0};
    for (int64_t i = 0; i < h->nparticles; i++) {
        float p[6];
        rs_convert(&u, s->x[i], s->y[i], s->z[i], s->vx[i], s->vy[i], s->vz[i], p);
        for (int c = 0; c < 3; c++) {
            if (p[c] < pmin[c]) pmin[c] = p[c];
            if (p[c] > pmax[c]) pmax[c] = p[c];
            v2[c] += (double)p[c + 3] * p[c + 3];
        }
    }
    for (int c = 0; c < 3; c++)
        printf("axis %d        x in [%.6f, %.6f] Mpc/h   v_rms %.2f km/s\n", c, pmin[c], pmax[c],
               sqrt(v2[c] / (double)h->nparticles));
}

/* splitmix64 + Fisher-Yates: a reproducible permutation of 0..n-1 */
static uint64_t sm64(uint64_t *s)
{
    uint64_t z = (*s += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

static int64_t *permutation(int64_t n, uint64_t seed)
{
    int64_t *p = malloc((size_t)n * sizeof(int64_t));
    if (!p) die("out of memory");
    for (int64_t i = 0; i < n; i++) p[i] = i;
    for (int64_t i = n - 1; i > 0; i--) {
        int64_t j = (int64_t)(sm64(&seed) % (uint64_t)(i + 1));
        int64_t t = p[i]; p[i] = p[j]; p[j] = t;
    }
    return p;
}

static void write_internal(const struct snap *s, const char *outdir, int nfiles, int shuffle,
                           uint64_t seed, double shift, double force_res, int strict)
{
    const struct pm_header *h = &s->h;
    const int64_t n = h->nparticles;
    struct rs_units u;
    rs_units_init(&u, h->box, h->ngrid, h->aexpn, h->astep);
    const float sgrid = (float)(shift * h->ngrid);       /* the shift in grid units */
    int64_t *perm = shuffle ? permutation(n, seed) : NULL;
    if (mkdir(outdir, 0755) != 0 && errno != EEXIST) die("cannot create OUTDIR");
    char name[4096];
    snprintf(name, sizeof name, "%s/out", outdir);
    if (mkdir(name, 0755) != 0 && errno != EEXIST) die("cannot create OUTDIR/out");
    struct rs_particle *buf = malloc((1 << 20) * sizeof(struct rs_particle));
    if (!buf) die("out of memory");
    for (int b = 0; b < nfiles; b++) {
        int64_t lo = n * b / nfiles, hi = n * (b + 1) / nfiles;
        struct binary_output_header bh;
        memset(&bh, 0, sizeof bh);
        bh.magic = ROCKSTAR_MAGIC;
        bh.snap = 0; bh.chunk = b;
        bh.scale = h->aexpn; bh.Om = h->om; bh.Ol = h->oml; bh.h0 = h->hubble;
        bh.num_particles = hi - lo;
        bh.box_size = (float)h->box;
        bh.particle_mass = (float)particle_mass(h);
        bh.particle_type = PARTICLE_TYPE_FULL;
        snprintf(bh.rockstar_version, VERSION_MAX_SIZE, "glamharness");
        snprintf(name, sizeof name, "%s/glam.%d.rbin", outdir, b);
        FILE *f = fopen(name, "wb");
        if (!f || fwrite(&bh, sizeof bh, 1, f) != 1) die("cannot write an INTERNAL file");
        for (int64_t i0 = lo; i0 < hi; i0 += (1 << 20)) {
            int64_t m = (hi - i0 < (1 << 20)) ? hi - i0 : (1 << 20);
            for (int64_t k = 0; k < m; k++) {
                int64_t i = perm ? perm[i0 + k] : i0 + k;
                buf[k].id = i;
                rs_convert(&u, s->x[i] + sgrid, s->y[i] + sgrid, s->z[i] + sgrid,
                           s->vx[i], s->vy[i], s->vz[i], buf[k].pos);
            }
            if (fwrite(buf, sizeof(struct rs_particle), (size_t)m, f) != (size_t)m) die("short write");
        }
        if (fclose(f) != 0) die("close failed");
    }
    free(buf);
    free(perm);
    snprintf(name, sizeof name, "%s/rockstar_ref.cfg", outdir);
    FILE *c = fopen(name, "w");
    if (!c) die("cannot write rockstar_ref.cfg");
    fprintf(c, "# standalone MPI-Rockstar reference for the in-line finder (written by rs_harness)\n"
               "# GLAM step %d, a = %.6f, velocity epoch a_v = %.6f (AEXPN - ASTEP/2)%s\n"
               "FILE_FORMAT = \"INTERNAL\"\nINBASE = \"%s\"\nFILENAME = \"glam.<block>.rbin\"\n"
               "NUM_BLOCKS = %d\nNUM_SNAPS = 1\nSTARTING_SNAP = 0\nOUTBASE = \"%s/out\"\n"
               "OUTPUT_FORMAT = \"BINARY\"\nTEMPORAL_HALO_FINDING = 0\nPERIODIC = 1\nEXTRA_PROFILING = 0\n"
               "BOX_SIZE = %.9g\nTOTAL_PARTICLES = %lld\nOm = %.9g\nOl = %.9g\nh0 = %.9g\nSCALE_NOW = %.9g\n"
               "FORCE_RES = %.9g\nMASS_DEFINITION = \"vir\"\nMIN_HALO_OUTPUT_SIZE = 20\nMIN_HALO_PARTICLES = 10\n"
               "FOF_LINKING_LENGTH = 0.28\nFOF_FRACTION = 0.7\nBOUND_PROPS = 1\nUNBOUND_THRESHOLD = 0.5\n"
               "STRICT_SO_MASSES = %d\n",
            h->istep, h->aexpn, u.a_vel, shift != 0 ? "; positions shifted" : "",
            outdir, nfiles, outdir, h->box, (long long)n, h->om, h->oml, h->hubble, h->aexpn,
            force_res, strict);
    fclose(c);
    printf("wrote %d INTERNAL files (%lld particles%s%s) and %s\n", nfiles, (long long)n,
           shuffle ? ", shuffled" : "", shift != 0 ? ", shifted" : "", name);
}

static void run_fof(const struct snap *s, double b)
{
    const struct pm_header *h = &s->h;
    const int64_t n = h->nparticles;
    struct rs_units u;
    rs_units_init(&u, h->box, h->ngrid, h->aexpn, h->astep);
    float *px = malloc((size_t)n * sizeof(float)), *py = malloc((size_t)n * sizeof(float)),
          *pz = malloc((size_t)n * sizeof(float));
    if (!px || !py || !pz) die("out of memory");
    for (int64_t i = 0; i < n; i++) {
        px[i] = rs_pos(&u, s->x[i]);
        py[i] = rs_pos(&u, s->y[i]);
        pz[i] = rs_pos(&u, s->z[i]);
    }
    /* r exactly as Rockstar: AVG_PARTICLE_SPACING = cbrt(m_p/(Om rho_c)) (double), r = spacing*b as float */
    double spacing = cbrt(particle_mass(h) / (h->om * CRITICAL_DENSITY));
    float r = (float)(spacing * b);
    struct fof_stats st;
    if (fof_reference(n, px, py, pz, (float)h->box, r, 10, NULL, &st) < 0) exit(1);
    printf("fof b=%.3f r=%.7f Mpc/h: groups(>=10) %lld  members %lld  f %.4f  largest %lld\n", b, r,
           (long long)st.ngroups, (long long)st.nmembers, st.fraction, (long long)st.largest);
    free(px); free(py); free(pz);
}

/* T2: the adapter's parallel FoF (rs_fof.c) against the serial reference, exact partition equality */
static int run_fofcmp(const struct snap *s, double b, int skip_reference)
{
    const struct pm_header *h = &s->h;
    const int64_t n = h->nparticles;
    struct rs_units u;
    rs_units_init(&u, h->box, h->ngrid, h->aexpn, h->astep);
    double spacing = cbrt(particle_mass(h) / (h->om * CRITICAL_DENSITY));
    struct rs_cellindex cells;
    struct rs_fof_groups g;
    double t0 = omp_get_wtime();
    if (rs_fof_run(n, s->x, s->y, s->z, &u, h->nrow, spacing, b, 10, 1, &cells, &g) != 0) return 1;
    double t1 = omp_get_wtime();
    printf("rs_fof (%d threads): %.2f s  coarse %lld^3 cells  groups %lld  members %lld  f %.4f  largest %lld\n",
           omp_get_max_threads(), t1 - t0, (long long)cells.nc, (long long)g.ngroups, (long long)g.nmembers,
           (double)g.nmembers / (double)n, (long long)g.largest);
    /* canonical order: sizes non-increasing, ties by smallest member; members sorted within each group */
    int bad_order = 0;
    for (int64_t gi = 0; gi < g.ngroups; gi++) {
        int64_t a = g.offset[gi], e = g.offset[gi + 1];
        for (int64_t k = a + 1; k < e; k++) if (g.member[k] <= g.member[k - 1]) bad_order++;
        if (gi > 0) {
            int64_t sp = g.offset[gi] - g.offset[gi - 1], sz = e - a;
            if (sz > sp || (sz == sp && g.member[a] <= g.member[g.offset[gi - 1]])) bad_order++;
        }
    }
    int64_t gsum = 0;                                       /* group_of_sorted agrees with the member lists */
    for (int64_t k = 0; k < n; k++) if (g.group_of_sorted[k] >= 0) gsum++;
    printf("canonical order violations %d; labelled sorted positions %lld (members %lld)\n", bad_order,
           (long long)gsum, (long long)g.nmembers);
    int rc = (bad_order != 0 || gsum != g.nmembers);
    if (!skip_reference) {
        float *px = malloc((size_t)n * sizeof(float)), *py = malloc((size_t)n * sizeof(float)),
              *pz = malloc((size_t)n * sizeof(float));
        int64_t *label = malloc((size_t)n * sizeof(int64_t));
        if (!px || !py || !pz || !label) die("out of memory");
        for (int64_t i = 0; i < n; i++) { px[i] = rs_pos(&u, s->x[i]); py[i] = rs_pos(&u, s->y[i]); pz[i] = rs_pos(&u, s->z[i]); }
        struct fof_stats st;
        double t2 = omp_get_wtime();
        if (fof_reference(n, px, py, pz, (float)h->box, rs_fof_linking_length(spacing, b), 10, label, &st) < 0) exit(1);
        double t3 = omp_get_wtime();
        int64_t mismatch = 0, labelled = 0;
        for (int64_t i = 0; i < n; i++) if (label[i] >= 0) labelled++;
        for (int64_t gi = 0; gi < g.ngroups; gi++) {
            int64_t root = g.member[g.offset[gi]];          /* smallest member = the reference's label */
            for (int64_t k = g.offset[gi]; k < g.offset[gi + 1]; k++)
                if (label[g.member[k]] != root) mismatch++;
        }
        printf("reference (serial): %.2f s  groups %lld  members %lld; T2 %s: group count %s, members %s, "
               "membership mismatches %lld\n", t3 - t2, (long long)st.ngroups, (long long)st.nmembers,
               (st.ngroups == g.ngroups && labelled == g.nmembers && mismatch == 0) ? "PASS" : "FAIL",
               st.ngroups == g.ngroups ? "equal" : "DIFFER", labelled == g.nmembers ? "equal" : "DIFFER",
               (long long)mismatch);
        rc |= !(st.ngroups == g.ngroups && labelled == g.nmembers && mismatch == 0);
        free(px); free(py); free(pz); free(label);
    }
    rs_fof_free_groups(&g);
    rs_fof_free_cells(&cells);
    return rc;
}

static int run_inline(struct snap *s, const char *config, const char *outdir, double shift)
{
    const struct pm_header *h = &s->h;
    if (shift != 0) {                               /* in grid units, kept in [1, Ngrid+1) */
        const float sg = (float)(shift * h->ngrid), top = (float)h->ngrid + 1.0f;
        float *a[3] = {s->x, s->y, s->z};
        for (int c = 0; c < 3; c++)
            for (int64_t i = 0; i < h->nparticles; i++) {
                float v = a[c][i] + sg;
                if (v >= top) v -= (float)h->ngrid;
                a[c][i] = v;
            }
    }
    struct rs_inline_info info;
    printf("%s\n", rs_inline_version());
    int rc = rs_inline_run(h->nparticles, s->x, s->y, s->z, s->vx, s->vy, s->vz, h->box, h->ngrid, h->aexpn,
                           h->astep, h->om, h->oml, h->hubble, h->istep, h->nrealization, config, outdir, &info);
    printf("status %d (%s); groups %lld members %lld largest %lld; halos %lld internal, %lld selected, %lld distinct written; "
           "FoF %.2f gather %.2f kernel %.2f strict %.2f parents %.2f write %.2f total %.2f s; "
           "VmRSS %.3f -> %.3f GiB; estimates FoF %.3f kernel %.3f GiB\n",
           rc, rc ? info.reason : "ok", (long long)info.ngroups, (long long)info.nmembers, (long long)info.largest,
           (long long)info.nhalos_internal, (long long)info.nhalos, (long long)info.nhosts, info.t_fof, info.t_gather,
           info.t_kernel, info.t_strict, info.t_parents, info.t_write, info.t_total, info.rss_entry_gib,
           info.rss_max_gib, info.need_fof_gib, info.need_kernel_gib);
    return rc;
}

int main(int argc, char **argv)
{
    if (argc < 4) {
        fprintf(stderr, "usage: rs_harness info DIR STEP\n"
                        "       rs_harness internal DIR STEP OUTDIR NFILES [--shuffle SEED] [--shift F]"
                        " [--force-res R] [--strict]\n");
        return 2;
    }
    const char *cmd = argv[1], *dir = argv[2];
    int step = atoi(argv[3]);
    struct snap s;
    if (!strcmp(cmd, "info")) {
        load(dir, step, &s);
        print_info(&s);
        return 0;
    }
    if (!strcmp(cmd, "inline") && argc >= 6) {
        double shift = 0.0;
        const char *fail = NULL;
        int repeat = 1;
        for (int i = 6; i + 1 < argc; i += 2) {
            if (!strcmp(argv[i], "--shift")) shift = atof(argv[i + 1]);
            else if (!strcmp(argv[i], "--fail")) fail = argv[i + 1];
            else if (!strcmp(argv[i], "--repeat")) repeat = atoi(argv[i + 1]);
        }
        load(dir, step, &s);
        if (repeat > 1) {   /* residue test: VmRSS before each of `repeat` calls in one process */
            int rc = 0;
            for (int k = 0; k < repeat; k++) {
                struct rs_inline_info info;
                const struct pm_header *h = &s.h;
                rc |= rs_inline_run(h->nparticles, s.x, s.y, s.z, s.vx, s.vy, s.vz, h->box, h->ngrid, h->aexpn,
                                    h->astep, h->om, h->oml, h->hubble, h->istep, h->nrealization, argv[4], argv[5],
                                    &info);
                printf("call %d: status %d, VmRSS at entry %.4f GiB, max sampled %.4f GiB\n", k + 1, info.status,
                       info.rss_entry_gib, info.rss_max_gib);
            }
            return rc;
        }
        if (fail) {     /* error-path test: a failing call, then a clean one in the same process */
            setenv("RS_INLINE_TEST_FAIL", fail, 1);
            int rc1 = run_inline(&s, argv[4], argv[5], shift);
            unsetenv("RS_INLINE_TEST_FAIL");
            printf("injected '%s': status %d (expected %d)\n", fail, rc1, RS_FAILED);
            int rc2 = run_inline(&s, argv[4], argv[5], shift);
            return (rc1 == RS_FAILED && rc2 == RS_OK) ? 0 : 1;
        }
        return run_inline(&s, argv[4], argv[5], shift);
    }
    if (!strcmp(cmd, "fofcmp")) {
        load(dir, step, &s);
        int skip = (argc > 4 && !strcmp(argv[4], "--no-reference"));
        return run_fofcmp(&s, 0.28, skip);
    }
    if (!strcmp(cmd, "fof")) {
        load(dir, step, &s);
        run_fof(&s, argc > 4 ? atof(argv[4]) : 0.28);
        return 0;
    }
    if (!strcmp(cmd, "internal") && argc >= 6) {
        int nfiles = atoi(argv[5]), shuffle = 0, strict = 0;
        uint64_t seed = 0;
        double shift = 0, force_res = 0.25;
        for (int i = 6; i < argc; i++) {
            if (!strcmp(argv[i], "--shuffle") && i + 1 < argc) { shuffle = 1; seed = strtoull(argv[++i], NULL, 10); }
            else if (!strcmp(argv[i], "--shift") && i + 1 < argc) shift = atof(argv[++i]);
            else if (!strcmp(argv[i], "--force-res") && i + 1 < argc) force_res = atof(argv[++i]);
            else if (!strcmp(argv[i], "--strict")) strict = 1;
            else { fprintf(stderr, "rs_harness: unknown option %s\n", argv[i]); return 2; }
        }
        if (nfiles < 1) die("NFILES must be >= 1");
        load(dir, step, &s);
        print_info(&s);
        write_internal(&s, argv[4], nfiles, shuffle, seed, shift, force_res, strict);
        return 0;
    }
    fprintf(stderr, "rs_harness: unknown command %s\n", cmd);
    return 2;
}
