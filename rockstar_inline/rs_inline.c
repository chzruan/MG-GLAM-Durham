/* rs_inline.c: Rockstar haloes of one MG-GLAM output, in process.
 * Part of the MG-GLAM in-line Rockstar adapter; GNU GPL v3 or later. The kernel
 * loop and merge follow MPI-Rockstar's do_halo_finding / copy_halos
 * (src/mpi_main.cpp, commit 0ae386f), the output filter its _should_print
 * (src/io/meta_io.c).
 *
 * One call (rs_inline_run, through the guard in rs_guard.cpp):
 *  1. configuration: Rockstar's parser on Run<box>/Rockstar.config, without
 *     setup_config (which takes the particle spacing from Rockstar's defaults
 *     and raises the core-dump limit); the run's values are set explicitly. The
 *     file is read once per process (at the start-up check or the first call)
 *     and its checksum goes into every .meta, so all outputs of a run use the
 *     same settings and a file system error at a later output cannot change them;
 *  2. memory pre-check, FoF (rs_fof.c), second pre-check;
 *  3. gather: the members of every group into Rockstar's array-of-structs
 *     buffer (the global p), each group sorted by particle id (the kernel's
 *     result depends on the member order) and unwrapped around its first member;
 *  4. kernel: find_subs on every group in an OpenMP loop (largest first), each
 *     call guarded; halos merged in group order, so the catalogue does not
 *     depend on the thread count; per-thread scratch released on every thread;
 *  5. halo centres wrapped, Rockstar's output filter, the distinct haloes
 *     (rs_parents.c: Rockstar's find_parents; subhaloes are not written, as GLAM's
 *     fixed resolution does not resolve substructure, but the kernel still finds
 *     them: a halo's properties include its substructure's particles), strict SO
 *     masses (rs_strict.c), publication (rs_write.c).
 * Failures and skips publish an empty catalogue whose .meta says why, unless a
 * complete catalogue of the step already exists (a restart), which is kept. */
#define _GNU_SOURCE
#include <errno.h>
#include <malloc.h>
#include <math.h>
#include <omp.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "bounds.h"
#include "config_vars.h"
#include "fof.h"
#include "groupies.h"
#include "halo.h"
#include "io/read_config.h"
#include "particle.h"
#include "universal_constants.h"

#include "rs_inline.h"
#include "rs_inline_internal.h"

#if defined(__has_include)
#if __has_include("rs_version.h")
#include "rs_version.h"                                      /* build/, written by the makefile */
#endif
#endif
#ifndef RS_ROCKSTAR_COMMIT
#define RS_ROCKSTAR_COMMIT "unknown"
#endif
#ifndef RS_BUILD_FLAGS
#define RS_BUILD_FLAGS "unknown"
#endif

#define GIB 1073741824.0

extern struct particle *p;                                   /* rs_glue.c */
int rs_guarded_find_subs(struct fof *f, struct FOFInfo *fofinfo, struct HaloInfo *haloinfo, int inject,
                         char *msg, int msglen);             /* rs_guard.cpp */

static double cfg_ceiling_gib = 900;    /* GLAM_MEMORY_CEILING_GIB */
static int    cfg_strict = 1;           /* GLAM_STRICT_SO */
static double cfg_force_res, cfg_force_res_phys_max;
static char     cfg_path[4096];         /* the configuration loaded in this process ("" = none) */
static uint64_t cfg_hash;               /* FNV-1a 64 of its bytes */

/* the output being written, for the sidecars of skipped and failed outputs too */
static struct { int set; double aexpn, astep, box; int64_t ngrid, n; } run_ctx;

int rs_inline_available(void) { return 1; }

const char *rs_inline_version(void)
{
    return "MPI-Rockstar " RS_ROCKSTAR_COMMIT "; " RS_BUILD_FLAGS;
}

static double status_gib(const char *key)                   /* a "Vm...:  N kB" line of /proc/self/status */
{
    FILE *f = fopen("/proc/self/status", "r");
    char line[256];
    double kb = 0;
    size_t kl = strlen(key);
    if (!f) return 0;
    while (fgets(line, sizeof line, f))
        if (!strncmp(line, key, kl)) {
            kb = atof(line + kl);
            break;
        }
    fclose(f);
    return kb / 1048576.0;
}

static void note_rss(struct rs_inline_info *info, double r)
{
    if (r > info->rss_max_gib) info->rss_max_gib = r;
}

static void sample_rss(struct rs_inline_info *info) { note_rss(info, rs_vmrss_gib()); }

/* the memory this process may use: the smallest limit of its cgroups (v2 memory.max, v1
 * memory.limit_in_bytes, the cgroup and every ancestor, e.g. a Slurm job's) and MemTotal, in bytes;
 * src says where it came from. Allocation failures are rare under Linux overcommit, so going past it
 * means the OOM killer: the pre-checks are the protection. */
static double read_number(const char *path)
{
    FILE *f = fopen(path, "r");
    char buf[64] = "";
    if (!f) return INFINITY;
    if (!fgets(buf, sizeof buf, f)) buf[0] = 0;
    fclose(f);
    char *end;
    double v = strtod(buf, &end);
    return (end == buf || v >= 1.0e18) ? INFINITY : v;     /* "max" and v1's "unlimited" */
}

static double status_meminfo_total(void)
{
    FILE *f = fopen("/proc/meminfo", "r");
    char line[256];
    double kb = 0;
    if (!f) return INFINITY;
    while (fgets(line, sizeof line, f))
        if (!strncmp(line, "MemTotal:", 9)) {
            kb = atof(line + 9);
            break;
        }
    fclose(f);
    return kb > 0 ? kb * 1024.0 : INFINITY;
}

static double memory_limit_bytes(char *src, size_t srclen)
{
    double best = INFINITY;
    snprintf(src, srclen, "none found");
    FILE *f = fopen("/proc/self/cgroup", "r");
    char line[4096];
    while (f && fgets(line, sizeof line, f)) {
        line[strcspn(line, "\n")] = 0;
        char *c1 = strchr(line, ':'), *c2 = c1 ? strchr(c1 + 1, ':') : NULL;
        if (!c2) continue;
        *c2 = 0;
        const char *ctrl = c1 + 1;
        char cur[4096];
        snprintf(cur, sizeof cur, "%s", c2 + 1);
        int v2 = ctrl[0] == 0, v1 = 0;
        for (const char *t = ctrl; t && *t; t = strchr(t, ',') ? strchr(t, ',') + 1 : NULL)
            if (!strncmp(t, "memory", 6) && (t[6] == ',' || t[6] == 0)) v1 = 1;
        if (!v1 && !v2) continue;
        for (;;) {
            char file[8300];
            if (v2) snprintf(file, sizeof file, "/sys/fs/cgroup%s/memory.max", cur);
            else    snprintf(file, sizeof file, "/sys/fs/cgroup/memory%s/memory.limit_in_bytes", cur);
            double v = read_number(file);
            if (v < best) {
                best = v;
                snprintf(src, srclen, "cgroup %s", file);
            }
            char *slash = strrchr(cur, '/');
            if (!slash || cur[1] == 0) break;
            if (slash == cur) cur[1] = 0; else *slash = 0;
        }
    }
    if (f) fclose(f);
    double total = status_meminfo_total();
    if (total < best) {
        best = total;
        snprintf(src, srclen, "MemTotal");
    }
    return best;
}

/* ------------------------------------------------------------------ configuration */
struct load_arg { const char *path; char *msg; int msglen; };

/* Rockstar's parser and its conversions (config_to_* may reach exit(): run through rs_guarded_call) */
static int load_body(void *v)
{
    struct load_arg *a = v;
    struct configfile c;
    memset(&c, 0, sizeof c);
    load_config(&c, (char *)a->path);
    if (c.num_entries == 0) {               /* read_config.c returns no entries when it cannot open the file */
        free_config(c);
        snprintf(a->msg, a->msglen, "%s: no entries read", a->path);
        return -1;
    }
#define string(a, b)  a = config_to_string(&c, #a, b)
#define real(a, b)    a = config_to_real(&c, #a, b)
#define real3(a, b)   config_to_real3(&c, #a, a, b)
#define integer(a, b) a = config_to_real(&c, #a, b)
#include "config.template.h"
#undef string
#undef real
#undef real3
#undef integer
    cfg_ceiling_gib = config_to_real(&c, "GLAM_MEMORY_CEILING_GIB", 900);
    cfg_strict = (int)config_to_real(&c, "GLAM_STRICT_SO", 1);
    syntax_check(&c, "[rockstar_inline warning]");
    free_config(c);
    return 0;
}

static int load_configuration(const char *path, char *msg, int msglen)
{
    if (!path || !*path) {
        snprintf(msg, msglen, "no Rockstar configuration given");
        return -1;
    }
    /* the bytes once, for the checksum and a clear message if the file cannot be read */
    FILE *f = fopen(path, "rb");
    if (!f) {
        snprintf(msg, msglen, "Rockstar configuration %s: %s", path, strerror(errno));
        return -1;
    }
    uint64_t hash = 1469598103934665603ULL;
    int ch;
    long nbytes = 0;
    while ((ch = getc(f)) != EOF) {
        hash = (hash ^ (uint64_t)(unsigned char)ch) * 1099511628211ULL;
        nbytes++;
    }
    int rerr = ferror(f);
    fclose(f);
    if (rerr || nbytes == 0) {
        snprintf(msg, msglen, "Rockstar configuration %s: %s", path, rerr ? "read error" : "empty");
        return -1;
    }
    cfg_path[0] = 0;
    struct load_arg a = {path, msg, msglen};
    if (rs_guarded_call(load_body, &a, path, msg, msglen) != 0) return -1;
    cfg_force_res = FORCE_RES;
    cfg_force_res_phys_max = FORCE_RES_PHYS_MAX ? FORCE_RES_PHYS_MAX : FORCE_RES;   /* config.c:35-36 */
    if (!(FORCE_RES > 0) || !(FOF_LINKING_LENGTH > 0 && FOF_LINKING_LENGTH <= 0.5) || MIN_HALO_PARTICLES < 2 ||
        !(FOF_FRACTION > 0 && FOF_FRACTION < 1) || MIN_HALO_OUTPUT_SIZE < 1 || !(cfg_ceiling_gib > 0) ||
        MIN_PART_PER_BIN < 1 || MIN_PART_PER_BIN > MIN_SCALE_PART || MAX_SCALE_BINS < 1 || MAX_SCALE_BINS > 1000) {
        snprintf(msg, msglen, "%s: FORCE_RES > 0, 0 < FOF_LINKING_LENGTH <= 0.5, 0 < FOF_FRACTION < 1, "
                              "MIN_HALO_PARTICLES >= 2, MIN_HALO_OUTPUT_SIZE >= 1, GLAM_MEMORY_CEILING_GIB > 0, "
                              "1 <= MIN_PART_PER_BIN <= MIN_SCALE_PART and 1 <= MAX_SCALE_BINS <= 1000 are required",
                 path);
        return -1;
    }
    if (STRICT_SO_MASSES || UNFILTERED_HALO_OUTPUT)
        fprintf(stdout, " [rockstar_inline warning] %s: STRICT_SO_MASSES and UNFILTERED_HALO_OUTPUT are not used "
                        "in-line (the strict pass is GLAM_STRICT_SO; the output filter always applies)\n", path);
    snprintf(cfg_path, sizeof cfg_path, "%s", path);
    cfg_hash = hash;
    return 0;
}

/* the start-up check: loads (and keeps) the configuration */
int rs_inline_check_config(const char *config_path, char *msg, int msglen)
{
    msg[0] = 0;
    return load_configuration(config_path, msg, msglen);
}

/* the effective Rockstar configuration (every key of config.template.h), for the .meta */
struct mbuf { char *s; size_t len, cap; int err; };
static void madd(struct mbuf *b, const char *fmt, ...);

static void meta_config(struct mbuf *mb)
{
#define string(a, b)  madd(mb, "config." #a " = %s\n", a)
#define real(a, b)    madd(mb, "config." #a " = %.9g\n", (double)(a))
#define real3(a, b)   madd(mb, "config." #a " = %.9g %.9g %.9g\n", a[0], a[1], a[2])
#define integer(a, b) madd(mb, "config." #a " = %lld\n", (long long)(a))
#include "config.template.h"
#undef string
#undef real
#undef real3
#undef integer
}

/* the values of this run (setup_config is not used, see above) */
static void set_run_globals(int64_t n, double box, double aexpn, double om, double oml, double hubble)
{
    BOX_SIZE = box;
    TOTAL_PARTICLES = n;
    Om = om;
    Ol = oml;
    h0 = hubble;
    SCALE_NOW = aexpn;
    W0 = -1;
    WA = 0;
    NON_COSMOLOGICAL = 0;
    TEMPORAL_HALO_FINDING = 0;
    LIGHTCONE = 0;
    OUTPUT_LEVELS = 0;
    PARALLEL_IO = 1;
    PERIODIC = 1;
    PARTICLE_MASS = CRITICAL_DENSITY * box * box * box * om / (double)n;      /* config.c:21-23 */
    AVG_PARTICLE_SPACING = cbrt(PARTICLE_MASS / (Om * CRITICAL_DENSITY));     /* config.c:24-25 */
    FORCE_RES = cfg_force_res;
    FORCE_RES_PHYS_MAX = cfg_force_res_phys_max;
    if (FORCE_RES * SCALE_NOW > FORCE_RES_PHYS_MAX)                           /* rockstar.c:60-61 */
        FORCE_RES = FORCE_RES_PHYS_MAX / SCALE_NOW;
}

/* ------------------------------------------------------------------ metadata */
static void madd(struct mbuf *b, const char *fmt, ...)      /* any failure sets b->err */
{
    va_list ap;
    if (b->err) return;
    for (;;) {
        va_start(ap, fmt);
        int w = vsnprintf(b->s ? b->s + b->len : NULL, b->s ? b->cap - b->len : 0, fmt, ap);
        va_end(ap);
        if (w < 0) {
            b->err = 1;
            return;
        }
        if (b->s && (size_t)w < b->cap - b->len) {
            b->len += (size_t)w;
            return;
        }
        size_t nc = (b->cap ? 2 * b->cap : 16384) + (size_t)w;
        char *t = realloc(b->s, nc);
        if (!t) {
            b->err = 1;
            return;
        }
        b->s = t;
        b->cap = nc;
    }
}

static void madd_columns(struct mbuf *b)
{
    char cols[32768];
    if (rs_columns_text(cols, sizeof cols, (int)BOUND_PROPS) > 0) madd(b, "%s", cols);
    else b->err = 1;
}

static void meta_common(struct mbuf *b, const char *status, const char *reason, int step, int box_index)
{
    madd(b, "format = raw little-endian records, %zu bytes per halo: Rockstar's struct halo (312 bytes, "
            "MPI-Rockstar halo.h with OUTPUT_RVMAX, OUTPUT_NFW_CHI2, OUTPUT_INERTIA_TENSOR) followed by the "
            "adapter's columns; no header; columns listed below\n", sizeof(struct rs_record));
    madd(b, "record_bytes = %zu\n", sizeof(struct rs_record));
    madd(b, "status = %s\n", status);
    madd(b, "reason = %s\n", reason && *reason ? reason : "-");
    madd(b, "step = %d\nbox_index = %d\n", step, box_index);
    madd(b, "finder = in-line Rockstar of MG-GLAM (rockstar_inline/), %s\n", rs_inline_version());
    if (run_ctx.set) {
        const double a_v = run_ctx.aexpn - 0.5 * run_ctx.astep;
        madd(b, "a = %.8f\nz = %.8f\n", run_ctx.aexpn, 1.0 / run_ctx.aexpn - 1.0);
        madd(b, "a_v = %.8f\nvelocity_epoch = a_v = AEXPN - ASTEP/2 (GLAM's momenta lag the positions by half "
                "a step); v [km/s] = V * 100 * (Box/Ngrid) / a_v; ASTEP = %.8f\n", a_v, run_ctx.astep);
        madd(b, "Box = %.9g\nNgrid = %lld\nNparticles = %lld\n", run_ctx.box, (long long)run_ctx.ngrid,
             (long long)run_ctx.n);
    }
    if (cfg_path[0]) madd(b, "config_file = %s\nconfig_fnv1a64 = %016llx\n", cfg_path, (unsigned long long)cfg_hash);
}

/* the sidecar of a skipped or failed output; an existing complete catalogue of the step (a restart
 * that repeats an output) is kept instead */
static int publish_status(const char *dir, int step, int box_index, const char *status, const char *reason)
{
    if (rs_catalogue_ok(dir, step, box_index)) {
        fprintf(stdout, " Rockstar step %d: %s; keeping the existing complete catalogue of this step\n", step,
                status);
        return 0;
    }
    struct mbuf b = {0};
    char err[512] = "metadata: out of memory";
    meta_common(&b, status, reason, step, box_index);
    madd(&b, "rows = 0\n");
    madd_columns(&b);
    int rc = (b.s && !b.err) ? rs_write_catalogue(dir, step, box_index, NULL, 0, b.s, err, sizeof err) : -1;
    if (rc != 0) fprintf(stdout, " Rockstar step %d: cannot publish the %s catalogue: %s\n", step, status, err);
    free(b.s);
    return rc;
}

int rs_publish_failure(const char *catalog_dir, int step, int box_index, const char *reason)
{
    return publish_status(catalog_dir, step, box_index, "failed", reason);
}

static int finish(struct rs_inline_info *info, const char *dir, int step, int box_index, int status,
                  const char *reason, double t_start)
{
    info->status = status;
    snprintf(info->reason, sizeof info->reason, "%s", reason);
    info->t_total = omp_get_wtime() - t_start;
    p = NULL;
    malloc_trim(0);
    publish_status(dir, step, box_index, status == RS_SKIPPED ? "skipped" : "failed", reason);
    fprintf(stdout, " Rockstar step %d: %s (%s); %.1f s\n", step, status == RS_SKIPPED ? "SKIPPED" : "FAILED",
            reason, info->t_total);
    fflush(stdout);
    return status;
}

/* ------------------------------------------------------------------ kernel helpers */
static int cmp_particle_id(const void *a, const void *b)
{
    int64_t c = ((const struct particle *)a)->id, d = ((const struct particle *)b)->id;
    return c < d ? -1 : (c > d);
}

/* meta_io.c _should_print without the bounds (the centres are wrapped into the box) */
static int should_print(const struct halo *h)
{
    if (h->flags & ALWAYS_PRINT_FLAG) return 1;
    if ((h->num_p < MIN_HALO_OUTPUT_SIZE) || (h->m * UNBOUND_THRESHOLD >= h->mgrav) ||
        ((h->mgrav < 1.5 * PARTICLE_MASS) && UNBOUND_THRESHOLD > 0))
        return 0;
    return 1;
}

/* the group of internal halo i: gstart is non-decreasing, group g holds [gstart[g], gstart[g+1]) */
static int64_t group_of_halo(const int64_t *gstart, int64_t ng, int64_t i)
{
    int64_t lo = 0, hi = ng;
    while (hi - lo > 1) {
        int64_t mid = (lo + hi) / 2;
        if (gstart[mid] <= i) lo = mid; else hi = mid;
    }
    return lo;
}

/* ------------------------------------------------------------------ the run */
int rs_inline_run_impl(int64_t n, const float *x, const float *y, const float *z,
                       const float *vx, const float *vy, const float *vz, double box, int64_t ngrid,
                       double aexpn, double astep, double om, double oml, double hubble, int step,
                       int box_index, const char *config_path, const char *catalog_dir,
                       struct rs_inline_info *info)
{
    struct rs_inline_info local;
    if (!info) info = &local;
    memset(info, 0, sizeof *info);
    const double t_start = omp_get_wtime();
    info->rss_entry_gib = info->rss_max_gib = rs_vmrss_gib();
    const double hwm_entry = status_gib("VmHWM:");
    run_ctx.set = 1;
    run_ctx.aexpn = aexpn;
    run_ctx.astep = astep;
    run_ctx.box = box;
    run_ctx.ngrid = ngrid;
    run_ctx.n = n;
    char msg[512] = "";
#ifdef RS_INLINE_TEST_HOOKS     /* the harness's build only (makefile): kernel | assert | alloc */
    const char *test = getenv("RS_INLINE_TEST_FAIL");
    const int inject = test ? (!strcmp(test, "kernel") ? 1 : !strcmp(test, "assert") ? 2 : 0) : 0;
    const int fail_alloc = test && !strcmp(test, "alloc");
#else
    const int inject = 0, fail_alloc = 0;
#endif

    if ((!cfg_path[0] || strcmp(cfg_path, config_path ? config_path : "")) &&
        load_configuration(config_path, msg, sizeof msg) != 0)
        return finish(info, catalog_dir, step, box_index, RS_FAILED, msg, t_start);
    set_run_globals(n, box, aexpn, om, oml, hubble);
    struct rs_units u;
    rs_units_init(&u, box, ngrid, aexpn, astep);
    const int64_t nrow = llround(cbrt((double)n));
    const int T = omp_get_max_threads();
    char limit_src[8400];
    const double limit_gib = memory_limit_bytes(limit_src, sizeof limit_src) / GIB;
    const double ceiling = cfg_ceiling_gib < 0.97 * limit_gib ? cfg_ceiling_gib : 0.97 * limit_gib;

    /* ---- FoF ---- */
    info->need_fof_gib = rs_fof_peak_bytes(n, nrow, T) / GIB;
    if (info->rss_entry_gib + info->need_fof_gib > ceiling) {
        snprintf(msg, sizeof msg, "memory pre-check: VmRSS %.1f GiB + FoF %.1f GiB > ceiling %.1f GiB "
                 "(GLAM_MEMORY_CEILING_GIB %.1f, 97 %% of the %.1f GiB limit)", info->rss_entry_gib,
                 info->need_fof_gib, ceiling, cfg_ceiling_gib, limit_gib);
        return finish(info, catalog_dir, step, box_index, RS_SKIPPED, msg, t_start);
    }
    double t0 = omp_get_wtime();
    struct rs_cellindex cells;
    struct rs_fof_groups groups;
    int frc = rs_fof_run(n, x, y, z, &u, nrow, AVG_PARTICLE_SPACING, FOF_LINKING_LENGTH, MIN_HALO_PARTICLES, 0,
                         &cells, &groups);
    if (frc != RS_FOF_OK) {
        if (frc == RS_FOF_ECOORD)
            snprintf(msg, sizeof msg, "FoF: %lld particles with coordinates not finite or outside GLAM's range",
                     (long long)groups.nbad);
        else
            snprintf(msg, sizeof msg, "FoF: %s", rs_fof_strerror(frc));
        return finish(info, catalog_dir, step, box_index, RS_FAILED, msg, t_start);
    }
    info->t_fof = omp_get_wtime() - t0;
    note_rss(info, groups.rss_link_gib);
    sample_rss(info);
    const int64_t ng = groups.ngroups, nm = groups.nmembers;
    info->ngroups = ng;
    info->nmembers = nm;
    info->largest = groups.largest;

    /* ---- second pre-check: the member buffer (32 B per member), per-thread kernel scratch (~134 B per
     *      particle of the largest groups in flight, from the sizes of MPI-Rockstar's per-particle arrays),
     *      and the halos, twice (the threads' lists with their extra_halo_info, and the merged copy): ~1
     *      internal halo per 68 members measured at GLAM's resolution, 1 per 50 assumed ---- */
    {
        double inflight = (double)(ng < T ? ng : T) * (double)groups.largest;
        if (inflight > (double)nm) inflight = (double)nm;
        double need = 32.0 * nm + 140.0 * inflight +
                      ((double)nm / 50.0 + 1.0) * 2.0 * (double)(sizeof(struct halo) + sizeof(struct extra_halo_info));
        info->need_kernel_gib = need / GIB;
        double rss = rs_vmrss_gib();
        if (rss + info->need_kernel_gib > ceiling) {
            snprintf(msg, sizeof msg, "memory pre-check: VmRSS %.1f GiB + kernel %.1f GiB > ceiling %.1f GiB "
                     "(GLAM_MEMORY_CEILING_GIB %.1f, 97 %% of the %.1f GiB limit)", rss, info->need_kernel_gib,
                     ceiling, cfg_ceiling_gib, limit_gib);
            rs_fof_free_groups(&groups);
            rs_fof_free_cells(&cells);
            return finish(info, catalog_dir, step, box_index, RS_SKIPPED, msg, t_start);
        }
    }

    /* ---- gather ---- */
    t0 = omp_get_wtime();
    struct particle *buf = fail_alloc ? NULL : malloc((size_t)(nm > 0 ? nm : 1) * sizeof(struct particle));
    int64_t *cursor = malloc((size_t)(ng > 0 ? ng : 1) * sizeof(int64_t));
    if (!buf || !cursor) {
        free(buf);
        free(cursor);
        rs_fof_free_groups(&groups);
        rs_fof_free_cells(&cells);
        return finish(info, catalog_dir, step, box_index, RS_FAILED, "gather: out of memory", t_start);
    }
    memcpy(cursor, groups.offset, (size_t)ng * sizeof(int64_t));
#pragma omp parallel for schedule(static)
    for (int64_t k = 0; k < n; k++) {
        int64_t gi = groups.group_of_sorted[k];
        if (gi < 0) continue;
        int64_t slot;
#pragma omp atomic capture
        slot = cursor[gi]++;
        int64_t i = cells.perm[k];
        buf[slot].id = i;
        rs_convert(&u, x[i], y[i], z[i], vx[i], vy[i], vz[i], buf[slot].pos);
    }
    free(cursor);
    sample_rss(info);                                   /* the call's peak: cell index, labels and buffer */
    rs_fof_release_labels(&groups);
    /* each group by particle id, unwrapped around its first member (minimum image): a group must span
     * less than half the box minus the linking length, or it would be cut */
    const float boxf = (float)box, half = 0.5f * (float)box;
    const float max_span = half - rs_fof_linking_length(AVG_PARTICLE_SPACING, FOF_LINKING_LENGTH);
    int64_t wide_group = -1;
#pragma omp parallel for schedule(dynamic, 16)
    for (int64_t gi = 0; gi < ng; gi++) {
        struct particle *a = buf + groups.offset[gi];
        int64_t na = groups.offset[gi + 1] - groups.offset[gi];
        qsort(a, (size_t)na, sizeof(struct particle), cmp_particle_id);
        float lo[3], hi[3];
        for (int d = 0; d < 3; d++) lo[d] = hi[d] = a[0].pos[d];
        for (int64_t j = 1; j < na; j++)
            for (int d = 0; d < 3; d++) {
                float dx = a[j].pos[d] - a[0].pos[d];
                if (dx > half) a[j].pos[d] -= boxf;
                else if (dx < -half) a[j].pos[d] += boxf;
                if (a[j].pos[d] < lo[d]) lo[d] = a[j].pos[d];
                if (a[j].pos[d] > hi[d]) hi[d] = a[j].pos[d];
            }
        if (hi[0] - lo[0] >= max_span || hi[1] - lo[1] >= max_span || hi[2] - lo[2] >= max_span) {
#pragma omp critical(rs_inline_wide)
            if (wide_group < 0 || gi < wide_group) wide_group = gi;
        }
    }
    info->t_gather = omp_get_wtime() - t0;
    if (wide_group >= 0) {
        snprintf(msg, sizeof msg, "FoF group %lld (%lld particles) spans half the box or more", (long long)wide_group,
                 (long long)(groups.offset[wide_group + 1] - groups.offset[wide_group]));
        free(buf);
        rs_fof_free_groups(&groups);
        rs_fof_free_cells(&cells);
        return finish(info, catalog_dir, step, box_index, RS_FAILED, msg, t_start);
    }

    /* ---- kernel ----
     * One failing group abandons the whole output, on purpose: after a throw inside find_subs, Rockstar's
     * capacity counters (num_alloc_gh, num_alloced_halo_ids, alloced_metrics, fast3tree's allocated sizes)
     * can be ahead of the memory really allocated, so a thread must not call find_subs again until its
     * per-thread scratch is freed, which happens at the end of the parallel region. */
    t0 = omp_get_wtime();
    p = buf;
    calc_mass_definition();
    int64_t *blk_start = malloc((size_t)(ng > 0 ? ng : 1) * sizeof(int64_t));
    int64_t *blk_count = malloc((size_t)(ng > 0 ? ng : 1) * sizeof(int64_t));
    int64_t *gstart = malloc((size_t)(ng + 1) * sizeof(int64_t));
    struct halo *all_h = NULL;
    int failed = !blk_start || !blk_count || !gstart;
    char fail_msg[512] = "kernel bookkeeping: out of memory";
    const int64_t inject_group = inject ? ng / 2 : -1;
    int64_t total = 0;
#pragma omp parallel copyin(particle_thresh_dens, particle_rvir_dens, particle_rvir_dens_z0, min_dens_index, \
                            dynamical_time, scale_dx)
    {
        struct FOFInfo fofinfo;
        struct HaloInfo hinfo;
        init_fofinfo(&fofinfo);
        init_haloinfo(&hinfo);
        int64_t *mine = NULL, nmine = 0, capmine = 0;
#pragma omp for schedule(dynamic, 1)
        for (int64_t gi = 0; gi < ng; gi++) {
            if (__atomic_load_n(&failed, __ATOMIC_RELAXED)) continue;
            struct fof f;
            f.num_p = groups.offset[gi + 1] - groups.offset[gi];
            f.particles = buf + groups.offset[gi];
            int64_t h0 = hinfo.num_halos;
            char m2[400];
            if (rs_guarded_find_subs(&f, &fofinfo, &hinfo, gi == inject_group ? inject : 0, m2, sizeof m2) != 0) {
#pragma omp critical(rs_inline_fail)
                {
                    if (!failed) snprintf(fail_msg, sizeof fail_msg, "FoF group %lld (%lld particles): %s",
                                          (long long)gi, (long long)f.num_p, m2);
                    __atomic_store_n(&failed, 1, __ATOMIC_RELAXED);
                }
                continue;
            }
            blk_start[gi] = h0;
            blk_count[gi] = hinfo.num_halos - h0;
            if (nmine == capmine) {
                int64_t c2 = capmine ? 2 * capmine : 1024;
                int64_t *t = realloc(mine, (size_t)c2 * sizeof(int64_t));
                if (!t) {
#pragma omp critical(rs_inline_fail)
                    {
                        if (!failed) snprintf(fail_msg, sizeof fail_msg, "kernel bookkeeping: out of memory");
                        __atomic_store_n(&failed, 1, __ATOMIC_RELAXED);
                    }
                    continue;
                }
                mine = t;
                capmine = c2;
            }
            mine[nmine++] = gi;
        }
#pragma omp single
        {
            if (!failed) {
                gstart[0] = 0;
                for (int64_t gi = 0; gi < ng; gi++) gstart[gi + 1] = gstart[gi] + blk_count[gi];
                total = gstart[ng];
                all_h = malloc((size_t)(total > 0 ? total : 1) * sizeof(struct halo));
                if (!all_h) {
                    snprintf(fail_msg, sizeof fail_msg, "merge: out of memory");
                    failed = 1;
                }
            }
        }
        if (!failed)
            for (int64_t q = 0; q < nmine; q++) {
                const int64_t gi = mine[q];
                memcpy(all_h + gstart[gi], hinfo.halos + blk_start[gi], (size_t)blk_count[gi] * sizeof(struct halo));
            }
#pragma omp barrier
#pragma omp master
        sample_rss(info);                   /* the per-thread lists and the merged copy */
#pragma omp barrier
        free(mine);
        free_haloinfo(&hinfo);              /* the kernel's per-thread copies, trees, halo lists */
        free_fofinfo(&fofinfo);
        rs_free_potential_scratch();        /* scratch MPI-Rockstar keeps until exit */
        rs_free_subhalo_scratch();
    }
    p = NULL;
    free(buf);
    free(blk_start);
    free(blk_count);
    info->t_kernel = omp_get_wtime() - t0;
    if (failed) {
        free(all_h);
        free(gstart);
        rs_fof_free_groups(&groups);
        rs_fof_free_cells(&cells);
        return finish(info, catalog_dir, step, box_index, RS_FAILED, fail_msg, t_start);
    }
    info->nhalos_internal = total;

    /* ---- selection: Rockstar's output filter, then the distinct haloes ---- */
    int64_t *row = malloc((size_t)(total > 0 ? total : 1) * sizeof(int64_t));
    if (!row) {
        free(all_h); free(gstart);
        rs_fof_free_groups(&groups);
        rs_fof_free_cells(&cells);
        return finish(info, catalog_dir, step, box_index, RS_FAILED, "selection: out of memory", t_start);
    }
    int64_t nsel = 0, outside = 0;
    for (int64_t i = 0; i < total; i++) {
        wrap_into_box(all_h[i].pos);
        for (int d = 0; d < 3; d++) {                   /* wrap_into_box keeps pos = Box: map it to 0 */
            if (all_h[i].pos[d] == boxf) all_h[i].pos[d] = 0.0f;
            if (!(all_h[i].pos[d] >= 0.0f && all_h[i].pos[d] < boxf)) outside++;
        }
        row[i] = should_print(all_h + i) ? nsel++ : -1;
    }
    struct halo *sel = outside ? NULL : malloc((size_t)(nsel > 0 ? nsel : 1) * sizeof(struct halo));
    int64_t *selgroup = malloc((size_t)(nsel > 0 ? nsel : 1) * sizeof(int64_t));
    int64_t *pid = malloc((size_t)(nsel > 0 ? nsel : 1) * sizeof(int64_t));
    if (!sel || !selgroup || !pid) {
        free(sel); free(selgroup); free(pid);
        free(row); free(all_h); free(gstart);
        rs_fof_free_groups(&groups);
        rs_fof_free_cells(&cells);
        if (outside) snprintf(msg, sizeof msg, "%lld halo centre coordinates not finite or outside the box",
                              (long long)outside);
        return finish(info, catalog_dir, step, box_index, RS_FAILED, outside ? msg : "selection: out of memory",
                      t_start);
    }
#pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < total; i++) {
        if (row[i] < 0) continue;
        sel[row[i]] = all_h[i];
        selgroup[row[i]] = group_of_halo(gstart, ng, i);
    }
    free(row);
    free(all_h);
    free(gstart);
    rs_fof_free_groups(&groups);
    info->nhalos = nsel;

    t0 = omp_get_wtime();
    if (rs_guarded_find_parents(nsel, sel, pid, msg, sizeof msg) != 0) {
        free(sel); free(selgroup); free(pid);
        rs_fof_free_cells(&cells);
        return finish(info, catalog_dir, step, box_index, RS_FAILED, msg, t_start);
    }
    int64_t nhost = 0;
    for (int64_t i = 0; i < nsel; i++) if (pid[i] < 0) nhost++;
    struct rs_record *rec = calloc((size_t)(nhost > 0 ? nhost : 1), sizeof(struct rs_record));
    if (!rec) {
        free(sel); free(selgroup); free(pid);
        rs_fof_free_cells(&cells);
        return finish(info, catalog_dir, step, box_index, RS_FAILED, "records: out of memory", t_start);
    }
    for (int64_t i = 0, k = 0; i < nsel; i++) {
        if (pid[i] >= 0) continue;                      /* a subhalo: not written */
        struct rs_record *r = rec + k;
        r->h = sel[i];
        r->h.id = k;
        r->h.desc = -1;
        r->h.p_start = -1;
        r->fof_group = selgroup[i];
        r->m_strict = r->r_strict = NAN;
        for (int q = 0; q < 4; q++) r->alt_m_strict[q] = NAN;
        k++;
    }
    free(sel);
    free(selgroup);
    free(pid);
    info->nhosts = nhost;
    info->t_parents = omp_get_wtime() - t0;

    /* ---- strict SO masses of the distinct haloes ---- */
    t0 = omp_get_wtime();
    int strict_rc = 0;
    if (cfg_strict) strict_rc = rs_strict_pass(nhost, rec, x, y, z, &u, &cells);
    rs_fof_free_cells(&cells);
    info->t_strict = omp_get_wtime() - t0;
    if (strict_rc != 0) {
        free(rec);
        return finish(info, catalog_dir, step, box_index, RS_FAILED, "strict SO pass: out of memory", t_start);
    }

    /* ---- publication ---- */
    t0 = omp_get_wtime();
    const double hwm_exit = status_gib("VmHWM:");
    struct mbuf b = {0};
    char err[512] = "metadata: out of memory";
    meta_common(&b, "ok", "", step, box_index);
    madd(&b, "rows = %lld\nhalos_selected = %lld (after Rockstar's output filter; the %lld subhaloes among them are "
             "not written)\n", (long long)nhost, (long long)nsel, (long long)(nsel - nhost));
    madd(&b, "positions = x [Mpc/h comoving] = (x_grid - 1) * Box/Ngrid, wrapped into [0, Box); particle id = "
             "GLAM array index - 1\n");
    madd(&b, "Om = %.9g\nOl = %.9g\nh0 = %.9g\n", Om, Ol, h0);
    madd(&b, "fof_r_link = %.9g\nfof = exact periodic friends-of-friends of the adapter (rs_fof.c), d^2 < r^2 in "
             "float; Rockstar's own FoF additionally merges groups within 2r of particles with > 16 neighbours\n",
         (double)rs_fof_linking_length(AVG_PARTICLE_SPACING, FOF_LINKING_LENGTH));
    madd(&b, "FORCE_RES_configured = %.9g (config.FORCE_RES is the value used, capped at FORCE_RES_PHYS_MAX / a)\n",
         cfg_force_res);
    madd(&b, "bound_props = with BOUND_PROPS = 1 Rockstar computes every property except m from the bound "
             "particles (r, the alt_m_*, vmax, vrms, J, spins, shapes, ...)\n");
    madd(&b, "strict_so = %d\nselection = distinct haloes only: no parent in Rockstar's find_parents with the kernel "
             "radius r (ties in vmax taken in row order). Rockstar's kernel still finds substructure, and a halo's "
             "masses and properties include the particles of its substructure (num_child_particles)\n", cfg_strict);
    madd(&b, "member_order = particle id (each FoF group is handed to find_subs sorted by particle id)\n");
    madd(&b, "fof_groups = %lld\nfof_members = %lld\nfof_fraction = %.6f\nfof_largest = %lld\n", (long long)ng,
         (long long)nm, n > 0 ? (double)nm / (double)n : 0.0, (long long)info->largest);
    madd(&b, "halos_internal = %lld\n", (long long)total);
    madd(&b, "threads = %d\n", T);
    madd(&b, "time_fof_s = %.2f\ntime_gather_s = %.2f\ntime_kernel_s = %.2f\ntime_strict_s = %.2f\n"
             "time_parents_s = %.2f\n", info->t_fof, info->t_gather, info->t_kernel, info->t_strict,
         info->t_parents);
    madd(&b, "vmrss_entry_gib = %.3f\nvmrss_max_sampled_gib = %.3f\nvmhwm_entry_gib = %.3f\nvmhwm_exit_gib = %.3f "
             "(above the entry value only if this call set the process's peak)\n",
         info->rss_entry_gib, info->rss_max_gib, hwm_entry, hwm_exit);
    madd(&b, "need_fof_gib = %.3f\nneed_kernel_gib = %.3f\nmemory_ceiling_gib = %.1f\n"
             "memory_ceiling_configured_gib = %.1f\nmemory_limit_gib = %.1f (%s)\n", info->need_fof_gib,
         info->need_kernel_gib, ceiling, cfg_ceiling_gib, limit_gib, limit_src);
    meta_config(&b);
    madd_columns(&b);
    int wrc = (b.s && !b.err) ? rs_write_catalogue(catalog_dir, step, box_index, rec, nhost, b.s, err, sizeof err)
                              : -1;
    free(b.s);
    free(rec);
    info->t_write = omp_get_wtime() - t0;
    if (wrc != 0) return finish(info, catalog_dir, step, box_index, RS_FAILED, err, t_start);
    malloc_trim(0);
    info->t_total = omp_get_wtime() - t_start;
    info->status = RS_OK;
    fprintf(stdout, " Rockstar step %d: %lld distinct haloes (of %lld selected) from %lld FoF groups (f = %.4f, "
                    "largest %lld); FoF %.1f s, gather %.1f s, kernel %.1f s, strict %.1f s, parents %.1f s, "
                    "write %.1f s, total %.1f s; VmRSS %.2f -> max %.2f GiB\n",
            step, (long long)nhost, (long long)nsel, (long long)ng, n > 0 ? (double)nm / (double)n : 0.0,
            (long long)info->largest, info->t_fof, info->t_gather, info->t_kernel, info->t_strict, info->t_parents,
            info->t_write, info->t_total, info->rss_entry_gib, info->rss_max_gib);
    fflush(stdout);
    return RS_OK;
}
