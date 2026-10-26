/* rs_inline.c: Rockstar haloes of one MG-GLAM output, in process.
 * Part of the MG-GLAM in-line Rockstar adapter; GNU GPL v3 or later. The kernel
 * loop and merge follow MPI-Rockstar's do_halo_finding / copy_halos
 * (src/mpi_main.cpp, commit 0ae386f), the output filter its _should_print
 * (src/io/meta_io.c).
 *
 * One call (rs_inline_run, through the guard in rs_guard.cpp):
 *  1. configuration: Rockstar's parser on Run<box>/Rockstar.config, without
 *     setup_config (which takes the particle spacing from Rockstar's defaults
 *     and raises the core-dump limit); the run's values are set explicitly;
 *  2. memory pre-check, FoF (rs_fof.c), second pre-check;
 *  3. gather: the members of every group into Rockstar's array-of-structs
 *     buffer (the global p), each group sorted by particle id (the kernel's
 *     result depends on the member order) and unwrapped around its first member;
 *  4. kernel: find_subs on every group in an OpenMP loop (largest first), each
 *     call guarded; halos merged in group order, so the catalogue does not
 *     depend on the thread count; per-thread scratch released on every thread;
 *  5. halo centres wrapped, Rockstar's output filter, strict SO masses
 *     (rs_strict.c), parents (rs_parents.c), publication (rs_write.c).
 * Failures and skips publish an empty catalogue whose .meta says why. */
#define _GNU_SOURCE
#include <errno.h>
#include <malloc.h>
#include <math.h>
#include <omp.h>
#include <stdarg.h>
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
int rs_guarded_load_config(struct configfile *c, const char *path, char *msg, int msglen);

static double cfg_ceiling_gib = 900;    /* GLAM_MEMORY_CEILING_GIB */
static int    cfg_strict = 1;           /* GLAM_STRICT_SO */
static double cfg_force_res, cfg_force_res_phys_max;

int rs_inline_available(void) { return 1; }

const char *rs_inline_version(void)
{
    return "MPI-Rockstar " RS_ROCKSTAR_COMMIT "; " RS_BUILD_FLAGS;
}

static double vmrss_gib(void)
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

static void sample_rss(struct rs_inline_info *info)
{
    double r = vmrss_gib();
    if (r > info->rss_max_gib) info->rss_max_gib = r;
}

/* ------------------------------------------------------------------ configuration */
static int load_configuration(const char *path, char *msg, int msglen)
{
    struct stat st;
    if (!path || !*path || stat(path, &st) != 0) {
        snprintf(msg, msglen, "Rockstar configuration %s not found", path ? path : "(none)");
        return -1;
    }
    struct configfile c;
    memset(&c, 0, sizeof c);
    if (rs_guarded_load_config(&c, path, msg, msglen) != 0) return -1;
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
    cfg_force_res = FORCE_RES;
    cfg_force_res_phys_max = FORCE_RES_PHYS_MAX ? FORCE_RES_PHYS_MAX : FORCE_RES;   /* config.c:35-36 */
    if (!(FORCE_RES > 0) || !(FOF_LINKING_LENGTH > 0 && FOF_LINKING_LENGTH < 1) || MIN_HALO_PARTICLES < 2 ||
        !(FOF_FRACTION > 0 && FOF_FRACTION < 1) || MIN_HALO_OUTPUT_SIZE < 1 || !(cfg_ceiling_gib > 0)) {
        snprintf(msg, msglen, "%s: FORCE_RES > 0, 0 < FOF_LINKING_LENGTH < 1, 0 < FOF_FRACTION < 1, "
                              "MIN_HALO_PARTICLES >= 2, MIN_HALO_OUTPUT_SIZE >= 1 and GLAM_MEMORY_CEILING_GIB > 0 "
                              "are required", path);
        return -1;
    }
    return 0;
}

int rs_inline_check_config(const char *config_path, char *msg, int msglen)
{
    msg[0] = 0;
    return load_configuration(config_path, msg, msglen);
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
struct mbuf { char *s; size_t len, cap; };

static void madd(struct mbuf *b, const char *fmt, ...)
{
    va_list ap;
    for (;;) {
        va_start(ap, fmt);
        int w = vsnprintf(b->s ? b->s + b->len : NULL, b->s ? b->cap - b->len : 0, fmt, ap);
        va_end(ap);
        if (w < 0) return;
        if (b->s && (size_t)w < b->cap - b->len) {
            b->len += (size_t)w;
            return;
        }
        size_t nc = (b->cap ? 2 * b->cap : 16384) + (size_t)w;
        char *t = realloc(b->s, nc);
        if (!t) return;
        b->s = t;
        b->cap = nc;
    }
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
}

static int publish_status(const char *dir, int step, int box_index, const char *status, const char *reason)
{
    struct mbuf b = {0};
    char cols[16384], err[512];
    meta_common(&b, status, reason, step, box_index);
    madd(&b, "rows = 0\n");
    if (rs_columns_text(cols, sizeof cols) > 0) madd(&b, "%s", cols);
    int rc = b.s ? rs_write_catalogue(dir, step, box_index, NULL, 0, b.s, err, sizeof err) : -1;
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

static int fix_link(int64_t *link, int64_t bstart, int64_t bcount, int64_t gstart)
{
    if (*link == -1) return 0;
    if (*link < bstart || *link >= bstart + bcount) return -1;    /* links stay inside a FoF group */
    *link = *link - bstart + gstart;
    return 0;
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
    info->rss_entry_gib = info->rss_max_gib = vmrss_gib();
    char msg[512] = "";
    const char *test = getenv("RS_INLINE_TEST_FAIL");    /* test hooks: kernel | assert | alloc */
    const int inject = test ? (!strcmp(test, "kernel") ? 1 : !strcmp(test, "assert") ? 2 : 0) : 0;
    const int fail_alloc = test && !strcmp(test, "alloc");

    if (load_configuration(config_path, msg, sizeof msg) != 0)
        return finish(info, catalog_dir, step, box_index, RS_FAILED, msg, t_start);
    set_run_globals(n, box, aexpn, om, oml, hubble);
    struct rs_units u;
    rs_units_init(&u, box, ngrid, aexpn, astep);
    const int64_t nrow = llround(cbrt((double)n));
    const int T = omp_get_max_threads();

    /* ---- FoF ---- */
    info->need_fof_gib = rs_fof_peak_bytes(n, nrow, T) / GIB;
    if (info->rss_entry_gib + info->need_fof_gib > cfg_ceiling_gib) {
        snprintf(msg, sizeof msg, "memory pre-check: VmRSS %.1f GiB + FoF %.1f GiB > GLAM_MEMORY_CEILING_GIB %.1f",
                 info->rss_entry_gib, info->need_fof_gib, cfg_ceiling_gib);
        return finish(info, catalog_dir, step, box_index, RS_SKIPPED, msg, t_start);
    }
    double t0 = omp_get_wtime();
    struct rs_cellindex cells;
    struct rs_fof_groups groups;
    if (rs_fof_run(n, x, y, z, &u, nrow, AVG_PARTICLE_SPACING, FOF_LINKING_LENGTH, MIN_HALO_PARTICLES, 0,
                   &cells, &groups) != 0)
        return finish(info, catalog_dir, step, box_index, RS_FAILED, "FoF: out of memory", t_start);
    info->t_fof = omp_get_wtime() - t0;
    sample_rss(info);
    const int64_t ng = groups.ngroups, nm = groups.nmembers;
    info->ngroups = ng;
    info->nmembers = nm;
    info->largest = groups.largest;

    /* ---- second pre-check: buffer, per-thread kernel scratch (~120 B per particle of the largest
     *      groups in flight, measured in MPI-Rockstar's allocations), halos (~1 per 100 members) ---- */
    {
        double inflight = (double)(ng < T ? ng : T) * (double)groups.largest;
        if (inflight > (double)nm) inflight = (double)nm;
        double need = 32.0 * nm + 120.0 * inflight + ((double)nm / 100.0 + 1.0) * (2.0 * 360.0 + 376.0);
        info->need_kernel_gib = need / GIB;
        double rss = vmrss_gib();
        if (rss + info->need_kernel_gib > cfg_ceiling_gib) {
            snprintf(msg, sizeof msg, "memory pre-check: VmRSS %.1f GiB + kernel %.1f GiB > "
                     "GLAM_MEMORY_CEILING_GIB %.1f", rss, info->need_kernel_gib, cfg_ceiling_gib);
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
    rs_fof_release_labels(&groups);
    const float boxf = (float)box, half = 0.5f * (float)box;
#pragma omp parallel for schedule(dynamic, 16)
    for (int64_t gi = 0; gi < ng; gi++) {
        struct particle *a = buf + groups.offset[gi];
        int64_t na = groups.offset[gi + 1] - groups.offset[gi];
        qsort(a, (size_t)na, sizeof(struct particle), cmp_particle_id);
        for (int64_t j = 1; j < na; j++)
            for (int d = 0; d < 3; d++) {
                float dx = a[j].pos[d] - a[0].pos[d];
                if (dx > half) a[j].pos[d] -= boxf;
                else if (dx < -half) a[j].pos[d] += boxf;
            }
    }
    info->t_gather = omp_get_wtime() - t0;
    sample_rss(info);

    /* ---- kernel ---- */
    t0 = omp_get_wtime();
    p = buf;
    calc_mass_definition();
    int64_t *blk_start = malloc((size_t)(ng > 0 ? ng : 1) * sizeof(int64_t));
    int64_t *blk_count = malloc((size_t)(ng > 0 ? ng : 1) * sizeof(int64_t));
    int64_t *gstart = malloc((size_t)(ng + 1) * sizeof(int64_t));
    struct halo *all_h = NULL;
    struct extra_halo_info *all_e = NULL;
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
            if (failed) continue;
            struct fof f;
            f.num_p = groups.offset[gi + 1] - groups.offset[gi];
            f.particles = buf + groups.offset[gi];
            int64_t h0 = hinfo.num_halos;
            char m2[300];
            if (rs_guarded_find_subs(&f, &fofinfo, &hinfo, gi == inject_group ? inject : 0, m2, sizeof m2) != 0) {
#pragma omp critical(rs_inline_fail)
                {
                    if (!failed) snprintf(fail_msg, sizeof fail_msg, "FoF group %lld (%lld particles): %s",
                                          (long long)gi, (long long)f.num_p, m2);
                    failed = 1;
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
                        failed = 1;
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
                all_e = malloc((size_t)(total > 0 ? total : 1) * sizeof(struct extra_halo_info));
                if (!all_h || !all_e) {
                    snprintf(fail_msg, sizeof fail_msg, "merge: out of memory");
                    failed = 1;
                }
            }
        }
        if (!failed) {
            for (int64_t q = 0; q < nmine; q++) {
                const int64_t gi = mine[q], bs = blk_start[gi], bc = blk_count[gi], gs = gstart[gi];
                memcpy(all_h + gs, hinfo.halos + bs, (size_t)bc * sizeof(struct halo));
                for (int64_t j = 0; j < bc; j++) {
                    struct extra_halo_info e = hinfo.extra_info[bs + j];
                    if (fix_link(&e.child, bs, bc, gs) || fix_link(&e.next_cochild, bs, bc, gs) ||
                        fix_link(&e.prev_cochild, bs, bc, gs) || fix_link(&e.sub_of, bs, bc, gs)) {
#pragma omp critical(rs_inline_fail)
                        {
                            if (!failed) snprintf(fail_msg, sizeof fail_msg,
                                                  "merge: a halo link leaves FoF group %lld", (long long)gi);
                            failed = 1;
                        }
                    }
                    all_e[gs + j] = e;
                }
            }
        }
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
    sample_rss(info);
    if (failed) {
        free(all_h);
        free(all_e);
        free(gstart);
        rs_fof_free_groups(&groups);
        rs_fof_free_cells(&cells);
        return finish(info, catalog_dir, step, box_index, RS_FAILED, fail_msg, t_start);
    }
    info->nhalos_internal = total;

    /* ---- selection, records ---- */
    int64_t *row = malloc((size_t)(total > 0 ? total : 1) * sizeof(int64_t));
    if (!row) {
        free(all_h); free(all_e); free(gstart);
        rs_fof_free_groups(&groups);
        rs_fof_free_cells(&cells);
        return finish(info, catalog_dir, step, box_index, RS_FAILED, "records: out of memory", t_start);
    }
    int64_t nsel = 0;
    for (int64_t i = 0; i < total; i++) {
        wrap_into_box(all_h[i].pos);
        row[i] = should_print(all_h + i) ? nsel++ : -1;
    }
    struct rs_record *rec = calloc((size_t)(nsel > 0 ? nsel : 1), sizeof(struct rs_record));
    if (!rec) {
        free(row); free(all_h); free(all_e); free(gstart);
        rs_fof_free_groups(&groups);
        rs_fof_free_cells(&cells);
        return finish(info, catalog_dir, step, box_index, RS_FAILED, "records: out of memory", t_start);
    }
#pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < total; i++) {
        if (row[i] < 0) continue;
        struct rs_record *r = rec + row[i];
        r->h = all_h[i];
        r->h.id = row[i];
        r->h.desc = -1;
        r->h.p_start = -1;
        r->ps_parent = all_e[i].sub_of >= 0 ? row[all_e[i].sub_of] : -1;
        r->fof_group = group_of_halo(gstart, ng, i);
        r->m_strict = r->r_strict = NAN;
        for (int k = 0; k < 4; k++) r->alt_m_strict[k] = NAN;
    }
    free(row);
    free(all_h);
    free(all_e);
    free(gstart);
    rs_fof_free_groups(&groups);
    info->nhalos = nsel;

    /* ---- strict SO masses ---- */
    t0 = omp_get_wtime();
    int strict_rc = 0;
    if (cfg_strict) strict_rc = rs_strict_pass(nsel, rec, x, y, z, &u, &cells);
    rs_fof_free_cells(&cells);
    info->t_strict = omp_get_wtime() - t0;
    if (strict_rc != 0) {
        free(rec);
        return finish(info, catalog_dir, step, box_index, RS_FAILED, "strict SO pass: out of memory", t_start);
    }

    /* ---- parents ---- */
    t0 = omp_get_wtime();
    if (rs_find_parents(nsel, rec, cfg_strict) != 0) {
        free(rec);
        return finish(info, catalog_dir, step, box_index, RS_FAILED, "parents: out of memory", t_start);
    }
    for (int64_t i = 0; i < nsel; i++) if (rec[i].pid < 0) info->nhosts++;
    info->t_parents = omp_get_wtime() - t0;

    /* ---- publication ---- */
    t0 = omp_get_wtime();
    struct mbuf b = {0};
    char cols[16384], err[512];
    meta_common(&b, "ok", "", step, box_index);
    madd(&b, "rows = %lld\n", (long long)nsel);
    madd(&b, "hosts = %lld\n", (long long)info->nhosts);
    madd(&b, "a = %.8f\nz = %.8f\n", aexpn, 1.0 / aexpn - 1.0);
    madd(&b, "a_v = %.8f\nvelocity_epoch = a_v = AEXPN - ASTEP/2 (GLAM's momenta lag the positions by half a "
             "step); v [km/s] = V * 100 * (Box/Ngrid) / a_v; ASTEP = %.8f\n", u.a_vel, astep);
    madd(&b, "positions = x [Mpc/h comoving] = (x_grid - 1) * Box/Ngrid, wrapped into [0, Box); particle id = "
             "GLAM array index - 1\n");
    madd(&b, "Box = %.9g\nNgrid = %lld\nNparticles = %lld\n", box, (long long)ngrid, (long long)n);
    madd(&b, "Om = %.9g\nOl = %.9g\nh0 = %.9g\n", Om, Ol, h0);
    madd(&b, "PARTICLE_MASS = %.9g\nAVG_PARTICLE_SPACING = %.9g\n", PARTICLE_MASS, AVG_PARTICLE_SPACING);
    madd(&b, "FOF_LINKING_LENGTH = %.9g\nfof_r_link = %.9g\nfof = exact periodic friends-of-friends of the "
             "adapter (rs_fof.c), d^2 < r^2 in float; Rockstar's own FoF additionally merges groups within 2r of "
             "particles with > 16 neighbours\n", FOF_LINKING_LENGTH,
         (double)rs_fof_linking_length(AVG_PARTICLE_SPACING, FOF_LINKING_LENGTH));
    madd(&b, "FORCE_RES = %.9g\nFORCE_RES_configured = %.9g\nFORCE_RES_PHYS_MAX = %.9g\n", FORCE_RES,
         cfg_force_res, FORCE_RES_PHYS_MAX);
    madd(&b, "MASS_DEFINITION = %s\nMASS_DEFINITION2 = %s\nMASS_DEFINITION3 = %s\nMASS_DEFINITION4 = %s\n"
             "MASS_DEFINITION5 = %s\n", MASS_DEFINITION, MASS_DEFINITION2, MASS_DEFINITION3, MASS_DEFINITION4,
         MASS_DEFINITION5);
    madd(&b, "MIN_HALO_PARTICLES = %lld\nMIN_HALO_OUTPUT_SIZE = %lld\nFOF_FRACTION = %.9g\nBOUND_PROPS = %lld\n"
             "UNBOUND_THRESHOLD = %.9g\nINITIAL_METRIC_SCALING = %.9g\nINCLUDE_HOST_POTENTIAL_RATIO = %.9g\n"
             "EXACT_LL_CALC = %lld\n", (long long)MIN_HALO_PARTICLES, (long long)MIN_HALO_OUTPUT_SIZE,
         FOF_FRACTION, (long long)BOUND_PROPS, UNBOUND_THRESHOLD, INITIAL_METRIC_SCALING,
         INCLUDE_HOST_POTENTIAL_RATIO, (long long)EXACT_LL_CALC);
    madd(&b, "strict_so = %d\npid_radius = %s\n", cfg_strict, cfg_strict ? "r_strict" : "r");
    madd(&b, "member_order = particle id (each FoF group is handed to find_subs sorted by particle id)\n");
    madd(&b, "config_file = %s\n", config_path);
    madd(&b, "fof_groups = %lld\nfof_members = %lld\nfof_fraction = %.6f\nfof_largest = %lld\n", (long long)ng,
         (long long)nm, n > 0 ? (double)nm / (double)n : 0.0, (long long)info->largest);
    madd(&b, "halos_internal = %lld\n", (long long)total);
    madd(&b, "threads = %d\n", T);
    madd(&b, "time_fof_s = %.2f\ntime_gather_s = %.2f\ntime_kernel_s = %.2f\ntime_strict_s = %.2f\n"
             "time_parents_s = %.2f\n", info->t_fof, info->t_gather, info->t_kernel, info->t_strict,
         info->t_parents);
    madd(&b, "vmrss_entry_gib = %.3f\nvmrss_max_sampled_gib = %.3f\nneed_fof_gib = %.3f\nneed_kernel_gib = %.3f\n"
             "memory_ceiling_gib = %.1f\n", info->rss_entry_gib, info->rss_max_gib, info->need_fof_gib,
         info->need_kernel_gib, cfg_ceiling_gib);
    if (rs_columns_text(cols, sizeof cols) > 0) madd(&b, "%s", cols);
    const int have_meta = b.s != NULL;
    int wrc = have_meta ? rs_write_catalogue(catalog_dir, step, box_index, rec, nsel, b.s, err, sizeof err) : -1;
    free(b.s);
    free(rec);
    info->t_write = omp_get_wtime() - t0;
    if (wrc != 0)
        return finish(info, catalog_dir, step, box_index, RS_FAILED, have_meta ? err : "metadata: out of memory",
                      t_start);
    malloc_trim(0);
    info->t_total = omp_get_wtime() - t_start;
    info->status = RS_OK;
    fprintf(stdout, " Rockstar step %d: %lld haloes (%lld hosts) from %lld FoF groups (f = %.4f, largest %lld);"
                    " FoF %.1f s, gather %.1f s, kernel %.1f s, strict %.1f s, parents %.1f s, write %.1f s,"
                    " total %.1f s; VmRSS %.2f -> max %.2f GiB\n",
            step, (long long)nsel, (long long)info->nhosts, (long long)ng, n > 0 ? (double)nm / (double)n : 0.0,
            (long long)info->largest, info->t_fof, info->t_gather, info->t_kernel, info->t_strict, info->t_parents,
            info->t_write, info->t_total, info->rss_entry_gib, info->rss_max_gib);
    fflush(stdout);
    return RS_OK;
}
