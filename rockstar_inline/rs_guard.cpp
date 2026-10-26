// rs_guard.cpp: the only C++ of the in-line Rockstar adapter. Every call into
// Rockstar that may end in exit() (which MPI-Rockstar's error.h turns into a
// thrown rockstar_error) or in a failed assertion (shadow/assert.h, also a
// throw) goes through a guard here, so no exception ever crosses an OpenMP
// region (that would call std::terminate) or reaches GLAM's Fortran.
//
// Part of the MG-GLAM in-line Rockstar adapter; GNU GPL v3 or later.
#include <cstdio>
#include <cstring>
#include <new>

#include "error.h"
extern "C" {
#include "fof.h"
#include "groupies.h"
#include "io/read_config.h"
#include "rs_inline.h"
}

extern "C" void rs_assert_fail(const char *expr, const char *file, int line);
extern "C" int rs_inline_run_impl(int64_t n, const float *x, const float *y, const float *z,
                                  const float *vx, const float *vy, const float *vz, double box,
                                  int64_t ngrid, double aexpn, double astep, double om, double oml,
                                  double hubble, int step, int box_index, const char *config_path,
                                  const char *catalog_dir, struct rs_inline_info *info);
extern "C" int rs_publish_failure(const char *catalog_dir, int step, int box_index, const char *reason);

static void describe(const rockstar_error &e, char *msg, int msglen)
{
    std::snprintf(msg, msglen, "Rockstar error %d at %s:%d", e.code, e.file ? e.file : "?", e.line);
}

// inject: 0 none; 1 a Rockstar exit(); 2 a failed assertion (test hooks of rs_inline.c)
extern "C" int rs_guarded_find_subs(struct fof *f, struct FOFInfo *fofinfo, struct HaloInfo *haloinfo,
                                    int inject, char *msg, int msglen)
{
    try {
        if (inject == 1) rockstar_abort_impl(1, __FILE__, __LINE__);
        if (inject == 2) rs_assert_fail("injected test failure", __FILE__, __LINE__);
        find_subs(f, fofinfo, haloinfo);
        return 0;
    } catch (const rockstar_error &e) {
        describe(e, msg, msglen);
    } catch (const std::bad_alloc &) {
        std::snprintf(msg, msglen, "out of memory in the per-group kernel");
    } catch (...) {
        std::snprintf(msg, msglen, "unknown exception in the per-group kernel");
    }
    return -1;
}

extern "C" int rs_guarded_load_config(struct configfile *c, const char *path, char *msg, int msglen)
{
    try {
        load_config(c, const_cast<char *>(path));
        return 0;
    } catch (const rockstar_error &e) {
        std::snprintf(msg, msglen, "%s: not a usable Rockstar configuration (FILE_FORMAT must be set)", path);
    } catch (...) {
        std::snprintf(msg, msglen, "%s: cannot be parsed", path);
    }
    return -1;
}

extern "C" int rs_inline_run(int64_t n, const float *x, const float *y, const float *z, const float *vx,
                             const float *vy, const float *vz, double box, int64_t ngrid, double aexpn,
                             double astep, double om, double oml, double hubble, int step, int box_index,
                             const char *config_path, const char *catalog_dir, struct rs_inline_info *info)
{
    char reason[512] = "";
    try {
        return rs_inline_run_impl(n, x, y, z, vx, vy, vz, box, ngrid, aexpn, astep, om, oml, hubble, step,
                                  box_index, config_path, catalog_dir, info);
    } catch (const rockstar_error &e) {
        describe(e, reason, sizeof reason);
    } catch (const std::bad_alloc &) {
        std::snprintf(reason, sizeof reason, "out of memory");
    } catch (...) {
        std::snprintf(reason, sizeof reason, "unknown exception");
    }
    // an error outside the guarded kernel calls (serial Rockstar code): the run's
    // temporaries leak, but the simulation continues and the output says why
    std::fprintf(stdout, " Rockstar step %d: FAILED (%s)\n", step, reason);
    std::fflush(stdout);
    if (info) {
        info->status = RS_FAILED;
        std::snprintf(info->reason, sizeof info->reason, "%s", reason);
    }
    try {
        rs_publish_failure(catalog_dir, step, box_index, reason);
    } catch (...) {
    }
    return RS_FAILED;
}
