/* rs_stub.c: the in-line Rockstar interface (rs_inline.h) of a GLAM build
 * without Rockstar. It contains no Rockstar code; GLAM refuses a run that asks
 * for Rockstar haloes with it (CheckRockstarConfig), so it is never reached at
 * an output. Build the adapter with `make ROCKSTAR_DIR=<MPI-Rockstar source>`. */
#include <stdio.h>
#include <string.h>

#include "rs_inline.h"

int rs_inline_available(void) { return 0; }

const char *rs_inline_version(void) { return "not built"; }

int rs_inline_check_config(const char *config_path, char *msg, int msglen)
{
    (void)config_path;
    snprintf(msg, msglen, "this GLAM build has no in-line Rockstar: rebuild with make ROCKSTAR_DIR=<MPI-Rockstar source>");
    return -1;
}

int rs_inline_run(int64_t n, const float *x, const float *y, const float *z,
                  const float *vx, const float *vy, const float *vz,
                  double box, int64_t ngrid, double aexpn, double astep,
                  double om, double oml, double hubble, int mg_flag, int mg_model, double fr0,
                  const char *header, int step, int box_index,
                  const char *config_path, const char *catalog_dir, struct rs_inline_info *info)
{
    (void)n; (void)x; (void)y; (void)z; (void)vx; (void)vy; (void)vz; (void)box; (void)ngrid; (void)aexpn;
    (void)astep; (void)om; (void)oml; (void)hubble; (void)mg_flag; (void)mg_model; (void)fr0; (void)header;
    (void)step; (void)box_index; (void)config_path; (void)catalog_dir;
    if (info) {
        memset(info, 0, sizeof *info);
        info->status = RS_NOT_BUILT;
        snprintf(info->reason, sizeof info->reason, "this GLAM build has no in-line Rockstar");
    }
    return RS_NOT_BUILT;
}
