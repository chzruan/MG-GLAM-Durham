/* rs_inline.h: the interface of the in-line Rockstar adapter, called by MG-GLAM
 * at primary outputs (PMP2mod_rockstar.f90) and by the test harness.
 *
 * Part of the MG-GLAM in-line Rockstar adapter. The adapter links code from
 * MPI-Rockstar (https://github.com/Tomoaki-Ishiyama/mpi-rockstar) and Rockstar
 * (Behroozi et al. 2013), both GPLv3, and is distributed under the GNU General
 * Public License, version 3 or later. A GLAM build without the adapter links
 * rs_stub.c instead, which contains no Rockstar code.
 *
 * Nothing in here ends the calling process: every Rockstar error (exit() inside
 * Rockstar throws, assertions throw), allocation failure or failed memory
 * pre-check comes back as a status, and a catalogue stating the status is
 * published so the run's post-processing knows what happened. GLAM's particle
 * arrays are only read. */
#ifndef RS_INLINE_H
#define RS_INLINE_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

enum rs_status { RS_OK = 0, RS_SKIPPED = 1, RS_FAILED = 2, RS_NOT_BUILT = 3 };

struct rs_inline_info {
    int     status;                     /* enum rs_status */
    char    reason[512];                /* why, when not RS_OK */
    int64_t ngroups, nmembers, largest; /* FoF groups (>= MIN_HALO_PARTICLES), their particles, the largest */
    int64_t nhalos_internal;            /* halos the kernel found */
    int64_t nhalos;                     /* rows written (MIN_HALO_OUTPUT_SIZE, UNBOUND_THRESHOLD) */
    int64_t nhosts;                     /* rows with pid = -1 */
    double  t_fof, t_gather, t_kernel, t_strict, t_parents, t_write, t_total;   /* seconds */
    double  rss_entry_gib, rss_max_gib; /* VmRSS at entry and the largest sampled during the call */
    double  need_fof_gib, need_kernel_gib;  /* the pre-check estimates */
};

/* 1 when the adapter is built in, 0 for the stub */
int rs_inline_available(void);

/* "MPI-Rockstar <commit>, flags ..." of this build */
const char *rs_inline_version(void);

/* Parses the configuration file (Rockstar keys + GLAM_* keys) without running
 * anything; 0 when usable, else -1 and a message. For the start-up check. */
int rs_inline_check_config(const char *config_path, char *msg, int msglen);

/* Finds the haloes of one output and publishes
 * <catalog_dir>/Rockstar.<step:04d>.<box_index:04d>.bin + .meta.
 *   x..vz         GLAM's particle arrays (grid units, code velocities), n entries
 *   box           Mpc/h;  ngrid  GLAM's mesh;  aexpn, astep  the step's a and da
 *   om, oml, hubble  cosmology
 * Returns the status (also in *info when info is not NULL). */
int rs_inline_run(int64_t n, const float *x, const float *y, const float *z,
                  const float *vx, const float *vy, const float *vz,
                  double box, int64_t ngrid, double aexpn, double astep,
                  double om, double oml, double hubble, int step, int box_index,
                  const char *config_path, const char *catalog_dir, struct rs_inline_info *info);

#ifdef __cplusplus
}
#endif

#endif
