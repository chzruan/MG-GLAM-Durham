/* pm_read.h: read a GLAM PM snapshot (PMcrd.<step>.DAT + PMcrs<i>.<step>.DAT)
 * into the six SoA arrays GLAM holds in memory (grid units, code velocities).
 *
 * Part of the MG-GLAM in-line Rockstar adapter (test harness); GNU GPL v3 or
 * later, see rs_convert.h.
 *
 * Format (PMP2mod_tools.f90, ReadDataPM / WriteDataPM; GLAM is compiled with
 * -convert big_endian):
 *   PMcrd: one Fortran sequential unformatted record (4-byte big-endian length
 *          markers) of 529 bytes: HEADER char(45); AEXPN, AEXP0, AMPLT, ASTEP
 *          real*4; ISTEP int*4; PARTW, TINTG, EKIN, EKIN1, EKIN2, AU0, AEU0
 *          real*4; NROW, NGRID, Nrealization, Nseed int*4; Om, OmL, hubble
 *          real*4; Nparticles int*8; extras(100) real*4 with extras(100) = Box.
 *   PMcrs<i>: direct access, records of 6 x 2^20 big-endian real*4 (X, Y, Z,
 *          VX, VY, VZ of 2^20 particles, no markers), 512 records per file;
 *          the last record is written full length (its tail is stale).
 * step < 0 reads the checkpoint names PMcrd.DAT / PMcrs<i>.DAT. */
#ifndef PM_READ_H
#define PM_READ_H

#include <stdint.h>

struct pm_header {
    char    title[46];
    float   aexpn, aexp0, amplt, astep;
    int32_t istep;
    float   partw, tintg, ekin, ekin1, ekin2, au0, aeu0;
    int32_t nrow, ngrid, nrealization, nseed;
    float   om, oml, hubble;
    int64_t nparticles;
    float   extras[100];
    double  box;            /* extras[99] */
};

/* 0 on success; on failure a message on stderr and -1 */
int pm_read_header(const char *dir, int step, struct pm_header *h);

/* fill x..vz (each h->nparticles floats, allocated by the caller); returns the
 * number of positions outside [1, Ngrid+1) (left unchanged), or -1 on error */
int64_t pm_read_particles(const char *dir, int step, const struct pm_header *h,
                          float *x, float *y, float *z,
                          float *vx, float *vy, float *vz);

#endif
