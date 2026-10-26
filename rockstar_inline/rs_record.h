/* rs_record.h: one row of the in-line Rockstar catalogue
 * CATALOGS/Rockstar.<step>.<box>.bin (raw little-endian records, no header; the
 * .meta sidecar lists every column with its type, byte offset and unit).
 *
 * Part of the MG-GLAM in-line Rockstar adapter; GNU GPL v3 or later.
 *
 * The record starts with Rockstar's own struct halo (halo.h of MPI-Rockstar
 * 0ae386f, built with -DOUTPUT_RVMAX -DOUTPUT_NFW_CHI2 -DOUTPUT_INERTIA_TENSOR:
 * 312 bytes, HALO_FORMAT_REVISION 2), as Rockstar's halos_*.bin would hold it,
 * with id = the row number, desc = -1 and p_start = -1 (no member particle list
 * is written). The adapter's own columns follow. */
#ifndef RS_RECORD_H
#define RS_RECORD_H

#include <stdint.h>

#include "halo.h"

struct rs_record {
    struct halo h;
    int64_t pid;              /* parent as Rockstar's find_parents (parents.c): -1 for hosts */
    int64_t upid;             /* the root of the pid chain: -1 for hosts */
    int64_t ps_parent;        /* Rockstar's phase-space parent (extra_halo_info.sub_of) as a row, -1 if none */
    int64_t fof_group;        /* canonical FoF group (b = FOF_LINKING_LENGTH) the halo was found in */
    float   m_strict;         /* strict SO mass of MASS_DEFINITION, all particles (Msun/h) */
    float   r_strict;         /* its radius (kpc/h comoving) */
    float   alt_m_strict[4];  /* strict SO masses of MASS_DEFINITION2..5 (Msun/h) */
    int32_t strict_crossed;   /* bit 0: MASS_DEFINITION reached its threshold, bits 1-4: the others */
    int32_t pad;
};

#if defined(__cplusplus)
static_assert(sizeof(struct halo) == 312, "struct halo must be the 312-byte OUTPUT_* layout");
static_assert(sizeof(struct rs_record) == 376, "rs_record layout");
#else
_Static_assert(sizeof(struct halo) == 312, "struct halo must be the 312-byte OUTPUT_* layout");
_Static_assert(sizeof(struct rs_record) == 376, "rs_record layout");
#endif

#endif
