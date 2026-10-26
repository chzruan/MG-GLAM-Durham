/* rs_inline_internal.h: declarations shared by the files of the in-line Rockstar
 * adapter (not part of its interface, rs_inline.h).
 * Part of the MG-GLAM in-line Rockstar adapter; GNU GPL v3 or later. */
#ifndef RS_INLINE_INTERNAL_H
#define RS_INLINE_INTERNAL_H

#include <stddef.h>
#include <stdint.h>

#include "rs_convert.h"
#include "rs_fof.h"
#include "rs_record.h"

#ifdef __cplusplus
extern "C" {
#endif

/* rs_write.c */
int rs_columns_text(char *buf, size_t cap);
int rs_write_catalogue(const char *dir, int step, int box, const struct rs_record *rec, int64_t nrows,
                       const char *meta, char *err, int errlen);

/* rs_parents.c: pid / upid of n records (positions wrapped into the box) */
int rs_find_parents(int64_t n, struct rs_record *rec, int use_strict_radius);

/* rs_strict.c: strict SO masses of n records from every particle (GLAM arrays via the cell index);
 * needs calc_mass_definition() on the calling thread (its thresholds are copied into the team) */
int rs_strict_pass(int64_t n, struct rs_record *rec, const float *x, const float *y, const float *z,
                   const struct rs_units *u, const struct rs_cellindex *cells);   /* 0, or -1: out of memory */

/* rs_wrap_potential.c, rs_wrap_subhalo_metric.c: per-thread scratch MPI-Rockstar never frees */
void rs_free_potential_scratch(void);
void rs_free_subhalo_scratch(void);

#ifdef __cplusplus
}
#endif

#endif
