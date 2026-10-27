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

/* rs_write.c: the column list of the .meta (descriptions follow BOUND_PROPS); -1 if cap is too small */
int rs_columns_text(char *buf, size_t cap, int bound_props);
int rs_write_catalogue(const char *dir, int step, int box, const struct rs_record *rec, int64_t nrows,
                       const char *meta, char *err, int errlen);
/* 1 when <dir>/Rockstar.<step>.<box>.meta exists, says status = ok and was written for this output: the same
 * a and N_p (as the .meta prints them) and, unless config_hash is NULL, the same config_fnv1a64 */
int rs_catalogue_ok(const char *dir, int step, int box, const char *a, const char *n_p, const char *config_hash);
/* removes the fixed-name temporaries .Rockstar.<step>.<box>.{bin,meta}.tmp an interrupted write left */
void rs_remove_temporaries(const char *dir, int step, int box);

/* rs_parents.c: the find_parents parent (a row, or -1: distinct) of n halos (positions wrapped into the
 * box) with the kernel radius r; 0, or -1: out of memory. A Rockstar error inside throws: call it
 * through rs_guarded_find_parents (rs_guard.cpp), which releases the call's memory (rs_parents_cleanup) */
int rs_find_parents(int64_t n, const struct halo *h, int64_t *pid);
void rs_parents_cleanup(void);

/* rs_guard.cpp */
int rs_guarded_call(int (*fn)(void *), void *arg, const char *what, char *msg, int msglen);
int rs_guarded_find_parents(int64_t n, const struct halo *h, int64_t *pid, char *msg, int msglen);

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
