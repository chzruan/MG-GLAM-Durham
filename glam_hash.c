/*-------------------------------------------------------------------
 *  glam_hash.c: ID hash for the DEGRACE particle subsample (item I9 of
 *  docs/degrace_companion_implementation_plan.md), called from the
 *  Fortran module Subsample (PMP2mod_subsample.f90) through
 *  iso_c_binding. Written in C because unsigned 64-bit wrap-around
 *  arithmetic is well defined here and not in Fortran.
 *
 *  Selection rule (never change; the Python reference in
 *  test/extras/check_subsample.py and every downstream reader
 *  depend on it):
 *      keep particle ID  <=>  splitmix64(ID ^ SALT) % M == 0
 *  with ID the 0-based lattice index, SALT = 0x5EED00DE6C0A2026 and
 *  M the subsampling factor (M = 1 keeps everything).
 *
 *  The block routines evaluate the rule for the contiguous ID range
 *  [id0, id0 + n) so that Fortran does not call into C per particle.
 *
 *  Build (see the makefile): icx -O2 -c glam_hash.c
 *-------------------------------------------------------------------*/
#include <stdint.h>

uint64_t splitmix64(uint64_t x)
{
    uint64_t z = x + 0x9E3779B97F4A7C15ULL;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

/* 1 if the particle with this ID is kept, 0 otherwise. */
int glam_select(uint64_t id, uint64_t salt, uint64_t m)
{
    return (splitmix64(id ^ salt) % m) == 0;
}

/* Number of kept IDs in [id0, id0 + n). */
int64_t glam_count_block(uint64_t id0, int64_t n, uint64_t salt, uint64_t m)
{
    int64_t  i, c = 0;
    if (m == 1) return n;               /* full mode: every ID is kept */
    for (i = 0; i < n; i++)
        c += (splitmix64((id0 + (uint64_t)i) ^ salt) % m) == 0;
    return c;
}

/* Offsets (0-based, relative to id0, increasing) of the kept IDs in
 * [id0, id0 + n), written to idx[0 .. count-1]; returns count.
 * idx must have room for n entries (n <= 2^31 - 1). */
int64_t glam_select_indices(uint64_t id0, int64_t n, uint64_t salt, uint64_t m, int32_t *idx)
{
    int64_t  i, c = 0;
    if (m == 1) {                       /* full mode: every ID is kept */
        for (i = 0; i < n; i++) idx[i] = (int32_t)i;
        return n;
    }
    for (i = 0; i < n; i++)
        if ((splitmix64((id0 + (uint64_t)i) ^ salt) % m) == 0)
            idx[c++] = (int32_t)i;
    return c;
}
