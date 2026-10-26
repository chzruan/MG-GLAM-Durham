/* fof_ref.h: serial reference periodic FoF (see fof_ref.c).
 * Part of the MG-GLAM in-line Rockstar adapter; GNU GPL v3 or later. */
#ifndef FOF_REF_H
#define FOF_REF_H

#include <stdint.h>

struct fof_stats {
    int64_t ngroups;     /* groups with >= nmin members */
    int64_t nmembers;    /* particles in those groups */
    int64_t largest;     /* members of the largest group */
    double  fraction;    /* nmembers / n */
};

/* positions in [0, box) as float; r the float linking length; label (optional,
 * n entries) gets the smallest member index of a kept group, or -1.
 * Returns the number of kept groups, or -1 on error. */
int64_t fof_reference(int64_t n, const float *x, const float *y, const float *z,
                      float box, float r, int64_t nmin, int64_t *label,
                      struct fof_stats *st);

#endif
