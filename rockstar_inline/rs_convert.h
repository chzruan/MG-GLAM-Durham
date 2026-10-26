/* rs_convert.h: GLAM particle units -> Rockstar particle units, the single place
 * where this conversion is done (the in-line finder and the test harness both
 * include it).
 *
 * Part of the MG-GLAM in-line Rockstar adapter. The adapter links code from
 * MPI-Rockstar (https://github.com/Tomoaki-Ishiyama/mpi-rockstar) and Rockstar
 * (Behroozi et al. 2013), both GPLv3, and is distributed under the GNU General
 * Public License, version 3 or later (https://www.gnu.org/licenses/gpl-3.0).
 *
 * GLAM in memory: x in grid units, [1, Ngrid+1); V the code momentum, half a
 * leapfrog step behind the positions (PMP2mod_subsample.f90 header).
 * Rockstar: comoving position in Mpc/h, [0, Box); physical peculiar velocity
 * in km/s; int64 particle id.
 *   x [Mpc/h] = (x_grid - 1) * Box / Ngrid, wrapped into [0, Box)
 *   v [km/s]  = V * 100 * (Box / Ngrid) / a_v,   a_v = AEXPN - ASTEP/2
 *   id        = array index (0-based) = GLAM index - 1 (= Gadget ID for ic2pm ICs)
 * The arithmetic is done in double and rounded once to float, which is what
 * Rockstar's trees and our FoF both use. A coordinate that rounds to Box is
 * set to 0. */
#ifndef RS_CONVERT_H
#define RS_CONVERT_H

#include <stdint.h>

struct rs_units {
    double box;        /* Mpc/h */
    double xscale;     /* Box / Ngrid: Mpc/h per grid unit */
    double vscale;     /* 100 * Box / Ngrid / a_v: km/s per code velocity unit */
    double a_pos;      /* AEXPN, the epoch of the positions (SCALE_NOW) */
    double a_vel;      /* AEXPN - ASTEP/2, the epoch of the velocities */
    float  box_f;      /* (float)box */
};

static inline void rs_units_init(struct rs_units *u, double box, int64_t ngrid,
                                 double aexpn, double astep)
{
    u->box    = box;
    u->xscale = box / (double)ngrid;
    u->a_pos  = aexpn;
    u->a_vel  = aexpn - 0.5 * astep;
    u->vscale = 100.0 * u->xscale / u->a_vel;
    u->box_f  = (float)box;
}

/* one coordinate: grid units -> [0, Box) Mpc/h as float */
static inline float rs_pos(const struct rs_units *u, float xg)
{
    double x = ((double)xg - 1.0) * u->xscale;
    if (x >= u->box) x -= u->box;      /* x_grid is in [1, Ngrid+1): at most one wrap */
    else if (x < 0.0) x += u->box;
    float xf = (float)x;
    if (!(xf < u->box_f) || xf < 0.0f) xf = 0.0f;   /* rounded up to Box (or NaN) */
    return xf;
}

/* 1 when rs_pos maps xg faithfully: finite and at most one box outside [1, Ngrid+1) */
static inline int rs_pos_ok(const struct rs_units *u, float xg)
{
    double x = ((double)xg - 1.0) * u->xscale;
    return x >= -u->box && x < 2.0 * u->box;           /* false for NaN */
}

static inline float rs_vel(const struct rs_units *u, float vg)
{
    return (float)((double)vg * u->vscale);
}

/* fill pos[6] (x, y, z, vx, vy, vz) of Rockstar's struct particle */
static inline void rs_convert(const struct rs_units *u, float xg, float yg, float zg,
                              float vxg, float vyg, float vzg, float pos[6])
{
    pos[0] = rs_pos(u, xg);
    pos[1] = rs_pos(u, yg);
    pos[2] = rs_pos(u, zg);
    pos[3] = rs_vel(u, vxg);
    pos[4] = rs_vel(u, vyg);
    pos[5] = rs_vel(u, vzg);
}

#endif
