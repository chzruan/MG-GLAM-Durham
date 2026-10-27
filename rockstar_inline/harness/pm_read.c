/* pm_read.c: GLAM PM snapshot reader for the Rockstar test harness; see pm_read.h.
 * Part of the MG-GLAM in-line Rockstar adapter; GNU GPL v3 or later. */
#define _FILE_OFFSET_BITS 64
#define _POSIX_C_SOURCE 200809L
#include "pm_read.h"

#include <sys/types.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PM_RECORD     (1L << 20)          /* particles per record */
#define PM_RECS_FILE  512                 /* records per PMcrs file */
#define PM_HEADER_LEN 529

static uint32_t be32(const unsigned char *b)
{
    return ((uint32_t)b[0] << 24) | ((uint32_t)b[1] << 16) | ((uint32_t)b[2] << 8) | b[3];
}

static float bef(const unsigned char *b)
{
    uint32_t u = be32(b);
    float f;
    memcpy(&f, &u, 4);
    return f;
}

static int32_t bei(const unsigned char *b) { return (int32_t)be32(b); }

static int64_t bei64(const unsigned char *b)
{
    return (int64_t)(((uint64_t)be32(b) << 32) | be32(b + 4));
}

static void name_crd(char *buf, size_t n, const char *dir, int step)
{
    if (step < 0) snprintf(buf, n, "%s/PMcrd.DAT", dir);
    else          snprintf(buf, n, "%s/PMcrd.%04d.DAT", dir, step);
}

static void name_crs(char *buf, size_t n, const char *dir, int step, int ifile)
{
    if (step < 0) snprintf(buf, n, "%s/PMcrs%d.DAT", dir, ifile);
    else          snprintf(buf, n, "%s/PMcrs%d.%04d.DAT", dir, ifile, step);
}

int pm_read_header(const char *dir, int step, struct pm_header *h)
{
    char name[4096];
    unsigned char rec[4 + PM_HEADER_LEN + 4];
    name_crd(name, sizeof name, dir, step);
    FILE *f = fopen(name, "rb");
    if (!f) { fprintf(stderr, "pm_read: cannot open %s: %s\n", name, strerror(errno)); return -1; }
    size_t got = fread(rec, 1, sizeof rec, f);
    fclose(f);
    if (got != sizeof rec || be32(rec) != PM_HEADER_LEN || be32(rec + 4 + PM_HEADER_LEN) != PM_HEADER_LEN) {
        fprintf(stderr, "pm_read: %s: not a big-endian PM header record of %d bytes\n", name, PM_HEADER_LEN);
        return -1;
    }
    const unsigned char *b = rec + 4;
    memset(h, 0, sizeof *h);
    memcpy(h->title, b, 45);
    h->aexpn = bef(b + 45);  h->aexp0 = bef(b + 49);  h->amplt = bef(b + 53);  h->astep = bef(b + 57);
    h->istep = bei(b + 61);
    h->partw = bef(b + 65);  h->tintg = bef(b + 69);  h->ekin = bef(b + 73);   h->ekin1 = bef(b + 77);
    h->ekin2 = bef(b + 81);  h->au0 = bef(b + 85);    h->aeu0 = bef(b + 89);
    h->nrow = bei(b + 93);   h->ngrid = bei(b + 97);  h->nrealization = bei(b + 101); h->nseed = bei(b + 105);
    h->om = bef(b + 109);    h->oml = bef(b + 113);   h->hubble = bef(b + 117);
    h->nparticles = bei64(b + 121);
    for (int i = 0; i < 100; i++) h->extras[i] = bef(b + 129 + 4 * i);
    h->box = h->extras[99];
    if (h->nparticles <= 0 || h->ngrid <= 0 || !(h->box > 0)) {
        fprintf(stderr, "pm_read: %s: implausible header (N=%lld Ngrid=%d Box=%g)\n", name,
                (long long)h->nparticles, h->ngrid, h->box);
        return -1;
    }
    return 0;
}

int64_t pm_read_particles(const char *dir, int step, const struct pm_header *h,
                          float *x, float *y, float *z, float *vx, float *vy, float *vz)
{
    const int64_t n = h->nparticles;
    const int64_t npages = (n - 1) / PM_RECORD + 1;
    const size_t reclen = 6 * PM_RECORD * sizeof(float);
    unsigned char *rec = malloc(reclen);
    if (!rec) { fprintf(stderr, "pm_read: out of memory\n"); return -1; }
    float *dst[6] = {x, y, z, vx, vy, vz};
    const float lo = 1.0f, hi = (float)h->ngrid + 1.0f;
    int64_t outside = 0;
    FILE *f = NULL;
    int ifile = -1;
    char name[4096];
    for (int64_t page = 0; page < npages; page++) {
        int want = (int)(page / PM_RECS_FILE);
        if (want != ifile) {
            if (f) fclose(f);
            ifile = want;
            name_crs(name, sizeof name, dir, step, ifile);
            f = fopen(name, "rb");
            if (!f) { fprintf(stderr, "pm_read: cannot open %s: %s\n", name, strerror(errno)); free(rec); return -1; }
        }
        if (fseeko(f, (off_t)(page % PM_RECS_FILE) * (off_t)reclen, SEEK_SET) != 0 ||
            fread(rec, 1, reclen, f) != reclen) {
            fprintf(stderr, "pm_read: %s: short read at record %lld\n", name, (long long)(page % PM_RECS_FILE));
            fclose(f); free(rec); return -1;
        }
        const int64_t off = page * PM_RECORD;
        const int64_t nin = (page == npages - 1) ? n - off : PM_RECORD;
        for (int c = 0; c < 6; c++) {
            const unsigned char *b = rec + (size_t)c * PM_RECORD * sizeof(float);
            float *d = dst[c] + off;
            for (int64_t i = 0; i < nin; i++) d[i] = bef(b + 4 * i);
        }
        for (int c = 0; c < 3; c++)
            for (int64_t i = 0; i < nin; i++) {
                float v = dst[c][off + i];
                if (!(v >= lo && v < hi)) outside++;
            }
    }
    if (f) fclose(f);
    free(rec);
    return outside;
}
