/* hash_ref.c: reference values of the glam_hash.c routines for the Python
 * check (test/extras/check_subsample.py, check a). Writes, little-endian:
 *   int64 n; uint64 ids[n]; uint64 splitmix64(ids ^ SALT)[n];
 *   uint8 keep_M100[n]; uint8 keep_M7[n]
 * IDs: 0 .. n/2-1 sequential (like real particle IDs), then xorshift64*
 * full-range random values, with the edge cases 0, 2^63-1, 2^63, 2^64-1
 * and SALT itself at the end. Also checks glam_count_block and
 * glam_select_indices against glam_select on the same IDs; exits 1 on
 * a mismatch.
 * Usage: hash_ref.exe <output file> [n]      (default n = 10^7) */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

uint64_t splitmix64(uint64_t x);
int glam_select(uint64_t id, uint64_t salt, uint64_t m);
int64_t glam_count_block(uint64_t id0, int64_t n, uint64_t salt, uint64_t m);
int64_t glam_select_indices(uint64_t id0, int64_t n, uint64_t salt, uint64_t m, int32_t *idx);

#define SALT 0x5EED00DE6C0A2026ULL

int main(int argc, char **argv)
{
    int64_t n = 10000000, i, c, k, j;
    uint64_t *ids, *h, x = 0x1234567887654321ULL;
    uint8_t *s100, *s7;
    int32_t *idx;
    FILE *f;
    if (argc < 2) { fprintf(stderr, "usage: hash_ref.exe <out> [n]\n"); return 2; }
    if (argc > 2) n = atoll(argv[2]);
    ids = malloc(n*sizeof(uint64_t)); h = malloc(n*sizeof(uint64_t));
    s100 = malloc(n); s7 = malloc(n);
    for (i = 0; i < n; i++) {
        if (i < n/2) ids[i] = (uint64_t)i;
        else { x ^= x >> 12; x ^= x << 25; x ^= x >> 27; ids[i] = x * 0x2545F4914F6CDD1DULL; }
    }
    if (n >= 5) {
        ids[n-5] = 0; ids[n-4] = 0x7FFFFFFFFFFFFFFFULL; ids[n-3] = 0x8000000000000000ULL;
        ids[n-2] = 0xFFFFFFFFFFFFFFFFULL; ids[n-1] = SALT;
    }
    for (i = 0; i < n; i++) {
        h[i] = splitmix64(ids[i] ^ SALT);
        s100[i] = (uint8_t)glam_select(ids[i], SALT, 100);
        s7[i]   = (uint8_t)glam_select(ids[i], SALT, 7);
    }
    /* block routines against the per-ID routine on the sequential range, in blocks of 1000 */
    idx = malloc(1000*sizeof(int32_t));
    for (i = 0; i < n/2; i += 1000) {
        int64_t nb = (n/2 - i < 1000) ? n/2 - i : 1000;
        c = 0;
        for (j = 0; j < nb; j++) c += s100[i+j];
        if (glam_count_block((uint64_t)i, nb, SALT, 100) != c) { fprintf(stderr, "count_block mismatch at %ld\n", (long)i); return 1; }
        k = glam_select_indices((uint64_t)i, nb, SALT, 100, idx);
        if (k != c) { fprintf(stderr, "select_indices count mismatch at %ld\n", (long)i); return 1; }
        for (j = 0, c = 0; j < nb; j++) if (s100[i+j]) { if (idx[c] != (int32_t)j) { fprintf(stderr, "select_indices offset mismatch at %ld\n", (long)(i+j)); return 1; } c++; }
    }
    f = fopen(argv[1], "wb");
    if (!f) { perror(argv[1]); return 2; }
    fwrite(&n, sizeof(int64_t), 1, f);
    fwrite(ids, sizeof(uint64_t), n, f);
    fwrite(h, sizeof(uint64_t), n, f);
    fwrite(s100, 1, n, f);
    fwrite(s7, 1, n, f);
    fclose(f);
    printf("hash_ref: %ld IDs written; block routines agree with glam_select\n", (long)n);
    return 0;
}
