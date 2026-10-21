/* glam_kmp_alloc.c: route Fortran ALLOCATE/DEALLOCATE to glibc malloc.
 *
 * With -qopenmp the Intel Fortran runtime (libifcoremt) resolves its weak
 * references kmp_aligned_malloc / kmp_realloc / kmp_free against libiomp5,
 * which serves them from an embedded scalable (TBB-style) allocator whose
 * per-thread and large-object caches keep freed blocks and never give them
 * back while the program does not allocate. Measured in the memory audit of
 * 2026-10-21 (docs/memory-audit-20261021.md, cut 1): every halo-finder call
 * left its freed temporaries resident through all later force steps
 * (41-47 MB per thread, ~7 GiB at 128 threads); glibc's malloc_trim cannot
 * touch that memory. Defining the three symbols here makes the executable's
 * definitions win over libiomp5's (standard ELF interposition), so every
 * Fortran allocation goes to glibc: large arrays are mmap'd and unmapped on
 * deallocation, and TrimHeap (malloc_trim) releases the rest. The OpenMP
 * runtime itself does not use this API internally. No numerical effect.
 */
#include <stdlib.h>
#include <string.h>

void *kmp_aligned_malloc(size_t size, size_t alignment)
{
    void *p = NULL;
    if (alignment < sizeof(void *) || (alignment & (alignment - 1)) != 0)
        alignment = sizeof(void *) > alignment ? sizeof(void *) : 64;
    if (posix_memalign(&p, alignment, size ? size : 1) != 0)
        return NULL;
    return p;
}

void *kmp_malloc(size_t size) { return malloc(size ? size : 1); }

void *kmp_calloc(size_t n, size_t size) { return calloc(n ? n : 1, size ? size : 1); }

void *kmp_realloc(void *p, size_t size) { return realloc(p, size ? size : 1); }

void kmp_free(void *p) { free(p); }
