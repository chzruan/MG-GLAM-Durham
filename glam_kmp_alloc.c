/* glam_kmp_alloc.c: Fortran ALLOCATE/DEALLOCATE through glibc malloc, and the
 * heap trim the halo finder calls afterwards.
 *
 * Why. With -qopenmp the Intel Fortran runtime (libifcoremt.so.5, oneAPI
 * 2024.2) resolves its weak references kmp_aligned_malloc / kmp_realloc /
 * kmp_free against libiomp5.so, which serves them from an embedded scalable
 * (TBB-style) allocator. That allocator's large-object cache and per-thread
 * caches keep freed blocks until further allocation activity: after every
 * halo-finder call its freed temporaries stayed resident through all later
 * force steps (41-47 MB per thread, ~7 GiB at 128 threads; memory audit of
 * 2026-10-21, docs/memory-audit-20261021.md section 5, cut 1). glibc's
 * malloc_trim, the MALLOC_* tunables and -heap-arrays cannot touch that
 * memory, and the allocator's own clean-up entry point is not exported.
 *
 * How. These five definitions live in the executable, so by ELF symbol
 * interposition (the executable is first in the lookup scope; libiomp5 has no
 * DT_SYMBOLIC) libifcoremt's references bind here and every ALLOCATE goes to
 * glibc: large arrays are mmap'd and unmapped on DEALLOCATE, the rest is
 * released by glam_trim_heap. libiomp5 does not call these entry points
 * internally (its exported kmp_* are stubs onto kmpc_*, which stay inside
 * libiomp5), and the Fortran runtime's other allocation paths (I/O buffers,
 * copy-in/out temporaries) already use glibc malloc directly. What does NOT
 * come here: arrays under an explicit `!$omp allocate` directive or with
 * FASTMEM, which libifcoremt sends to __kmpc_aligned_alloc / __kmpc_free
 * (MG-GLAM uses neither). The runtime's documented FOR_DISABLE_KMP_MALLOC
 * switch is read but, per the 2024.2 disassembly, does not stop the kmp_*
 * path, so it is no substitute for this file.
 *
 * Checks. `make check-alloc-imports` lists the allocator entry points the
 * installed libifcoremt imports and fails if any name outside the known set
 * appears (a newer runtime could allocate through a name not defined here,
 * splitting allocate and free between two allocators: the first DEALLOCATE
 * would then die with "free(): invalid pointer"). Static Intel linking
 * (-static-intel / -qopenmp-link=static) fails loudly with duplicate
 * definitions; the makefile links the Intel runtime shared (-shared-intel).
 * No numerical effect: results were shown independent of allocation
 * placement (section 5 of the audit note).
 */
#include <stdlib.h>
#if defined(__INTEL_LLVM_COMPILER) || defined(__INTEL_COMPILER)
#include <omp.h>    /* Intel's omp.h declares the kmp_* prototypes: a signature mismatch fails to compile */
#endif
#if defined(__GLIBC__)
#include <malloc.h>
#endif

void *kmp_aligned_malloc(size_t size, size_t alignment)
{
    void *p = NULL;
    if (alignment < sizeof(void *) || (alignment & (alignment - 1)) != 0) {
        /* posix_memalign needs a power of two that is a multiple of sizeof(void*);
           the Fortran runtime asks for 32 or more, so this branch is defensive */
        size_t a = sizeof(void *);
        while (a < alignment) a <<= 1;
        alignment = a;
    }
    if (posix_memalign(&p, alignment, size ? size : 1) != 0)
        return NULL;    /* the Fortran runtime reports "insufficient virtual memory" / stat= 41 */
    return p;
}

void *kmp_malloc(size_t size) { return malloc(size ? size : 1); }

void *kmp_calloc(size_t n, size_t size) { return calloc(n ? n : 1, size ? size : 1); }

void *kmp_realloc(void *p, size_t size) { return realloc(p, size ? size : 1); }

void kmp_free(void *p) { free(p); }

/* Called by the halo finder (PMP2linker.f90 TrimHeap) after it has freed its
   work arrays: returns every free page of every glibc arena to the kernel.
   A no-op on a non-glibc libc, so the Fortran side stays portable. */
void glam_trim_heap(void)
{
#if defined(__GLIBC__)
    (void) malloc_trim(0);
#endif
}
