# Independent review of the first halo-property repair

Reviewed integrated commit `be03816a13170f5787897ea8b2b5e6725b092122`, including
halo-property repair `8b85c9fca71e3e7e8a5cbc54a5e1d1b264c6107b`. These are
findings against that reviewed revision; subsequent repairs must be recorded
separately. Root checkout was read only throughout this review.

## Findings

1. **P2: the documented outermost SO crossing depends on the initial Cell
   bracket.** `GetHalo:1044–1068` stops gathering when a sampled outer radius
   first becomes underdense, then picks the outermost root only inside that
   bracket. A nonmonotonic enclosed-density profile can rise again beyond it.
   For ten cold particles at radius 0.1 and 200 at radius 1.2 Mpc/h, equal mass
   1e12 Msun/h, Om0=0.3, and overdensity 200, the exact downward roots are
   0.525271249 and 1.449183548 Mpc/h. Keeping the particles and all physical
   parameters fixed, Cell=1 returns the smaller root and ten bound particles;
   Cell=2 returns the larger root and 210 bound particles. Both checked and
   optimized binaries return status zero. Define the outermost root over the
   verified physical search domain, rather than a Cell-dependent subinterval.
   This is a contract failure on a controlled profile, not an estimate of its
   incidence in cosmological catalogues.

2. **P2 at large supported populations: both heapsorts can overflow their
   default-integer child indices.** `BdmHaloSortRadii:1282` and
   `BdmHaloSortIds:1312` execute `child=2*child` after reaching a leaf. With
   n=1,200,000,000 and parent=1,100,000,000, the next child wraps to
   -2,094,967,296 and still satisfies `child<=last`. Subsequent array access
   can therefore use a negative subscript. These counts are below 1200^3.
   The bounded reproduction executes the index arithmetic, not a billion-row
   allocation; the complete large-array failure was not run. Use int64 index
   arithmetic or guard a leaf before doubling its index.

## Checks without a finding

Twenty-eight bounded probes use checked/FP-trapping and optimized GNU builds.
They exercise zero/one/nine/ten initial particles, a single exactly central
survivor after hot contaminants leave, zero/one/coincident-centre spherical
potential cases, repeated shell radii and distinct-pair energy, and a second
empty call that must clear previous output and membership. Artificial IDs
larger than 2^31 also exercised storage width at this reviewed revision; these
were not valid production original-row mappings. The follow-up regression
uses actual row identities and tests large integer values directly in the
sort helper. All these edge checks passed. Initial populations below ten
return insufficient-population status; a singleton central survivor is finite
and has the unresolved-Vmax status.

Static inspection also confirms that survivor compaction writes only to rows
already visited, that a successful final pass uses the same membership and
bulk velocity for its binding and property calculation, and that the spherical
potential excludes each particle's self term and counts each pair once.
This does not certify the spherical approximation as exact aspherical binding.

Reproduction:

```sh
micromamba run -n cosemu python3 -B BDM-refine/repairs/20260906/tests/halo_review.py
```

[halo_review_results.json](halo_review_results.json) records the original
review's source/test hashes and complete evidence. The original review driver
records findings as data, rather than presenting these findings as passing
physics assertions.
