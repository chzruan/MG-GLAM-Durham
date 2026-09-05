# Follow-up to the two P2 review findings

Both findings were reproduced against the pre-fix Python tools on branch
`fix/bdm-duplicate-hosts-20260905` (baseline commit `133d220`). The repair remains
on this branch; `cz` and the Fortran finder are unchanged by this follow-up.

## Float32 HDF5 arithmetic

For stored float32 total masses near `1e12` and `1.0115794e12`, evaluating their
logarithms in float32 gives a difference of **0.0050001144**, incorrectly failing
the strict `<0.005` condition. The same stored values evaluated in float64 give
**0.00499998736**, correctly passing. The original cleaner consequently removed
one row from equivalent ASCII input and zero from HDF5.

`catalogue_core.py` now promotes arithmetic operands before periodic wrapping,
position/velocity differences, norms, logarithms, mass selections, mass-bin
counts and radial-velocity statistics. Exact mass/count and defining-field
equality still use the original arrays. Output writers continue to copy the
original values and preserve HDF5 storage dtypes; no rounding tolerance or
scientific threshold was changed.

Regression cases cover the mass-difference example, separations and relative
speeds just inside their thresholds, mass selection/bin edges and radial
velocity statistics. Equivalent float32 and float64 values give the same
diagnostics, and ASCII/HDF5 masks agree. Float32 data and integer ID/count
datasets remain unchanged in the retained HDF5 rows.

## Catalogue published without a receipt

`clean_catalogue.py` now handles an existing catalogue whose sidecar is absent:

1. Re-read and hash the source, recompute the mask, and construct a freshly
   validated expected output in a unique temporary file.
2. Compare the orphan with that expected output. ASCII requires an identical
   complete-file SHA256. HDF5 requires every dataset's values, dtype and shape,
   every file/dataset attribute, and non-row metadata to match. HDF5 container
   layout can differ while these contents remain identical.
3. Reuse the verified catalogue without replacing its inode or bytes, and
   publish its missing sidecar. The receipt records the actual existing-file
   SHA256 and that recovery occurred. Its tool revision identifies this
   verification; the unavailable original producer commit is recorded as null.

For new products, the complete sidecar is prepared and closed before the first
publication link. The catalogue and then the sidecar are published with
no-clobber links. A process interruption between those links can still leave an
orphan; the verified recovery path handles it on the next campaign resume.
Handled failures remove only that attempt's uniquely named staging files.
Neither a mismatching orphan nor a completed catalogue/receipt is overwritten.

Fault-injection tests exercise failure during sidecar writing, failure after
catalogue publication but before receipt publication, and subsequent campaign
resume for ASCII and HDF5. Recovered catalogues retain their original SHA256
and inode. A reserialized but equivalent HDF5 file is accepted; altered retained
values, scalar metadata, attributes or dtypes are rejected without modification.

## Validation and scope

All **10 Python tests** pass, including the existing periodic-chain, merger-
control, original-row-mask and byte-preservation checks plus the new precision
and publication tests. The new regression draft failed on the original tools.
The command, logs and exact updated source hashes are saved together in
`provenance/review-fixes-validation.json`.

The completed 740-product campaign read ASCII as float64 and published all
receipts successfully. These fixes do not require rewriting those products;
their original execution hashes, receipts and scientific results are retained.
Any float32 HDF5 outputs produced by the old arithmetic should instead be
regenerated to a new destination and compared, since ordinary resume protects
completed receipts. Casting to float64 cannot recover precision lost in storage.

No Fortran source was changed or full production finder replay attempted in
this follow-up. The previously recorded Fortran regressions and snapshot
availability limitations remain as documented in `RESULTS.md`.
