/* rs_wrap_subhalo_metric.c: MPI-Rockstar's subhalo_metric.c, plus the release of
 * its per-thread search results (free_subtree releases the tree and the metric
 * array, not the results).
 * Part of the MG-GLAM in-line Rockstar adapter; GNU GPL v3 or later. */
#include "subhalo_metric.c"

/* call on every OpenMP thread (the pointers are threadprivate) */
void rs_free_subhalo_scratch(void)
{
    free_subtree();
    if (subtree_res) {
        fast3tree_results_free(subtree_res);
        subtree_res = NULL;
    }
}
