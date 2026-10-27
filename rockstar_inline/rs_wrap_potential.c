/* rs_wrap_potential.c: MPI-Rockstar's potential.c, plus the release of its
 * per-thread tree, which MPI-Rockstar never frees (it lives until the process
 * ends; inside GLAM it would stay resident through every later force step).
 * Part of the MG-GLAM in-line Rockstar adapter; GNU GPL v3 or later. */
#include "potential.c"

/* call on every OpenMP thread (the tree pointers are threadprivate) */
void rs_free_potential_scratch(void)
{
    if (p_tree) fast3tree_free(&p_tree);
    if (p_res) {
        fast3tree_results_free(p_res);
        p_res = NULL;
    }
}
