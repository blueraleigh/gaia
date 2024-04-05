#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <tskit.h>


SEXP
C_treeseq_intervals(SEXP tr)
{
    tsk_tree_t *tree = (tsk_tree_t *)R_ExternalPtrAddr(tr);
    tsk_size_t n = tsk_treeseq_get_num_trees(tree->tree_sequence);
    SEXP ans = PROTECT(Rf_allocMatrix(REALSXP, n, 6));
    int k;
    int ret;
    tsk_id_t u;
    double age;    
    for (ret = tsk_tree_first(tree), k = 0;
         ret == TSK_TREE_OK;
         ret = tsk_tree_next(tree), ++k)
    {
        REAL(ans)[k + 0*n] = tree->interval.left;
        REAL(ans)[k + 1*n] = tree->interval.right;
        REAL(ans)[k + 2*n] = tree->interval.right - tree->interval.left;
        REAL(ans)[k + 3*n] = tree->num_edges;
        REAL(ans)[k + 4*n] = tsk_tree_get_num_roots(tree);
        REAL(ans)[k + 5*n] = -1;
        for (u = tree->left_child[tree->virtual_root]; 
             u != TSK_NULL; 
             u = tree->right_sib[u])
        {
            tsk_tree_get_time(tree, u, &age);
            if (age > REAL(ans)[k + 5*n])
                REAL(ans)[k + 5*n] = age;
        }
    }
    UNPROTECT(1);
    return ans;
}
