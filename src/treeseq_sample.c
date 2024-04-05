#include <R.h>
#include <Rinternals.h>

#include <tskit.h>
#include "error.h"


/* Randomly sample a tree from a tree sequence. Trees are sampled using 
** weights proportional to the length of the genome each tree covers */
static int
tsx_tree_sample(tsk_tree_t *tree, int *tree_index)
{
    int ret;
    int k = 1;
    double L = tree->tree_sequence->tables->sequence_length;
    double w = unif_rand() * L;
    ret = tsk_tree_first(tree);
    check_tsk_error(ret);
    w -= tree->interval.right - tree->interval.left;
    while (ret == TSK_TREE_OK && w > 0)
    {
        ++k;
        ret = tsk_tree_next(tree);
        check_tsk_error(ret);
        w -= tree->interval.right - tree->interval.left;
    }
    *tree_index = k;
    return ret;
}


// here we use the convention that the first tree has tree_index = 1
static int
tsx_tree_at(tsk_tree_t *tree, int tree_index)
{
    int ret;
    int k = 0;
    ret = tsk_tree_first(tree);
    check_tsk_error(ret);
    while (ret == TSK_TREE_OK && ++k < tree_index)
    {
        ret = tsk_tree_next(tree);
        check_tsk_error(ret);
    }
    return ret;
}


SEXP
C_treeseq_sample(SEXP tr, SEXP tree_index)
{
    int ret;
    int k = *INTEGER(tree_index);
    tsk_tree_t *tree = (tsk_tree_t *)R_ExternalPtrAddr(tr);
    GetRNGstate();
    if (k <= 0)
        ret = tsx_tree_sample(tree, &k);
    else
        ret = tsx_tree_at(tree, k);
    PutRNGstate();
    check_tsk_error(ret);
    return Rf_ScalarInteger(k);
}
