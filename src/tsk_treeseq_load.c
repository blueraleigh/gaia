#include <R.h>
#include <Rinternals.h>

#include <tskit.h>
#include "error.h"

void
tsx_treeseq_exptr_free(SEXP treeseq)
{
    tsk_treeseq_t *ts = (tsk_treeseq_t *)R_ExternalPtrAddr(treeseq);
    tsk_treeseq_free(ts);
    free(ts);
    ts = NULL;
}


void
tsx_tree_exptr_free(SEXP tr)
{
    tsk_tree_t *tree = (tsk_tree_t *)R_ExternalPtrAddr(tr);
    tsk_tree_free(tree);
    free(tree);
    tree = NULL;
}


SEXP
C_treeseq_load(SEXP filename)
{
    tsk_treeseq_t *ts = malloc(sizeof(*ts));
    tsk_tree_t *tree = malloc(sizeof(*tree));
    if (!ts || !tree)
        TSX_ERROR("memory allocation failed");
    int rc = tsk_treeseq_load(ts, CHAR(STRING_ELT(filename, 0)), 0);
    check_tsk_error(rc);
    rc = tsk_tree_init(tree, ts, TSK_SAMPLE_LISTS);
    check_tsk_error(rc);
    SEXP exptr1 = PROTECT(R_MakeExternalPtr(ts, NULL, NULL));
    R_RegisterCFinalizer(exptr1, tsx_treeseq_exptr_free);
    SEXP exptr2 = PROTECT(R_MakeExternalPtr(tree, NULL, NULL));
    R_RegisterCFinalizer(exptr2, tsx_tree_exptr_free);
    UNPROTECT(2);
    return Rf_list2(exptr1, exptr2);
}
