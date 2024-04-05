#include <R.h>
#include <Rinternals.h>

#include <tskit.h>
#include "error.h"


SEXP
C_treeseq_indexes(SEXP treeseq)
{
    tsk_treeseq_t *ts = (tsk_treeseq_t *)R_ExternalPtrAddr(treeseq);

    int num_edges = ts->tables->indexes.num_edges;
    const tsk_id_t *restrict I = ts->tables->indexes.edge_insertion_order;
    const tsk_id_t *restrict O = ts->tables->indexes.edge_removal_order;

    SEXP ans = PROTECT(Rf_allocMatrix(INTSXP, num_edges, 2));

    memcpy(INTEGER(ans), (int *)I, num_edges * sizeof(int));
    memcpy(INTEGER(ans) + num_edges, (int *)O, num_edges * sizeof(int));

    UNPROTECT(1);
    return ans;
}
