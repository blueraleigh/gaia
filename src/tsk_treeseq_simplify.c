#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>

#include <tskit.h>
#include "error.h"

// defined in tsx_treeseq_load.c
void
tsx_treeseq_exptr_free(SEXP treeseq);

void
tsx_tree_exptr_free(SEXP treeseq);


SEXP
C_treeseq_simplify(SEXP treeseq, SEXP sample_ids, SEXP filter_flags)
{
    int ret;
    tsk_size_t num_samples = (tsk_size_t)Rf_length(sample_ids);
    tsk_treeseq_t *ts = (tsk_treeseq_t *)R_ExternalPtrAddr(treeseq);
    tsk_size_t N = tsk_treeseq_get_num_nodes(ts);
    const tsk_id_t *samples = (tsk_id_t *)INTEGER(sample_ids);
    /*
    tsk_flags_t options = TSK_SIMPLIFY_FILTER_SITES |
                          TSK_SIMPLIFY_FILTER_POPULATIONS |
                          TSK_SIMPLIFY_FILTER_INDIVIDUALS;*/
    tsk_treeseq_t *tss = malloc(sizeof(*tss));
    tsk_flags_t options = (tsk_flags_t)(*INTEGER(filter_flags));
    tsk_tree_t *tree = malloc(sizeof(*tree));
    tsk_id_t *node_map = (tsk_id_t *) R_alloc(N, sizeof(*node_map));
    if (!tss || !tree)
        TSX_ERROR("memory allocation failed");
    ret = tsk_treeseq_simplify(ts,samples,num_samples,options,tss,node_map);
    check_tsk_error(ret);
    ret = tsk_tree_init(tree, tss, TSK_SAMPLE_LISTS);
    check_tsk_error(ret);
    SEXP exptr1 = PROTECT(R_MakeExternalPtr(tss, NULL, NULL));
    R_RegisterCFinalizer(exptr1, tsx_treeseq_exptr_free);
    SEXP exptr2 = PROTECT(R_MakeExternalPtr(tree, NULL, NULL));
    R_RegisterCFinalizer(exptr2, tsx_tree_exptr_free);
    SEXP nodemap = PROTECT(Rf_allocVector(INTSXP, N));
    memcpy(INTEGER(nodemap), (int *)node_map, N * sizeof(int));
    UNPROTECT(3);
    return Rf_list3(exptr1, exptr2, nodemap);
}


/* The edge table postfixed by 1 is simplified from the edge table post fixed
** by 2. We want to fix the edge id's in the simplified table so that we can
** match them to the edges in the full table. This is only valid when the
** tree sequence was simplified using TSK_SIMPLIFY_KEEP_UNARY */
static void
fixup_edge_ids(
    int num_edges1,
    int *edge_id1,
    double *left1,
    double *right1,
    tsk_id_t *parent_id1,
    tsk_id_t *child_id1,
    int num_edges2,
    double *left2,
    double *right2,
    tsk_id_t *parent_id2,
    tsk_id_t *child_id2,
    tsk_id_t *node_map)
{
    int i;
    int j;
    j = 0;
    for (i = 0; i < num_edges1; ++i)
    {
        while (!(
            parent_id1[i] == node_map[parent_id2[j]] && 
            child_id1[i] == node_map[child_id2[j]] &&
            (left1[i] >= left2[j] && right1[i] <= right2[j])))
        {
            ++j;
            if (j >= num_edges2) Rf_error("invalid edge id");
        }
        edge_id1[i] = j;
    }
}


SEXP
C_treeseq_fixup_edge_ids(
    SEXP treeseq_simplified, SEXP treeseq_full, SEXP node_map)
{

    tsk_treeseq_t *tss = (tsk_treeseq_t *)R_ExternalPtrAddr(treeseq_simplified);
    tsk_treeseq_t *ts = (tsk_treeseq_t *)R_ExternalPtrAddr(treeseq_full);

    int num_edges = tss->tables->edges.num_rows;

    SEXP edge_id = PROTECT(Rf_allocVector(INTSXP, num_edges));

    fixup_edge_ids(
        num_edges,
        INTEGER(edge_id),
        tss->tables->edges.left,
        tss->tables->edges.right,
        tss->tables->edges.parent,
        tss->tables->edges.child,
        ts->tables->edges.num_rows,
        ts->tables->edges.left,
        ts->tables->edges.right,
        ts->tables->edges.parent,
        ts->tables->edges.child,
        INTEGER(node_map)
    );

    UNPROTECT(1);
    return edge_id;
}
