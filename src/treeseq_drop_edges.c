#include <R.h>
#include <Rinternals.h>
#include <assert.h>
#include <tskit.h>

#include "error.h"

// defined in tsx_treeseq_load.c
void
tsx_treeseq_exptr_free(SEXP treeseq);

void
tsx_tree_exptr_free(SEXP treeseq);

// drop edges in the tree sequence that leave the given parents
// or enter the given children
static int
tsx_treeseq_drop_edges(
    tsk_treeseq_t *dest,
    const tsk_treeseq_t *src,
    int *restrict drop_parent,
    int *restrict drop_child)
{
    int ret;
    tsk_table_collection_t *tables = malloc(sizeof(*tables));

    if (!tables)
        TSX_ERROR("memory allocation failed");
    
    ret = tsk_table_collection_init(tables, 0);
    if (ret != 0)
    {
        tsk_table_collection_free(tables);
        free(tables);
        goto out;
    }
    const double *restrict edge_left = src->tables->edges.left;
    const double *restrict edge_right = src->tables->edges.right;
    const tsk_id_t *restrict edge_parent = src->tables->edges.parent;
    const tsk_id_t *restrict edge_child = src->tables->edges.child;

    int i;
    tsk_id_t u;
    tsk_id_t v;
    int num_edges = tsk_treeseq_get_num_edges(src);

    tsk_edge_table_t *edges = &tables->edges;

    tsk_individual_table_copy(
        &src->tables->individuals, &tables->individuals, TSK_NO_INIT);
    tsk_node_table_copy(
        &src->tables->nodes, &tables->nodes, TSK_NO_INIT);
    tsk_site_table_copy(
        &src->tables->sites, &tables->sites, TSK_NO_INIT);
    tsk_mutation_table_copy(
        &src->tables->mutations, &tables->mutations, TSK_NO_INIT);
    tsk_migration_table_copy(
        &src->tables->migrations, &tables->migrations, TSK_NO_INIT);
    tsk_population_table_copy(
        &src->tables->populations, &tables->populations, TSK_NO_INIT);

    tables->sequence_length = src->tables->sequence_length;

    tsk_table_collection_set_time_units(
        tables, src->tables->time_units, src->tables->time_units_length);

    for (i = 0; i < num_edges; ++i)
    {
        u = edge_parent[i];
        v = edge_child[i];
        if (drop_parent[u] || drop_child[v])
            continue;
        tsk_edge_table_add_row(edges, edge_left[i], edge_right[i], u, v,
            NULL, 0);
    }
    tables->indexes.num_edges = edges->num_rows;
    tables->indexes.edge_insertion_order = NULL;
    tables->indexes.edge_removal_order = NULL;
    
    tsk_table_sorter_t sorter;

    tsk_table_sorter_init(&sorter, tables, 0);
    tsk_table_sorter_run(&sorter, NULL);
    tsk_table_sorter_free(&sorter);
    tsk_table_collection_build_index(tables, 0);

    tsk_flags_t options = TSK_SIMPLIFY_FILTER_SITES |
                          TSK_SIMPLIFY_FILTER_POPULATIONS |
                          TSK_SIMPLIFY_FILTER_INDIVIDUALS;

    tsk_table_collection_simplify(
        tables, src->samples, src->num_samples, options, NULL);

    ret = tsk_treeseq_init(dest, tables, 
        TSK_TS_INIT_BUILD_INDEXES | TSK_TAKE_OWNERSHIP);

out:
    return ret;
}


SEXP
C_treeseq_drop_edges(SEXP treeseq, SEXP drop_parent, SEXP drop_child)
{
    tsk_treeseq_t *ts = malloc(sizeof(*ts));
    tsk_tree_t *tree = malloc(sizeof(*tree));
    if (!ts || !tree)
        TSX_ERROR("memory allocation failed");
    int rc = tsx_treeseq_drop_edges(
        ts,
        (tsk_treeseq_t *)R_ExternalPtrAddr(treeseq),
        INTEGER(drop_parent),
        INTEGER(drop_child)
    );
    if (rc < 0)
    {
        free(ts);
        free(tree);
        TSX_ERROR(tsk_strerror(rc));
    }
    rc = tsk_tree_init(tree, ts, TSK_SAMPLE_LISTS);
    if (rc < 0)
    {
        free(ts);
        free(tree);
        TSX_ERROR(tsk_strerror(rc));
    }
    SEXP exptr1 = PROTECT(R_MakeExternalPtr(ts, NULL, NULL));
    R_RegisterCFinalizer(exptr1, tsx_treeseq_exptr_free);
    SEXP exptr2 = PROTECT(R_MakeExternalPtr(tree, NULL, NULL));
    R_RegisterCFinalizer(exptr2, tsx_tree_exptr_free);
    UNPROTECT(2);
    return Rf_list2(exptr1, exptr2);
}
