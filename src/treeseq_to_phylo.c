#include <R.h>
#include <Rinternals.h>

#include <tskit.h>
#include "error.h"


static SEXP
to_phylo1(tsk_tree_t *tree, tsk_id_t root, tsk_id_t *nodes, int *node_map)
{
    int ret;
    tsk_id_t u;
    tsk_size_t num_nodes, j;
    tsk_id_t *parent = tree->parent;
    tsk_id_t N = tree->virtual_root;
    ret = tsk_tree_preorder_from(tree, root, nodes, &num_nodes);
    check_tsk_error(ret);
    int i, k = 0;
    int ntip = 0;
    int nnode = 0;
    double len;
    for (j = 0; j < num_nodes; ++j)
    {
        u = nodes[j];
        if (tree->num_children[u] == 0)
            node_map[u] = ++ntip;
        else
            node_map[u] = ++nnode;
        if (parent[u] == TSK_NULL || parent[u] == N)
            continue;
        ++k;
    }
    for (j = 0; j < num_nodes; ++j)
    {
        u = nodes[j];
        if (tree->num_children[u] == 0)
            continue;
        node_map[u] += ntip;
    }
    SEXP phy = PROTECT(Rf_allocVector(VECSXP, 6));
    SEXP phylo_names = PROTECT(Rf_allocVector(STRSXP, 6));
    SET_VECTOR_ELT(phy, 0, Rf_allocMatrix(INTSXP, k, 2));       // $edge
    SET_VECTOR_ELT(phy, 3, Rf_allocVector(REALSXP, k));         // $edge.length
    SET_VECTOR_ELT(phy, 1, Rf_allocVector(STRSXP, ntip));       // $tip.label
    SET_VECTOR_ELT(phy, 2, Rf_ScalarInteger(nnode));            // $Nnode
    SET_VECTOR_ELT(phy, 4, Rf_allocVector(INTSXP, nnode+ntip)); // $node.id
    SET_VECTOR_ELT(phy, 5, Rf_allocVector(INTSXP, k));          // $edge.id
    SET_STRING_ELT(phylo_names, 0, Rf_mkChar("edge"));
    SET_STRING_ELT(phylo_names, 1, Rf_mkChar("tip.label"));
    SET_STRING_ELT(phylo_names, 2, Rf_mkChar("Nnode"));
    SET_STRING_ELT(phylo_names, 3, Rf_mkChar("edge.length"));
    SET_STRING_ELT(phylo_names, 4, Rf_mkChar("node.id"));
    SET_STRING_ELT(phylo_names, 5, Rf_mkChar("edge.id"));
    char node_label[1024];
    for (j = 0, i = 0; j < num_nodes; ++j)
    {
        u = nodes[j];
        snprintf(node_label, 1023, "%d", (int)u);
        SET_INTEGER_ELT(VECTOR_ELT(phy, 4), node_map[u] - 1, (int)u);
        if (tree->num_children[u] == 0)
        {
            SET_STRING_ELT(VECTOR_ELT(phy, 1), node_map[u] - 1, 
                Rf_mkChar(node_label));
        }
        if (parent[u] == TSK_NULL || parent[u] == N)
            continue;
        tsk_tree_get_branch_length(tree, u, &len);
        INTEGER(VECTOR_ELT(phy, 0))[i + 0*k] = node_map[parent[u]];
        INTEGER(VECTOR_ELT(phy, 0))[i + 1*k] = node_map[u];
        REAL(VECTOR_ELT(phy, 3))[i] = len;
        INTEGER(VECTOR_ELT(phy, 5))[i] = tree->edge[u];
        ++i;
    }
    Rf_setAttrib(phy, R_ClassSymbol, Rf_mkString("phylo"));
    Rf_setAttrib(phy, R_NamesSymbol, phylo_names);
    Rf_setAttrib(phy, Rf_install("order"), Rf_mkString("cladewise"));
    UNPROTECT(2);
    return phy;
}


static SEXP
to_phylo(tsk_tree_t *tree, tsk_id_t *nodes, int *node_map)
{
    tsk_id_t u;
    const tsk_id_t *left_child = tree->left_child;
    const tsk_id_t *right_sib = tree->right_sib;
    const tsk_id_t N = tree->virtual_root;
    int i;
    int num_roots = (int)tsk_tree_get_num_roots(tree);
    if (num_roots > 1)
    {
        SEXP ans = PROTECT(Rf_allocVector(VECSXP, num_roots));
        for (u = left_child[N], i = 0; u != TSK_NULL; u = right_sib[u], ++i)
            SET_VECTOR_ELT(ans, i, to_phylo1(tree, u, nodes, node_map));
        Rf_setAttrib(ans, R_ClassSymbol, Rf_mkString("multiPhylo"));
        UNPROTECT(1);
        return ans;
    }
    return to_phylo1(tree, tree->left_child[N], nodes, node_map);
}


SEXP
C_treeseq_to_phylo(SEXP tr)
{
    tsk_tree_t *tree = (tsk_tree_t *)R_ExternalPtrAddr(tr);   
    tsk_size_t N = tsk_treeseq_get_num_nodes(tree->tree_sequence);
    void *buffer = malloc(N * sizeof(tsk_id_t) + N * sizeof(int));
    tsk_id_t *nodes = (tsk_id_t *)buffer;
    int *node_map = (int *)(nodes + N);
    SEXP ans = PROTECT(to_phylo(tree, nodes, node_map));
    free(buffer);
    UNPROTECT(1);
    return ans;
}
