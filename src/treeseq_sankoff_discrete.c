#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <assert.h>
#include <tskit.h>

#include "treeseq_sankoff.h"
#include "error.h"
#include "fequals.h"

/* Sample MPR histories of geographic locations of genetic ancestors
** using generalized (Sankoff) parsimony.
**
** Geographic locations are stored as the integers 0, 1, ..., n-1, 
** where n is the number of locations.
**
** Requires a user-defined cost matrix giving the cost of transitioning
** from each state to every other state.
**
** For more information on the basic algorithm for a single tree see
**
**   Sankoff, D. and P. Rousseau 1975. Math. Programming 9: 240-246.
*/

typedef struct mpr {
    int num_states;
    double *g;
    double *h;
    double *f;
    const double *cost;
} mpr_t;


typedef struct mpr_summary {
    double *F;
    double *node_weight;
    double *tree_length;
    double mean_tree_length;
    double tree_weight;
} mpr_summary_t;


static void
calc_stem_cost(tsk_id_t u, tsk_id_t v, tsx_tree_t *tree)
{
    assert(u != TSK_NULL);
    assert(v != TSK_NULL);
    mpr_t *mpr = (mpr_t *)(tree->mpr);
    int i;
    int j;
    int n = mpr->num_states;
    int vn = v*n;
    double *restrict g_v = mpr->g + vn;
    double *restrict h_v = mpr->h + vn;
    const double *restrict cost = mpr->cost;
    double t, min_t, scale = 1;
    if (tree->time)
    {
        scale = 1 / (tree->time[u] - tree->time[v]);
    }
    for (i = 0; i < n; ++i)
    {
        min_t = R_PosInf;
        for (j = 0; j < n; ++j)
        {
            t = scale * cost[i+j*n] + g_v[j];
            if (t < min_t)
                min_t = t;
        }
        h_v[i] = min_t;
    }
}


static void
calc_final_cost(tsk_id_t v, tsx_tree_t *tree)
{
    mpr_t *mpr = (mpr_t *)(tree->mpr);
    int i;
    int j;
    tsk_id_t u = tree->parent[v];
    int n = mpr->num_states;
    int vn = v*n;
    if (u == TSK_NULL)
    {
        memcpy(mpr->f + vn, mpr->g + vn, n*sizeof(double));
        return;
    }
    int un = u*n;
    double *restrict g_v = mpr->g + vn;
    double *restrict h_v = mpr->h + vn;
    double *restrict f_v = mpr->f + vn;
    double *restrict f_u = mpr->f + un;
    const double *restrict cost = mpr->cost;
    double t, min_t, scale = 1;
    if (tree->time)
    {
        scale = 1 / (tree->time[u] - tree->time[v]);
    }
    for (i = 0; i < n; ++i)
    {
        min_t = R_PosInf;
        for (j = 0; j < n; ++j)
        {
            t = f_u[j] - h_v[j] + scale * cost[j+i*n] + g_v[i];
            if (t < min_t)
                min_t = t;
        }
        f_v[i] = min_t;
    }
}


static void
increment_node_cost(tsk_id_t u, tsk_id_t v, tsx_tree_t *tree, int sign)
{
    mpr_t *mpr = (mpr_t *)(tree->mpr);
    int i;
    int n = mpr->num_states;
    double *restrict g_u = mpr->g + u*n;
    double *restrict h_v = mpr->h + v*n;
    for (i = 0; i < n; ++i)
        g_u[i] += sign*h_v[i];
}


static void
mean_tree_length(
    tsx_tree_t *tree, int t_index, double t_left, double t_right, void *params)
{
    int i;
    tsk_id_t u;
    tsk_id_t virtual_root = tree->virtual_root;
    tsk_id_t *left_child = tree->left_child;
    tsk_id_t *right_sib = tree->right_sib;
    mpr_t *mpr = (mpr_t *)(tree->mpr);
    double min_s;
    double *g;
    double len = 0;
    int n = mpr->num_states;
    for (u = left_child[virtual_root]; u != TSK_NULL; u = right_sib[u])
    {
        min_s = R_PosInf;
        g = mpr->g + u*n;
        for (i = 0; i < n; ++i)
        {
            if (g[i] < min_s)
                min_s = g[i];
        }
        len += min_s;
    }
    len /= tree->num_edges;
    
    double w = t_right - t_left;

    mpr_summary_t *mpr_summary = (mpr_summary_t *)params;

    mpr_summary->tree_weight += w;
    mpr_summary->tree_length[t_index] = len;
    mpr_summary->mean_tree_length += 
        w * (len - mpr_summary->mean_tree_length) / mpr_summary->tree_weight;
}


static void
update_mean_costs(tsx_tree_t *tree, int t_index, double t_left, double t_right, 
    tsk_id_t node_id, void *params)
{
    mpr_t *mpr = (mpr_t *)(tree->mpr);
    mpr_summary_t *mpr_summary = (mpr_summary_t *)params;

    int i;
    int n = mpr->num_states;

    double w = t_right - t_left;

    mpr_summary->node_weight[node_id] += w;

    double total_weight = mpr_summary->node_weight[node_id];

    const double *restrict f = mpr->f + n * node_id;
    double *restrict F = mpr_summary->F + n * node_id;

    for (i = 0; i < n; ++i)
    {
        if (R_FINITE(F[i]))
            F[i] += w * (f[i] - F[i]) / total_weight;
    }
}


SEXP C_treeseq_discrete_mpr(
    SEXP treeseq,
    SEXP use_brlen, 
    SEXP num_states,
    SEXP g,
    SEXP cost)
{
    tsx_tree_t tree;
    mpr_t mpr;
    mpr_summary_t mpr_summary;
    tsk_treeseq_t *ts = (tsk_treeseq_t *) R_ExternalPtrAddr(treeseq);
    int num_nodes = (int) tsk_treeseq_get_num_nodes(ts);
    int num_trees = (int) tsk_treeseq_get_num_trees(ts);

    int n = *INTEGER(num_states);
    
    SEXP h = PROTECT(Rf_allocMatrix(REALSXP, n, num_nodes));
    SEXP f = PROTECT(Rf_allocMatrix(REALSXP, n, num_nodes));

    SEXP F = PROTECT(Rf_allocMatrix(REALSXP, n, num_nodes));
    SEXP node_weight = PROTECT(Rf_allocVector(REALSXP, num_nodes));
    SEXP tree_length = PROTECT(Rf_allocVector(REALSXP, num_trees));

    memset(REAL(h), 0, n * num_nodes * sizeof(double));
    memset(REAL(f), 0, n * num_nodes * sizeof(double));
    memset(REAL(F), 0, n * num_nodes * sizeof(double));
    memset(REAL(node_weight), 0, num_nodes * sizeof(double));
    memset(REAL(tree_length), 0, num_trees * sizeof(double));

    mpr_summary.F = REAL(F);
    mpr_summary.node_weight = REAL(node_weight);

    mpr_summary.tree_weight = 0;
    mpr_summary.mean_tree_length = 0;
    mpr_summary.tree_length = REAL(tree_length);

    mpr.g = REAL(g);
    mpr.h = REAL(h);
    mpr.f = REAL(f);
    mpr.num_states = n;
    mpr.cost = REAL(cost);

    int rc = tsx_tree_init(&tree, ts, *INTEGER(use_brlen));

    if (rc)
        goto out;

    tree.mpr = (void *)(&mpr);
    tree.calc_stem_cost = &calc_stem_cost;
    tree.calc_final_cost = &calc_final_cost;
    tree.increment_node_cost = &increment_node_cost;

    tsx_treeseq_sankoff(
        ts,
        &tree,
        (void *)(&mpr_summary),
        &mean_tree_length,
        (void *)(&mpr_summary),
        &update_mean_costs
    );

out:
    tsx_tree_free(&tree);
    UNPROTECT(5);
    return Rf_list3(
        Rf_ScalarReal(mpr_summary.mean_tree_length), tree_length, F);
}


SEXP C_treeseq_discrete_mpr_minimize(SEXP F, SEXP index1)
{
    int i;
    int j;
    int k;
    int indx1 = *INTEGER(index1);
    int num_nodes = INTEGER(Rf_getAttrib(F, R_DimSymbol))[1];
    int num_states = *INTEGER(Rf_getAttrib(F, R_DimSymbol));

    SEXP ans = PROTECT(Rf_allocVector(INTSXP, num_nodes));
    int *x = INTEGER(ans);
    const double *p = REAL(F);
    const double *q;
    double min_score;
    GetRNGstate();
    for (i = 0; i < num_nodes; ++i)
    {
        k = 0;
        min_score = R_PosInf;
        q = p + i*num_states;
        for (j = 0; j < num_states; ++j)
        {
            if (q[j] < min_score)
                min_score = q[j];
        }
        for (j = 0; j < num_states; ++j)
        {
            if (fequals(q[j], min_score))
            {
                ++k;
                if (unif_rand() < (1 / (double)k))
                {
                    x[i] = j + indx1;
                }
            }
        }
    }
    PutRNGstate();
    UNPROTECT(1);
    return ans;
}
