#include <R.h>
#include <Rinternals.h>

#include <assert.h>
#include <tskit.h>

#include "treeseq_sankoff.h"
#include "error.h"
#include "fequals.h"

/* Sample MPR histories of geographic locations of genetic ancestors
** using linear parsimony.
**
** For each non-sample node we calculate the parameters of a convex
** piecewise linear function 


       { b + a[0] * x                    if x <= x[0]
       { s[0] + a[1] * (x - x[0])        if x[0] < x <= x[1] 
F(x) = { ...
       { s[k-1] + a[k] * (x - x[k-1])    if x[k-1] < x <= x[k]
       { s[k]   + a[k+1] * (x - x[k])    if x > x[k]

where s[0] = b + a[0]*x[0] and s[i+1] = s[i] + a[i+1]*(x[i+1] - x[i])
for 0 <= i < k

** Which gives the minimum sum of distances between all
** ancestor-descendant pairs needed to explain the spatial
** distribution of samples assuming the non-sample node is in
** location x.
**
** For more information on the basic algorithm for a single tree see
**
**  Miklós Csurös. 2008. Ancestral Reconstruction by Asymmetric Wagner
**  Parsimony over Continuous Characters and Squared Parsimony over
**  Distributions. Pp 72-86. In: Nelson, C.E., Vialette, S. (eds)
**  Comparative Genomics.
*/


// convex piecewise linear function.
// slopes and breakpoints are stored in non-decreasing order
typedef struct plf {
    double *slope;
    double *breakpoint;
    double intercept;
    int num_breaks;
    int max_num_breaks;
} plf_t;


typedef struct mpr {
    plf_t **g;
    plf_t **h;
    plf_t **f;
    // workspace
    plf_t *tmp;
    // observed states
    const double *x;
    int num_dims;
} mpr_t;


typedef struct mpr_summary {
    plf_t **F;
    plf_t *tmp;
    double *node_weight;
    double *tree_length;
    double mean_tree_length;
    double tree_weight;
    double tol;
} mpr_summary_t;


static int
colinear(double a, double b, double tol)
{
    int sign_a = (0 < a) - (a < 0);
    int sign_b = (0 < b) - (b < 0);
    if (sign_a != sign_b) return 0;
    if (fabs(a - b) > tol) return 0;
    return 1;
}


static int
plf_init(plf_t *f, int n)
{
    f->slope = calloc(n+1, sizeof(double));
    f->breakpoint = calloc(n, sizeof(double));
    f->intercept = 0;
    f->num_breaks = 0;
    f->max_num_breaks = n;
    if (f->slope == NULL || f->breakpoint == NULL)
    {
        TSX_WARN("memory allocation failed");
        return 1;
    }
    return 0;
}


static void
plf_free(plf_t *f)
{
    free(f->slope);
    free(f->breakpoint);
}


static double
plf_min_score(plf_t *f)
{
    int i;
    double shift;
    double *restrict fs = f->slope;
    double *restrict fb = f->breakpoint;
    i = 0;
    shift = f->intercept + fs[0]*fb[0];
    while (fs[i+1] < 0)
    {
        ++i;
        shift += fs[i] * (fb[i] - fb[i-1]);
    }
    return shift;
}


static int
plf_min_index(plf_t *f)
{
    int i;
    double *restrict fs = f->slope;
    i = 0;
    while (fs[i+1] < 0)
    {
        ++i;
    }
    return i;
}


// ret(x) = min_z { c(x,z) + f(z) } where the cost c(x,z)
// of going from x to z is b * abs(x - z)
static void
plf_min(plf_t *f, double b, plf_t *ret)
{
    assert(b > 0);
    assert(f != ret);

    int i;
    int x_left;
    int num_breaks = 0;
    double shift;
    double minus_b = -b;

    double *restrict fs = f->slope;
    double *restrict fb = f->breakpoint;

    double *restrict rs = ret->slope;
    double *restrict rb = ret->breakpoint;

    if (minus_b <= fs[0])
    {
        x_left = 0;
        ret->intercept = f->intercept;
    }
    else
    {
        i = 0;
        shift = f->intercept + fs[0]*fb[0];
        while (minus_b >= fs[i+1])
        {
            ++i;
            shift += fs[i] * (fb[i] - fb[i-1]);
        }
        x_left = i + 1;
        ret->intercept = shift + b * fb[i];
        rs[0] = minus_b;
        rb[0] = fb[i];
        ++num_breaks;
        assert (num_breaks <= ret->max_num_breaks);
    }
    if (b >= fs[f->num_breaks])
    {
        for (i = x_left; i < f->num_breaks; ++i)
        {
            rs[num_breaks] = fs[i];
            rb[num_breaks] = fb[i];
            ++num_breaks;
            assert (num_breaks <= ret->max_num_breaks);
        }
        rs[num_breaks] = fs[f->num_breaks];
    }
    else
    {
        i = x_left;
        while (b >= fs[i])
        {
            assert(i < f->num_breaks);
            rs[num_breaks] = fs[i];
            rb[num_breaks] = fb[i];
            ++num_breaks;
            ++i;
            assert (num_breaks <= ret->max_num_breaks);
        }
        // check for colinearity to avoid superfluous breakpoints
        if (fequals(b, rs[num_breaks - 1]))
            --num_breaks;
        else
            rs[num_breaks] = b;
    }

    assert (num_breaks >= 0);
    
#ifndef NDEBUG
    if (num_breaks == 0) assert (fequals(rs[0], 0));
#endif

    ret->num_breaks = num_breaks;
}


// ret = ascale*a + sign*bscale*b
static void
plf_add(plf_t *a, plf_t *b, int sign, double ascale, double bscale, plf_t *ret)
{
    assert(a != ret);
    assert(b != ret);
    assert(a != b);
    
    int i;
    int j;
    int num_breaks;

    double prev_slope;
    double next_slope;
    double next_break;

    double *restrict breakpoint = ret->breakpoint;
    double *restrict slope = ret->slope;

    double *restrict ab = a->breakpoint;
    double *restrict as = a->slope;
    double *restrict bb = b->breakpoint;
    double *restrict bs = b->slope;

    i = 0;
    j = 0;
    num_breaks = 0;
    prev_slope = R_NegInf;

    while (i < a->num_breaks && j < b->num_breaks)
    {
        if (ab[i] < bb[j])
        {
            next_break = ab[i];
            next_slope = ascale * as[i] + sign * bscale * bs[j];
            ++i;
        }
        else if (ab[i] > bb[j])
        {
            next_break = bb[j];
            next_slope = ascale * as[i] + sign * bscale * bs[j];
            ++j;
        }
        else
        {
            next_break = ab[i];
            next_slope = ascale * as[i] + sign * bscale * bs[j];
            ++i;
            ++j;
        }
        if (fequals(prev_slope, next_slope))
        {
            assert (num_breaks > 0);
            breakpoint[num_breaks-1] = next_break;
        }
        else
        {
            breakpoint[num_breaks] = next_break;
            slope[num_breaks] = next_slope;
            prev_slope = next_slope;
            ++num_breaks;
            assert (num_breaks <= ret->max_num_breaks);
        }
    }
    
    while (i < a->num_breaks)
    {
        next_break = ab[i];
        next_slope = ascale * as[i] + sign * bscale * bs[j];
        ++i;
        if (fequals(prev_slope, next_slope))
        {
            assert (num_breaks > 0);
            breakpoint[num_breaks-1] = next_break;
        }
        else
        {
            breakpoint[num_breaks] = next_break;
            slope[num_breaks] = next_slope;
            prev_slope = next_slope;
            ++num_breaks;
            assert (num_breaks <= ret->max_num_breaks);
        }
    }
    while (j < b->num_breaks)
    {
        next_break = bb[j];
        next_slope = ascale * as[i] + sign * bscale * bs[j];
        ++j;
        if (fequals(prev_slope, next_slope))
        {
            assert (num_breaks > 0);
            breakpoint[num_breaks-1] = next_break;
        }
        else
        {
            breakpoint[num_breaks] = next_break;
            slope[num_breaks] = next_slope;
            prev_slope = next_slope;
            ++num_breaks;
            assert (num_breaks <= ret->max_num_breaks);
        }
    }
    
    // don't forget the final slope after the last breakpoint
    next_slope = ascale * as[i] + sign * bscale * bs[j];

    //assert (num_breaks > 0);

    if (fequals(prev_slope, next_slope))
        --num_breaks;
    else
        slope[num_breaks] = next_slope;

    assert (num_breaks >= 0);

#ifndef NDEBUG
    if (num_breaks == 0) assert (fequals(slope[0], 0));
#endif

    ret->num_breaks = num_breaks;
    ret->intercept = ascale * a->intercept + sign * bscale * b->intercept;
}


static void
plf_shift(plf_t *a, double h)
{
    int i;
    a->intercept += a->slope[0]*h;
    for (i = 0; i < a->num_breaks; ++i)
        a->breakpoint[i] += h;
}


// ret = (1-t)*a + t*b
static void
plf_lerp(plf_t *a, plf_t *b, double t, plf_t *ret)
{
    assert (t >= 0 && t <= 1);
    plf_add(a, b, +1, 1-t, t, ret);
}


static void
plf_defrag(plf_t *f, plf_t *tmp, double tol)
{
    int i;
    int num_breaks = 0;
    double *restrict s = f->slope;
    double *restrict b = f->breakpoint;
    double *restrict ts = tmp->slope;
    double *restrict tb = tmp->breakpoint;
    ts[0] = s[0];
    tb[0] = b[0];
    ++num_breaks;
    for (i = 1; i < f->num_breaks; ++i)
    {
        if (colinear(s[i-1], s[i], tol))
        {
            assert (num_breaks > 0);
            tb[num_breaks-1] = b[i];
        }
        else
        {
            ts[num_breaks] = s[i];
            tb[num_breaks] = b[i];
            ++num_breaks;
        }
    }
    if (colinear(s[i], s[i-1], tol))
    {
        assert (num_breaks > 0);
        --num_breaks;
    }
    else
    {
        ts[num_breaks] = s[i];
    }
    f->num_breaks = num_breaks;
    memcpy(b, tb, num_breaks * sizeof(*b));
    memcpy(s, ts, (num_breaks + 1) * sizeof(*s));
}


static void
calc_stem_cost(tsk_id_t u, tsk_id_t v, tsx_tree_t *tree)
{
    assert(u != TSK_NULL);
    assert(v != TSK_NULL);
    
    mpr_t *mpr = (mpr_t *)(tree->mpr);

    plf_t *g;
    plf_t *h;

    int i;
    int num_dims = mpr->num_dims;
    double b = 1;
    double minus_b;

    if (tree->time)
        b = 1 / (tree->time[u] - tree->time[v]);

    minus_b = -b;

    if (!(tree->node_flags[v] & TSK_NODE_IS_SAMPLE))
    {
        for (i = 0; i < num_dims; ++i)
        {
            g = mpr->g[i] + v;
            h = mpr->h[i] + v;
            plf_min(g, b, h);
        }
    }
    else
    {
        const double *restrict x = mpr->x + v*num_dims;   
        for (i = 0; i < num_dims; ++i)
        {
            h = mpr->h[i] + v;
            h->breakpoint[0] = x[i];
            h->slope[0] = minus_b;
            h->slope[1] = b;
            h->intercept = x[i] < 0 ? minus_b*x[i] : b*x[i];
            h->num_breaks = 1;
        }
    }
}


static void
increment_node_cost(tsk_id_t u, tsk_id_t v, tsx_tree_t *tree, int sign)
{
    mpr_t *mpr = (mpr_t *)(tree->mpr);

    int i;
    int num_dims = mpr->num_dims;

    plf_t *g;
    plf_t *h;
    plf_t *tmp = mpr->tmp;
    
    for (i = 0; i < num_dims; ++i)
    {
        g = mpr->g[i] + u;
        h = mpr->h[i] + v;

        plf_add(g, h, sign, 1, 1, tmp);

        g->num_breaks = tmp->num_breaks;
        g->intercept = tmp->intercept;
        
        memcpy(g->breakpoint, tmp->breakpoint, tmp->num_breaks*sizeof(double));
        memcpy(g->slope, tmp->slope, (tmp->num_breaks+1)*sizeof(double));
    }
}


static void
calc_final_cost(tsk_id_t v, tsx_tree_t *tree)
{
    if (tree->node_flags[v] & TSK_NODE_IS_SAMPLE) return;

    mpr_t *mpr = (mpr_t *)(tree->mpr);

    int i;
    int num_dims = mpr->num_dims;

    tsk_id_t u = tree->parent[v];

    plf_t *g;
    plf_t *h;
    plf_t *f;

    if (u == TSK_NULL)
    {
        for (i = 0; i < num_dims; ++i)
        {
            g = mpr->g[i] + v;
            f = mpr->f[i] + v;
            f->num_breaks = g->num_breaks;
            f->intercept = g->intercept;
            memcpy(f->slope, g->slope,
                (g->num_breaks+1)*sizeof(double));
            memcpy(f->breakpoint, g->breakpoint,
                g->num_breaks*sizeof(double));
        }
        return;
    }

    plf_t *fu;
    plf_t *tmp = mpr->tmp;

    double b = 1;

    if (tree->time)
        b = 1 / (tree->time[u] - tree->time[v]);

    if (tree->node_flags[v] & TSK_NODE_IS_SAMPLE)
    {
        for (i = 0; i < num_dims; ++i)
        {
            h = mpr->h[i] + v;
            f = mpr->f[i] + v;
            fu = mpr->f[i] + u;
            plf_add(fu, h, -1, 1, 1, tmp);
            plf_min(tmp, b, f);
        }
    }
    else
    {
        for (i = 0; i < num_dims; ++i)
        {
            g = mpr->g[i] + v;
            h = mpr->h[i] + v;
            f = mpr->f[i] + v;
            fu = mpr->f[i] + u;
            plf_add(fu, h, -1, 1, 1, f);
            plf_min(f, b, tmp);
            plf_add(tmp, g, +1, 1, 1, f);
        }
    }
}


static void
mean_tree_length(
    tsx_tree_t *tree, int t_index, double t_left, double t_right, void *params)
{
    tsk_id_t u;
    tsk_id_t virtual_root = tree->virtual_root;
    tsk_id_t *left_child = tree->left_child;
    tsk_id_t *right_sib = tree->right_sib;
    mpr_t *mpr = (mpr_t *)(tree->mpr);
    int i;
    int num_dims = mpr->num_dims;
    double len = 0;
    for (u = left_child[virtual_root]; u != TSK_NULL; u = right_sib[u])
    {
        if (!(tree->node_flags[u] & TSK_NODE_IS_SAMPLE))
        {
            for (i = 0; i < num_dims; ++i)
                len += plf_min_score(mpr->g[i] + u);
        }
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
    if (tree->node_flags[node_id] & TSK_NODE_IS_SAMPLE) return;

    mpr_t *mpr = (mpr_t *)(tree->mpr);
    mpr_summary_t *mpr_summary = (mpr_summary_t *)params;

    int i;
    int num_dims = mpr->num_dims;

    double w = t_right - t_left;

    mpr_summary->node_weight[node_id] += w;

    double total_weight = mpr_summary->node_weight[node_id];

    plf_t *f;
    plf_t *F;

    for (i = 0; i < num_dims; ++i)
    {
        f = mpr->f[i] + node_id;
        F = mpr_summary->F[i] + node_id;

        plf_add(f, F, -1, 1, 1, mpr->tmp);
        plf_add(F, mpr->tmp, +1, 1, w / total_weight, mpr_summary->tmp);

        F->num_breaks = mpr_summary->tmp->num_breaks;
        F->intercept = mpr_summary->tmp->intercept;
        
        memcpy(F->slope, mpr_summary->tmp->slope,
            (mpr_summary->tmp->num_breaks+1)*sizeof(double));
        memcpy(F->breakpoint, mpr_summary->tmp->breakpoint,
            (mpr_summary->tmp->num_breaks)*sizeof(double));

        plf_defrag(F, mpr_summary->tmp, mpr_summary->tol);
    }
}


SEXP C_treeseq_linear_mpr(
    SEXP treeseq,
    SEXP use_brlen, 
    SEXP x,
    SEXP nx,
    SEXP tol)
{
    tsx_tree_t tree;
    mpr_t mpr;
    mpr_summary_t mpr_summary;
    tsk_treeseq_t *ts = (tsk_treeseq_t *) R_ExternalPtrAddr(treeseq);
    int i;
    int j;
    int num_nodes = (int) tsk_treeseq_get_num_nodes(ts);
    int num_trees = (int) tsk_treeseq_get_num_trees(ts);

    int num_dims = *INTEGER(Rf_getAttrib(x, R_DimSymbol));

    int *num_sites = INTEGER(nx);
    int max_num_sites = 0;
    
    plf_t **g = calloc(num_dims, sizeof(plf_t *));
    plf_t **h = calloc(num_dims, sizeof(plf_t *));
    plf_t **f = calloc(num_dims, sizeof(plf_t *));
    plf_t **F = calloc(num_dims, sizeof(plf_t *));

    for (j = 0; j < num_dims; ++j)
    {
        g[j] = calloc(num_nodes, sizeof(plf_t));
        h[j] = calloc(num_nodes, sizeof(plf_t));
        f[j] = calloc(num_nodes, sizeof(plf_t));
        F[j] = calloc(num_nodes, sizeof(plf_t));
        for (i = 0; i < num_nodes; ++i)
        {
            if (plf_init(g[j] + i, num_sites[j])) goto out;
            if (plf_init(h[j] + i, num_sites[j])) goto out;
            if (plf_init(f[j] + i, num_sites[j])) goto out;
            if (plf_init(F[j] + i, num_sites[j])) goto out;
        }
        if (num_sites[j] > max_num_sites)
            max_num_sites = num_sites[j];
    }

    plf_t tmp1;
    plf_t tmp2;

    if (plf_init(&tmp1, max_num_sites)) goto out;
    if (plf_init(&tmp2, max_num_sites)) goto out;

    SEXP node_weight = PROTECT(Rf_allocVector(REALSXP, num_nodes));
    SEXP tree_length = PROTECT(Rf_allocVector(REALSXP, num_trees));

    memset(REAL(node_weight), 0, num_nodes * sizeof(double));
    memset(REAL(tree_length), 0, num_trees * sizeof(double));

    mpr_summary.F = F;
    mpr_summary.tmp = &tmp2;
    mpr_summary.node_weight = REAL(node_weight);
    mpr_summary.tree_weight = 0;
    mpr_summary.mean_tree_length = 0;
    mpr_summary.tree_length = REAL(tree_length);
    mpr_summary.tol = *REAL(tol);

    mpr.g = g;
    mpr.h = h;
    mpr.f = f;
    mpr.tmp = &tmp1;
    mpr.x = REAL(x);
    mpr.num_dims = num_dims;

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

    SEXP ans = PROTECT(Rf_allocVector(VECSXP, num_dims));

    for (j = 0; j < num_dims; ++j)
    {
        SET_VECTOR_ELT(
            ans
            , j
            , Rf_allocVector(VECSXP, num_nodes)
        );
        for (i = 0; i < num_nodes; ++i)
        {
            SET_VECTOR_ELT(
                VECTOR_ELT(ans, j)
                , i
                , Rf_allocVector(VECSXP, 3)
            );
            SET_VECTOR_ELT(
                VECTOR_ELT(VECTOR_ELT(ans, j), i)
                , 0
                , Rf_ScalarReal((F[j] + i)->intercept)
            );
            SET_VECTOR_ELT(
                VECTOR_ELT(VECTOR_ELT(ans, j), i)
                , 1
                , Rf_allocVector(REALSXP, (F[j] + i)->num_breaks + 1)
            );
            SET_VECTOR_ELT(
                VECTOR_ELT(VECTOR_ELT(ans, j), i)
                , 2
                , Rf_allocVector(REALSXP, (F[j] + i)->num_breaks)
            );
            memcpy(
                REAL(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(ans, j), i), 1))
                , (F[j] + i)->slope
                , ((F[j] + i)->num_breaks + 1) * sizeof(double)
            );
            memcpy(
                REAL(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(ans, j), i), 2))
                , (F[j] + i)->breakpoint
                , (F[j] + i)->num_breaks * sizeof(double)
            );
        }
    }

out:
    
    for (j = 0; j < num_dims; ++j)
    {
        for (i = 0; i < num_nodes; ++i)
        {
            if (g[j]) plf_free(g[j] + i);
            if (h[j]) plf_free(h[j] + i);
            if (f[j]) plf_free(f[j] + i);
            if (F[j]) plf_free(F[j] + i);
        }
        free(g[j]);
        free(h[j]);
        free(f[j]);
        free(F[j]);
    }
    plf_free(&tmp1);
    plf_free(&tmp2);
    free(g);
    free(h);
    free(f);
    free(F);
    tsx_tree_free(&tree);
    UNPROTECT(3);
    return Rf_list3(
        Rf_ScalarReal(mpr_summary.mean_tree_length), tree_length, ans);
}


SEXP C_treeseq_linear_mpr_minimize(SEXP F)
{
    int i;
    int j;
    int k;
    int l;
    int num_dims = Rf_length(F);
    int n = Rf_length(VECTOR_ELT(F, 0));

    SEXP ans = PROTECT(Rf_allocMatrix(REALSXP, n, num_dims));
    double *restrict x = REAL(ans);
    double *restrict a;
    double *restrict b;
    GetRNGstate();
    for (j = 0; j < num_dims; ++j)
    {
        for (i = 0; i < n; ++i)
        {
            a = REAL(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(F, j), i), 1));
            b = REAL(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(F, j), i), 2));
            k = 0;
            while (a[k+1] < 0)
                ++k;
            x[i+j*n] = b[k];
            l = 1;
            while (a[k+1] == 0)
            {
                ++l;
                ++k;
                if (unif_rand() < (1 / (double)l))
                    x[i+j*n] = b[k];
            }
        }
    }
    PutRNGstate();
    UNPROTECT(1);
    return ans;
}


static double
compute_score(double x, double intercept, double *slope, double *breakpoint,
    int num_breaks)
{
    int j;
    double score = 0;
    if (x < breakpoint[0])
    {
        score = intercept + x*slope[0];
    }
    else
    {
        j = 0;
        score = intercept + slope[0]*breakpoint[0];
        while (j < (num_breaks-1) && x > breakpoint[j+1])
        {
            ++j;
            score += slope[j] * (breakpoint[j] - breakpoint[j-1]);
        }
        score += slope[j+1] * (x - breakpoint[j]);
    }
    return score;
}


SEXP C_treeseq_linear_mpr_minimize_discrete(SEXP F, SEXP sites)
{
    int i;
    int j;
    int k;
    int num_dims = Rf_length(F);
    int n = Rf_length(VECTOR_ELT(F, 0));
    int num_breaks;
    int num_sites = *INTEGER(Rf_getAttrib(sites, R_DimSymbol));
    
    SEXP ans = PROTECT(Rf_allocVector(INTSXP, n));
    
    int *restrict x = INTEGER(ans);
    double *restrict coords = REAL(sites);

    double intercept;
    double *restrict a;
    double *restrict b;

    double score;
    double min_score;

    for (i = 0; i < n; ++i)
    {
        min_score = R_PosInf;
        for (k = 0; k < num_sites; ++k)
        {
            score = 0;
            for (j = 0; j < num_dims; ++j)
            {
                intercept = *REAL(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(F, j), i), 0));
                a = REAL(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(F, j), i), 1));
                b = REAL(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(F, j), i), 2));
                num_breaks = Rf_length(
                    VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(F, j), i), 2));
                score += compute_score(
                    coords[k+j*num_sites], intercept, a, b, num_breaks);
            }
            if (score < min_score)
            {
                min_score = score;
                x[i] = k + 1; // add 1 for R index
            }
        }
    }
    UNPROTECT(1);
    return ans;
}
