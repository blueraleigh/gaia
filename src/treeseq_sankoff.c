#include <R.h>
#include <Rinternals.h>

#include <assert.h>

#include "error.h"
#include "treeseq_sankoff.h"


//
//   (toward root) +u+----------+v+ (toward tips)
//
static void
update_node_and_stem_cost(tsk_id_t u, tsk_id_t v, tsx_tree_t *tree, int sign)
{
    if (tree->parent[u] != TSK_NULL)
    {
        // subtract u's stem cost from its parent's node cost.
        // in the next call to this function u replaces v and
        // we will update the node cost of u's parent with u's
        // new stem cost.
        tree->increment_node_cost(tree->parent[u], u, tree, -1);
    }
    if (sign == +1)
        tree->calc_stem_cost(u, v, tree);
    tree->increment_node_cost(u, v, tree, sign);
}


static void
remove_branch(tsk_id_t p, tsk_id_t c, tsx_tree_t *tree)
{
    tsk_id_t *restrict left_child = tree->left_child;
    tsk_id_t *restrict right_child = tree->right_child;
    tsk_id_t *restrict left_sib = tree->left_sib;
    tsk_id_t *restrict right_sib = tree->right_sib;
    tsk_id_t *restrict parent = tree->parent;
    tsk_id_t lsib = left_sib[c];
    tsk_id_t rsib = right_sib[c];
    if (lsib == TSK_NULL)
        left_child[p] = rsib;
    else
        right_sib[lsib] = rsib;
    if (rsib == TSK_NULL)
        right_child[p] = lsib;
    else
        left_sib[rsib] = lsib;
    parent[c] = TSK_NULL;
    left_sib[c] = TSK_NULL;
    right_sib[c] = TSK_NULL;
}


static void
insert_branch(tsk_id_t p, tsk_id_t c, tsx_tree_t *tree)
{
    tsk_id_t *restrict left_child = tree->left_child;
    tsk_id_t *restrict right_child = tree->right_child;
    tsk_id_t *restrict left_sib = tree->left_sib;
    tsk_id_t *restrict right_sib = tree->right_sib;
    tsk_id_t *restrict parent = tree->parent;
    parent[c] = p;
    tsk_id_t u = right_child[p];
    if (u == TSK_NULL)
    {
        left_child[p] = c;
        left_sib[c] = TSK_NULL;
        right_sib[c] = TSK_NULL;
    }
    else
    {
        right_sib[u] = c;
        left_sib[c] = u;
        right_sib[c] = TSK_NULL;
    }
    right_child[p] = c;
}


static void
remove_root(tsk_id_t root, tsx_tree_t *tree)
{
    remove_branch(tree->virtual_root, root, tree);
}


static void
insert_root(tsk_id_t root, tsx_tree_t *tree)
{
    insert_branch(tree->virtual_root, root, tree);
    tree->parent[root] = TSK_NULL;
}


static void
remove_edge(tsk_id_t p, tsk_id_t c, tsx_tree_t *tree)
{
    assert(p != TSK_NULL);
    assert(c != TSK_NULL);
    assert(p != tree->virtual_root);
    assert(tree->parent[c] == p);
    remove_branch(p, c, tree);
    tree->num_edges -= 1;
    int path_end_was_root;
    tsk_id_t path_end, u = p, v = c;
    tsk_id_t *restrict parent = tree->parent;
    int *restrict num_samples = tree->num_samples;
    update_node_and_stem_cost(u, v, tree, -1);
    path_end = u;
    path_end_was_root = num_samples[u] > 0 ? 1 : 0;
    num_samples[u] -= num_samples[c];
    v = u;
    u = parent[u];
    while (u != TSK_NULL)
    {
        update_node_and_stem_cost(u, v, tree, +1);
        path_end = u;
        path_end_was_root = num_samples[u] > 0 ? 1 : 0;
        num_samples[u] -= num_samples[c];
        v = u;
        u = parent[u];
    }
    if (path_end_was_root && (num_samples[path_end] == 0))
        remove_root(path_end, tree);
    if (num_samples[c] > 0)
        insert_root(c, tree);
}


static void
insert_edge(tsk_id_t p, tsk_id_t c, tsx_tree_t *tree)
{
    assert(p != TSK_NULL);
    assert(c != TSK_NULL);
    assert(tree->parent[c] == TSK_NULL);
    int path_end_was_root;
    tsk_id_t path_end, u = p, v = c;
    tsk_id_t *restrict parent = tree->parent;
    int *restrict num_samples = tree->num_samples;
    update_node_and_stem_cost(u, v, tree, +1);
    path_end = u;
    path_end_was_root = num_samples[u] > 0 ? 1 : 0;
    num_samples[u] += num_samples[c];
    v = u;
    u = parent[u];
    while (u != TSK_NULL)
    {
        update_node_and_stem_cost(u, v, tree, +1);
        path_end = u;
        path_end_was_root = num_samples[u] > 0 ? 1 : 0;
        num_samples[u] += num_samples[c];
        v = u;
        u = parent[u];
    }
    if (num_samples[c] > 0)
        remove_root(c, tree);
    if ((num_samples[path_end] > 0) && !path_end_was_root)
        insert_root(path_end, tree);
    insert_branch(p, c, tree);
    tree->num_edges += 1;
    assert(tree->parent[c] == p);
}


void
tsx_treeseq_sankoff(
    const tsk_treeseq_t *ts,
    tsx_tree_t *tree,
    void *tree_params,
    void (*TREE_FUN)(tsx_tree_t *, int, double, double, void *),
    void *node_params,
    void (*NODE_FUN)(tsx_tree_t *, int, double, double, tsk_id_t, void *))
{
    tsk_id_t u;
    tsk_id_t v;
    tsk_id_t h;
    tsk_id_t tj;
    tsk_id_t tk;
    int t_index;
    int stack_top;
    double t_left;
    double t_right;

    const tsk_size_t num_nodes = ts->tables->nodes.num_rows;
    const tsk_id_t num_edges = (tsk_id_t) ts->tables->edges.num_rows;
    const double sequence_length = ts->tables->sequence_length;
    const tsk_id_t *restrict I = ts->tables->indexes.edge_insertion_order;
    const tsk_id_t *restrict O = ts->tables->indexes.edge_removal_order;
    const double *restrict edge_left = ts->tables->edges.left;
    const double *restrict edge_right = ts->tables->edges.right;
    const tsk_id_t *restrict edge_parent = ts->tables->edges.parent;
    const tsk_id_t *restrict edge_child = ts->tables->edges.child;

    tsk_id_t virtual_root = tree->virtual_root;
    tsk_id_t *restrict right_child = tree->right_child;
    tsk_id_t *restrict left_sib = tree->left_sib;

    tsk_id_t *restrict stack = malloc((num_nodes + 1) * sizeof(*stack));

    if (!stack)
    {
        TSX_WARN("memory allocation failed");
        goto out;
    }

    // iterate over the trees
    tj = 0;
    tk = 0;
    t_index = 0;
    t_left = 0;
    while (tj < num_edges || t_left < sequence_length)
    {
        while (tk < num_edges && edge_right[O[tk]] == t_left)
        {
            h = O[tk];
            tk++;
            v = edge_child[h];
            u = edge_parent[h];
            remove_edge(u, v, tree);
        }
        
        while (tj < num_edges && edge_left[I[tj]] == t_left)
        {
            h = I[tj];
            tj++;
            v = edge_child[h];
            u = edge_parent[h];
            insert_edge(u, v, tree);
        }

        t_right = sequence_length;
        
        if (tj < num_edges)
            t_right = TSK_MIN(t_right, edge_left[I[tj]]);

        if (tk < num_edges)
            t_right = TSK_MIN(t_right, edge_right[O[tk]]);

        if (TREE_FUN)
            TREE_FUN(tree, t_index, t_left, t_right, tree_params);

        if (NODE_FUN)
        {
            stack_top = -1;
            for (u = right_child[virtual_root]; u != TSK_NULL; u = left_sib[u])
            {
                stack_top++;
                stack[stack_top] = u;
            }
            while (stack_top >= 0)
            {
                u = stack[stack_top];
                stack_top--;
                for (v = right_child[u]; v != TSK_NULL; v = left_sib[v])
                {
                    stack_top++;
                    stack[stack_top] = v;
                }
                tree->calc_final_cost(u, tree);
                NODE_FUN(tree, t_index, t_left, t_right, u, node_params);
            }
        }

        t_left = t_right;
        ++t_index;
    }
out:
    free(stack);
}


int
tsx_tree_init(tsx_tree_t *tree, const tsk_treeseq_t *ts, int time)
{
    tsk_id_t u;
    tsk_size_t j;
    tsk_size_t num_nodes = ts->tables->nodes.num_rows;
    tsk_size_t num_nodes1 = num_nodes + 1;
    tsk_id_t *parent = malloc(num_nodes1 * sizeof(*parent));
    tsk_id_t *left_child = malloc(num_nodes1 * sizeof(*left_child));
    tsk_id_t *right_child = malloc(num_nodes1 * sizeof(*right_child));
    tsk_id_t *left_sib = malloc(num_nodes1 * sizeof(*left_sib));
    tsk_id_t *right_sib = malloc(num_nodes1 * sizeof(*right_sib));
    int *num_samples = calloc(num_nodes1, sizeof(*num_samples));
    if (   parent == NULL 
        || left_child == NULL
        || right_child == NULL
        || left_sib == NULL
        || right_sib == NULL 
        || num_samples == NULL)
    {
        TSX_WARN("memory allocation failed");
        goto out;
    }
    memset(left_child, 0xff, num_nodes1 * sizeof(*left_child));
    memset(right_child, 0xff, num_nodes1 * sizeof(*right_child));
    memset(left_sib, 0xff, num_nodes1 * sizeof(*left_sib));
    memset(right_sib, 0xff, num_nodes1 * sizeof(*right_sib));
    memset(parent, 0xff, num_nodes1 * sizeof(*parent));
    tree->virtual_root = (tsk_id_t)num_nodes;
    tree->left_child = left_child;
    tree->right_child = right_child;
    tree->left_sib = left_sib;
    tree->right_sib = right_sib;
    tree->parent = parent;
    tree->time = time ? ts->tables->nodes.time : NULL;
    tree->node_flags = ts->tables->nodes.flags;
    assert (tree->node_flags);
    tree->num_samples = num_samples;
    tree->num_edges = 0;
    tree->mpr = NULL;
    tree->calc_stem_cost = NULL;
    tree->calc_final_cost = NULL;
    tree->increment_node_cost = NULL;
    for (j = 0; j < ts->num_samples; ++j)
    {
        u = ts->samples[j];
        num_samples[u] = 1;
        insert_root(u, tree);
    }
    return 0;
out:
    free(parent);
    free(left_child);
    free(right_child);
    free(left_sib);
    free(right_sib);
    free(num_samples);
    return 1;
}


void
tsx_tree_free(tsx_tree_t *tree)
{
    free(tree->left_child);
    free(tree->right_child);
    free(tree->left_sib);
    free(tree->right_sib);
    free(tree->parent);
    free(tree->num_samples);
}
