#ifndef TSX_SANKOFF_H
#define TSX_SANKOFF_H

#include <tskit.h>

typedef struct tsx_tree_t {
    // quintuply linked tree structure
    tsk_id_t virtual_root;
    tsk_id_t *left_child;
    tsk_id_t *right_child;
    tsk_id_t *left_sib;
    tsk_id_t *right_sib;
    tsk_id_t *parent;
    const double *time;
    const tsk_flags_t *node_flags;
    int *num_samples; // Number of samples in the tree
    int num_edges; // Number of edges in the tree
    void *mpr;  // Parsimony data structure
    void (*calc_stem_cost)(tsk_id_t, tsk_id_t, struct tsx_tree_t *); // Stem cost function
    void (*calc_final_cost)(tsk_id_t, struct tsx_tree_t *); // Final cost function
    void (*increment_node_cost)(tsk_id_t, tsk_id_t, struct tsx_tree_t *, int); // Node cost function
} tsx_tree_t;

int
tsx_tree_init(tsx_tree_t *, const tsk_treeseq_t *, int);

void
tsx_tree_free(tsx_tree_t *);

// This function performs the two-pass generalized parsimony algorithm.
// TREE_FUN is a pointer to a function that is called once on each tree
// after the initial downpass phase of the algorithm completes when the 
// node cost and stem cost functions have been constructed. NODE_FUN is a
// pointer to a function that is called once for each node in a tree during
// the uppass phase of the alorithm after the final cost function for the
// node has been constructed.
void
tsx_treeseq_sankoff(
    const tsk_treeseq_t *ts,
    tsx_tree_t *tree,
    void *tree_params,
    void (*TREE_FUN)(tsx_tree_t *, int, double, double, void *),
    void *node_params,
    void (*NODE_FUN)(tsx_tree_t *, int, double, double, tsk_id_t, void *)
);

#endif
