#ifndef TSX_GRAPH_H
#define TSX_GRAPH_H


// Note that if there is an edge connecting states (i, j) we require
// there to be an edge connecting states (j, i), although the edge
// weights may differ. this is verified from the calling function in R
typedef struct graph {
    int num_states;
    // shortest path distances between all pairs of states
    double *distances;
    // adjacency matrix in compressed sparse column format:
    // p[j] gives the index of the weight that begins column j.
    // i.e., defining k := p[j], A[i[k], j] = weights[k]
    int *p;
    // row index of each non-zero weight
    int *i;
    // non-zero weights
    double *weights;
} graph_t;


void
graph_sample_path(graph_t *, int, int, int *, int *);

#endif
