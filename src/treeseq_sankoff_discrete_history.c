#include <R.h>
#include <Rinternals.h>
#include <tskit.h>

#include "graph.h"


/* Each migration history comprises two vectors giving the visited states
** and the times those states were entered into. A migration history always
** begins and ends with the two endpoint states and times of the lineage.
** For example, if a lineage that began at time 10 and ended at time 0 
** remained in state 1 for its entire duration, its migration history would 
** look like
**
** path_state: [  1   1 ]
** path_time:  [  0  10 ]
**
** On the hand, if it ended in state 2 and visited states 3, 4, and 5 to get
** there it would look like
**
** path_state: [  2   2   5   4   3   1 ]
** path_time:  [  0   2   4   6   8  10 ]
**
** We distribute the lineage duration equally among the visited states in a
** migration history.
*/
SEXP C_treeseq_discrete_mpr_edge_history(
    SEXP treeseq,
    SEXP node_states,
    SEXP cost_matrix,
    SEXP adjacency_matrix,
    SEXP index1)
{
    graph_t graph;

    /* indx = 0: node states are indexed from 0 */
    /* indx = 1: node states are indexed from 1 */
    int indx1 = *INTEGER(index1);
    int num_states = *INTEGER(Rf_getAttrib(cost_matrix, R_DimSymbol));
    const int *restrict node_state = INTEGER(node_states);

    graph.num_states = num_states;
    graph.distances = REAL(cost_matrix);
    graph.i = INTEGER(R_do_slot(adjacency_matrix, Rf_mkString("i")));
    graph.p = INTEGER(R_do_slot(adjacency_matrix, Rf_mkString("p")));
    graph.weights = REAL(R_do_slot(adjacency_matrix, Rf_mkString("x")));

    tsk_treeseq_t *ts = (tsk_treeseq_t *)R_ExternalPtrAddr(treeseq);

    int num_edges = ts->tables->edges.num_rows;

    int path_length;
    int *path = (int *) R_alloc(num_states, sizeof(*path));

    const tsk_id_t *restrict edge_parent = ts->tables->edges.parent;
    const tsk_id_t *restrict edge_child = ts->tables->edges.child;
    const double *time = ts->tables->nodes.time;

    SEXP paths = PROTECT(Rf_allocVector(VECSXP, num_edges));
    SEXP times = PROTECT(Rf_allocVector(VECSXP, num_edges));
    SEXP path_lengths = PROTECT(Rf_allocVector(INTSXP, num_edges));

    int i;
    int j;
    int from;
    int to;
    double event_time;
    double waiting_time;
    tsk_id_t u;
    tsk_id_t v;
    GetRNGstate();
    for (i = 0; i < num_edges; ++i)
    {
        u = edge_parent[i];
        v = edge_child[i];
        from = node_state[u] - indx1;
        to = node_state[v] - indx1;
        if (from != to)
        {
            graph_sample_path(&graph, from, to, &path_length, path);
        }
        else
        {
            path_length = 1;
            path[0] = from;
        }
        waiting_time = (time[u] - time[v]) / path_length;
        event_time = time[v] + waiting_time;
        SET_VECTOR_ELT(paths, i, Rf_allocVector(INTSXP, path_length + 1));
        SET_VECTOR_ELT(times, i, Rf_allocVector(REALSXP, path_length + 1));
        SET_INTEGER_ELT(path_lengths, i, path_length + 1);
        SET_INTEGER_ELT(VECTOR_ELT(paths, i), 0, to + indx1);
        SET_REAL_ELT(VECTOR_ELT(times, i), 0, time[v]);
        for (j = 0; j < path_length; ++j)
        {
            // ordered from younger to older (from child node to parent node)
            SET_INTEGER_ELT(VECTOR_ELT(paths, i), j + 1, path[j] + indx1);
            SET_REAL_ELT(VECTOR_ELT(times, i), j + 1, event_time);
            event_time += waiting_time;
        }
    }
    PutRNGstate();
    
    UNPROTECT(3);
    return Rf_list3(paths, times, path_lengths);
}
