#include <R.h>
#include <Rinternals.h>
#include <tskit.h>

#include "error.h"

/* Test if a migration history contains a transition from state i to
** state j during the time interval [interval_start, interval_end).
** Each migration history comprises two vectors giving the visited states
** and the times those states were entered into. A migration history
** always begins and ends with the two endpoint states and times of the lineage.
** For example, if a lineage that began at time 10 and ended at time 0 
** remained in state 1 for its entire duration, its migration history would 
** look like
**
** path_state: [  1   1 ]
** path_time:  [  0  10 ]
**
** On the hand, if it ended in state 2 and visited states 3, 4, and 5 it would
** look like
**
** path_state: [  2   2   5   4   3   1 ]
** path_time:  [  0   2   4   6   8  10 ]
*/
static int
path_contains(
    int path_length,
    int *path_state,
    double *path_time,
    int i,
    int j,
    double interval_start,
    double interval_end)
{
    int k;
    int to;
    int from;
    double et;
    to = path_state[1];
    for (k = 2; k < path_length; ++k)
    {
        et = path_time[k-1];
        from = path_state[k];
        if (from == i && to == j)
        {
            if (et >= interval_start && et < interval_end)
                return 1;
        }
        to = from;
    }
    return 0;
}


/* Count the number of samples in each set descended from u and multiply that
** number by the ancestral segment length. Note that in the calling context, 
** the node u is always the root of a local tree. These counts are then
** distributed into each of the time bins spanned by the lineages on the path
** from a sample to the root.
*/
static void
count_sample_segments_descended_from_node(
    tsk_id_t u,
    int num_times,
    int num_sample_sets,
    double segment_length,
    double *restrict times,
    double *restrict node_time,
    int *restrict sample_set,
    tsk_id_t *restrict left_sample,
    tsk_id_t *restrict right_sample,
    tsk_id_t *restrict next_sample,
    double *restrict counts)
{
    int k;
    int it;
    int jt;
    int mflag;
    int sample_set_id;
    tsk_id_t sample_id;
    it = findInterval2(times, num_times, node_time[u], 0, 0, 0, 1, &mflag);
    if (mflag == -1) return;
    if (mflag == 1) --it;
    --it;
    sample_id = left_sample[u];
    while (sample_id != TSK_NULL)
    {
        sample_set_id = sample_set[sample_id];
        if (sample_set_id >= 0)
        {
            jt = findInterval2(
                times, num_times, node_time[sample_id], 0, 0, 0, 1, &mflag);
            if (mflag != 1)
            {
                if (jt > 0) --jt;
                for (k = jt; k <= it; ++k)
                {
                    counts[sample_set_id + k * num_sample_sets] += 
                        segment_length;
                }
            }
        }
        if (sample_id == right_sample[u])
            break;
        sample_id = next_sample[sample_id];
    }
}


/* Count the number of samples in each set descended from each migration on
** the edge (u, v) and multiply that number by the ancestral segment length.
** Distribute those counts into the time bins during which the migrations
** occurred.
*/
static void
count_sample_segments_descended_from_edge_migrations(
    tsk_id_t u,
    tsk_id_t v,
    int num_times,
    int num_state_sets,
    int num_sample_sets,
    double segment_length,
    int path_length,
    int *restrict path,
    double *restrict path_time,
    int *restrict path_offset,
    int *restrict path_states,
    double *restrict path_times,
    double *restrict times,
    double *restrict node_time,
    int *restrict state_set,
    int *restrict sample_set,
    tsk_id_t *restrict edge,
    tsk_id_t *restrict parent,
    tsk_id_t *restrict left_sample,
    tsk_id_t *restrict right_sample,
    tsk_id_t *restrict next_sample,
    double *restrict counts,
    double *restrict max_age)
{
    int i;
    int j;
    int k;
    int t;
    int to;
    int from;
    int mflag;
    int visited;
    int m_index;
    int s_index;
    int t_index;
    int sample_set_id;
    tsk_id_t sample_id;
    int num_state_sets_squared = num_state_sets * num_state_sets;
    int offset = num_state_sets_squared * num_sample_sets;
    double event_time;
    t = 1;
    j = path[1];
    for (k = 2; k < path_length; ++k) // iterating into the past
    {
        i = path[k];
        event_time = path_time[k-1];
        t = findInterval2(
            times, num_times, event_time, 0, 0, 0, t, &mflag);
        if (mflag == 1) break;
        if (mflag == 0)
        {
            visited = 0;
            while (parent[u] != TSK_NULL && node_time[u] < times[t])
            {
                if (path_contains(
                    path_offset[edge[u]+1] - path_offset[edge[u]],
                    path_states + path_offset[edge[u]],
                    path_times + path_offset[edge[u]],
                    i,
                    j,
                    times[t-1],
                    times[t]))
                {
                    visited = 1;
                    break;
                }
                u = parent[u];
            }
            // if we haven't already seen an i to j transition
            // from a parental edge in the same time bin
            if (!visited)
            {
                from = state_set[i];
                to = state_set[j];
                m_index = from + to * num_state_sets;
                t_index = (t-1) * offset;
                if (event_time > max_age[m_index])
                    max_age[m_index] = event_time;
                sample_id = left_sample[v];
                while (sample_id != TSK_NULL)
                {
                    sample_set_id = sample_set[sample_id];
                    if (sample_set_id >= 0)
                    {
                        s_index = num_state_sets_squared * sample_set_id;
                        counts[t_index + s_index + m_index] += segment_length;                            
                    }
                    if (sample_id == right_sample[v])
                        break;
                    sample_id = next_sample[sample_id];
                }
            }
        }
        j = i;
    }
}


SEXP C_treeseq_discrete_mpr_ancestry_flux(
    SEXP tr,
    SEXP r_path_offsets,
    SEXP r_path_states,
    SEXP r_path_times,
    SEXP r_num_state_sets,
    SEXP r_state_sets,
    SEXP r_num_sample_sets,
    SEXP r_sample_sets,
    SEXP r_times)
{
    tsk_tree_t *tree = (tsk_tree_t *) R_ExternalPtrAddr(tr);

    if (!tsk_tree_has_sample_lists(tree))
    {
        TSX_ERROR("sample list tracking not enabled for this tree sequence");
        return R_NilValue;
    }

    const tsk_treeseq_t *ts = tree->tree_sequence;

    tsk_id_t *restrict edge = tree->edge;
    tsk_id_t *restrict parent = tree->parent;
    tsk_id_t *restrict right_child = tree->right_child;
    tsk_id_t *restrict left_sib = tree->left_sib;
    tsk_id_t *restrict left_sample = tree->left_sample;
    tsk_id_t *restrict right_sample = tree->right_sample;
    tsk_id_t *restrict next_sample = tree->next_sample;
    double *restrict node_time = ts->tables->nodes.time;
    
    int num_times = Rf_length(r_times);
    int num_time_bins = num_times - 1;
    int num_state_sets = *INTEGER(r_num_state_sets);
    int num_state_sets_squared = num_state_sets * num_state_sets;
    int num_sample_sets = *INTEGER(r_num_sample_sets);
    int *restrict state_set = INTEGER(r_state_sets);
    int *restrict sample_set = INTEGER(r_sample_sets);
    int *restrict path_offset = INTEGER(r_path_offsets);
    int *restrict path_states = INTEGER(r_path_states);
    double *restrict times = REAL(r_times);
    double *restrict path_times = REAL(r_path_times);
    
    SEXP dims = PROTECT(Rf_allocVector(INTSXP, 4));
    SET_INTEGER_ELT(dims, 0, num_state_sets);
    SET_INTEGER_ELT(dims, 1, num_state_sets);
    SET_INTEGER_ELT(dims, 2, num_sample_sets);
    SET_INTEGER_ELT(dims, 3, num_time_bins);
    SEXP flux = PROTECT(Rf_allocArray(REALSXP, dims));
    SEXP flux_scale = PROTECT(
        Rf_allocMatrix(REALSXP, num_sample_sets, num_time_bins));
    SEXP flux_max_age = PROTECT(
        Rf_allocMatrix(REALSXP, num_state_sets, num_state_sets));

    Rf_setAttrib(flux, Rf_install("scale"), flux_scale);
    Rf_setAttrib(flux, Rf_install("max.age"), flux_max_age);

    double *numer = REAL(flux);
    double *denom = REAL(flux_scale);
    double *max_age = REAL(flux_max_age);
    
    Memzero(denom, num_sample_sets * num_time_bins);
    Memzero(numer, num_state_sets_squared * num_sample_sets * num_time_bins);
    Memzero(max_age, num_state_sets_squared);

    int rc;
    int stack_top;
    
    int num_nodes = tsk_treeseq_get_num_nodes(ts);
    int virtual_root = num_nodes;
    
    int *path;
    int path_length;
    double *path_time;

    tsk_id_t *restrict stack = (tsk_id_t *) R_alloc(
        num_nodes + 1, sizeof(*stack));

    tsk_id_t h;
    tsk_id_t u;
    tsk_id_t v;

    double segment_length;
    double sequence_length = ts->tables->sequence_length;

    for (
        rc = tsk_tree_first(tree);
        rc == TSK_TREE_OK;
        rc = tsk_tree_next(tree))
    {
        segment_length = tree->interval.right - tree->interval.left;
        segment_length /= sequence_length;
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
            
            h = edge[u];
            
            // true for roots (including isolated sample nodes)
            if (h == TSK_NULL)
            {
                count_sample_segments_descended_from_node(
                    u,
                    num_times,
                    num_sample_sets,
                    segment_length,
                    times,
                    node_time,
                    sample_set,
                    left_sample,
                    right_sample,
                    next_sample,
                    denom
                );
                continue;
            }
            
            v = u;
            u = parent[v];

            if (node_time[u] < times[0] || node_time[v] > times[num_time_bins])
                continue;

            path = path_states + path_offset[h];
            path_length = path_offset[h+1] - path_offset[h];
            path_time = path_times + path_offset[h];

            count_sample_segments_descended_from_edge_migrations(
                u,
                v,
                num_times,
                num_state_sets,
                num_sample_sets,
                segment_length,
                path_length,
                path,
                path_time,
                path_offset,
                path_states,
                path_times,
                times,
                node_time,
                state_set,
                sample_set,
                edge,
                parent,
                left_sample,
                right_sample,
                next_sample,
                numer,
                max_age
            );
        }
    }
    for (int i = 0; i < num_time_bins; ++i)
    {
        for (int j = 0; j < num_sample_sets; ++j)
        {
            if (denom[j] > 0)
            {
                for (int k = 0; k < num_state_sets_squared; ++k)
                    numer[k] /= denom[j];
            }
            numer += num_state_sets_squared;
        }
        denom += num_sample_sets;
    }
    UNPROTECT(4);
    return flux;
}
