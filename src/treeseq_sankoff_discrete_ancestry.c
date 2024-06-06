#include <R.h>
#include <Rinternals.h>
#include <assert.h>
#include <tskit.h>

#include "error.h"


SEXP C_treeseq_discrete_mpr_ancestry(
    SEXP tr,
    SEXP r_path_offsets,
    SEXP r_path_states,
    SEXP r_path_times,
    SEXP r_node_states,
    SEXP r_num_state_sets,
    SEXP r_state_sets,
    SEXP r_num_sample_sets,
    SEXP r_sample_sets,
    SEXP r_times)
{
    tsk_tree_t *tree = (tsk_tree_t *)R_ExternalPtrAddr(tr);
    const tsk_treeseq_t *ts = tree->tree_sequence;
    int num_nodes = tsk_treeseq_get_num_nodes(ts);
    int virtual_root = num_nodes;

    if (!tsk_tree_has_sample_lists(tree))
    {
        TSX_ERROR("sample list tracking not enabled for this tree sequence");
        return R_NilValue;
    }

    tsk_id_t *restrict parent = tree->parent;
    tsk_id_t *restrict right_child = tree->right_child;
    tsk_id_t *restrict left_sib = tree->left_sib;
    tsk_id_t *restrict left_sample = tree->left_sample;
    tsk_id_t *restrict right_sample = tree->right_sample;
    tsk_id_t *restrict next_sample = tree->next_sample;
    double *restrict node_time = ts->tables->nodes.time;

    tsk_id_t u;
    tsk_id_t v;
    tsk_id_t sample_id;

    int num_times = Rf_length(r_times);
    int num_state_sets = *INTEGER(r_num_state_sets);
    int num_sample_sets = *INTEGER(r_num_sample_sets);
    int *restrict state_set = INTEGER(r_state_sets);
    int *restrict sample_set = INTEGER(r_sample_sets);
    int *restrict node_state = INTEGER(r_node_states);
    int *restrict path_offset = INTEGER(r_path_offsets);
    int *restrict path_states = INTEGER(r_path_states);
    double *restrict times = REAL(r_times);
    double *restrict path_times = REAL(r_path_times);

    SEXP dims = PROTECT(Rf_allocVector(INTSXP, 3));
    INTEGER(dims)[0] = num_state_sets;
    INTEGER(dims)[1] = num_sample_sets;
    INTEGER(dims)[2] = num_times;
    SEXP ans = PROTECT(Rf_allocArray(REALSXP, dims));
    SEXP ans_attr = PROTECT(
        Rf_allocMatrix(REALSXP, num_sample_sets, num_times));

    Rf_setAttrib(ans, Rf_install("scale"), ans_attr);

    int h;
    int k;
    int j;
    int t;
    int stack_top;
    int rc;
    int it;
    int jt;
    int mflag;
    int state_id;
    int state_set_id;
    int sample_set_id;
    int path_length;
    int *path;
    double *path_time;
    int offset = num_state_sets * num_sample_sets;

    double segment_length;
    double sequence_length = ts->tables->sequence_length;

    double denom;
    double *numer;
    double *ret = REAL(ans);
    Memzero(ret, offset*num_times);

    tsk_flags_t *node_flags = ts->tables->nodes.flags;

    tsk_id_t *restrict stack = (tsk_id_t *) R_alloc(
        num_nodes + 1, sizeof(*stack));

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
            
            h = tree->edge[u];
            
            if (h == TSK_NULL)
            {
                if (node_flags && (node_flags[u] & TSK_NODE_IS_SAMPLE))
                {
                    jt = findInterval2(
                        times, num_times, node_time[u], 0, 0, 0, 1, &mflag);
                    if (mflag == 0)
                    {
                        if (jt > 0) --jt;
                        numer = ret + offset*jt;
                        state_set_id = state_set[node_state[u]];
                        sample_set_id = sample_set[u];
                        if (sample_set_id >= 0)
                        {
                            numer[state_set_id+sample_set_id*num_state_sets] += 
                                segment_length;
                        }
                    }
                }
                continue;
            }
            
            v = u;
            u = parent[v];

            if (node_time[u] < times[0] || node_time[v] > times[num_times-1])
                continue;

            jt = findInterval2(
                times, num_times, node_time[v], 0, 0, 0, 1, &mflag);

            if (mflag == 1) continue;

            it = findInterval2(
                times, num_times, node_time[u], 0, 0, 0, jt, &mflag);

            if (mflag == -1) continue;
            assert (it >= jt);
            --it;
            if (jt > 0) --jt;

            path = path_states + path_offset[h];
            path_length = path_offset[h+1] - path_offset[h];
            path_time = path_times + path_offset[h];

            t = 1;
            for (k = jt; k <= it; ++k)
            {
                t = findInterval2(
                    path_time, path_length, times[k], 0, 0, 0, t, &mflag);
                if (mflag == 0)
                {
                    numer = ret + offset*k;
                    state_id = path[t];
                    state_set_id = state_set[state_id];
                    sample_id = left_sample[v];
                    while (sample_id != TSK_NULL)
                    {
                        sample_set_id = sample_set[sample_id];
                        if (sample_set_id >= 0)
                        {
                            numer[state_set_id+sample_set_id*num_state_sets] += 
                                segment_length;
                        }
                        if (sample_id == right_sample[v])
                            break;
                        sample_id = next_sample[sample_id];
                    }
                }
            }
        }
    }

    numer = ret;
    for (k = 0; k < num_times; ++k)
    {
        for (h = 0; h < num_sample_sets; ++h)
        {
            denom = 0;
            for (j = 0; j < num_state_sets; ++j)
                denom += numer[j + h * num_state_sets];
            SET_REAL_ELT(ans_attr, h + k * num_sample_sets, denom);
            if (denom > 0)
            {
                for (j = 0; j < num_state_sets; ++j)
                    numer[j + h * num_state_sets] /= denom;
            }
        }
        numer += offset;
    }
    UNPROTECT(3);
    return ans;
}
