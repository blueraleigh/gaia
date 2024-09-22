#include <R.h>
#include <Rmath.h>
#include <assert.h>

#include "graph.h"
#include "fequals.h"


// called by sample_path. sample the previous state on a shortest path 
// from `source` to `target` given the `current`, which is initially set
// to `target`
static int
prev_state(graph_t *g, int current, int source)
{
    int i;
    int k;
    int s = -1;
    double d;
    double D;
    
    int num_states = g->num_states;

    // neighbors that enter into the current state and their edge weights
    int num_neighbors = g->p[current+1] - g->p[current];
    const int *restrict neighbors = g->i + g->p[current];
    const double *restrict weights = g->weights + g->p[current];
    const double *restrict distances = g->distances;

    k = 0;

    // shortest graph distance from source state to current state
    D = distances[source + current*num_states];

    for (i = 0; i < num_neighbors; ++i)
    {
        if (neighbors[i] == source)
            return source;
        // shortest graph distance from source state to
        // current state that goes through this neighbor
        d = distances[source + neighbors[i]*num_states] + weights[i];
        
        if (fequals(d, D))
        {
            ++k;
            if (unif_rand() < (1 / (double)k))
            {
                s = neighbors[i];
            }
        }
    }
    assert (s >= 0);
    return s;
}


// Samples a shortest path from `source` to `target` in reverse iteration
// order. That is, to iterate from `source` to `target` use the following
// idiom:
//   
//   graph_sample_path(graph, source, target, &path_len, path);
//
//   for (i = path_len - 1; i >= 0; --i)
//   {
//        state = path[i];
//   }
//
void
graph_sample_path(graph_t *g, int source, int target, int *ret_path_length, 
    int *ret_path)
{
    int path_length = 0;
    int current = target;
    ret_path[path_length++] = current;
    while (current != source)
    {
        current = prev_state(g, current, source);
        ret_path[path_length++] = current;
    }
    *ret_path_length = path_length;
}
