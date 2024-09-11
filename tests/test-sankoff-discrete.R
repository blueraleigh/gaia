// This file compares the output of gaia to a more basic implementation
// programmed directly in R.

library(gaia)

stem.cost = function(node, state, parent, time, node_cost, stem_cost, cost) {
    u = node + 1
    if (parent[u] == -1) return (stem_cost)
    b = if (!is.null(time)) (time[parent[u]+1] - time[u]) else 1
    num_states = nrow(cost)
    for (i in 1:num_states) {
        min_cost = Inf
        for (j in 1:num_states) {
            tmp = (1/b)*cost[i,j] + node_cost[u,j]
            if (tmp < min_cost)
                min_cost = tmp
        }
        stem_cost[u,i] = min_cost
    }
    return (stem_cost)
}

node.cost = function(node, state, parent, time, node_cost, stem_cost, cost) {
    u = node + 1
    children = which(parent == node)
    num_children = length(children)
    if (num_children == 0) return (node_cost)
    node_cost[u,] = colSums(stem_cost[children,,drop=FALSE])
    return (node_cost)
}

final.cost = function(node, state, parent, time, node_cost, stem_cost,
    final_cost, cost)
{
    u = node + 1
    if (parent[u] == -1) {
        final_cost[u,] = node_cost[u,]
        return (final_cost)
    }
    v = u
    u = parent[v] + 1
    b = if (!is.null(time)) (time[u] - time[v]) else 1
    num_states = nrow(cost)
    for (i in 1:num_states) {
        min_cost = Inf
        for (j in 1:num_states) {
            tmp = final_cost[u,j] - stem_cost[v,j] + (1/b)*cost[j,i] + node_cost[v,i]
            if (tmp < min_cost)
                min_cost = tmp
        }
        final_cost[v, i] = min_cost
    }
    return (final_cost)
}

mpr = function(state, parent, time, postorder, num_nodes, cost) {
    node_cost = matrix(0, num_nodes, nrow(cost))
    stem_cost = matrix(0, num_nodes, nrow(cost))
    final_cost = matrix(0, num_nodes, nrow(cost))
    node_cost[1:3, ] = Inf
    node_cost[cbind(1:3, state)] = 0
    for (node in postorder)
    {
        node_cost = node.cost(node,state,parent,time,node_cost,stem_cost,cost)
        stem_cost = stem.cost(node,state,parent,time,node_cost,stem_cost,cost)
    }
    for (node in rev(postorder))
    {
        final_cost = final.cost(
            node,state,parent,time,node_cost,stem_cost,final_cost,cost)
    }
    return (structure(final_cost, node_cost=node_cost, stem_cost=stem_cost))
}


# Example tree sequence:
#
#      [0, 20)                     [20, 80)                  [80, 100)
#
#    +----6----+
#    |         |                                           +----5----+
#    |         |                                           |         |
#    |     +---4---+              +---4---+                |     +---4---+
#    |     |       |              |       |                |     |       |
#    |     |       |           +--3--+    |                |     |       |  
#    |     |       |           |     |    |                |     |       |
#    0     1       2           0     2    1                0     1       2
#
#
# node_id    time        state
#  0          0          3
#  1          0          1
#  2          0          1
#  3          0.15       -
#  4          0.60       -
#  5          0.80       -
#  6          1.00       -
#

ts = treeseq_load("./test.trees")

cost = matrix(0, 3, 3)
cost[1,2] = cost[2,1] = cost[2,3] = cost[3,2] = 1
cost[1,3] = cost[3,1] = 2

t1 = c(6,4,4,-1,6,-1,-1)
t2 = c(3,4,3,4,-1,-1,-1)
t3 = c(5,4,4,-1,5,-1,-1)
time = c(0,0,0,0.15,0.6,0.8,1.0)
state = c(3L,1L,1L)

S = treeseq_discrete_mpr(ts, cbind(node_id=0:2, state_id=state), cost, TRUE)$mpr

s1 = mpr(state, t1, time, c(1,2,4,0,6), 7, cost)
s2 = mpr(state, t2, time, c(0,2,3,1,4), 7, cost)
s3 = mpr(state, t3, time, c(1,2,4,0,5), 7, cost)
s4 = .2*s1[5,] + .6*s2[5,] + .2*s3[5,]


stopifnot(all.equal(S[, 4], s2[4, ]))
stopifnot(all.equal(S[, 5], s4))
stopifnot(all.equal(S[, 6], s3[6, ]))
stopifnot(all.equal(S[, 7], s1[7, ]))


S = treeseq_discrete_mpr(ts, cbind(node_id=0:2, state_id=state), cost, FALSE)$mpr

s1 = mpr(state, t1, NULL, c(1,2,4,0,6), 7, cost)
s2 = mpr(state, t2, NULL, c(0,2,3,1,4), 7, cost)
s3 = mpr(state, t3, NULL, c(1,2,4,0,5), 7, cost)
s4 = .2*s1[5,] + .6*s2[5,] + .2*s3[5,]


stopifnot(all.equal(S[, 4], s2[4, ]))
stopifnot(all.equal(S[, 5], s4))
stopifnot(all.equal(S[, 6], s3[6, ]))
stopifnot(all.equal(S[, 7], s1[7, ]))

