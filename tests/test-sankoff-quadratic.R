library(gaia)

stem.cost = function(node, x, y, parent, time, node_cost, stem_cost) {
    u = node + 1
    if (parent[u] == -1) return (stem_cost)
    b = if (!is.null(time)) time[parent[u]+1] - time[u] else 1
    num_children = length(which(parent == node))
    if (num_children > 0)
    {
        p1 = node_cost[u, 1] / (b*node_cost[u, 1] + 1)
        p2 = node_cost[u, 2] / (b*node_cost[u, 1] + 1)
        p3 = node_cost[u, 3] / (b*node_cost[u, 1] + 1)
        p4 = node_cost[u, 4] - 
            ((node_cost[u, 2])^2 / (4*(node_cost[u, 1]+ 1/b))) -
            ((node_cost[u, 3])^2 / (4*(node_cost[u, 1]+ 1/b)))
    }
    else
    {
        p1 = 1 / b
        p2 = -(2*x[u]) / b
        p3 = -(2*y[u]) / b
        p4 = (x[u]^2 + y[u]^2) / b
    }
    stem_cost[u,] = c(p1,p2,p3,p4)
    return (stem_cost)
}

node.cost = function(node, x, y, parent, time, node_cost, stem_cost) {
    u = node + 1
    children = which(parent == node)
    num_children = length(children)
    if (num_children == 0) return (node_cost)
    node_cost[u,] = colSums(stem_cost[children,,drop=FALSE])
    return (node_cost)
}

final.cost = function(node, x, y, parent, time, node_cost, stem_cost,
    final_cost)
{
    u = node + 1
    if (parent[u] == -1) {
        final_cost[u,] = node_cost[u,]
        return (final_cost)
    }
    v = u
    u = parent[v] + 1
    b = if (!is.null(time)) time[u] - time[v] else 1
    p = final_cost[u, ] - stem_cost[v, ]
    final_cost[v,1] = p[1] / (b*p[1]+1)
    final_cost[v,2] = p[2] / (b*p[1]+1)
    final_cost[v,3] = p[3] / (b*p[1]+1)
    final_cost[v,4] = p[4] - 
            ((p[2])^2 / (4*(p[1]+ 1/b))) -
            ((p[3])^2 / (4*(p[1]+ 1/b)))
    final_cost[v,] = final_cost[v,] + node_cost[v,]
    return (final_cost)
}

mpr = function(x, y, parent, time, postorder, num_nodes) {
    node_cost = matrix(0, num_nodes, 4)
    stem_cost = matrix(0, num_nodes, 4)
    final_cost = matrix(0, num_nodes, 4)
    for (node in postorder)
    {
        node_cost = node.cost(node,x,y,parent,time,node_cost,stem_cost)
        stem_cost = stem.cost(node,x,y,parent,time,node_cost,stem_cost)
    }
    for (node in rev(postorder))
    {
        final_cost = final.cost(
            node,x,y,parent,time,node_cost,stem_cost,final_cost)
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
# node_id    time        x      y
#  0          0          0.50   2.70
#  1          0          1.25   0.41
#  2          0          2.00   1.50
#  3          0.15       -      -
#  4          0.60       -      -
#  5          0.80       -      -
#  6          1.00       -      -
#

ts = treeseq_load("./test.trees")

t1 = c(6,4,4,-1,6,-1,-1)
t2 = c(3,4,3,4,-1,-1,-1)
t3 = c(5,4,4,-1,5,-1,-1)
time = c(0,0,0,0.15,0.6,0.8,1.0)
x = c(0.5, 1.25, 2)
y = c(2.7, .41, 1.5)

S = treeseq_quadratic_mpr(ts, cbind(node_id=0:2, x=x, y=y), TRUE)$mpr

s1 = mpr(x, y, t1, time, c(1,2,4,0,6), 7)
s2 = mpr(x, y, t2, time, c(0,2,3,1,4), 7)
s3 = mpr(x, y, t3, time, c(1,2,4,0,5), 7)
s4 = .2*s1[5,] + .6*s2[5,] + .2*s3[5,]


stopifnot(all.equal(S[, 4], s2[4, ]))
stopifnot(all.equal(S[, 5], s4))
stopifnot(all.equal(S[, 6], s3[6, ]))
stopifnot(all.equal(S[, 7], s1[7, ]))


S = treeseq_quadratic_mpr(ts, cbind(node_id=0:2, x=x, y=y), FALSE)$mpr

s1 = mpr(x, y, t1, NULL, c(1,2,4,0,6), 7)
s2 = mpr(x, y, t2, NULL, c(0,2,3,1,4), 7)
s3 = mpr(x, y, t3, NULL, c(1,2,4,0,5), 7)
s4 = .2*s1[5,] + .6*s2[5,] + .2*s3[5,]


stopifnot(all.equal(S[, 4], s2[4, ]))
stopifnot(all.equal(S[, 5], s4))
stopifnot(all.equal(S[, 6], s3[6, ]))
stopifnot(all.equal(S[, 7], s1[7, ]))

