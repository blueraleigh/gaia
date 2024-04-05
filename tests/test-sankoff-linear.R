library(gaia)

new_plf = function(max_num_breaks) {
    structure(local({
        max_num_breaks = max_num_breaks
        slope = numeric(max_num_breaks)
        breakpoint = numeric(max_num_breaks)
        intercept = 0
        num_breaks = 0L
        function(x) {
            if (num_breaks == 0L) return (0)
            if (x < breakpoint[1]) return (intercept + x*slope[1])
            j = 1L
            score = intercept + slope[1]*breakpoint[1]
            while (j < num_breaks && x > breakpoint[j+1])
            {
                j = j + 1L
                score = score + slope[j] * (breakpoint[j] - breakpoint[j-1])
            }
            score = score + slope[j+1] * (x - breakpoint[j])
            return (score)
        }
    }), class=c("plf", "function"))
}


# ret(x) = min_z { c(x,z) + f(z) } where the cost c(x,z)
# of going from x to z is b * abs(x - z)
plf_min = function(f, b) {
    num_breaks = 0L
    minus_b = -b

    fs = environment(f)$slope
    fb = environment(f)$breakpoint
    fi = environment(f)$intercept
    fn = environment(f)$num_breaks

    ret = new_plf(environment(f)$max_num_breaks)
    rs = environment(ret)$slope
    rb = environment(ret)$breakpoint

    if (minus_b <= fs[1])
    {
        x_left = 1L
        ri = fi
    }
    else
    {
        i = 1L
        shift = fi + fs[1]*fb[1]
        while (minus_b >= fs[i+1])
        {
            i = i + 1L
            shift = shift + fs[i] * (fb[i] - fb[i-1])
        }
        x_left = i + 1L
        ri = shift + b * fb[i]
        rs[1] = minus_b
        rb[1] = fb[i]
        num_breaks = num_breaks + 1L
    }
    if (b >= fs[fn+1])
    {
        for (i in x_left:fn)
        {
            num_breaks = num_breaks + 1L
            rs[num_breaks] = fs[i]
            rb[num_breaks] = fb[i]
        }
        rs[num_breaks + 1] = fs[fn + 1]
    }
    else
    {
        i = x_left
        while (b >= fs[i])
        {
            num_breaks = num_breaks + 1L
            rs[num_breaks] = fs[i]
            rb[num_breaks] = fb[i]
            i = i + 1L
        }
        # check for colinearity to avoid superfluous breakpoints
        if (isTRUE(all.equal(b, rs[num_breaks])))
            num_breaks = num_breaks - 1L
        else
            rs[num_breaks + 1] = b
    }

    assign("slope", rs[1:(num_breaks+1)], envir=environment(ret))
    assign("breakpoint", rb[1:num_breaks], envir=environment(ret))
    assign("num_breaks", num_breaks, envir=environment(ret))
    assign("intercept", ri, envir=environment(ret))

    return (ret)
}


plf_add = function(a, b, ascale, bscale, sign) {
    ab = environment(a)$breakpoint
    bb = environment(b)$breakpoint
    as = environment(a)$slope
    bs = environment(b)$slope
    ai = environment(a)$intercept
    bi = environment(b)$intercept
    an = environment(a)$num_breaks
    bn = environment(b)$num_breaks

    obj = new_plf(environment(a)$max_num_breaks)
    slope = environment(obj)$slope
    breakpoint = environment(obj)$breakpoint

    i = 1L
    j = 1L
    num_breaks = 0L
    prev_slope = -Inf

    while (i <= an && j <= bn)
    {
        if (ab[i] < bb[j])
        {
            next_break = ab[i]
            next_slope = ascale * as[i] + sign * bscale * bs[j]
            i = i + 1L
        }
        else if (ab[i] > bb[j])
        {
            next_break = bb[j]
            next_slope = ascale * as[i] + sign * bscale * bs[j]
            j = j + 1L
        }
        else
        {
            next_break = ab[i]
            next_slope = ascale * as[i] + sign * bscale * bs[j]
            i = i + 1L
            j = j + 1L
        }
        if (isTRUE(all.equal(prev_slope, next_slope)))
        {
            breakpoint[num_breaks] = next_break
        }
        else
        {
            num_breaks = num_breaks + 1L
            breakpoint[num_breaks] = next_break
            slope[num_breaks] = next_slope
            prev_slope = next_slope
        }
    }
    
    while (i <= an)
    {
        next_break = ab[i]
        next_slope = ascale * as[i] + sign * bscale * bs[j]
        i = i + 1L
        if (isTRUE(all.equal(prev_slope, next_slope)))
        {
            breakpoint[num_breaks] = next_break
        }
        else
        {
            num_breaks = num_breaks + 1L
            breakpoint[num_breaks] = next_break
            slope[num_breaks] = next_slope
            prev_slope = next_slope
        }
    }
    while (j <= bn)
    {
        next_break = bb[j]
        next_slope = ascale * as[i] + sign * bscale * bs[j]
        j = j + 1L
        if (isTRUE(all.equal(prev_slope, next_slope)))
        {
            breakpoint[num_breaks] = next_break
        }
        else
        {
            num_breaks = num_breaks + 1L
            breakpoint[num_breaks] = next_break
            slope[num_breaks] = next_slope
            prev_slope = next_slope
        }
    }
    
    # don't forget the final slope after the last breakpoint
    next_slope = ascale * as[i] + sign * bscale * bs[j]

    if (isTRUE(all.equal(prev_slope, next_slope)))
        num_breaks = num_breaks - 1L
    else
        slope[num_breaks+1L] = next_slope

    assign("slope", slope[1:(num_breaks+1)], envir=environment(obj))
    assign("breakpoint", breakpoint[1:num_breaks], envir=environment(obj))
    assign("num_breaks", num_breaks, envir=environment(obj))
    assign("intercept", ascale * ai + sign * bscale * bi, envir=environment(obj))

    return (obj)
}


stem.cost = function(node, x, parent, time, node_cost, stem_cost) {
    u = node + 1
    if (parent[u] == -1) return (stem_cost)
    b = if (!is.null(time)) 1 / (time[parent[u]+1] - time[u]) else 1
    num_children = length(which(parent == node))
    if (num_children > 0)
    {
        h = plf_min(node_cost[[u]], b)   
    }
    else
    {
        h = new_plf(3L)
        environment(h)$breakpoint = x[u]
        environment(h)$slope = c(-b, b)
        environment(h)$intercept = if (x[u] < 0) -b*x[u] else b*x[u]
        environment(h)$num_breaks = 1L
    }
    stem_cost[[u]] = h
    return (stem_cost)
}

node.cost = function(node, x, parent, time, node_cost, stem_cost) {
    u = node + 1
    children = which(parent == node)
    num_children = length(children)
    if (num_children == 0) return (node_cost)
    tmp = stem_cost[[children[1]]]
    for (child in children[2:num_children]) {
        tmp = plf_add(tmp, stem_cost[[child]], 1, 1, 1)
    }
    node_cost[[u]] = tmp
    return (node_cost)
}

final.cost = function(node, x, parent, time, node_cost, stem_cost,
    final_cost)
{
    u = node + 1
    children = which(parent == node)
    num_children = length(children)
    if (num_children == 0) return (final_cost)
    if (parent[u] == -1) {
        final_cost[[u]] = node_cost[[u]]
        return (final_cost)
    }
    v = u
    u = parent[v] + 1
    b = if (!is.null(time)) 1 / (time[u] - time[v]) else 1
    tmp = plf_add(final_cost[[u]], stem_cost[[v]], 1, 1, -1)
    tmp = plf_min(tmp, b)
    final_cost[[v]] = plf_add(tmp, node_cost[[v]], 1, 1, 1)
    return (final_cost)
}

mpr = function(x, parent, time, postorder, num_nodes) {
    node_cost = vector("list", num_nodes)
    stem_cost = vector("list", num_nodes)
    final_cost = vector("list", num_nodes)
    for (node in postorder)
    {
        node_cost = node.cost(node,x,parent,time,node_cost,stem_cost)
        stem_cost = stem.cost(node,x,parent,time,node_cost,stem_cost)
    }
    for (node in rev(postorder))
    {
        final_cost = final.cost(
            node,x,parent,time,node_cost,stem_cost,final_cost)
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

S = treeseq_linear_mpr(ts, cbind(node_id=0:2, x=x), TRUE)$mpr

s1 = mpr(x, t1, time, c(1,2,4,0,6), 7)
s2 = mpr(x, t2, time, c(0,2,3,1,4), 7)
s3 = mpr(x, t3, time, c(1,2,4,0,5), 7)
s4 = plf_add(plf_add(s1[[5]], s2[[5]], 0.2, 0.6, 1), s3[[5]], 1, 0.2, 1)



stopifnot(all.equal(S[[1]][[4]], as.list(environment(s2[[4]]))[c(2,4,3)], 
    use.names=FALSE, check.attributes=FALSE))
stopifnot(all.equal(S[[1]][[5]], as.list(environment(s4))[c(2,4,3)], 
    use.names=FALSE, check.attributes=FALSE))
stopifnot(all.equal(S[[1]][[6]], as.list(environment(s3[[6]]))[c(2,4,3)], 
    use.names=FALSE, check.attributes=FALSE))
stopifnot(all.equal(S[[1]][[7]], as.list(environment(s1[[7]]))[c(2,4,3)], 
    use.names=FALSE, check.attributes=FALSE))


S = treeseq_linear_mpr(ts, cbind(node_id=0:2, x=x), FALSE)$mpr

s1 = mpr(x, t1, NULL, c(1,2,4,0,6), 7)
s2 = mpr(x, t2, NULL, c(0,2,3,1,4), 7)
s3 = mpr(x, t3, NULL, c(1,2,4,0,5), 7)
s4 = plf_add(plf_add(s1[[5]], s2[[5]], 0.2, 0.6, 1), s3[[5]], 1, 0.2, 1)


stopifnot(all.equal(S[[1]][[4]], as.list(environment(s2[[4]]))[c(2,4,3)], 
    use.names=FALSE, check.attributes=FALSE))
stopifnot(all.equal(S[[1]][[5]], as.list(environment(s4))[c(2,4,3)], 
    use.names=FALSE, check.attributes=FALSE))
stopifnot(all.equal(S[[1]][[6]], as.list(environment(s3[[6]]))[c(2,4,3)], 
    use.names=FALSE, check.attributes=FALSE))
stopifnot(all.equal(S[[1]][[7]], as.list(environment(s1[[7]]))[c(2,4,3)], 
    use.names=FALSE, check.attributes=FALSE))

