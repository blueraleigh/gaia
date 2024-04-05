#' Return the tree sequence edge table
#'
#' @param ts A \code{treeseq} object.
#' @return The tree sequence edge table as a \code{data.frame} object.
treeseq_edges = function(ts)
{
    stopifnot(inherits(ts, "treeseq"))
    if (is.null(edges <- attr(ts, "edges")))
    {
        edges = .Call(C_treeseq_edges, ts@treeseq)
        attr(ts, "edges") = edges
    }
    edges
}


#' Return the number of edges in a tree sequence
#'
#' @param ts A \code{treeseq} object.
treeseq_num_edges = function(ts)
{
    stopifnot(inherits(ts, "treeseq"))
    if (is.null(num_edges <- attr(ts, "num_edges")))
    {
        num_edges = nrow(treeseq_edges(ts))
        attr(ts, "num_edges") = num_edges
    }
    num_edges
}


#' Drop edges from the tree sequence
#'
#' @param ts A \code{treeseq} object.
#' @param parent An integer vector of node id's. Edges that emanate from nodes
#' whose id is in \code{parent} will be dropped.
#' @param child An integer vector of node id's. Edges that enter nodes whose 
#' id is in \code{child} will be dropped.
#' @return A new \code{treeseq} object with the specified edges removed.
treeseq_drop_edges = function(ts, parent, child)
{
    stopifnot(inherits(ts, "treeseq"))
    stopifnot(!missing(parent) || !missing(child))
    drop.parent = integer(treeseq_num_nodes(ts))
    if (!missing(parent))
        drop.parent[parent+1L] = 1L
    drop.child = integer(treeseq_num_nodes(ts))
    if (!missing(child))
        drop.child[child+1L] = 1L
    handle = .Call(C_treeseq_drop_edges, ts@treeseq, drop.parent, drop.child)
    ts_new = treeseq()
    ts_new@treeseq = handle[[1L]]
    ts_new@tree = handle[[2L]]
    ts_new
}
