#' Return the tree sequence node table
#'
#' @param ts A \code{treeseq} object.
#' @return The tree sequence node table as a \code{data.frame} object.
treeseq_nodes = function(ts)
{
    stopifnot(inherits(ts, "treeseq"))
    if (is.null(nodes <- attr(ts, "nodes")))
    {
        nodes = .Call(C_treeseq_nodes, ts@treeseq)
        attr(ts, "nodes") = nodes
    }
    nodes
}

#' Return the number of nodes in a tree sequence
#'
#' @param ts A \code{treeseq} object.
treeseq_num_nodes = function(ts)
{
    stopifnot(inherits(ts, "treeseq"))
    if (is.null(num_nodes <- attr(ts, "num_nodes")))
    {
        num_nodes = nrow(treeseq_nodes(ts))
        attr(ts, "num_nodes") = num_nodes
    }
    num_nodes
}

#' Return the number of samples in a tree sequence
#'
#' @param ts A \code{treeseq} object.
treeseq_num_samples = function(ts)
{
    stopifnot(inherits(ts, "treeseq"))
    if (is.null(num_samples <- attr(ts, "num_samples")))
    {
        num_samples = sum(treeseq_nodes(ts)[, "is_sample"])
        attr(ts, "num_samples") = num_samples
    }
    num_samples
}