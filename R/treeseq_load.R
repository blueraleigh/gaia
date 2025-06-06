#' Load a tree sequence from a file
#'
#' @description
#' Creates a \code{treeseq} object by loading a tree sequence from a .trees file. 
#' The tree sequence must be in the tskit format. This function initializes both
#' the full tree sequence and its first local tree.
#'
#' @param filename Character string giving the path to a .trees file
#'
#' @return A \code{treeseq} object containing:
#'   \describe{
#'     \item{treeseq}{External pointer to the full tree sequence}
#'     \item{tree}{External pointer to the current local tree, initially the first tree}
#'   }
#'
#' @details
#' The .trees file must be in the tskit format. Tree sequences can be created and 
#' exported to this format using tools like msprime, SLiM, or tskit's Python API.
#'
#' The returned object provides access to both the complete tree sequence and its 
#' constituent local trees. The tree slot initially points to the first local tree, 
#' but can be updated to other trees using \code{treeseq_sample}.
#'
#' @seealso
#' \code{\link{treeseq_write}} for saving tree sequences
#' \code{\link{treeseq_sample}} for accessing different local trees
#'
#' @examples
#' # Load tree sequence
#' ts = treeseq_load(system.file("extdata", "test.trees", package="gaia"))
#' 
#' # Access node table showing sample nodes 0-2 and internal nodes 3-6
#' nodes = treeseq_nodes(ts)
#' 
#' # Access edge table showing topology changes across three intervals
#' edges = treeseq_edges(ts)
#' 
#' # View first local tree (interval [0,20))
#' tree = treeseq_to_phylo(ts)
treeseq_load = function(filename)
{
    handle = .Call(C_treeseq_load, normalizePath(filename))
    ts = treeseq()
    ts@treeseq = handle[[1L]]
    ts@tree = handle[[2L]]
    ts
}
