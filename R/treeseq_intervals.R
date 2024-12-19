#' Get information about local trees in a tree sequence
#'
#' @description
#' Returns details about each local tree in a tree sequence, including their genomic 
#' spans, number of edges, number of roots, and maximum root age.
#'
#' @param ts A \code{treeseq} object
#'
#' @return A matrix with columns:
#'   \describe{
#'     \item{left}{Left endpoint of tree's genomic interval (inclusive)}
#'     \item{right}{Right endpoint of tree's genomic interval (exclusive)}
#'     \item{length}{Length of genomic interval (right - left)}
#'     \item{num_edges}{Number of edges in the tree}
#'     \item{num_roots}{Number of root nodes in the tree}
#'     \item{max_root_age}{Time of the oldest root in the tree}
#'   }
#'
#' @details
#' A tree sequence consists of a sequence of local trees along a genome, where each 
#' tree represents the genealogical relationships for a specific genomic interval. 
#' This function summarizes key properties of each local tree.
#'
#' Genomic intervals use the standard [left, right) convention where left is inclusive 
#' and right is exclusive. The units of these intervals (e.g., base pairs, genetic 
#' map positions) depend on how the tree sequence was created.
#'
#' A tree can have multiple roots when there are lineages that have not yet coalesced 
#' at the top of the tree. The age of a root is its time value in the node table.
#'
#' @seealso
#' \code{\link{treeseq_sample}} for accessing specific local trees
#' \code{\link{treeseq_to_phylo}} for converting local trees to phylo objects
#'
#' @examples
#' # Load example tree sequence
#' ts <- treeseq_load(system.file("extdata", "example.trees", package="gaia"))
#'
#' # Get information about all local trees
#' trees <- treeseq_intervals(ts)
#'
#' # Find trees with multiple roots
#' multi_root <- which(trees[,"num_roots"] > 1)
#'
#' # Find oldest tree
#' oldest <- which.max(trees[,"max_root_age"])
#'
#' @export
treeseq_intervals = function(ts)
{
    stopifnot(inherits(ts, "treeseq"))
    structure(
        .Call(C_treeseq_intervals, ts@tree)
        , dimnames=list(
            NULL
            , c("left","right","length","num_edges","num_roots","max_root_age")
        )
    )
}
