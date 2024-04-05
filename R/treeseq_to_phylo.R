#' Extract a local tree from a tree sequence
#'
#' @param ts A \code{treeseq} object.
#' @return A \code{phlyo} object representing the current local tree.
#' @details Two additional components are added to the returned \code{phylo}
#' object
#' \describe{
#' \item{$node.id}{A vector that maps the node indices in the \code{phylo}
#' object to their corresponding identifiers in the tree sequence. E.g.,
#' \code{$node.id[i]} will return the tree sequence node id of the node with
#' index \code{i} in the \code{phylo} object.} 
#' \item{$edge.id}{A vector that maps the edge indices in the \code{phylo}
#' object to their corresponding identifiers in the tree sequence.}
#'}
#' NB: both of these components are invalidated if the \code{phylo} object is
#' reordered by a call to \code{ape::reorder.phylo}.
#' @note A call to this function should be preceded by a call to 
#' \code{treeseq_sample}.
treeseq_to_phylo = function(ts)
{
    stopifnot(inherits(ts, "treeseq"))
    .Call(C_treeseq_to_phylo, ts@tree)
}
