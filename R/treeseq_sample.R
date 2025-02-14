#' Select a local tree from a tree sequence
#'
#' @description
#' Updates the current local tree pointer in a tree sequence object to reference 
#' either a specific tree or a randomly chosen tree. The selected tree can then be 
#' accessed using other functions.
#'
#' @param ts A \code{treeseq} object
#' @param at Integer specifying which local tree to select (1-based indexing). If 
#'   less than or equal to 0 (default), a random tree is selected with probability 
#'   proportional to its genomic span.
#'
#' @return Integer giving the index of the selected tree (1-based)
#'
#' @details
#' A tree sequence consists of a sequence of local trees along a genome. This 
#' function updates which local tree is currently being referenced by the tree 
#' sequence object.
#'
#' When at ≤ 0, a tree is chosen randomly with probability proportional to its 
#' genomic span. This means trees covering longer genomic intervals are more likely 
#' to be selected.
#'
#' The selected tree becomes available for:
#' - Converting to phylo format via \code{treeseq_to_phylo}
#' - Exploring its structure and properties
#' - Other tree-specific operations
#'
#' @seealso
#' \code{\link{treeseq_to_phylo}} for converting the selected tree
#' \code{\link{treeseq_intervals}} for information about all local trees
#'
#' @examples
#' # Load tree sequence
#' ts = treeseq_load(system.file("extdata", "test.trees", package = "gaia"))
#' 
#' # Select first tree (spanning [0,20))
#' treeseq_sample(ts, at=1)
#' phy_1 = treeseq_to_phylo(ts)
#'
#' # Select second tree (spanning [20,80))
#' treeseq_sample(ts, at=2)
#' phy_2 = treeseq_to_phylo(ts)
#'
#' # Select a random tree, weighted by genomic span
#' # Most likely to select second tree as it spans 60% of sequence
#' treeseq_sample(ts)
#' phy_rand = treeseq_to_phylo(ts)
#'
treeseq_sample = function(ts, at=-1L)
{
    stopifnot(inherits(ts, "treeseq"))
    storage.mode(at) = "integer"
    
    .Call(C_treeseq_sample, ts@tree, at)
}
