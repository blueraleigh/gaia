#' Convert current local tree to phylo object
#'
#' @description
#' Creates an ape-compatible phylogenetic tree object (\code{phylo}) representing 
#' the current local tree in a tree sequence. The resulting object can be used 
#' with ape's plotting and analysis functions.
#'
#' @param ts A \code{treeseq} object
#'
#' @return A \code{phylo} object with additional components:
#'   \describe{
#'     \item{node.id}{Integer vector mapping phylo node indices to tree sequence 
#'       node IDs}
#'     \item{edge.id}{Integer vector mapping phylo edge indices to tree sequence 
#'       edge IDs}
#'   }
#'
#' @details
#' The function converts the current local tree (as indicated by \code{ts@tree}) 
#' to ape's phylo format. The local tree can be changed using \code{treeseq_sample}.
#'
#' Two additional components are added to help map between tree sequence and phylo 
#' objects:
#' - node.id maps from phylo node indices to tree sequence node IDs
#' - edge.id maps from phylo edge indices to tree sequence edge IDs
#'
#' Note: These mappings become invalid if the phylo object is reordered (e.g., via 
#' \code{ape::reorder.phylo}).
#'
#' @section Warning:
#' The function assumes \code{treeseq_sample} has been called to select a local 
#' tree. Without this, it will convert the first tree by default.
#'
#' @seealso
#' \code{\link{treeseq_sample}} for selecting local trees
#'
#' @examples
#' # Load tree sequence
#' ts <- treeseq_load(system.file("extdata", "test.trees", package = "gaia"))
#' 
#' # Select first tree ([0,20) interval)
#' # Sample functionality not currently working, but it should work like this:
#' # treeseq_sample(ts, at=1)
#' 
#' # Convert to phylo object - should show:
#' # - Samples 0-2 at time 0
#' # - Internal nodes 4 and 6 at times 0.6 and 1.0
#' phy <- treeseq_to_phylo(ts)
#' 
#' # Select second tree ([20,80) interval)
#' # Sample functionality not currently working, but it should work like this:
#' # treeseq_sample(ts, at=2)
#' 
#' # Convert to phylo object - should show:
#' # - Samples 0-2 at time 0
#' # - Internal nodes 3 and 4 at times 0.15 and 0.6
#' phy2 <- treeseq_to_phylo(ts)
#'
#' @export
treeseq_to_phylo = function(ts)
{
    stopifnot(inherits(ts, "treeseq"))
    .Call(C_treeseq_to_phylo, ts@tree)
}
