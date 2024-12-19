#' S4 class representing a tree sequence
#'
#' @description
#' A class that encapsulates a tree sequence - a data structure recording the 
#' genealogical relationships between genetic sequences. The class maintains both 
#' the complete tree sequence and a pointer to the current local tree being 
#' examined.
#'
#' @slot treeseq External pointer to the complete tree sequence
#' @slot tree External pointer to the current local tree
#'
#' @details
#' Tree sequences efficiently encode the evolutionary relationships between sequences 
#' by storing a sequence of local trees along a genome. Each local tree represents 
#' the genealogical relationships for a specific genomic interval.
#'
#' The class maintains two key components:
#' 1. The complete tree sequence containing all genealogical information
#' 2. A pointer to one local tree, which can be updated to examine different 
#'    genomic positions
#'
#' The current local tree can be changed using \code{treeseq_sample} and viewed 
#' using \code{treeseq_to_phylo}.
#'
#' @seealso
#' \code{\link{treeseq_load}} for creating tree sequence objects
#' \code{\link{treeseq_sample}} for accessing different local trees
#'
#' @examples
#' # Create a new empty tree sequence object
#' ts <- new("treeseq")
#'
#' # More commonly, load from file
#' ts <- treeseq_load(system.file("extdata", "example.trees", package="gaia"))
#'
#' @export
treeseq = setClass(
    "treeseq"
    , slots=list(
          treeseq="externalptr"
        , tree="externalptr")
)
