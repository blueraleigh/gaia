#' Write a tree sequence to file
#'
#' @description
#' Saves a tree sequence to a .trees file in tskit format. The saved file can be 
#' loaded back into R or used with other software that supports the tskit format.
#'
#' @param ts A \code{treeseq} object
#' @param filename Character string giving the path where the tree sequence should 
#'   be saved. The path will be created if it doesn't exist.
#'
#' @return None (called for side effect)
#'
#' @details
#' The tree sequence is saved in tskit's native .trees format. This format 
#' efficiently encodes all tables (nodes, edges, individuals, etc.) and their 
#' relationships.
#'
#' The saved file can be:
#' - Loaded back into R using \code{treeseq_load}
#' - Used with tskit's Python API
#' - Used with other software supporting the tskit format
#'
#' @seealso
#' \code{\link{treeseq_load}} for loading saved tree sequences
#'
#' @examples
#' # Load example tree sequence
#' ts <- treeseq_load(system.file("extdata", "example.trees", package="gaia"))
#'
#' # Save to new file
#' treeseq_write(ts, "output.trees")
#'
#' @export
treeseq_write = function(ts, filename)
{
    stopifnot(inherits(ts, "treeseq"))
    .Call(C_treeseq_write, ts@treeseq, normalizePath(filename, mustWork=FALSE))
    invisible()
}
