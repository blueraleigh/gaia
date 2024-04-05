#' Return the tree sequence individuals table
#'
#' @param ts A \code{treeseq} object.
#' @return The tree sequence individuals table as a \code{matrix} object.
treeseq_individuals = function(ts)
{
    stopifnot(inherits(ts, "treeseq"))
    .Call(C_treeseq_individuals, ts@treeseq)
}