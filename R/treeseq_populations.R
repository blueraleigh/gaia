#' Return the tree sequence populations table
#'
#' @param ts A \code{treeseq} object.
#' @return The tree sequence populations table as a \code{matrix} object.
treeseq_populations = function(ts)
{
    stopifnot(inherits(ts, "treeseq"))
    .Call(C_treeseq_populations, ts@treeseq)
}
