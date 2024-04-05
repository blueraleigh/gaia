#' Write a tree sequence to a file
#'
#' @param ts A \code{treeseq} object.
#' @param filename The name of the file to write to.
treeseq_write = function(ts, filename)
{
    stopifnot(inherits(ts, "treeseq"))
    .Call(C_treeseq_write, ts@treeseq, normalizePath(filename, mustWork=FALSE))
    invisible()
}
