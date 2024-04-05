#' Seek to a local tree in a tree sequence
#'
#' @param ts A \code{treeseq} object.
#' @param at The index of the local tree to seek to. If less than or equal to 
#' zero (default), a random index is used.
#' @details While this function does not return anything, it modifies the
#' underlying representation of the tree sequence so that \code{ts@tree}
#' points to a specific local tree as determined by the value of the \code{at}
#' parameter.
#' @return The index of the local tree that \code{ts@tree} points to.
treeseq_sample = function(ts, at=-1L)
{
    stopifnot(inherits(ts, "treeseq"))
    storage.mode(at) = "integer"
    
    .Call(C_treeseq_sample, ts@tree, at)
}
