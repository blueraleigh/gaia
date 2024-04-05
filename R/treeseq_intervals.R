#' Return the genomic spans of local trees in a tree sequence
#'
#' @param ts A \code{treeseq} object.
#' @return A matrix containing a row for each local tree in the tree sequence
#' specifying the genomic interval over which it relates the samples.
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
