#' An S4 class to wrap a tskit tree sequence
#'
#' @slot treeseq The tree sequence object.
#' @slot tree A single marginal tree in the tree sequence.
treeseq = setClass(
    "treeseq"
    , slots=list(
          treeseq="externalptr"
        , tree="externalptr")
)
