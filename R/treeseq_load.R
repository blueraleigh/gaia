#' Load a tree sequence from a file
#'
#' @param filename The name of the file with the tree sequence.
#' @return An object of class \code{treeseq}.
treeseq_load <- function(filename) {
  handle <- .Call(C_treeseq_load, normalizePath(filename))
  ts <- treeseq()
  ts@treeseq <- handle[[1L]]
  ts@tree <- handle[[2L]]
  ts
}
