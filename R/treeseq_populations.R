#' Access the populations table of a tree sequence
#'
#' @description
#' Returns the populations table from a tree sequence, which records information 
#' about distinct populations referenced by nodes. Each population can have 
#' associated metadata.
#'
#' @param ts A \code{treeseq} object
#'
#' @return A matrix with columns:
#'   \describe{
#'     \item{population_id}{Unique identifier for each population (0-based indexing)}
#'     \item{metadata}{Raw vectors containing population-specific metadata}
#'   }
#'
#' @details
#' The populations table stores information about distinct populations in the tree 
#' sequence. Each population represents a group of nodes that can be referenced in 
#' the nodes table via population IDs.
#'
#' Metadata can store arbitrary information about each population in a raw vector 
#' format. The interpretation of metadata depends on how the tree sequence was 
#' created.
#'
#' Population IDs use 0-based indexing and are referenced by the population_id 
#' column in the nodes table.
#'
#' @seealso
#' \code{\link{treeseq_nodes}} for accessing node information
#'
#' @examples
#' # Load tree sequence
#' ts <- treeseq_load(system.file("extdata", "test.trees", package = "gaia"))
#' 
#' # Get populations table
#' pops <- treeseq_populations(ts)
#'
#' @export
treeseq_populations = function(ts)
{
    stopifnot(inherits(ts, "treeseq"))
    .Call(C_treeseq_populations, ts@treeseq)
}
