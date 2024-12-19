#' Access the individuals table of a tree sequence
#'
#' @description
#' Returns the individuals table from a tree sequence, which records information about 
#' distinct individuals (organisms) represented in the tree sequence. Each individual 
#' can have spatial locations, parent-child relationships, and metadata.
#'
#' @param ts A \code{treeseq} object
#'
#' @return A matrix with columns:
#'   \describe{
#'     \item{individual_id}{Unique identifier for each individual (0-based indexing)}
#'     \item{location}{Numeric vector of spatial coordinates}
#'     \item{parents}{Integer vector of parent individual IDs (0-based indexing)}
#'     \item{metadata}{Raw vector containing individual-specific metadata}
#'   }
#'
#' @details
#' The individuals table stores information about distinct organisms in the tree 
#' sequence. Each individual can:
#' - Have multiple nodes associated with it (e.g., maternal and paternal genomes)
#' - Be located in space via coordinates
#' - Have recorded parents
#' - Carry additional metadata
#'
#' Location coordinates and parent IDs are stored as variable-length vectors for each 
#' individual. The interpretation of location coordinates (e.g., as geographic 
#' coordinates, plot coordinates) depends on how the tree sequence was created.
#'
#' Individual IDs use 0-based indexing and are referenced by the individual_id 
#' column in the nodes table.
#'
#' @seealso
#' \code{\link{treeseq_nodes}} for accessing node information
#' \code{\link{treeseq_populations}} for accessing population information
#'
#' @examples
#' # Load tree sequence
#' ts <- treeseq_load(system.file("extdata", "test.trees", package = "gaia"))
#' 
#' # Get individuals table
#' inds <- treeseq_individuals(ts)
#'
#' @export
treeseq_individuals = function(ts)
{
    stopifnot(inherits(ts, "treeseq"))
    .Call(C_treeseq_individuals, ts@treeseq)
}