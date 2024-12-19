#' Access the nodes table of a tree sequence
#'
#' @description
#' Returns the nodes table from a tree sequence, which contains core information 
#' about all nodes including their times, population assignments, sample flags, and 
#' individual associations.
#'
#' @param ts A \code{treeseq} object
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{node_id}{Unique identifier for each node (0-based indexing)}
#'     \item{is_sample}{Logical indicating whether node represents a sample}
#'     \item{time}{Time ago when node existed (present = 0, increasing into past)}
#'     \item{population_id}{ID of population node belongs to (0-based indexing)}
#'     \item{individual_id}{ID of individual node belongs to (0-based indexing)}
#'     \item{metadata}{Raw vector containing node-specific metadata}
#'   }
#'
#' @details
#' The nodes table is a fundamental component of tree sequences, storing information 
#' about every node that appears in the genealogical trees. Each node represents a 
#' haploid genome that existed at a particular time.
#'
#' Nodes can be:
#' - Samples (observed data) or non-samples (inferred ancestors)
#' - Associated with specific populations
#' - Associated with specific individuals (e.g., maternal/paternal genomes)
#' - Annotated with arbitrary metadata
#'
#' Times are measured backwards from the present (0), increasing into the past. The 
#' time units (e.g., generations, years) depend on how the tree sequence was created.
#'
#' @seealso
#' \code{\link{treeseq_edges}} for accessing edge information
#' \code{\link{treeseq_individuals}} for accessing individual information
#' \code{\link{treeseq_populations}} for accessing population information
#'
#' @examples
#' # Load tree sequence
#' ts <- treeseq_load(system.file("extdata", "test.trees", package = "gaia"))
#' 
#' # Get nodes table
#' nodes <- treeseq_nodes(ts)
#' 
#' # Find sample nodes (0-2)
#' samples <- nodes[nodes$is_sample, ]
#' 
#' # Find oldest nodes (nodes 5 and 6 at times 0.8 and 1.0)
#' oldest <- nodes[nodes$time > 0.7, ]
#'
#' @export
treeseq_nodes = function(ts)
{
    stopifnot(inherits(ts, "treeseq"))
    if (is.null(nodes <- attr(ts, "nodes")))
    {
        nodes = .Call(C_treeseq_nodes, ts@treeseq)
        attr(ts, "nodes") = nodes
    }
    nodes
}

#' Count the number of nodes in a tree sequence
#'
#' @description
#' Returns the total number of nodes in a tree sequence. Includes all nodes 
#' (samples and non-samples) that appear in any local tree.
#'
#' @param ts A \code{treeseq} object
#'
#' @return Integer giving the number of nodes
#'
#' @details
#' The number of nodes corresponds to the number of rows in the node table. Each 
#' node represents a haploid genome that existed at some point in time and appears 
#' in at least one local tree.
#'
#' This count includes:
#' - Sample nodes (observed data)
#' - Non-sample nodes (inferred ancestors)
#' - Root nodes
#' - Internal nodes
#'
#' @seealso
#' \code{\link{treeseq_nodes}} for accessing the full node table
#' \code{\link{treeseq_num_samples}} for counting just sample nodes
#'
#' @examples
#' # Load tree sequence
#' ts <- treeseq_load(system.file("extdata", "test.trees", package = "gaia"))
#' 
#' # Count nodes (should be 7: samples 0-2 and internal nodes 3-6)
#' n_nodes <- treeseq_num_nodes(ts)
#'
#' @export
treeseq_num_nodes = function(ts)
{
    stopifnot(inherits(ts, "treeseq"))
    if (is.null(num_nodes <- attr(ts, "num_nodes")))
    {
        num_nodes = nrow(treeseq_nodes(ts))
        attr(ts, "num_nodes") = num_nodes
    }
    num_nodes
}

#' Count the number of samples in a tree sequence
#'
#' @description
#' Returns the number of sample nodes in a tree sequence. Samples represent the 
#' observed data (e.g., sequenced genomes) from which the genealogies were inferred.
#'
#' @param ts A \code{treeseq} object
#'
#' @return Integer giving the number of samples
#'
#' @details
#' Samples are nodes that represent directly observed data, as opposed to inferred 
#' ancestors. They are marked with the is_sample flag in the nodes table.
#'
#' The count includes all nodes flagged as samples, regardless of:
#' - Their time (though samples typically exist at time 0)
#' - Their population assignments
#' - Which individual they belong to
#'
#' @seealso
#' \code{\link{treeseq_nodes}} for accessing the full node table
#' \code{\link{treeseq_num_nodes}} for counting all nodes
#'
#' @examples
#' # Load tree sequence
#' ts <- treeseq_load(system.file("extdata", "test.trees", package = "gaia"))
#' 
#' # Count samples (should be 3: nodes 0-2)
#' n_samples <- treeseq_num_samples(ts)
#'
#' @export
treeseq_num_samples = function(ts)
{
    stopifnot(inherits(ts, "treeseq"))
    if (is.null(num_samples <- attr(ts, "num_samples")))
    {
        num_samples = sum(treeseq_nodes(ts)[, "is_sample"])
        attr(ts, "num_samples") = num_samples
    }
    num_samples
}