#' Access the edge table of a tree sequence
#'
#' @description
#' Returns the edge table from a tree sequence, which defines the genealogical 
#' relationships between nodes. Each edge records a parent-child relationship along 
#' with the genomic interval over which that relationship holds.
#'
#' @param ts A \code{treeseq} object
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{edge_id}{Unique identifier for each edge (0-based indexing)}
#'     \item{left}{Left endpoint of the genomic interval (inclusive)}
#'     \item{right}{Right endpoint of the genomic interval (exclusive)}
#'     \item{parent_id}{Node ID of parent (0-based indexing)}
#'     \item{child_id}{Node ID of child (0-based indexing)}
#'   }
#'
#' @details
#' The edge table is one of the core components of a tree sequence, defining how 
#' genetic material is inherited between nodes. Each row describes a parent-child 
#' relationship that exists over a specific genomic interval.
#'
#' The genomic interval [left, right) uses the standard convention where left is 
#' inclusive and right is exclusive. Intervals are in units defined by the tree 
#' sequence (typically base pairs or genetic map positions).
#'
#' Node IDs reference entries in the node table, which can be accessed via 
#' \code{treeseq_nodes}. Both edge_id and node IDs use 0-based indexing.
#'
#' @seealso
#' \code{\link{treeseq_nodes}} for accessing node information
#' \code{\link{treeseq_num_edges}} for counting edges
#' \code{\link{treeseq_drop_edges}} for removing edges
#'
#' @examples
#' # Load tree sequence
#' ts = treeseq_load(system.file("extdata", "test.trees", package="gaia"))
#' 
#' # Get edge table
#' edges = treeseq_edges(ts)
#' 
#' # Find edges from node 4 (appears in all three trees)
#' parent_edges = edges[edges$parent_id == 4, ]
#' 
#' # Find edges to node 1 (changes parent across trees)
#' child_edges = edges[edges$child_id == 1, ]
treeseq_edges = function(ts)
{
    stopifnot(inherits(ts, "treeseq"))
    if (is.null(edges <- attr(ts, "edges")))
    {
        edges = .Call(C_treeseq_edges, ts@treeseq)
        attr(ts, "edges") = edges
    }
    edges
}


#' Count the number of edges in a tree sequence
#'
#' @description
#' Returns the total number of edges in a tree sequence. Each edge represents a 
#' parent-child relationship over a specific genomic interval.
#'
#' @param ts A \code{treeseq} object
#'
#' @return Integer giving the number of edges
#'
#' @details
#' The number of edges corresponds to the number of rows in the edge table. Each 
#' edge defines an ancestor-descendant relationship that exists over some portion 
#' of the genome.
#'
#' This count includes all edges, regardless of their genomic span. The same 
#' parent-child pair may be represented by multiple edges if their relationship 
#' exists over disconnected genomic intervals.
#'
#' @seealso
#' \code{\link{treeseq_edges}} for accessing the full edge table
#' \code{\link{treeseq_drop_edges}} for removing edges
#'
#' @examples
#' # Load tree sequence
#' ts = treeseq_load(system.file("extdata", "test.trees", package="gaia"))
#' 
#' # Count edges
#' n_edges = treeseq_num_edges(ts)
treeseq_num_edges = function(ts)
{
    stopifnot(inherits(ts, "treeseq"))
    if (is.null(num_edges <- attr(ts, "num_edges")))
    {
        num_edges = nrow(treeseq_edges(ts))
        attr(ts, "num_edges") = num_edges
    }
    num_edges
}


#' Remove edges from a tree sequence
#'
#' @description
#' Creates a new tree sequence with specified edges removed. Edges can be removed 
#' based on their parent nodes, child nodes, or both. The resulting tree sequence 
#' is automatically simplified to maintain a valid topology.
#'
#' @param ts A \code{treeseq} object
#' @param parent Integer vector of node IDs. All edges emanating from these parent 
#'   nodes will be removed. Uses 0-based indexing.
#' @param child Integer vector of node IDs. All edges entering these child nodes 
#'   will be removed. Uses 0-based indexing.
#'
#' @return A new \code{treeseq} object with the specified edges removed and the 
#'   topology simplified
#'
#' @details
#' This function allows selective pruning of relationships from the tree sequence. 
#' Edges can be removed based on their parent nodes (all edges coming from specific 
#' ancestors), their child nodes (all edges going to specific descendants), or both.
#'
#' After removing edges, the tree sequence is automatically simplified to:
#' 1. Remove any nodes that become disconnected
#' 2. Remove redundant nodes (those with only one child)
#' 3. Merge adjacent edges with identical parent-child relationships
#' 4. Update sample node flags
#'
#' At least one of parent or child must be specified. If both are specified, edges 
#' matching either criterion are removed.
#'
#' @seealso
#' \code{\link{treeseq_edges}} for viewing the edge table
#' \code{\link{treeseq_simplify}} for general tree sequence simplification
#'
#' @examples
#' # Load tree sequence 
#' ts = treeseq_load(system.file("extdata", "test.trees", package="gaia"))
#' 
#' # Remove all edges from node 4 (internal node present in all trees)
#' ts2 = treeseq_drop_edges(ts, parent=4)
#' 
#' # Remove edges to sample nodes 0 and 1
#' ts3 = treeseq_drop_edges(ts, child=c(0,1))
#' 
#' # Remove both sets of edges
#' ts4 = treeseq_drop_edges(ts, parent=4, child=c(0,1))
treeseq_drop_edges = function(ts, parent, child)
{
    stopifnot(inherits(ts, "treeseq"))
    stopifnot(!missing(parent) || !missing(child))
    drop.parent = integer(treeseq_num_nodes(ts))
    if (!missing(parent))
        drop.parent[parent+1L] = 1L
    drop.child = integer(treeseq_num_nodes(ts))
    if (!missing(child))
        drop.child[child+1L] = 1L
    handle = .Call(C_treeseq_drop_edges, ts@treeseq, drop.parent, drop.child)
    ts_new = treeseq()
    ts_new@treeseq = handle[[1L]]
    ts_new@tree = handle[[2L]]
    ts_new
}
