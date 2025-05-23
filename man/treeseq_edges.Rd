% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/treeseq_edges.R
\name{treeseq_edges}
\alias{treeseq_edges}
\title{Access the edge table of a tree sequence}
\usage{
treeseq_edges(ts)
}
\arguments{
\item{ts}{A \code{treeseq} object}
}
\value{
A data frame with columns:
  \describe{
    \item{edge_id}{Unique identifier for each edge (0-based indexing)}
    \item{left}{Left endpoint of the genomic interval (inclusive)}
    \item{right}{Right endpoint of the genomic interval (exclusive)}
    \item{parent_id}{Node ID of parent (0-based indexing)}
    \item{child_id}{Node ID of child (0-based indexing)}
  }
}
\description{
Returns the edge table from a tree sequence, which defines the genealogical 
relationships between nodes. Each edge records a parent-child relationship along 
with the genomic interval over which that relationship holds.
}
\details{
The edge table is one of the core components of a tree sequence, defining how 
genetic material is inherited between nodes. Each row describes a parent-child 
relationship that exists over a specific genomic interval.

The genomic interval [left, right) uses the standard convention where left is 
inclusive and right is exclusive. Intervals are in units defined by the tree 
sequence (typically base pairs or genetic map positions).

Node IDs reference entries in the node table, which can be accessed via 
\code{treeseq_nodes}. Both edge_id and node IDs use 0-based indexing.
}
\examples{
# Load tree sequence
ts = treeseq_load(system.file("extdata", "test.trees", package="gaia"))

# Get edge table
edges = treeseq_edges(ts)

# Find edges from node 4 (appears in all three trees)
parent_edges = edges[edges$parent_id == 4, ]

# Find edges to node 1 (changes parent across trees)
child_edges = edges[edges$child_id == 1, ]
}
\seealso{
\code{\link{treeseq_nodes}} for accessing node information
\code{\link{treeseq_num_edges}} for counting edges
\code{\link{treeseq_drop_edges}} for removing edges
}
