% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/treeseq_nodes.R
\name{treeseq_nodes}
\alias{treeseq_nodes}
\title{Access the nodes table of a tree sequence}
\usage{
treeseq_nodes(ts)
}
\arguments{
\item{ts}{A \code{treeseq} object}
}
\value{
A data frame with columns:
  \describe{
    \item{node_id}{Unique identifier for each node (0-based indexing)}
    \item{is_sample}{Logical indicating whether node represents a sample}
    \item{time}{Time ago when node existed (present = 0, increasing into past)}
    \item{population_id}{ID of population node belongs to (0-based indexing)}
    \item{individual_id}{ID of individual node belongs to (0-based indexing)}
    \item{metadata}{Raw vector containing node-specific metadata}
  }
}
\description{
Returns the nodes table from a tree sequence, which contains core information 
about all nodes including their times, population assignments, sample flags, and 
individual associations.
}
\details{
The nodes table is a fundamental component of tree sequences, storing information 
about every node that appears in the genealogical trees. Each node represents a 
haploid genome that existed at a particular time.

Nodes can be:
- Samples (observed data) or non-samples (inferred ancestors)
- Associated with specific populations
- Associated with specific individuals (e.g., maternal/paternal genomes)
- Annotated with arbitrary metadata

Times are measured backwards from the present (0), increasing into the past. The 
time units (e.g., generations, years) depend on how the tree sequence was created.
}
\examples{
# Load tree sequence
ts = treeseq_load(system.file("extdata", "test.trees", package="gaia"))

# Get nodes table
nodes = treeseq_nodes(ts)

# Find sample nodes (0-2)
samples = nodes[nodes$is_sample, ]

# Find oldest nodes (nodes 5 and 6 at times 0.8 and 1.0)
oldest = nodes[nodes$time > 0.7, ]
}
\seealso{
\code{\link{treeseq_edges}} for accessing edge information
\code{\link{treeseq_individuals}} for accessing individual information
\code{\link{treeseq_populations}} for accessing population information
}
