% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/treeseq_sample.R
\name{treeseq_sample}
\alias{treeseq_sample}
\title{Seek to a local tree in a tree sequence}
\usage{
treeseq_sample(ts, at = -1L)
}
\arguments{
\item{ts}{A \code{treeseq} object.}

\item{at}{The index of the local tree to seek to. If less than or equal to 
zero (default), a random index is used.}
}
\value{
The index of the local tree that \code{ts@tree} points to.
}
\description{
Seek to a local tree in a tree sequence
}
\details{
While this function does not return anything, it modifies the
underlying representation of the tree sequence so that \code{ts@tree}
points to a specific local tree as determined by the value of the \code{at}
parameter.
}
