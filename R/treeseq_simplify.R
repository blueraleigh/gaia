#' Simplify a tree sequence
#'
#' @description
#' Creates a simplified tree sequence containing only the genealogical relationships 
#' relevant to a specified subset of samples. The simplified sequence maintains all 
#' necessary ancestry while removing extraneous nodes and edges.
#'
#' @param ts A \code{treeseq} object
#' @param samples Integer vector of node IDs specifying which samples to retain 
#'   (0-based indexing)
#' @param node.map Logical indicating whether to return a mapping from old to new 
#'   node IDs (default: FALSE)
#' @param filter.sites Logical indicating whether to remove sites not ancestral to 
#'   retained samples (default: TRUE)
#' @param filter.populations Logical indicating whether to remove unused populations 
#'   (default: TRUE)
#' @param filter.individuals Logical indicating whether to remove unused individuals 
#'   (default: TRUE)
#' @param no.filter.nodes Logical indicating whether to retain nodes not ancestral 
#'   to samples (default: FALSE)
#' @param no.update.sample.flags Logical indicating whether to preserve original 
#'   sample flags (default: FALSE)
#' @param reduce.to.site.topology Logical indicating whether to remove edges not 
#'   supporting variants (default: FALSE)
#' @param keep.unary Logical indicating whether to retain nodes with single child 
#'   (default: FALSE)
#' @param keep.input.roots Logical indicating whether to retain original root nodes 
#'   (default: FALSE)
#' @param keep.unary.in.individuals Logical indicating whether to retain unary nodes 
#'   in individuals (default: FALSE)
#'
#' @return If node.map=FALSE (default), returns a new simplified \code{treeseq} 
#'   object. If node.map=TRUE, returns a \code{treeseq} object with an additional 
#'   attribute 'node.map' containing integer vector mapping original to new node IDs.
#'
#' @details
#' Simplification reduces a tree sequence to only the genealogical relationships 
#' necessary to describe the ancestry of a specified set of samples. This process:
#' 1. Removes nodes and edges not ancestral to retained samples
#' 2. Optionally removes sites, populations, and individuals no longer referenced
#' 3. Merges redundant edges
#' 4. Updates sample flags
#' 
#' The process preserves all genealogical information relevant to the retained 
#' samples while potentially greatly reducing the size of the tree sequence.
#'
#' Multiple filtering options control what information is retained or removed during 
#' simplification.
#'
#' @seealso
#' \code{\link{treeseq_drop_edges}} for selective edge removal
#'
#' @examples
#' # Load tree sequence 
#' ts <- treeseq_load(system.file("extdata", "test.trees", package = "gaia"))
#' 
#' # Simplify to just nodes 0 and 1
#' ts2 <- treeseq_simplify(ts, samples=c(0,1))
#' 
#' # Get mapping between old and new node IDs for all nodes 0-6
#' ts3 <- treeseq_simplify(ts, samples=c(0,1), node.map=TRUE)
#' node_map <- attr(ts3, "node.map")
#'
#' @export
treeseq_simplify = function(
    ts,
    samples,
    node.map=FALSE,
    filter.sites=TRUE,
    filter.populations=TRUE,
    filter.individuals=TRUE,
    no.filter.nodes=FALSE,
    no.update.sample.flags=FALSE,
    reduce.to.site.topology=FALSE,
    keep.unary=FALSE,
    keep.input.roots=FALSE,
    keep.unary.in.individuals=FALSE)
{
    stopifnot(inherits(ts, "treeseq"))
    storage.mode(samples) = "integer"
    TSK_SIMPLIFY_FILTER_SITES = bitwShiftL(1L, 0L)
    TSK_SIMPLIFY_FILTER_POPULATIONS = bitwShiftL(1L, 1L)
    TSK_SIMPLIFY_FILTER_INDIVIDUALS = bitwShiftL(1L, 2L)
    TSK_SIMPLIFY_NO_FILTER_NODES = bitwShiftL(1L, 7L)
    TSK_SIMPLIFY_NO_UPDATE_SAMPLE_FLAGS = bitwShiftL(1L, 8L)
    TSK_SIMPLIFY_REDUCE_TO_SITE_TOPOLOGY = bitwShiftL(1L, 3L)
    TSK_SIMPLIFY_KEEP_UNARY = bitwShiftL(1L, 4L)
    TSK_SIMPLIFY_KEEP_INPUT_ROOTS = bitwShiftL(1L, 5L)
    TSK_SIMPLIFY_KEEP_UNARY_IN_INDIVIDUALS = bitwShiftL(1L, 6L)
    flags = 0L
    if (filter.sites)
        flags = bitwOr(flags, TSK_SIMPLIFY_FILTER_SITES)
    if (filter.populations)
        flags = bitwOr(flags, TSK_SIMPLIFY_FILTER_POPULATIONS)
    if (filter.individuals)
        flags = bitwOr(flags, TSK_SIMPLIFY_FILTER_INDIVIDUALS)
    if (no.filter.nodes)
        flags = bitwOr(flags, TSK_SIMPLIFY_NO_FILTER_NODES)
    if (no.update.sample.flags)
        flags = bitwOr(flags, TSK_SIMPLIFY_NO_UPDATE_SAMPLE_FLAGS)
    if (reduce.to.site.topology)
        flags = bitwOr(flags, TSK_SIMPLIFY_REDUCE_TO_SITE_TOPOLOGY)
    if (keep.unary)
        flags = bitwOr(flags, TSK_SIMPLIFY_KEEP_UNARY)
    if (keep.input.roots)
        flags = bitwOr(flags, TSK_SIMPLIFY_KEEP_INPUT_ROOTS)
    if (keep.unary.in.individuals)
        flags = bitwOr(flags, TSK_SIMPLIFY_KEEP_UNARY_IN_INDIVIDUALS)
    handle = .Call(C_treeseq_simplify, ts@treeseq, samples, flags)
    ts = treeseq()
    ts@treeseq = handle[[1L]]
    ts@tree = handle[[2L]]
    if (!node.map)
        return (ts)
    return (structure(ts, node.map=handle[[3]]))
}
