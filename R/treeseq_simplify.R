#' Simplify a tree sequence
#'
#' @param ts A \code{treeseq} object.
#' @param samples A vector of sample id's.
#' @return A new tree sequence pruned to only include the samples listed in
#' the \code{samples} parameter.
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
