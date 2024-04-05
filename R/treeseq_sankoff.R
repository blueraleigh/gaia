#' Compute minimum migration costs for genetic ancestors in a tree sequence
#'
#' @param ts A \code{treeseq} object.
#' @param sample_locations An integer matrix giving the geographic locations
#' of sampled genomes. The 'node_id' column specifies the node identifier of
#' of each sample in the tree sequence, and the 'state_id' columns specifies the
#' geographic state where each sample is found. Note that node identifiers use
#' a 0-based index while state identifiers use a 1-based index. 
#' @param cost_matrix A numeric matrix giving the state-to-state migration
#' costs that will be used by the optimization. All costs should be non-negative
#' and any non-zero costs on the diagonal will be ignored.
#' @param use_brlen A logical specifying whether or not to scale the migration
#' costs by inverse branch lengths. The default (FALSE) ignores branch lengths
#' in the optimization.
#' @details Returns a numeric vector \code{f} for each ancestral node in the
#' tree sequence such that \code{f[i]} gives the average (taken over all base
#' pairs in the ancestor) minimum migration cost needed to explain the
#' geographic distribution of \emph{all} sample nodes when the ancestor is
#' located in geographic state \code{i}.
treeseq_discrete_mpr = function(ts, sample_locations, cost_matrix, 
    use_brlen=FALSE)
{
    stopifnot(inherits(ts, "treeseq"))
    stopifnot(is.matrix(sample_locations))
    stopifnot(is.matrix(cost_matrix))
    stopifnot(!is.null(colnames(sample_locations)))
    stopifnot(all(colnames(sample_locations) %in% c("node_id","state_id")))
    storage.mode(sample_locations) = "integer"
    stopifnot(nrow(cost_matrix) == ncol(cost_matrix))
    stopifnot(all(cost_matrix >= 0L))
    diag(cost_matrix) = 0
    storage.mode(cost_matrix) = "double"
    num_states = nrow(cost_matrix)
    stopifnot(all(sample_locations[,"state_id"] > 0L))
    stopifnot(all(sample_locations[,"state_id"] <= num_states))
    N = nrow(treeseq_nodes(ts))
    G = matrix(0, num_states, N)
    sample_ids = sample_locations[, "node_id"] + 1L
    state_ids = sample_locations[, "state_id"]
    G[, sample_ids] = Inf
    G[cbind(state_ids, sample_ids)] = 0
    structure(.Call(
        C_treeseq_discrete_mpr
        , ts@treeseq
        , as.integer(use_brlen)
        , num_states
        , G
        , cost_matrix),
        names=c("mean_tree_length", "tree_length", "mpr"),
        class=c("discrete", "mpr")
    )
}


#' Assign genetic ancestors to a set of discrete geographic locations
#'
#' @param obj The result of \code{treeseq_discrete_mpr}.
#' @param index1 A logical indicating whether or not the returned state
#' assignments should be indexed from 1 (TRUE) or 0 (FALSE).
#' @details Assigns genetic ancestors in a tree sequence to the geographic
#' location that minimizes the overall average migration cost.
treeseq_discrete_mpr_minimize = function(obj, index1=TRUE)
{
    stopifnot(inherits(obj, "discrete") && inherits(obj, "mpr"))
    stopifnot(is.logical(index1))
    .Call(C_treeseq_discrete_mpr_minimize, obj$mpr, as.integer(index1))
}


#' Sample a migration history for each edge in a tree sequence
#'
#' @param ts An object of class \code{treeseq}.
#' @param obj The result of \code{treeseq_discrete_mpr}.
#' @param cost_matrix A numeric matrix giving the state-to-state migration
#' costs. All costs should be non-negative and any non-zero costs on the 
#' diagonal will be ignored.
#' @param adjacency_matrix A binary matrix specifying geographic connections.
#' Geographic states \code{i} and \code{j} are connected if the corresponding
#' entry in \code{adjacency_matrix} is 1.
#' @param index1 A logical indicating whether or not the returned state
#' assignments should be indexed from 1 (TRUE) or 0 (FALSE).
#' @details Samples a minimum cost migration history for each edge in a tree
#' sequence and returns the history as a data.frame. Each history is a sequence 
#' of visited geographic states ordered from most recent (present) to least 
#' recent (past). The first element in the history is the state at the end
#' of the edge and the last element in the history is the state at the
#' beginning of the edge. Intermediate elements are the states visited along
#' the way. If there are no intermediate elements, the edge remained in the
#' same state for its entire duration.
treeseq_discrete_mpr_edge_history = function(ts, obj, cost_matrix,
    adjacency_matrix, index1=TRUE)
{
    stopifnot(inherits(ts, "treeseq"))
    stopifnot(inherits(obj, "discrete") && inherits(obj, "mpr"))
    stopifnot(is.matrix(cost_matrix))
    stopifnot(is.matrix(adjacency_matrix))
    stopifnot(isSymmetric(adjacency_matrix))
    stopifnot(nrow(cost_matrix) == ncol(cost_matrix))
    stopifnot(all(cost_matrix >= 0L))
    stopifnot(is.logical(index1))
    diag(cost_matrix) = 0
    diag(adjacency_matrix) = 0
    stopifnot(all(rowSums(adjacency_matrix) > 0))
    stopifnot(all(colSums(adjacency_matrix) > 0))
    storage.mode(cost_matrix) = "double"
    if (!inherits(adjacency_matrix, "dgCMatrix"))
    {
        adjacency_matrix = Matrix::Matrix(adjacency_matrix, sparse=TRUE)
        adjacency_matrix = methods::as(adjacency_matrix, "generalMatrix")
        stopifnot(inherits(adjacency_matrix, "dgCMatrix"))
    }

    node_states = treeseq_discrete_mpr_minimize(obj, index1)

    h = .Call(
        C_treeseq_discrete_mpr_edge_history
        , ts@treeseq
        , node_states
        , cost_matrix
        , adjacency_matrix
        , as.integer(index1)
    )
    H = data.frame(
        edge_id=rep.int(0:(length(h[[1]])-1L), times=h[[3]]),
        state_id=do.call(c, h[[1]]),
        time=do.call(c, h[[2]])
    )
    structure(
        H
        , node.state=node_states
        , path.offset=c(cumsum(h[[3]]) - h[[3]], nrow(H)) + as.integer(index1)
    )
}


#' Spatiotemporal ancestry coefficients
#'
#' @param ts An object of class \code{treeseq}.
#' @param obj The result of \code{treeseq_discrete_mpr}.
#' @param cost_matrix A numeric matrix giving the state-to-state migration
#' costs. All costs should be non-negative and any non-zero costs on the 
#' diagonal will be ignored.
#' @param adjacency_matrix A binary matrix specifying geographic connections.
#' Geographic states \code{i} and \code{j} are connected if the corresponding
#' entry in \code{adjacency_matrix} is 1.
#' @param times A set of times ordered from present to past at which to compute
#' ancestry coefficients.
#' @param state_sets An integer vector with one entry for each geographic
#' state specifying the set to which that state belongs. Sets must be numbered
#' beginning from 1.
#' @param sample_sets An integer vector with one entry for each sample
#' specifying the set to which that sample belongs. Sets must be numbered
#' beginning from 1. A value of 0 may be used to indicate that the sample
#' belongs to none of the sets and hence is to be excluded from the calculation.
treeseq_discrete_mpr_ancestry = function(ts, obj, cost_matrix,
    adjacency_matrix, times, state_sets, sample_sets)
{
    stopifnot(inherits(ts, "treeseq"))
    stopifnot(inherits(obj, "discrete") && inherits(obj, "mpr"))
    stopifnot(!is.unsorted(times))
    stopifnot(all(times >= 0))
    num_states = nrow(cost_matrix)
    num_samples = treeseq_num_samples(ts)
    if (missing(state_sets))
        state_sets = 1:num_states
    if (missing(sample_sets))
        sample_sets = 1:num_samples
    stopifnot(is.integer(state_sets))
    stopifnot(is.integer(sample_sets))
    stopifnot(all(state_sets > 0))
    stopifnot(all(sample_sets >= 0))
    stopifnot(all(tabulate(state_sets) > 0))
    stopifnot(all(tabulate(sample_sets) > 0))
    stopifnot(length(state_sets) == num_states)
    stopifnot(length(sample_sets) == num_samples)
    e = treeseq_discrete_mpr_edge_history(
        ts, obj, cost_matrix, adjacency_matrix, FALSE)
    .Call(
        C_treeseq_discrete_mpr_ancestry
        , ts@tree
        , attr(e, "path.offset")
        , e$state_id
        , e$time
        , attr(e, "node.state")
        , max(state_sets)
        , state_sets - 1L
        , max(sample_sets)
        , sample_sets - 1L
        , times
    )
}


#' Spatiotemporal ancestry flux coefficients
#'
#' @param ts An object of class \code{treeseq}.
#' @param obj The result of \code{treeseq_discrete_mpr}.
#' @param cost_matrix A numeric matrix giving the state-to-state migration
#' costs. All costs should be non-negative and any non-zero costs on the 
#' diagonal will be ignored.
#' @param adjacency_matrix A binary matrix specifying geographic connections.
#' Geographic states \code{i} and \code{j} are connected if the corresponding
#' entry in \code{adjacency_matrix} is 1.
#' @param times A set of times ordered from present to past giving the intervals
#' used for calculating flux coefficients.
#' @param state_sets An integer vector with one entry for each geographic
#' state specifying the set to which that state belongs. Sets must be numbered
#' beginning from 1.
#' @param sample_sets An integer vector with one entry for each sample
#' specifying the set to which that sample belongs. Sets must be numbered
#' beginning from 1. A value of 0 may be used to indicate that the sample
#' belongs to none of the sets and hence is to be excluded from the calculation.
treeseq_discrete_mpr_ancestry_flux = function(ts, obj, cost_matrix,
    adjacency_matrix, times, state_sets, sample_sets)
{
    stopifnot(inherits(ts, "treeseq"))
    stopifnot(inherits(obj, "discrete") && inherits(obj, "mpr"))
    stopifnot(!is.unsorted(times))
    stopifnot(all(times >= 0))
    num_states = nrow(cost_matrix)
    num_samples = treeseq_num_samples(ts)
    if (missing(state_sets))
        state_sets = 1:num_states
    if (missing(sample_sets))
        sample_sets = rep(1L, num_samples)
    stopifnot(is.integer(state_sets))
    stopifnot(is.integer(sample_sets))
    stopifnot(all(state_sets > 0))
    stopifnot(all(sample_sets >= 0))
    stopifnot(all(tabulate(state_sets) > 0))
    stopifnot(all(tabulate(sample_sets) > 0))
    stopifnot(length(state_sets) == num_states)
    stopifnot(length(sample_sets) == num_samples)
    num_state_sets = as.numeric(max(state_sets))
    num_sample_sets = as.numeric(max(sample_sets))
    num_time_bins = length(times) - 1
    storage = num_state_sets * num_state_sets * num_sample_sets * num_time_bins
    if (storage > .Machine$integer.max)
        stop("storage requirements too large")
    e = treeseq_discrete_mpr_edge_history(
        ts, obj, cost_matrix, adjacency_matrix, FALSE)
    .Call(
        C_treeseq_discrete_mpr_ancestry_flux
        , ts@tree
        , attr(e, "path.offset")
        , e$state_id
        , e$time
        , as.integer(num_state_sets)
        , state_sets - 1L
        , as.integer(num_sample_sets)
        , sample_sets - 1L
        , times
    )
}


treeseq_quadratic_mpr = function(ts, sample_locations, use_brlen=FALSE)
{
    stopifnot(inherits(ts, "treeseq"))
    stopifnot(is.matrix(sample_locations))
    stopifnot(!is.null(colnames(sample_locations)))
    stopifnot(colnames(sample_locations)[1] == "node_id")
    storage.mode(sample_locations) = "double"
    N = nrow(treeseq_nodes(ts))
    num_samples = nrow(sample_locations)
    x = matrix(0, ncol(sample_locations)-1L, N)
    sample_ids = sample_locations[, 1] + 1L
    for (i in 1:num_samples)
    {
        sample_id = sample_locations[i, 1]
        x[, sample_id+1L] = sample_locations[i, -1]
    }
    L = .Call(
        C_treeseq_quadratic_mpr
        , ts@treeseq
        , as.integer(use_brlen)
        , x
    )
    names(L) = c("mean_tree_length", "tree_length", "mpr")
    structure(L, class=c("quadratic", "mpr"))
}


treeseq_quadratic_mpr_minimize = function(obj)
{
    stopifnot(inherits(obj, "quadratic") && inherits(obj, "mpr"))
    .Call(C_treeseq_quadratic_mpr_minimize, obj$mpr)
}


treeseq_quadratic_mpr_minimize_discrete = function(obj, sites)
{
    stopifnot(inherits(obj, "quadratic") && inherits(obj, "mpr"))
    stopifnot(is.matrix(sites))
    storage.mode(sites) = "double"
    .Call(
        C_treeseq_quadratic_mpr_minimize_discrete
        , obj$mpr
        , sites)
}


treeseq_linear_mpr = function(ts, sample_locations, use_brlen=FALSE,
    tol=0.01)
{
    stopifnot(inherits(ts, "treeseq"))
    stopifnot(is.matrix(sample_locations))
    stopifnot(!is.null(colnames(sample_locations)))
    stopifnot(tol > 0 && tol < 1)
    stopifnot(colnames(sample_locations)[1] == "node_id")
    storage.mode(sample_locations) = "double"
    N = nrow(treeseq_nodes(ts))
    num_samples = nrow(sample_locations)
    x = matrix(0, ncol(sample_locations)-1L, N)
    sample_ids = sample_locations[, "node_id"] + 1L
    for (i in 1:num_samples)
    {
        sample_id = sample_locations[i, 1]
        x[, sample_id+1L] = sample_locations[i, -1]
    }
    nx = apply(sample_locations[,-1,drop=FALSE], 2, function(s) length(unique(s)))
    L = .Call(
        C_treeseq_linear_mpr
        , ts@treeseq
        , as.integer(use_brlen)
        , x
        , nx
        , tol
    )
    names(L) = c("mean_tree_length", "tree_length", "mpr")
    structure(L, class=c("linear", "mpr"), max_num_breakpoints=nx)
}


treeseq_linear_mpr_minimize = function(obj)
{
    stopifnot(inherits(obj, "linear") && inherits(obj, "mpr"))
    .Call(C_treeseq_linear_mpr_minimize, obj$mpr)
}


treeseq_linear_mpr_minimize_discrete = function(obj, sites)
{
    stopifnot(inherits(obj, "linear") && inherits(obj, "mpr"))
    .Call(
        C_treeseq_linear_mpr_minimize_discrete
        , obj$mpr
        , sites)
}
