#' Compute minimum migration costs for geographic ancestry reconstruction
#'
#' @description
#' Uses generalized (Sankoff) parsimony to compute the minimum migration costs needed 
#' to explain sampled geographic locations under different possible ancestral state 
#' assignments. For each ancestral node and potential geographic state, calculates 
#' the total migration cost required to explain all sampled descendant locations when 
#' that ancestor is assigned to that state.
#'
#' @param ts A \code{treeseq} object, typically loaded via \code{\link{treeseq_load}}
#' @param sample_locations An integer matrix with two columns:
#'   \describe{
#'     \item{node_id}{Node identifiers for sampled genomes (0-based indexing)}
#'     \item{state_id}{Geographic state assignments for samples (1-based indexing)}
#'   }
#' @param cost_matrix A symmetric numeric matrix where entry [i,j] gives the migration
#'   cost between states i and j. Must have non-negative values. Diagonal elements 
#'   (representing costs of remaining in the same state) are ignored.
#' @param use_brlen Logical indicating whether to scale migration costs by inverse
#'   branch lengths (TRUE) or treat all branches equally (FALSE, default)
#'
#' @return A list with class \code{c("discrete", "mpr")} containing:
#'   \describe{
#'     \item{mean_tree_length}{Genome-wide average migration cost}
#'     \item{tree_length}{Vector giving the migration cost for each local tree}
#'     \item{mpr}{Matrix where entry [i,j] gives the minimum migration cost needed
#'       to explain sample locations when node i is placed in state j}
#'   }
#'
#' @details
#' This function implements generalized parsimony on tree sequences using the 
#' algorithm of Sankoff and Rousseau (1975). It calculates minimum migration cost 
#' functions for each node by proceeding up the tree sequence from samples to 
#' ancestral nodes. These costs are then used to determine optimal geographic 
#' locations using \code{treeseq_discrete_mpr_minimize}.
#'
#' For each non-sample node, the function computes the minimum total migration cost 
#' needed to explain the geographic locations of all sampled descendants under each 
#' possible geographic state assignment for that node. These costs are averaged 
#' across all base pairs in the genome that inherit from the node.
#'
#' Migration costs are specified through the cost_matrix parameter. The matrix must
#' be symmetric and have non-negative values. Entry [i,j] gives the cost of migrating
#' from state i to state j. Diagonal elements are ignored since they represent the
#' cost of remaining in the same state, which is always 0.
#'
#' @seealso
#' \code{\link{treeseq_discrete_mpr_minimize}} for determining optimal state assignments
#' \code{\link{treeseq_discrete_mpr_edge_history}} for reconstructing migration histories
#' \code{\link{treeseq_quadratic_mpr}} for continuous-space reconstruction
#' \code{\link{treeseq_linear_mpr}} for linear cost reconstruction
#'
#' @references
#' Sankoff, D. and Rousseau, P. (1975) Locating the vertices of a Steiner tree in 
#' arbitrary space. \emph{Mathematical Programming}, 9:240-246.
#'
#' @examples
#' # Load example tree sequence
#' ts <- treeseq_load(system.file("extdata", "example.trees", package="gaia"))
#'
#' # Define sample locations (3 samples in 2 states)
#' samples <- data.frame(
#'   node_id = 0:2,        # Node IDs use 0-based indexing
#'   state_id = c(1,2,1)   # States use 1-based indexing
#' )
#'
#' # Create cost matrix for 2 states
#' costs <- matrix(c(
#'   0, 1,  # Cost of 1 to move between states
#'   1, 0
#' ), 2, 2)
#'
#' # Compute migration costs
#' mpr_costs <- treeseq_discrete_mpr(ts, samples, costs)
#'
#' # Find optimal states using the costs
#' states <- treeseq_discrete_mpr_minimize(mpr_costs)
#'
#' @export
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


#' Determine optimal geographic states from minimum migration costs
#'
#' @description
#' Uses the migration costs computed by \code{treeseq_discrete_mpr} to identify the 
#' optimal geographic state for each ancestral node - the state that minimizes the 
#' total migration cost needed to explain sample locations. When multiple states 
#' achieve the minimum cost, one is chosen randomly.
#'
#' @param obj Result object from \code{\link{treeseq_discrete_mpr}}
#' @param index1 Logical indicating whether returned state assignments should use
#'   1-based indexing (TRUE, default) or 0-based indexing (FALSE)
#'
#' @return An integer vector giving the optimal geographic state assignment for
#'   each node in the tree sequence. States are indexed from 1 (if index1=TRUE)
#'   or 0 (if index1=FALSE).
#'
#' @details
#' For each node, this function examines the migration costs computed by 
#' \code{treeseq_discrete_mpr} to identify the geographic state(s) that achieve
#' the minimum total cost. When multiple states achieve the minimum cost (equally 
#' parsimonious reconstructions), one state is selected randomly using uniform 
#' probability.
#'
#' The function maintains consistent indexing with the input data by default
#' (index1=TRUE), matching the 1-based state indexing used in sample_locations
#' when calling \code{treeseq_discrete_mpr}. Set index1=FALSE to use 0-based
#' indexing, which may be more convenient when working with tree sequence node IDs.
#'
#' @seealso
#' \code{\link{treeseq_discrete_mpr}} for computing migration costs
#' \code{\link{treeseq_discrete_mpr_edge_history}} for detailed migration histories
#'
#' @examples
#' # Load example tree sequence
#' ts <- treeseq_load(system.file("extdata", "example.trees", package="gaia"))
#'
#' # Define sample locations (3 samples in 2 states)
#' samples <- data.frame(
#'   node_id = 0:2,
#'   state_id = c(1,2,1)
#' )
#'
#' # Create cost matrix for 2 states
#' costs <- matrix(c(
#'   0, 1,
#'   1, 0
#' ), 2, 2)
#'
#' # Compute migration costs
#' mpr_costs <- treeseq_discrete_mpr(ts, samples, costs)
#'
#' # Find optimal states (1-based indexing)
#' states <- treeseq_discrete_mpr_minimize(mpr_costs)
#'
#' # Find optimal states (0-based indexing)
#' states0 <- treeseq_discrete_mpr_minimize(mpr_costs, index1=FALSE)
#'
#' @export
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

#' Compute quadratic parsimony costs for continuous geographic reconstruction
#'
#' @description
#' Uses squared change parsimony to compute the minimum squared distance costs needed
#' to explain sampled geographic coordinates under different possible ancestral 
#' locations. For each ancestral node, calculates parameters of a quadratic function 
#' that gives the total squared distance cost required to explain all sampled 
#' descendant locations when that ancestor is placed at any point in continuous space.
#'
#' @param ts A \code{treeseq} object, typically loaded via \code{\link{treeseq_load}}
#' @param sample_locations A numeric matrix with columns:
#'   \describe{
#'     \item{node_id}{Node identifiers for sampled genomes (0-based indexing)}
#'     \item{Additional columns}{Coordinate values for each spatial dimension (x, y, etc.)}
#'   }
#' @param use_brlen Logical indicating whether to scale distances by inverse
#'   branch lengths (TRUE) or treat all branches equally (FALSE, default)
#'
#' @return A list with class \code{c("quadratic", "mpr")} containing:
#'   \describe{
#'     \item{mean_tree_length}{Genome-wide average squared distance cost}
#'     \item{tree_length}{Vector giving the squared distance cost for each local tree}
#'     \item{mpr}{Array containing quadratic function parameters for each node that
#'       define the squared distance cost as a function of ancestor location}
#'   }
#'
#' @details
#' This function implements squared change parsimony reconstruction on tree sequences
#' following Maddison (1991). For each non-sample node, it computes parameters of a
#' quadratic function that gives the minimum total squared distance cost needed to
#' explain the geographic coordinates of all sampled descendants when that ancestor
#' is placed at any location in continuous space. These parameters are then used to
#' determine optimal geographic locations using \code{treeseq_quadratic_mpr_minimize}.
#'
#' The method assumes ancestral locations can be anywhere in continuous space, rather
#' than restricted to discrete states like \code{treeseq_discrete_mpr}. It minimizes
#' the sum of squared Euclidean distances between connected nodes in the tree.
#'
#' @seealso
#' \code{\link{treeseq_quadratic_mpr_minimize}} for determining optimal locations
#' \code{\link{treeseq_discrete_mpr}} for discrete state reconstruction
#' \code{\link{treeseq_linear_mpr}} for linear distance reconstruction
#'
#' @references
#' Maddison, W.P. (1991) Square-Change Parsimony Reconstructions of Ancestral
#' States for Continuous-Valued Characters on a Phylogenetic Tree. 
#' \emph{Systematic Zoology}, 40(3):304-314.
#'
#' @examples
#' # Load example tree sequence
#' ts <- treeseq_load(system.file("extdata", "example.trees", package="gaia"))
#'
#' # Define sample locations (3 samples in 2D space)
#' samples <- data.frame(
#'   node_id = 0:2,
#'   x = c(0.5, 1.25, 2.0),
#'   y = c(2.7, 0.41, 1.5)
#' )
#'
#' # Compute squared distance costs
#' mpr_costs <- treeseq_quadratic_mpr(ts, samples)
#'
#' # Find optimal locations
#' locations <- treeseq_quadratic_mpr_minimize(mpr_costs)
#'
#' @export
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

#' Determine optimal continuous locations from quadratic parsimony costs
#'
#' @description
#' Uses the squared distance costs computed by \code{treeseq_quadratic_mpr} to 
#' identify the optimal geographic location for each ancestral node - the location that 
#' minimizes the total squared distance cost needed to explain sample locations. For 
#' each node, the optimal location is found by minimizing the node's quadratic cost 
#' function.
#'
#' @param obj Result object from \code{\link{treeseq_quadratic_mpr}}
#'
#' @return A numeric matrix where each row gives the optimal geographic coordinates
#'   for a node. Columns correspond to spatial dimensions (x, y, etc.) in the same
#'   order as the input sample_locations.
#'
#' @details
#' For each node, this function minimizes the quadratic cost function computed by
#' \code{treeseq_quadratic_mpr} to find the location in continuous space that
#' achieves the minimum total squared distance cost. The quadratic form ensures
#' a unique minimum for each node.
#'
#' The optimal location for a node represents the spatial position that minimizes
#' the sum of squared Euclidean distances to its connected nodes (parent and
#' children) in the tree, weighted by the proportion of the genome inherited
#' through each connection.
#'
#' @seealso
#' \code{\link{treeseq_quadratic_mpr}} for computing squared distance costs
#' \code{\link{treeseq_quadratic_mpr_minimize_discrete}} for discrete location assignments
#'
#' @examples
#' # Load example tree sequence
#' ts <- treeseq_load(system.file("extdata", "example.trees", package="gaia"))
#'
#' # Define sample locations (3 samples in 2D space)
#' samples <- data.frame(
#'   node_id = 0:2,
#'   x = c(0.5, 1.25, 2.0),
#'   y = c(2.7, 0.41, 1.5)
#' )
#'
#' # Compute squared distance costs
#' mpr_costs <- treeseq_quadratic_mpr(ts, samples)
#'
#' # Find optimal continuous locations
#' locations <- treeseq_quadratic_mpr_minimize(mpr_costs)
#'
#' @export
treeseq_quadratic_mpr_minimize = function(obj)
{
    stopifnot(inherits(obj, "quadratic") && inherits(obj, "mpr"))
    .Call(C_treeseq_quadratic_mpr_minimize, obj$mpr)
}

#' Assign ancestral nodes to discrete locations using quadratic parsimony costs
#'
#' @description
#' Uses the squared distance costs computed by \code{treeseq_quadratic_mpr} to assign 
#' each ancestral node to the discrete location (from a provided set) that minimizes 
#' the total squared distance cost needed to explain sample locations. Unlike 
#' \code{treeseq_quadratic_mpr_minimize}, which allows locations anywhere in continuous 
#' space, this function restricts assignments to a specified set of discrete locations.
#'
#' @param obj Result object from \code{\link{treeseq_quadratic_mpr}}
#' @param sites A numeric matrix where each row represents a possible location,
#'   with columns corresponding to spatial dimensions (x, y, etc.) in the same
#'   order as used in the original sample_locations
#'
#' @return An integer vector giving the index of the optimal site from the sites
#'   matrix for each node in the tree sequence. Indices are 1-based.
#'
#' @details
#' For each node, this function evaluates the quadratic cost function computed by
#' \code{treeseq_quadratic_mpr} at each candidate location provided in the sites
#' matrix. The location achieving the minimum cost is selected as optimal for
#' that node.
#'
#' This function is useful when ancestral locations must be restricted to a
#' discrete set of possibilities (e.g., habitat patches, sampling locations)
#' despite using a continuous-space reconstruction method.
#'
#' @seealso
#' \code{\link{treeseq_quadratic_mpr}} for computing squared distance costs
#' \code{\link{treeseq_quadratic_mpr_minimize}} for continuous location assignments
#'
#' @examples
#' # Load example tree sequence
#' ts <- treeseq_load(system.file("extdata", "example.trees", package="gaia"))
#'
#' # Define sample locations (3 samples in 2D space)
#' samples <- data.frame(
#'   node_id = 0:2,
#'   x = c(0.5, 1.25, 2.0),
#'   y = c(2.7, 0.41, 1.5)
#' )
#'
#' # Define possible ancestral locations
#' sites <- matrix(c(
#'   0.5, 2.7,  # Site 1
#'   1.25, 0.41,  # Site 2
#'   2.0, 1.5  # Site 3
#' ), ncol=2, byrow=TRUE)
#'
#' # Compute squared distance costs
#' mpr_costs <- treeseq_quadratic_mpr(ts, samples)
#'
#' # Find optimal discrete locations
#' location_indices <- treeseq_quadratic_mpr_minimize_discrete(mpr_costs, sites)
#'
#' @export
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

#' Compute linear parsimony costs for continuous geographic reconstruction
#'
#' @description
#' Uses linear parsimony to compute the minimum absolute distance costs needed
#' to explain sampled geographic coordinates under different possible ancestral 
#' locations. For each ancestral node, calculates parameters of a piecewise linear 
#' function that gives the total absolute distance cost required to explain all 
#' sampled descendant locations when that ancestor is placed at any point in 
#' continuous space.
#'
#' @param ts A \code{treeseq} object, typically loaded via \code{\link{treeseq_load}}
#' @param sample_locations A numeric matrix with columns:
#'   \describe{
#'     \item{node_id}{Node identifiers for sampled genomes (0-based indexing)}
#'     \item{Additional columns}{Coordinate values for each spatial dimension (x, y, etc.)}
#'   }
#' @param use_brlen Logical indicating whether to scale distances by inverse
#'   branch lengths (TRUE) or treat all branches equally (FALSE, default)
#' @param tol Numeric tolerance for merging breakpoints in piecewise linear functions
#'   (default: 0.01)
#'
#' @return A list with class \code{c("linear", "mpr")} containing:
#'   \describe{
#'     \item{mean_tree_length}{Genome-wide average absolute distance cost}
#'     \item{tree_length}{Vector giving the absolute distance cost for each local tree}
#'     \item{mpr}{Array containing piecewise linear function parameters for each node
#'       that define the absolute distance cost as a function of ancestor location}
#'   }
#'
#' @details
#' This function implements linear (Manhattan distance) parsimony reconstruction on
#' tree sequences following the approach of Csűrös (2008). For each non-sample node,
#' it computes parameters of a piecewise linear function that gives the minimum total
#' absolute distance cost needed to explain the geographic coordinates of all sampled
#' descendants when that ancestor is placed at any location in continuous space.
#' These parameters are then used to determine optimal geographic locations using
#' \code{treeseq_linear_mpr_minimize}.
#'
#' The method assumes ancestral locations can be anywhere in continuous space, rather
#' than restricted to discrete states like \code{treeseq_discrete_mpr}. It minimizes
#' the sum of absolute (L1) distances between connected nodes in the tree.
#'
#' @seealso
#' \code{\link{treeseq_linear_mpr_minimize}} for determining optimal locations
#' \code{\link{treeseq_discrete_mpr}} for discrete state reconstruction
#' \code{\link{treeseq_quadratic_mpr}} for squared distance reconstruction
#'
#' @references
#' Csűrös, M. (2008) Ancestral Reconstruction by Asymmetric Wagner Parsimony over
#' Continuous Characters and Squared Parsimony over Distributions. In: Nelson, C.E.,
#' Vialette, S. (eds) Comparative Genomics. RECOMB-CG 2008. Lecture Notes in
#' Computer Science, vol 5267. Springer, Berlin, Heidelberg.
#'
#' @examples
#' # Load example tree sequence
#' ts <- treeseq_load(system.file("extdata", "example.trees", package="gaia"))
#'
#' # Define sample locations (3 samples in 2D space)
#' samples <- data.frame(
#'   node_id = 0:2,
#'   x = c(0.5, 1.25, 2.0),
#'   y = c(2.7, 0.41, 1.5)
#' )
#'
#' # Compute absolute distance costs
#' mpr_costs <- treeseq_linear_mpr(ts, samples)
#'
#' # Find optimal locations
#' locations <- treeseq_linear_mpr_minimize(mpr_costs)
#'
#' @export
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

#' Determine optimal continuous locations from linear parsimony costs
#'
#' @description
#' Uses the absolute distance costs computed by \code{treeseq_linear_mpr} to identify 
#' the optimal geographic location for each ancestral node - the location that 
#' minimizes the total absolute distance cost needed to explain sample locations. For 
#' each node, the optimal location is found by minimizing the node's piecewise linear 
#' cost function.
#'
#' @param obj Result object from \code{\link{treeseq_linear_mpr}}
#'
#' @return A numeric matrix where each row gives the optimal geographic coordinates
#'   for a node. Columns correspond to spatial dimensions (x, y, etc.) in the same
#'   order as the input sample_locations.
#'
#' @details
#' For each node, this function minimizes the piecewise linear cost function 
#' computed by \code{treeseq_linear_mpr} to find the location in continuous space 
#' that achieves the minimum total absolute distance cost. Due to the L1 (Manhattan) 
#' distance metric used, optimal locations are often at breakpoints of the piecewise 
#' linear function.
#'
#' The optimal location for a node represents the spatial position that minimizes
#' the sum of absolute distances to its connected nodes (parent and children) in 
#' the tree, weighted by the proportion of the genome inherited through each 
#' connection.
#'
#' @seealso
#' \code{\link{treeseq_linear_mpr}} for computing absolute distance costs
#' \code{\link{treeseq_linear_mpr_minimize_discrete}} for discrete location assignments
#'
#' @examples
#' # Load example tree sequence
#' ts <- treeseq_load(system.file("extdata", "example.trees", package="gaia"))
#'
#' # Define sample locations (3 samples in 2D space)
#' samples <- data.frame(
#'   node_id = 0:2,
#'   x = c(0.5, 1.25, 2.0),
#'   y = c(2.7, 0.41, 1.5)
#' )
#'
#' # Compute absolute distance costs
#' mpr_costs <- treeseq_linear_mpr(ts, samples)
#'
#' # Find optimal continuous locations
#' locations <- treeseq_linear_mpr_minimize(mpr_costs)
#'
#' @export
treeseq_linear_mpr_minimize = function(obj)
{
    stopifnot(inherits(obj, "linear") && inherits(obj, "mpr"))
    .Call(C_treeseq_linear_mpr_minimize, obj$mpr)
}

#' Assign ancestral nodes to discrete locations using linear parsimony costs
#'
#' @description
#' Uses the absolute distance costs computed by \code{treeseq_linear_mpr} to assign 
#' each ancestral node to the discrete location (from a provided set) that minimizes 
#' the total absolute distance cost needed to explain sample locations. Unlike 
#' \code{treeseq_linear_mpr_minimize}, which allows locations anywhere in continuous 
#' space, this function restricts assignments to a specified set of discrete locations.
#'
#' @param obj Result object from \code{\link{treeseq_linear_mpr}}
#' @param sites A numeric matrix where each row represents a possible location,
#'   with columns corresponding to spatial dimensions (x, y, etc.) in the same
#'   order as used in the original sample_locations
#'
#' @return An integer vector giving the index of the optimal site from the sites
#'   matrix for each node in the tree sequence. Indices are 1-based.
#'
#' @details
#' For each node, this function evaluates the piecewise linear cost function 
#' computed by \code{treeseq_linear_mpr} at each candidate location provided in 
#' the sites matrix. The location achieving the minimum cost is selected as optimal 
#' for that node.
#'
#' This function is useful when ancestral locations must be restricted to a
#' discrete set of possibilities (e.g., habitat patches, sampling locations)
#' despite using a continuous-space reconstruction method.
#'
#' @seealso
#' \code{\link{treeseq_linear_mpr}} for computing absolute distance costs
#' \code{\link{treeseq_linear_mpr_minimize}} for continuous location assignments
#'
#' @examples
#' # Load example tree sequence
#' ts <- treeseq_load(system.file("extdata", "example.trees", package="gaia"))
#'
#' # Define sample locations (3 samples in 2D space)
#' samples <- data.frame(
#'   node_id = 0:2,
#'   x = c(0.5, 1.25, 2.0),
#'   y = c(2.7, 0.41, 1.5)
#' )
#'
#' # Define possible ancestral locations
#' sites <- matrix(c(
#'   0.5, 2.7,  # Site 1
#'   1.25, 0.41,  # Site 2
#'   2.0, 1.5  # Site 3
#' ), ncol=2, byrow=TRUE)
#'
#' # Compute absolute distance costs
#' mpr_costs <- treeseq_linear_mpr(ts, samples)
#'
#' # Find optimal discrete locations
#' location_indices <- treeseq_linear_mpr_minimize_discrete(mpr_costs, sites)
#'
#' @export
treeseq_linear_mpr_minimize_discrete = function(obj, sites)
{
    stopifnot(inherits(obj, "linear") && inherits(obj, "mpr"))
    .Call(
        C_treeseq_linear_mpr_minimize_discrete
        , obj$mpr
        , sites)
}
