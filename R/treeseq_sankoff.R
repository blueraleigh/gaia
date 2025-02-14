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
#' # Load tree sequence
#' ts = treeseq_load(system.file("extdata", "test.trees", package="gaia"))
#' 
#' # Define sample locations - nodes 0 and 2 in one state, node 1 in another
#' state = c(2L,1L,1L)
#' samples = cbind(node_id=0:2, state_id=state)
#' 
#' # Create cost matrix for 2 states
#' # Migration between states has cost 1
#' costs = matrix(c(0, 1, 1, 0), 2, 2)
#' 
#' # Compute migration costs
#' # Will calculate costs at nodes 3-6 based on samples 0-2
#' mpr_costs = treeseq_discrete_mpr(ts, samples, costs)
#' 
#' # Can also use branch lengths to scale costs
#' mpr_costs_bl = treeseq_discrete_mpr(ts, samples, costs, use_brlen=TRUE)
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
#' # Load tree sequence
#' ts = treeseq_load(system.file("extdata", "test.trees", package="gaia"))
#' 
#' # Set up discrete state reconstruction
#' state = c(2L,1L,1L)
#' samples = cbind(node_id=0:2, state_id=state)
#' costs = matrix(c(0, 1, 1, 0), 2, 2)
#' mpr_costs = treeseq_discrete_mpr(ts, samples, costs)
#' 
#' # Get optimal states for all nodes (0-6)
#' states = treeseq_discrete_mpr_minimize(mpr_costs)
#' 
#' # Same with 0-based state indexing
#' states0 = treeseq_discrete_mpr_minimize(mpr_costs, index1=FALSE)
treeseq_discrete_mpr_minimize = function(obj, index1=TRUE)
{
    stopifnot(inherits(obj, "discrete") && inherits(obj, "mpr"))
    stopifnot(is.logical(index1))
    .Call(C_treeseq_discrete_mpr_minimize, obj$mpr, as.integer(index1))
}


#' Sample migration paths for each edge in a tree sequence
#'
#' @description
#' For each edge in a tree sequence, samples a minimum-cost migration path between 
#' the geographic states of parent and child nodes. The path represents the sequence 
#' of state transitions that could have occurred along the branch, consistent with 
#' the parsimony reconstruction.
#'
#' @param ts A \code{treeseq} object
#' @param obj Result object from \code{treeseq_discrete_mpr}
#' @param cost_matrix A symmetric numeric matrix where entry [i,j] gives the migration
#'   cost between states i and j. Must have non-negative values. Diagonal elements
#'   are ignored.
#' @param adjacency_matrix A binary matrix specifying allowed transitions between states.
#'   Entry [i,j] should be 1 if direct transitions are allowed between states i and j,
#'   0 otherwise. Must be symmetric.
#' @param index1 Logical indicating whether state IDs should use 1-based indexing
#'   (TRUE, default) or 0-based indexing (FALSE)
#'
#' @return A data frame where each row represents a state in a migration path:
#'   \describe{
#'     \item{edge_id}{ID of the edge this state belongs to}
#'     \item{state_id}{Geographic state at this point in the path}
#'     \item{time}{Time at which this state occurred}
#'   }
#' Also contains attributes:
#'   \describe{
#'     \item{node.state}{Vector of state assignments for all nodes}
#'     \item{path.offset}{Index vector for finding states belonging to each edge}
#'   }
#'
#' @details
#' For each edge in the tree sequence, this function samples one possible migration
#' path between the states assigned to the parent and child nodes by the parsimony
#' reconstruction. The path represents a sequence of state transitions that could
#' have occurred along that branch while achieving the minimum total migration cost.
#'
#' When multiple minimum-cost paths exist between two states, one is chosen randomly.
#' The path always includes at least two points: the states at the beginning and end
#' of the edge. Additional points are added when state transitions occur along the edge.
#'
#' The adjacency matrix controls which direct transitions between states are allowed.
#' This can be used to enforce geographic constraints (e.g., requiring migration
#' through intermediate states).
#'
#' @seealso
#' \code{\link{treeseq_discrete_mpr}} for computing the parsimony reconstruction
#' \code{\link{treeseq_discrete_mpr_minimize}} for assigning states to nodes
#'
#' @examples
#' # Load tree sequence
#' ts = treeseq_load(system.file("extdata", "test.trees", package="gaia"))
#'
#' # Define sample states
#' state = c(2L,1L,1L)
#' samples = cbind(node_id=0:2, state_id=state)
#'
#' # Create cost matrix for 2 states
#' costs = matrix(c(0, 1, 1, 0), 2, 2)
#'
#' # Create adjacency matrix allowing transitions between states
#' adjacency = matrix(1, 2, 2)
#' # No self-transitions allowed in adjacency
#' diag(adjacency) = 0
#'
#' # Compute migration costs
#' mpr_costs = treeseq_discrete_mpr(ts, samples, costs)
#'
#' # Get migration histories for each edge
#' # This will show state changes along branches in all three trees
#' history = treeseq_discrete_mpr_edge_history(ts, mpr_costs, costs, adjacency)
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


#' Compute geographic ancestry coefficients through time
#'
#' @description
#' For each combination of sample subset and geographic region, calculates the 
#' proportion of genetic material inherited from ancestors in that region at 
#' different points in time. This allows tracking how the geographic distribution 
#' of ancestry changes as we move backwards in time.
#'
#' @param ts A \code{treeseq} object
#' @param obj Result object from \code{treeseq_discrete_mpr}
#' @param cost_matrix Symmetric numeric matrix of migration costs between states
#' @param adjacency_matrix Binary matrix specifying allowed transitions between states
#' @param times Numeric vector of time points at which to calculate ancestry,
#'   must be ordered from present (0) to past
#' @param state_sets Integer vector grouping states into regions. Length must match
#'   number of states, values indicate region membership (1-based)
#' @param sample_sets Integer vector grouping samples into subsets. Length must match
#'   number of samples, values indicate subset membership (1-based). Use 0 to exclude
#'   samples.
#'
#' @return A 3-dimensional array with dimensions:
#'   [region, sample_subset, time_point]
#' Values represent the proportion of genetic material that sample_subset inherits
#' from ancestors in region at each time_point.
#'
#' @details
#' This function traces genetic material backwards in time to determine its
#' geographic location at different time points. For each sample subset, time point,
#' and geographic region, it calculates the fraction of the genome that was located
#' in that region at that time.
#'
#' States can be grouped into regions using state_sets. For example, multiple
#' states might represent different locations within the same continent. Similarly,
#' samples can be grouped into subsets using sample_sets, allowing calculation
#' of ancestry coefficients for different populations or sampling locations.
#'
#' The function uses the migration paths sampled by edge_history to determine
#' locations of genetic material through time. Different random samplings of
#' migration paths may give slightly different results when multiple equally
#' parsimonious paths exist.
#'
#' @seealso
#' \code{\link{treeseq_discrete_mpr_edge_history}} for the underlying migration paths
#' \code{\link{treeseq_discrete_mpr_ancestry_flux}} for tracking migration between regions
#'
#' @examples
#' # Load tree sequence
#' ts = treeseq_load(system.file("extdata", "test.trees", package="gaia"))
#'
#' # Set up states, costs, and adjacency
#' state = c(2L,1L,1L)
#' samples = cbind(node_id=0:2, state_id=state)
#' costs = matrix(c(0,1,1,0), 2, 2)
#' adjacency = matrix(1, 2, 2)
#' diag(adjacency) = 0
#'
#' # Compute base MPR
#' mpr_costs = treeseq_discrete_mpr(ts, samples, costs)
#'
#' # Define timepoints to examine ancestry
#' # Must cover node times from 0 to 1.0
#' times = seq(0, 1.0, by=0.2)
#'
#' # Compute ancestry coefficients
#' ancestry = treeseq_discrete_mpr_ancestry(
#'   ts, mpr_costs, costs, adjacency, times, 
#'   state_sets=1:2,     # Keep states separate
#'   sample_sets=1:3     # Keep samples separate
#' )
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


#' Compute migration flux between geographic regions through time
#'
#' @description
#' Tracks the movement of genetic material between geographic regions over time by
#' computing the proportion of each sample subset's genome that migrated between
#' each pair of regions during different time intervals.
#'
#' @param ts A \code{treeseq} object
#' @param obj Result object from \code{treeseq_discrete_mpr}
#' @param cost_matrix Symmetric numeric matrix of migration costs between states
#' @param adjacency_matrix Binary matrix specifying allowed transitions between states
#' @param times Numeric vector of time points defining intervals for flux calculation,
#'   must be ordered from present (0) to past
#' @param state_sets Integer vector grouping states into regions. Length must match
#'   number of states, values indicate region membership (1-based)
#' @param sample_sets Integer vector grouping samples into subsets. Length must match
#'   number of samples, values indicate subset membership (1-based). Use 0 to exclude
#'   samples.
#'
#' @return A 4-dimensional array with dimensions:
#'   [source_region, dest_region, sample_subset, time_interval]
#' Values represent the proportion of sample_subset's genome that moved from
#' source_region to dest_region during each time_interval.
#'
#' @details
#' This function extends ancestry coefficient calculations by explicitly tracking
#' migrations between regions. For each time interval, it identifies genetic material
#' that changed regions during that interval and records the source and destination
#' regions.
#'
#' States can be grouped into regions using state_sets. For example, multiple
#' states might represent different locations within the same continent. Similarly,
#' samples can be grouped into subsets using sample_sets, allowing calculation
#' of migration flux for different populations or sampling locations.
#'
#' The function uses the migration paths sampled by edge_history to determine
#' when migrations occurred. Different random samplings of migration paths may
#' give slightly different results when multiple equally parsimonious paths exist.
#'
#' Time intervals are defined by consecutive pairs of values in the times parameter.
#' The flux value for an interval represents all migrations that occurred between
#' the start and end of that interval.
#'
#' @seealso
#' \code{\link{treeseq_discrete_mpr_edge_history}} for the underlying migration paths
#' \code{\link{treeseq_discrete_mpr_ancestry}} for static ancestry proportions
#'
#' @examples
#' # Load tree sequence
#' ts = treeseq_load(system.file("extdata", "test.trees", package="gaia"))
#'
#' # Set up states, costs, and adjacency 
#' state = c(2L,1L,1L)
#' samples = cbind(node_id=0:2, state_id=state)
#' costs = matrix(c(0,1,1,0), 2, 2)
#' adjacency = matrix(1, 2, 2)
#' diag(adjacency) = 0
#'
#' # Compute base MPR
#' mpr_costs = treeseq_discrete_mpr(ts, samples, costs)
#'
#' # Define timepoints for flux intervals
#' # Times should span node times (0-1.0)
#' times = seq(0, 1.0, by=0.2)
#'
#' # Compute migration flux between states over time
#' flux = treeseq_discrete_mpr_ancestry_flux(
#'   ts, mpr_costs, costs, adjacency, times,
#'   state_sets=1:2,     # Keep states separate
#'   sample_sets=1:3     # Keep samples separate
#' )
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
#' # Load tree sequence
#' ts = treeseq_load(system.file("extdata", "test.trees", package="gaia"))
#' 
#' # Define the x and y coordinates (locations)
#' x = c(0, 2, 1)
#' y = c(0, 0, 1)
#' 
#' samples = cbind(node_id=0:2, x=x, y=y)
#' 
#' # Compute squared distance costs for nodes 0-6
#' mpr_costs = treeseq_quadratic_mpr(ts, samples)
#' 
#' # Can use branch lengths to scale distances by time
#' mpr_costs_bl = treeseq_quadratic_mpr(ts, samples, use_brlen=TRUE)
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
#' # Load tree sequence
#' ts = treeseq_load(system.file("extdata", "test.trees", package="gaia"))
#' 
#' # Define the x and y coordinates (locations)
#' x = c(0, 2, 1)
#' y = c(0, 0, 1)
#' 
#' samples = cbind(node_id=0:2, x=x, y=y)
#' 
#' mpr_costs = treeseq_quadratic_mpr(ts, samples)
#' 
#' # Find optimal continuous locations for all nodes 0-6
#' locations = treeseq_quadratic_mpr_minimize(mpr_costs)
#' # Results should place internal nodes 3-6 at weighted averages
#' # of their descendants' positions
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
#' # Load tree sequence
#' ts = treeseq_load(system.file("extdata", "test.trees", package = "gaia"))
#' 
#' # Define the x and y coordinates (locations)
#' x = c(0, 2, 1)
#' y = c(0, 0, 1)
#' 
#' samples = cbind(node_id=0:2, x=x, y=y)
#' 
#' mpr_costs = treeseq_quadratic_mpr(ts, samples)
#' 
#' # Define possible locations matching sample coordinates
#' sites = matrix(c(
#'   0, 0,    # Site 1: matches node 0
#'   2, 0,    # Site 2: matches node 1
#'   1, 1     # Site 3: matches node 2
#' ), ncol=2, byrow=TRUE)
#' 
#' # Get optimal site assignment for all nodes 0-6
#' location_indices = treeseq_quadratic_mpr_minimize_discrete(mpr_costs, sites)
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
#' # Load tree sequence
#' ts = treeseq_load(system.file("extdata", "test.trees", package="gaia"))
#' 
#' # Define the x and y coordinates (locations)
#' x = c(0, 2, 1)
#' y = c(0, 0, 1)
#' 
#' samples = cbind(node_id=0:2, x=x, y=y)
#' 
#' # Compute absolute distance costs for nodes 0-6
#' mpr_costs = treeseq_linear_mpr(ts, samples)
#' 
#' # Use branch lengths to scale distances
#' mpr_costs_bl = treeseq_linear_mpr(ts, samples, use_brlen=TRUE)
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
#' # Load tree sequence
#' ts = treeseq_load(system.file("extdata", "test.trees", package="gaia"))
#' 
#' # Define the x and y coordinates (locations)
#' x = c(0, 2, 1)
#' y = c(0, 0, 1)
#' 
#' samples = cbind(node_id=0:2, x=x, y=y)
#' 
#' mpr_costs = treeseq_linear_mpr(ts, samples)
#' 
#' # Find optimal continuous locations for all nodes 0-6
#' # Unlike quadratic MPR, optimal locations often at sample coordinates
#' locations = treeseq_linear_mpr_minimize(mpr_costs)
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
#' # Load tree sequence
#' ts = treeseq_load(system.file("extdata", "test.trees", package="gaia"))
#' 
#' # Define the x and y coordinates (locations)
#' x = c(0, 2, 1)
#' y = c(0, 0, 1)
#' 
#' samples = cbind(node_id=0:2, x=x, y=y)
#' 
#' mpr_costs = treeseq_linear_mpr(ts, samples)
#' 
#' # Define possible locations matching sample coordinates
#' sites = matrix(c(
#'   0, 0,    # Site 1: matches node 0
#'   2, 0,    # Site 2: matches node 1 
#'   1, 1     # Site 3: matches node 2
#' ), ncol=2, byrow=TRUE)
#' 
#' # Get optimal site assignment for all nodes 0-6
#' location_indices = treeseq_linear_mpr_minimize_discrete(mpr_costs, sites)
treeseq_linear_mpr_minimize_discrete = function(obj, sites)
{
    stopifnot(inherits(obj, "linear") && inherits(obj, "mpr"))
    .Call(
        C_treeseq_linear_mpr_minimize_discrete
        , obj$mpr
        , sites)
}
