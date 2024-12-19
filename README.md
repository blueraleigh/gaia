# gaia: Geographic Ancestry Inference Algorithm

`gaia` is an R package for inferring the geographic locations of genetic ancestors using tree sequences. It implements three approaches to ancestral location reconstruction:

1. Discrete parsimony - for ancestors restricted to a finite set of locations
2. Squared change parsimony - for ancestors in continuous space, minimizing squared distances
3. Linear parsimony - for ancestors in continuous space, minimizing absolute distances

## Installation

Install directly from GitHub:
```r
remotes::install_github("blueraleigh/gaia")
```

Or clone and install locally:
```bash
git clone https://github.com/blueraleigh/gaia
cd gaia
R CMD BUILD .
R CMD INSTALL gaia_*.tar.gz
```

## Quick Start

### Working with Discrete Locations

```r
library(gaia)

# Load your tree sequence (.trees format)
ts <- treeseq_load("path/to/treesequence.trees")

# Define sample locations - each sample must be assigned to a discrete state
# node_id: Tree sequence node IDs (0-based)
# state_id: Location state IDs (1-based)
samples <- data.frame(
  node_id = c(0:10),               # Your sample node IDs
  state_id = c(1,1,2,2,3,1,2,3,2,1,3)  # Your state assignments
)

# Create cost matrix for migrations between states
# - Must be symmetric
# - Diagonal elements (self-transitions) are ignored
# - All costs must be non-negative
num_states <- max(samples$state_id)
costs <- matrix(1, num_states, num_states)  # Default cost of 1 between states
diag(costs) <- 0                           # No cost to stay in same state

# Optional: Define which states can transition to each other
adjacency <- matrix(1, num_states, num_states)  # Default: all transitions allowed
diag(adjacency) <- 0                           # Exclude self-transitions

# Compute migration costs
mpr <- treeseq_discrete_mpr(ts, samples, costs)

# Get optimal state assignments for ancestors
states <- treeseq_discrete_mpr_minimize(mpr)

# Get detailed migration histories
history <- treeseq_discrete_mpr_edge_history(ts, mpr, costs, adjacency)
```

### Working with Continuous Space

```r
# For ancestors in continuous space, provide sample coordinates
samples <- data.frame(
  node_id = c(0:10),          # Sample node IDs
  x = runif(11, 0, 10),       # X coordinates
  y = runif(11, 0, 10)        # Y coordinates
)

# Using squared distance (minimizes sum of squared Euclidean distances)
mpr_quad <- treeseq_quadratic_mpr(ts, samples)
locations_quad <- treeseq_quadratic_mpr_minimize(mpr_quad)

# Using absolute distance (minimizes sum of Manhattan distances)
mpr_lin <- treeseq_linear_mpr(ts, samples)
locations_lin <- treeseq_linear_mpr_minimize(mpr_lin)

# If ancestors must be restricted to specific locations despite using 
# continuous-space reconstruction:
possible_sites <- matrix(c(
  1.5, 2.0,  # Site 1 coordinates
  4.2, 3.1,  # Site 2 coordinates
  6.7, 5.5   # Site 3 coordinates
), ncol=2, byrow=TRUE)

locations_discrete <- treeseq_quadratic_mpr_minimize_discrete(mpr_quad, possible_sites)
```

## Key Functions

- `treeseq_discrete_mpr()` - Discrete state reconstruction
- `treeseq_quadratic_mpr()` - Continuous space reconstruction using squared distances
- `treeseq_linear_mpr()` - Continuous space reconstruction using absolute distances
- `treeseq_discrete_mpr_edge_history()` - Detailed migration histories
- `treeseq_discrete_mpr_ancestry()` - Ancestry coefficients through time
- `treeseq_discrete_mpr_ancestry_flux()` - Migration flux between regions

## Documentation

Access complete function documentation through R's help system:
```r
?treeseq_discrete_mpr
?treeseq_quadratic_mpr
?treeseq_linear_mpr
```

## References

Grundler et al. (2024) A geographic history of human genetic ancestry. bioRxiv doi: 10.1101/2024.03.27.586858

## License

CC-BY 4.0 International