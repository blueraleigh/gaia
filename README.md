# About
The gaia package performs generalized parsimony analysis on the local gene 
trees comprising an ancestral recombination graph (ARG) to infer the geographic 
locations of genetic ancestors. It relies on the [succinct tree sequence format](https://tskit.dev).

# Installation
Install from within an R session:

```r
remotes::install_github("blueraleigh/gaia")
```

Install from the command line:

```bash
git clone https://github.com/blueraleigh/gaia
cd gaia
R CMD BUILD .
R CMD INSTALL <name of tarball>
```

# Key functions

- `treeseq_discrete_mpr`: Discrete state space reconstruction using arbitrary costs
- `treeseq_quadratic_mpr`: Continuous state space reconstruction using squared distances  
- `treeseq_linear_mpr`: Continuous state space reconstruction using absolute distances  
- `treeseq_discrete_mpr_edge_history`: Detailed migration histories  
- `treeseq_discrete_mpr_ancestry`: Ancestry coefficients through time  
- `treeseq_discrete_mpr_ancestry_flux`: Migration flux between regions  

# Documentation

Documentation is a work in progress and is available through the R help system (e.g., `?treeseq_discrete_mpr`).
