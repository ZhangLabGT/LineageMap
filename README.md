# LineageMap

Reconstructing Spatially Resolved Cell Lineage Trees from Single-Cell Multi-Omics Data

## Overview

**LineageMap** is an R package for reconstructing cell lineage trees by integrating:

- **Lineage barcodes** (CRISPR/Cas9-based or similar)
- **Spatial coordinates** (spatial transcriptomics)
- **Gene expression states** (cell type or differentiation state)

The method uses a hierarchical approach: first building a backbone tree (NJ) via clustering, then refining local subtrees through likelihood-based optimization with spatial and state information.

## Installation

```r
# Install from GitHub
devtools::install_github("ZhangLabGT/LineageMap")
```

## Quick Start

```r
library(LineageMap)

# Build lineage tree
tree <- Build_LineageMap(
  muts = barcode_matrix,          # Cell x barcode matrix
  meta = metadata,                 # Data frame with spatial (x, y) and state columns
  state_lineages = state_network,  # State transition network
  max_Iter = 200,                  # Maximum iterations for local search
  threshold = 0.2,                 # Clustering threshold; could be automatically calculated with otsu function
  backbone_type = "majority",      # Backbone reconstruction method
  lambda1 = 0.05,                  # Weight for gene expression state
  lambda2 = 0.1,                   # Weight for spatial location
  alpha = 1                        # Hyperparameter for spatial likelihood
)

# Visualize the tree
plot(tree)
```

## Input Data Format

### 1. Barcode Matrix (`muts`)

A matrix where:

- **Rows**: cells (with row names as cell IDs, e.g., "cell_1", "cell_2")
- **Columns**: barcode positions
- **Values**: character states (0, 1, 2, ..., or "-" for missing)

```r
#           site_1 site_2 site_3 ...
# cell_1    "0"    "1"    "0"    ...
# cell_2    "0"    "1"    "2"    ...
# cell_3    "1"    "0"    "-"    ...
```

### 2. Metadata (`meta`)

A data frame with spatial and state information:

- Must contain columns: `x`, `y` (spatial coordinates), `state` (cell type/state)
- Row names should match the barcode matrix

```r
#        x      y     state
# cell_1 10.5   20.3  "A"
# cell_2 15.2   18.7  "B"
# cell_3 12.1   25.4  "A"
```

### 3. State Lineages (`state_lineages`)

A data structure defining the state transition network (from root state to leaf states).

## Key Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `threshold` | Clustering threshold for grouping similar cells | 0.2 |
| `backbone_type` | Method for backbone reconstruction: "majority", "NJ_median", or "NJ_mean" | "majority" |
| `max_Iter` | Maximum iterations for local search optimization | 200 |
| `lambda1` | Weight for gene expression state likelihood | 0.05 |
| `lambda2` | Weight for spatial location likelihood | 0.1 |
| `alpha` | Hyperparameter for spatial likelihood | 1 |

## Advanced Usage

For large datasets, use the parallelized version:

```r
tree <- Build_LineageMap_parallel(
  muts = barcode_matrix,
  meta = metadata,
  state_lineages = state_network,
  outer_cores = 12,    # Cores for parallelizing across subgroups
  inner_cores = 5      # Cores for parallelizing within each subgroup
)
```

## Method Overview

1. **Clustering**: Group cells based on barcode similarity using Louvain clustering
2. **Backbone Construction**: Build a coarse tree structure connecting cell groups
3. **Local Refinement**: For each group, perform local search to optimize subtree topology using:
   - Barcode likelihood
   - Spatial proximity (cells closer in space are more likely to be related)
   - State transitions (following the provided state lineage network)
4. **Tree Assembly**: Combine the backbone and refined subtrees into the final lineage tree


## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

