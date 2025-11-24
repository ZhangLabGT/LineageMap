#' Calculate hamming distance with dropouts
#'
#' @param a barcode a
#' @param b barcode b
#' @param dropout dropout pattern
#' @param w_mut weight between 2 mutated barcodes
#' @param w_ref weight between an unmutated state and a mutated state
#' @export
WeightedDropoutHamming <- function(a, b, dropout = "-", w_mut = 1.0, w_ref = 0.5) {
  # a, b are vectors of equal length (character or numeric states)

  valid <- (a != dropout & b != dropout)
  if (sum(valid) == 0) return(0)

  d <- 0
  for (i in which(valid)) {
    if (a[i] == b[i]) {
      next
    } else {
      if (a[i] != "0" & b[i] != "0") {
        d <- d + w_mut      # mutation vs different mutation
      } else {
        d <- d + w_ref      # mutation vs reference
      }
    }
  }
  return(d / sum(valid))  # normalize by valid positions
}

#' Calculate pairwise hamming distance matrix with dropouts
#'
#' @param X barcode matrix
#' @param dropout dropout pattern
#' @param w_mut weight between 2 mutated barcodes
#' @param w_ref weight between an unmutated state and a mutated state
#' @export
DropoutDistMatrix <- function(X, dropout = "-", w_mut = 1.0, w_ref = 0.5) {
  n <- nrow(X)
  dmat <- matrix(0, n, n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      d <- WeightedDropoutHamming(as.character(X[i,]), as.character(X[j,]),
                                  dropout = dropout, w_mut = w_mut, w_ref = w_ref)
      dmat[i,j] <- d
      dmat[j,i] <- d
    }
  }
  return(dmat)
}

#' Build similarity graph + Louvain
#'
#' @param X barcode matrix
#' @param dropout dropout pattern
#' @param threshold sparsity parameter for similarity graph
#' @param resolution clustering resolution for louvain
#' @param dmat pairwise distance matrix
#' @import igraph
#' @export
LouvainCluster <- function(X, dropout = "-", threshold = 0.2,resolution = 1,dmat = NULL) {
  if (is.null(dmat)){
    dmat <- DropoutDistMatrix(X, dropout = dropout)
    n <- nrow(X)
  }
  else{
    n <- nrow(dmat)
  }

  # similarity = 1 - distance
  sim <- 1 - dmat
  sim[is.na(sim)] <- 0

  edges <- which(sim > (1 - threshold) & upper.tri(sim), arr.ind = TRUE)
  default_membership  <- function(n){
    x <- 1:n
    names(x) <- 1:n
    x
  }
  if (nrow(edges) == 0) return(default_membership(n))

  g <- graph_from_data_frame(
    data.frame(
      from = edges[,1],
      to   = edges[,2],
      weight = sim[edges]
    ),
    directed = FALSE,
    vertices = data.frame(name = 1:n)
  )

  # Louvain community detection
  cl <- cluster_louvain(g, weights = E(g)$weight,resolution = resolution)
  return(membership(cl))
}

#' Build similarity graph + Louvain
#'
#' @param X barcode matrix
#' @param groups subgroup ids of cells
#' @param dmat pairwise distance matrix
#' @import data.table
#' @import igraph
#' @export
SubgroupTree <- function(X, groups,dmat = NULL) {
  group_id <- unique(groups)
  dt <- rbindlist(lapply(group_id,function(g){
    cells <- as.numeric(names(groups)[groups == g])
    data.table(group_name = paste0("cells",as.character(min(cells)),"-",as.character(max(cells))),id = g,cells = list(cells))
  }))
  X_unique <- sapply(group_id, function(g) {
    rows <- X[groups == g, , drop = FALSE]
    apply(rows, 2, function(col) {
      # majority vote consensus
      #names(sort(table(col), decreasing = TRUE))[1]

      col <- col[col != "-"]
      if (length(col) == 0) return("-")
      names(which.max(table(col)))
    })
  })
  X_unique <- t(X_unique)
  rownames(X_unique) <- dt$group_name
  dist_h <- DropoutDistMatrix(X_unique)
  Treenj <- nj(dist_h)
  tree_backbone <- multi2di(Treenj)
  tree_backbone$tip.label <- dt$group_name
  #tree_backbone$edge.length <- rep(1, length(tree_backbone$edge.length))
  return(list(tree_backbone = tree_backbone, backbone_meta = dt))
}


#' compute inter-group distances using pairwise distances
#'
#' @param X barcode matrix
#' @param groups subgroup ids of cells
#' @param dropout the dropout character
#' @param w_mut the weight between two mutations
#' @param w_ref weight between an unmutated state and a mutated state
#' @param summary_fn the function used to represent pairwise distance
#' @import data.table
#' @import igraph
#' @export
GroupDistanceMatrix <- function(X, groups, dropout = "-", w_mut = 1, w_ref = 0.5, summary_fn = median) {
  group_ids <- unique(groups)
  G <- length(group_ids)
  D <- matrix(0, G, G, dimnames = list(group_ids, group_ids))
  for (i in seq_len(G)) {
    for (j in i:G) {
      gi <- group_ids[i]; gj <- group_ids[j]
      rows_i <- which(groups == gi)
      rows_j <- which(groups == gj)
      dvals <- c()
      for (r in rows_i) for (s in rows_j) {
        if (r == s) next
        dvals <- c(dvals, WeightedDropoutHamming(as.character(X[r,]), as.character(X[s,]),
                                                 dropout = dropout, w_mut = w_mut, w_ref = w_ref))
      }
      if (length(dvals) == 0) dsum <- 0 else dsum <- summary_fn(dvals)
      D[i,j] <- D[j,i] <- dsum
    }
  }
  as.dist(D)
}

#' Build the backbone tree using NJ on pairwise distances between subgroups
#'
#' @param X barcode matrix
#' @param groups subgroup ids of cells
#' @import data.table
#' @import igraph
#' @export
SubgroupTree_NJ <- function(X, groups, dropout = "-", w_mut = 1, w_ref = 0.5, summary_fn = median) {
  group_id <- unique(groups)
  dt <- rbindlist(lapply(group_id,function(g){
    cells <- as.numeric(names(groups)[groups == g])
    data.table(
      group_name = paste0("cells", min(cells), "-", max(cells)),
      id = g,
      cells = list(cells)
    )
  }))

  # Compute inter-group distances
  dist_h <- GroupDistanceMatrix(X, groups, dropout = dropout, w_mut = w_mut, w_ref = w_ref, summary_fn = summary_fn)

  # Perform NJ on merged nodes
  Treenj <- nj(dist_h)
  Treenj <- multi2di(Treenj)
  #reenj$edge.length <- rep(1, length(Treenj$edge.length))
  Treenj$tip.label <- dt$group_name

  return(list(tree_backbone = Treenj, backbone_meta = dt))
}


LouvainClusterDynamic <- function(
    X,
    dropout = "-",
    threshold = NULL,           # if NULL → determine dynamically
    resolution = 1,
    dmat = NULL,
    method = c("otsu", "quantile", "dropout_rule")
) {
  method <- match.arg(method)

  # Compute distance matrix if needed
  if (is.null(dmat)) {
    dmat <- DropoutDistMatrix(X, dropout = dropout)
    n <- nrow(X)
  } else {
    n <- nrow(dmat)
  }

  sim <- 1 - dmat
  sim[is.na(sim)] <- 0
  sim_vals <- sim[upper.tri(sim)]

  # ------------------------
  # Dynamic thresholding
  # ------------------------
  if (is.null(threshold)) {

    # Estimate dropout level
    dropout_rate <- mean(X == dropout)
    # Ensures numerical safety
    dropout_rate <- max(min(dropout_rate, 0.99), 0.01)

    if (method == "otsu") {
      # Otsu thresholding
      threshold <- otsu_threshold(sim_vals)

    } else if (method == "quantile") {
      # Higher dropout -> higher quantile cutoff (be stricter)
      # Low dropout -> allow lower similarity edges
      q <- 0.90 + 0.05 * dropout_rate
      q <- min(q, 0.98)
      thr_sim <- quantile(sim_vals, probs = q, na.rm = TRUE)
      threshold <- 1 - thr_sim

    } else if (method == "dropout_rule") {
      # Simple dropout-based threshold rule
      # High dropout (>60%) → threshold ~ 0.35
      # Low dropout (<10%) → threshold ~ 0.15
      threshold <- 0.12 + 0.4 * dropout_rate
    }
  }

  # ------------------------
  # Build similarity graph
  # ------------------------
  edges <- which(sim > (1 - threshold) & upper.tri(sim), arr.ind = TRUE)

  default_membership <- function(n) {
    x <- 1:n
    names(x) <- 1:n
    x
  }

  if (nrow(edges) == 0) {
    return(default_membership(n))
  }

  g <- graph_from_data_frame(
    data.frame(
      from = edges[,1],
      to   = edges[,2],
      weight = sim[edges]
    ),
    directed = FALSE,
    vertices = data.frame(name = 1:n)
  )

  cl <- cluster_louvain(g, weights = E(g)$weight, resolution = resolution)
  return(membership(cl))
}


# Helper: Otsu thresholding for continuous values
otsu_threshold <- function(x, nbins = 200) {
  x <- x[is.finite(x)]
  h <- hist(x, breaks = nbins, plot = FALSE)
  counts <- h$counts
  mids <- h$mids
  total <- sum(counts)
  sumB <- 0
  wB <- 0
  maximum <- 0
  sum1 <- sum(mids * counts)
  for (i in 1:length(counts)) {
    wB <- wB + counts[i]
    if (wB == 0) next
    wF <- total - wB
    if (wF == 0) break
    sumB <- sumB + mids[i] * counts[i]
    mB <- sumB / wB
    mF <- (sum1 - sumB) / wF
    between <- wB * wF * (mB - mF)^2
    if (between > maximum) {
      threshold <- mids[i]
      maximum <- between
    }
  }
  return(1 - threshold)  # convert similarity threshold to distance threshold
}
