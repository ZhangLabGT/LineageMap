#' LineageMap main function: build tree with barcode, spatial coordinates and gene expression states
#'
#' @param muts lineage barcode matrix
#' @param meta spatial and state information of cells
#' @param state_lineages the lineages that makes the state network from the root state to leaf states
#' @param max_Iter the maximum iterations for local search
#' @param threshold the threshold for building the similarity map
#' @param backbone_type the method used to reconstruct backbone
#' @param lambda1 weight for gene expression state likelihood
#' @param lambda2 weight for spatial location likelihood
#' @param alpha hyperparameter for spatial likelihood
#' @import data.table
#' @import phangorn
#' @export
Build_LineageMap <- function(muts,meta, state_lineages, max_Iter = 200,threshold = 0.2,backbone_type = c("majority","NJ_median","NJ_mean"),lambda1 = 0.05,lambda2 = 0.1,alpha = 1) {
  labels <- rownames(muts)
  groups <- LouvainCluster(muts,threshold = threshold)
  if (length(unique(groups)) <= 2){
    groups <- LouvainCluster(muts,threshold = 1/(ncol(muts_leaves)))
  }

  if (backbone_type == "majority"){
    returnList <- SubgroupTree(muts,groups)
  }
  else if (backbone_type == "NJ_median"){
    returnList <- SubgroupTree_NJ(muts,groups,summary_fn = median)
  }

  tree_backbone <- returnList[[1]]
  dt <- returnList[[2]]

  subtree_list <- list()
  for (i in 1:nrow(dt)) {
    cellids <- unlist(dt$cells[i])
    #print(labels[cellids])
    muts_sub <- muts[cellids,]
    #rownames(muts_sub) <- cellids
    meta_sub <- meta[cellids,]
    labels_sub <- labels[cellids]
    if (length(cellids) > 1) {
      res <- FindBestTree(muts_sub,meta_sub, labels = labels_sub, state_lineages, maxIter = max_Iter,lambda1 = lambda1,lambda2 = lambda2,alpha = alpha)
      #subtree_opt <- res$best_tree
      subtree_opt <- res
      subtree_opt$name <- dt$group_name[i]
      subtree_list[[length(subtree_list) + 1]] <- subtree_opt
    } else {
      tree_backbone$tip.label[tree_backbone$tip.label == dt$group_name[i]] <- substr(labels_sub, 6, nchar(labels_sub))
    }
  }
  #browser()
  tree_final <- ConstructTree(tree_backbone, subtree_list)
}


#' Finding the best tree structure based on local search
#'
#' @param muts the lineage barcodes of cells
#' @param meta spatial coordinates and state information of cells
#' @param labels the labels of cells
#' @param state_lineages the lineages that makes the state network from the root state to leaf states
#' @param newick_lookup lookup table for internal nodes
#' @param maxIter maximum number of iterations of local search
#' @param n_chains the number of chains to run
#' @param lambda1 weight for gene expression state likelihood
#' @param lambda2 weight for spatial location likelihood
#' @param alpha hyperparameter for spatial likelihood
#' @import castor
FindBestTree <- function(muts,meta, labels, state_lineages, newick_lookup = NULL, maxIter = 200,n_chains = 5,lambda1 = 0.05,lambda2 = 0.1,alpha = 1) {

  #browser()
  seeds <- runif(10000, 1, 99999)
  leafs <- c()
  for (i in 1:length(labels)) {
    label <- labels[i]
    n <- nchar(label)
    if (substr(label, n, n) == ';') {
      labels[i] <- substr(label, 1, n - 1)
    }else {
      if (startsWith(label,"cell_")){
        leafs <- c(leafs, substr(label, 6, n))
      } else{
        leafs <- c(leafs,label)
      }
    }
  }

  if (length(labels) == 2){
    best_tree <- rtree(2)
    best_tree$tip.label <- rownames(muts)
    best_tree$edge.length <- rep(1, length(best_tree$edge.length))

    return(best_tree)
  }



  tree_init <- NJ_from_barcode(muts)


  lambda_values <- seq(0.5, 1, length.out = n_chains)
  results <- mclapply(
    1:n_chains,
    function(chain_id) run_chain(chain_id, seeds, muts, meta, state_lineages, maxIter = maxIter,lambda_restart = lambda_values[chain_id],lambda1 = lambda1,lambda2 = lambda2,alpha = alpha),
    mc.cores = n_chains
  )
  #browser()
  # Find best chain
  best_idx <- which.max(sapply(results, function(x) x$max_likelihood))
  best_result <- results[[best_idx]]
  max_likelihood <- best_result
  best_tree <- best_result$best_tree
  loc_infer <- best_result$best_loc

  return(best_tree)
}


#' LineageMap main function: build tree with barcode, spatial coordinates and gene expression states, using parallelized subtree searches
#'
#' @param muts lineage barcode matrix
#' @param meta spatial and state information of cells
#' @param state_lineages the lineages that makes the state network from the root state to leaf states
#' @param max_Iter the maximum iterations for local search
#' @param threshold the threshold for building the similarity map
#' @param backbone_type the method used to reconstruct backbone
#' @param outer_cores the number of cores per outer loop
#' @param inner_cores the number of cores per inner loop (within FindBestTree)
#' @param lambda1 weight for gene expression state likelihood
#' @param lambda2 weight for spatial location likelihood
#' @param alpha hyperparameter for spatial likelihood
#' @import data.table
#' @import phangorn
#' @export
Build_LineageMap_parallel <- function(muts, meta, state_lineages,
                             max_Iter = 200, threshold = 0.2,
                             backbone_type = c("majority", "NJ_median", "NJ_mean"),
                             outer_cores = 12, inner_cores = 5,lambda1 = 0.05,lambda2 = 0.1,alpha = 1) {
  labels <- rownames(muts)
  groups <- LouvainCluster(muts, threshold = threshold)
  if (length(unique(groups)) <= 2){
    groups <- LouvainCluster(muts, threshold = 1 / ncol(muts))
  }

  if (backbone_type == "majority") {
    returnList <- SubgroupTree(muts, groups)
  } else if (backbone_type == "NJ_median") {
    returnList <- SubgroupTree_NJ(muts, groups, summary_fn = median)
  }

  tree_backbone <- returnList[[1]]
  dt <- returnList[[2]]

  # Parallelise across subgroups
  subtree_list <- mclapply(1:nrow(dt), function(i) {
    cellids <- unlist(dt$cells[i])
    muts_sub <- muts[cellids, ]
    meta_sub <- meta[cellids, ]
    labels_sub <- labels[cellids]

    if (length(cellids) > 1) {
      res <- FindBestTree(muts_sub, meta_sub, labels = labels_sub,
                          state_lineages, maxIter = max_Iter,
                          n_chains = inner_cores,lambda1 = lambda1,lambda2 = lambda2,alpha = alpha)  # pass fewer chains
      res$tip.label <- sub("cell_","",labels[as.numeric(res$tip.label)])
      res$name <- dt$group_name[i]
      return(res)
    } else {
      return(list(single = TRUE, name = dt$group_name[i], label = substr(labels_sub, 6, nchar(labels_sub))))
    }
  }, mc.cores = outer_cores)

  # Update singletons
  for (sub in subtree_list) {
    if (!is.null(sub$single) && sub$single) {
      tree_backbone$tip.label[tree_backbone$tip.label == sub$name] <- sub$label
    }
  }

  # Merge all subtrees
  subtree_objs <- Filter(function(x) is.null(x$single), subtree_list)
  tree_final <- ConstructTree(tree_backbone, subtree_objs)
  return(tree_final)
}

#' LineageMap main function: build tree with pairwise distances, spatial coordinates and gene expression states, using parallelized subtree searches
#'
#' @param muts lineage barcode matrix
#' @param meta spatial and state information of cells
#' @param state_lineages the lineages that makes the state network from the root state to leaf states
#' @param max_Iter the maximum iterations for local search
#' @param threshold the threshold for building the similarity map
#' @param backbone_type the method used to reconstruct backbone
#' @param outer_cores the number of cores per outer loop
#' @param inner_cores the number of cores per inner loop (within FindBestTree)
#' @param lambda1 weight for gene expression state likelihood
#' @param lambda2 weight for spatial location likelihood
#' @param alpha hyperparameter for spatial likelihood
#' @import data.table
#' @import phangorn
#' @export
Build_LineageMap_pd <- function(muts, meta, state_lineages,
                                      max_Iter = 200, threshold = 0.2,
                                      backbone_type = c("majority", "NJ_median", "NJ_mean"),
                                      outer_cores = 12, inner_cores = 5,lambda1 = 0.05,lambda2 = 0.1,alpha = 1) {
  labels <- rownames(muts)
  rownames(meta) <- labels
  groups <- LouvainCluster(muts, threshold = threshold)
  #names(groups) <- sub("cell_", "",labels)
  if (length(unique(groups)) <= 2){
    groups <- LouvainCluster(muts, threshold = 1 / ncol(muts))
  }

  if (backbone_type == "majority") {
    returnList <- SubgroupTree(muts, groups)
  } else if (backbone_type == "NJ_median") {
    returnList <- SubgroupTree_NJ(muts, groups, summary_fn = median)
  }

  tree_backbone <- returnList[[1]]
  dt <- returnList[[2]]
  #browser()
  # Parallelise across subgroups
  subtree_list <- mclapply(1:nrow(dt), function(i) {
    cellids <- unlist(dt$cells[i])
    muts_sub <- muts[cellids, ]
    meta_sub <- meta[cellids, ]
    labels_sub <- labels[cellids]
    #browser()
    if (length(cellids) > 1) {
      res <- FindBestTree(muts_sub, meta_sub, labels = labels_sub,
                          state_lineages, maxIter = max_Iter,
                          n_chains = inner_cores,lambda1 = lambda1,lambda2 = lambda2,alpha = alpha)  # pass fewer chains
      #res$tip.label <- sub("cell_","",labels[as.numeric(res$tip.label)])
      res$name <- dt$group_name[i]
      return(res)
    } else {
      return(list(single = TRUE, name = dt$group_name[i], label = labels_sub))
    }
  }, mc.cores = outer_cores)

  # Update singletons
  for (sub in subtree_list) {
    if (!is.null(sub$single) && sub$single) {
      tree_backbone$tip.label[tree_backbone$tip.label == sub$name] <- sub$label
    }
  }

  # Merge all subtrees
  subtree_objs <- Filter(function(x) is.null(x$single), subtree_list)
  tree_final <- ConstructTree(tree_backbone, subtree_objs)
  return(tree_final)
}


run_chain <- function(chain_id, seeds, muts_leaves, df_spatial_continuous, phyla, maxIter = 1000,lambda_restart = 0.1,lambda1 = 0.05,lambda2 = 0.1,alpha = 1) {
  # --- Initialization ---
  set.seed(seeds[chain_id])

  muts <- as.data.frame(muts_leaves)
  loc_data <- df_spatial_continuous%>%dplyr::select(x,y)

  # Starting tree
  tree <- starting_tree_from_data(muts, loc_data, lambda = lambda_restart)
  #tree$edge.length <- rep(1, length(tree$edge.length))
  tree_init <- tree

  # Initial likelihood evaluation
  returnlist <- LikelihoodCal_ST(
    tree, muts, df_spatial_continuous$state, phyla, loc_data,
    lambda1 = 0, lambda2 = 0.1, alpha = 1
  )

  max_likelihood <- returnlist$likelihood
  current_likelihood <- max_likelihood
  loc_infer <- returnlist$loc_data
  best_loc <- loc_infer
  best_tree <- tree
  muts_internal <- returnlist$muts_internal
  #lambda_restart <- 0
  restart_count <- 0

  # Tracking containers
  likelihood_curve <- numeric()
  best_tree <- tree
  best_tree_list <- list(tree)
  best_loc_list <- list(loc_infer)
  maxl_list <- c(max_likelihood)

  LikelihoodCal_time <- 0
  TreeLC_time <- 0

  # --- Local Search Loop ---
  for (i in 1:maxIter) {
    set.seed(seeds[chain_id] + i)
    # Tree proposal
    c_time <- proc.time()
    tree_new <- TreeLC2(tree, muts_internal, loc_infer)
    TreeLC_time <- TreeLC_time + proc.time() - c_time

    # Likelihood
    c_time <- proc.time()
    returnlist <- LikelihoodCal_ST(
      tree_new, muts, df_spatial_continuous$state, phyla, loc_data,
      lambda1 = lambda1, lambda2 = lambda2, alpha = alpha
    )
    likelihood_new <- returnlist$likelihood
    loc_infer <- returnlist$loc_data
    muts_internal <- returnlist$muts_internal
    LikelihoodCal_time <- LikelihoodCal_time + proc.time() - c_time

    # Record likelihoods
    current_likelihood <- likelihood_new
    likelihood_curve <- c(likelihood_curve, max_likelihood)

    # Update best
    if (likelihood_new > max_likelihood) {
      max_likelihood <- likelihood_new
      tree <- tree_new
      best_tree <- tree
      best_loc <- loc_infer
    }

    # Progress
    if (i %% 100 == 0) {
      msg <- sprintf(
        "Chain %d | Iter %d | Current likelihood = %.4f | Best likelihood = %.4f",
        chain_id, i, current_likelihood, max_likelihood
      )
      print(msg)
    }

    # Local optima detection (same as your original)
    if (i > 100) {
      recent_vals <- likelihood_curve[(i - 20):i]
      if (length(unique(recent_vals)) == 1) {
        best_tree_list[[length(best_tree_list) + 1]] <- best_tree
        best_loc_list[[length(best_loc_list) + 1]] <- best_loc
        maxl_list <- c(maxl_list, max_likelihood)

        set.seed(seeds[chain_id] + i)
        #lambda_restart <- 1 - i / (2 * maxIter)
        tree <- tree_init

        returnlist <- LikelihoodCal_ST(
          tree, muts, df_spatial_continuous$state, phyla, loc_data,
          lambda1 = lambda1, lambda2 = lambda2, alpha = alpha
        )
        max_likelihood <- returnlist$likelihood
        loc_infer <- returnlist$loc_data
        muts_internal <- returnlist$muts_internal
      }
    }
  }

  return(list(
    best_tree = best_tree,
    best_loc = best_loc,
    likelihood_curve = likelihood_curve,
    max_likelihood = max_likelihood,
    all_trees = best_tree_list,
    all_locs = best_loc_list,
    all_maxl = maxl_list
  ))
}


ConstructTree <- function(tree_backbone,subtrees){
  tree <- tree_backbone
  #print(tree_backbone$tip.label)
  for (subtree in subtrees){
    bind_tip <- subtree$name
    #print(bind_tip)
    tree <- bind.tree(tree, subtree, where = which(tree$tip.label==bind_tip))
  }
  return(tree)
}


starting_tree_from_data <- function(X,loc,lambda = 0.9) {
  l <- dim(X)[1]
  #rownames(X) <- 1:l

  # ---- Barcode distances ----
  dist_barcode <-DropoutDistMatrix(X)

  # ---- Spatial distances ----
  dist_spatial <- as.matrix(dist(loc, method = "euclidean"))

  # ---- Normalize both to [0,1] ----
  if (sum(dist_barcode)>0){
    dist_barcode <- dist_barcode / max(dist_barcode)
  }
  if (sum(dist_spatial)>0){
    dist_spatial <- dist_spatial / max(dist_spatial)
  }

  # ---- Weighted combination ----
  dist_combined <- lambda * dist_barcode + (1 - lambda) * dist_spatial

  # ---- NJ tree ----
  Treenj <- nj(as.dist(dist_combined))
  tree <- multi2di(Treenj)
  #tree$edge.length <- rep(1, length(tree$edge.length))

  return(tree)
}
