# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   https://r-pkgs.org
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

hello <- function() {
  print("Hello, world!")
}



#' LinRace main function(Dif): asymmetric division based Neighbor Joining
#'
#' @param muts lineage barcode matrix
#' @param states cell states of single cells
#' @param counts gene expression data of single cells
#' @param state_lineages the lineages that makes the state network from the root state to leaf states
#' @param max_Iter the maximum iterations for local search
#' @param lambda1 hyperparameter for asymmetric division likelihood
#' @param lambda2 hyperparameter for neighbor distance likelihood
#' @import data.table
#' @import phangorn
#' @export
load_spatedsim <- function(header, ncells, mr, rd, rm, run) {
  content <- c("spatial_meta", "spatial_counts", "spatial_map", "meta", "expr", "cm")
  params <- paste0("_ncells_", ncells, "_mr_", mr, "_rd_", rd, "_rm_", rm, "_run_", run)

  items <- list()

  for (prefix in content) {
    suffix <- ifelse(prefix == "cm", ".txt", ".csv")
    path <- file.path(header, paste0(prefix, params, suffix))

    sep <- ifelse(prefix %in% c("spatial_meta", "meta"), ",",
                  ifelse(prefix == "spatial_map", "\t", " "))

    items[[length(items) + 1]] <- read.csv(path, sep = sep)
  }

  items[[length(items) + 1]] <- params
  names(items) <- content
  return(items)
}


#' Neighbor Joining on unique lineage barcodes
#'
#' @param X lineage barcode matrix
#' @import data.table
#' @import ape
#' @export
DivideMut <- function(X) {
  mut_table_str <- c()
  for (i in 1:nrow(X)) {
    barcode <- X[i,]
    cs <- paste(X[i,], collapse = '|')
    mut_table_str <- c(mut_table_str, cs)
  }
  dt <- as.data.table(mut_table_str)[, list(list(.I)), by = mut_table_str]
  colnames(dt) <- c("barcode", "cellids")

  #obtain tree backbone using nj()
  X_unique <- X[!duplicated(X),]
  l <- dim(X_unique)[1]
  rownames(X_unique) <- 1:l
  sim_tree <- list()

  for (j in 1:l) {
    sim_tree[[j]] <- as.character(X_unique[j,])
  }

  names(sim_tree) <- rownames(X_unique)
  #nmstrings <- c('0','-',c(1:100))
  nmstrings <- unique(unlist(sim_tree))
  sim_data <- phyDat(sim_tree, type = 'USER', levels = nmstrings)

  #dist_wh2 <- WH(sim_data, InfoW, dropout=TRUE)
  dist_h <- dist.hamming(sim_data)
  Treenj <- nj(dist_h)
  tree_backbone <- multi2di(Treenj)
  tree_backbone$edge.length <- rep(1, length(tree_backbone$edge.length))
  #tree_backbone$tip.label <- 1

  return(list(tree_backbone, dt))
}

is_bifurcating <- function(tree) {
  all(table(tree$edge[,1]) <= 2)
}

#' Rebalance a tree
#'
#' @param tree lineage tree
#' @import data.table
#' @import ape
#' @export
reroot_balance <- function(tree) {
  ncells <- length(tree$tip.label)
  node <- tree$edge[1,1]
  while (count_tips(tree,node) > ncells / 2) {
    ncells_bound <- 0
    children_nodes <- tree$edge[tree$edge[,1] == node, 2]
    child_max <- children_nodes[1]
    for (child in children_nodes) {
      #cells_child <- getDescendants(tree,child)
      ncells_child <- count_tips(tree,child)
      if (ncells_child >= ncells_bound) {
        child_max <- child
        ncells_bound <- ncells_child
      }
    }

    if (ncells_bound < ncells/2){
      break
    } else {
      node <- child_max
    }

  }

  tree <- root(tree,node = node)
  root <- length(tree$tip.label)+1
  new_clades <- tree$edge[tree$edge[,1] == root, 2]

  reg <- 0
  biggest_clade_index <- NULL
  biggest_clade_length <- NULL

  for (i in seq_along(new_clades)) {
    num_leaves <- count_tips(tree,i)
    if (num_leaves > reg) {
      reg <- num_leaves
      biggest_clade_index <- i
      biggest_clade_length <- tree$edge.length[i]
    }
  }
  new_clade_children_index <- setdiff(seq_along(new_clades), biggest_clade_index)
  new_clade_children <- new_clades[new_clade_children_index]

  subtree1 <- extract.clade(tree, new_clade_children[1])
  subtree2 <- extract.clade(tree, new_clade_children[2])
  subtree3 <- extract.clade(tree, new_clades[biggest_clade_index])

  # Merge the subtrees into a new single subtree
  new_root_branch <- read.tree(text="(A:1,B:1);")
  new_subtree1 <- bind.tree(new_root_branch, subtree1, where = which(new_root_branch$tip.label=="A"))
  new_subtree <- bind.tree(new_subtree1, subtree2, where = which(new_subtree1$tip.label=="B"))

  # Merge the new subtree with subtree3
  new_root_branch <- read.tree(text="(A:1,B:1):1;")
  new_subtree2 <- bind.tree(new_root_branch, new_subtree, where = which(new_root_branch$tip.label=="A"))
  final_tree <- bind.tree(new_subtree2, subtree3, where = which(new_subtree2$tip.label=="B"))

  #tree$edge <- rbind(tree$edge, new_clade)

  return(final_tree)
}

#' Neighbor Joining tree from lineage barcodes
#'
#' @param X lineage barcode matrix
#' @import data.table
#' @import ape
#' @export
NJ_from_barcode <- function(X) {
  l <- dim(X)[1]
  #rownames(X) <- 1:l
  sim_tree <- list()

  for (j in 1:l) {
    sim_tree[[j]] <- as.character(X[j,])
  }

  names(sim_tree) <- rownames(X)
  #nmstrings <- c('0','-',c(1:100))
  nmstrings <- unique(unlist(sim_tree))
  sim_data <- phyDat(sim_tree, type = 'USER', levels = nmstrings)

  #dist_wh2 <- WH(sim_data, InfoW, dropout=TRUE)
  dist_h <- dist.hamming(sim_data)
  Treenj <- nj(dist_h)
  tree <- multi2di(Treenj)
  tree$edge.length <- rep(1, length(tree$edge.length))
  #tree$tip.label <- 1

  return(tree)
}

#' Getting clone ids from tree
#'
#' @param tree input lineage tree
#' @param num_clones number of desired clones
#' @import ape
#' @export
Define_clones <- function(tree, num_clones = 4) {
  if (!is_bifurcating(tree)) {
    stop("Input tree is not bifurcating.")
  }

  if (length(tree$edge[tree$edge[,1] == tree$edge[1,1], 2]) != 2) {
    cat("Input tree root has more than 2 children. Will run reroot and binarize the tree.\n")
  }
  tree <- reroot_balance(tree)

  if (num_clones %% 2 != 0) {
    stop("Input number of clones should be a power of 2.")
  }

  Num_generation <- log2(num_clones)
  node_list <- list(tree$edge[1,1])

  for (i in seq_len(Num_generation)) {
    old_node_list <- node_list
    node_list <- list()
    for (node in old_node_list) {
      children <- tree$edge[tree$edge[,1] == node, 2]
      node_list <- append(node_list, as.list(children))
    }
  }

  cloneid <- rep(0, length(tree$tip.label))

  for (i in seq_len(num_clones)) {
    sub_tree <- extract.clade(tree, node = node_list[[i]])

    clone_leaf_list <- as.numeric(gsub("[^0-9]", "", sub_tree$tip.label)) - 1
    cloneid[clone_leaf_list + 1] <- i
  }

  #Z_clone <- Mask_clone(cloneid)

  return(cloneid)
}

count_tips <- function(tree, node) {
  if (node >= length(tree$tip.label)){
    subtree <- extract.clade(tree, node)
    return(length(subtree$tip.label))
  } else{
    return(1)
  }
}


NJ_ih <- function(X,locs){
  if (is.matrix(X)) X <- as.dist(X)
  if (anyNA(X))
    stop("missing values are not allowed in the distance matrix\nConsider using njs()")
  if (any(is.infinite(X)))
    stop("infinite values are not allowed in the distance matrix")
  N <- as.integer(attr(X, "Size"))
  if (N < 3) stop("cannot build an NJ tree with less than 3 observations")
  labels <- attr(X, "Labels")
  if (is.null(labels)) labels <- as.character(1:N)

  node_definition <- NULL

  while(N > 3){
    #if (N == 5){browser()}
    r <- Diver_cal(X,N)
    Q <- compute_Q(X,r,N)
    Q_min <- which(Q == min(Q))

    # If multiple candidates have the same Q value
    if (length(Q_min) == 1){
      coord_merge <- rowcol(Q_min[1],N)
    } else{

      coord_merge_candidates <- mapply(rowcol, Q_min, N)
      min_dist <- 99999
      min_index <- 0
      for (i in 1:ncol(coord_merge_candidates)){
        loc_1 <- locs[coord_merge_candidates[1,i],]
        loc_2 <- locs[coord_merge_candidates[2,i],]
        dist_candidate <- sqrt((loc_1$x - loc_2$x)^2 + (loc_1$y - loc_2$y)^2)

        if (dist_candidate < min_dist){
          min_dist <- dist_candidate
          min_index <- i
          coord_merge <- c(coord_merge_candidates[1,i],coord_merge_candidates[2,i])

        }
      }
      Q_min <- Q_min[min_index]
    }

    #browser()
    merge_member_1 = coord_merge[1]
    merge_member_2 = coord_merge[2]

    #compute distances between the merged node and the new node
    dist_merged_1 <- X[Q_min]/2 - (r[merge_member_1] - r[merge_member_2])/(2*(N-2))
    if (dist_merged_1 <= 0){
      dist_merged_1 <- 0
    } else if (dist_merged_1 > X[Q_min]){
      dist_merged_1 <- X[Q_min]
    }
    dist_merged_2 <- X[Q_min] - dist_merged_1

    node_definition <- sprintf('(%s:%f, %s:%f)',
                               labels[merge_member_1], dist_merged_1,labels[merge_member_2], dist_merged_2)
    #print(node_definition)
    X <- update_dist(X,merge_member_1,merge_member_2,labels,node_definition)
    locs <- update_locs(merge_member_1,merge_member_2,locs)

    N <- as.integer(attr(X, "Size"))
    labels <- attr(X, "Labels")
  }
  merge_member_1 = labels[1]
  merge_member_2 = labels[2]
  r <- Diver_cal(X,N)

  dist_merged_1 <- X[1]/2 - (r[1] - r[2])/(2*(N-2))
  if (dist_merged_1 <= 0){
    dist_merged_1 <- 0
  } else if (dist_merged_1 > X[1]){
    dist_merged_1 <- X[1]
  }
  dist_merged_2 <- X[1] - dist_merged_1

  # ...then determine their distance to the other remaining node
  node_definition = labels[3]
  internal_len = max(0.5*(X[2]+X[3]-X[1]),0)
  # ...and finally create the newick string describing the whole tree.
  newick = sprintf("(%s:%f, %s:%f, %s:%f);",merge_member_1, dist_merged_1,
                   node_definition, internal_len,
                   merge_member_2, dist_merged_2)

  # return the phylo object transformed from newick.
  tree_nj <- read.tree(text = newick)
  return(tree_nj)
}

update_locs <- function(merge_member_1,merge_member_2,locs){
  loc_1 <- locs[merge_member_1,]
  loc_2 <- locs[merge_member_2,]
  new_loc <- data.frame(x = (loc_1$x + loc_2$x)/2,y = (loc_1$y + loc_2$y)/2)
  results <- rbind(new_loc,locs)

  results <- results[-c(merge_member_1+1,merge_member_2+1),]

  return(results)
}

distdex<-function(i,j,n){
  n*(i-1) - i*(i-1)/2 + j-i
} #given row, column, and n, return index

# rowcol <- function(index, N) {
#   i <- ceiling((2*N - 1 - sqrt((2*N - 1)^2 - 8*index))/2)
#   j <- index - (i - 1) * (i - 2) / 2
#   print(paste("rowcol index:", index, "returns:", i, j))
#   return(c(i, j))
# }

rowcol<-function(ix,n) { #given index, return row and column
  nr=ceiling(n-(1+sqrt(1+4*(n^2-n-2*ix)))/2)
  nc=n-(2*n-nr+1)*nr/2+ix+nr
  return(c(nr,nc))
}

Diver_cal <- function(X,N){
  r <- rep(0,N)
  for (i in 1:length(X)){
    rc <- rowcol(i,N)
    r[rc] <- r[rc] + X[i]
  }
  return(r)
}

compute_Q <- function(X,r,N){
  Q <- X
  for (i in 1:length(X)){
    coord <- rowcol(i,N)
    Q[i] <- X[i] - sum(r[coord])/(N-2)
  }
  return(Q)
}

update_dist <- function(X,index1,index2,labels,new_node_id){
  X <- as.matrix(X)
  N <- length(labels)

  labels_new = c(new_node_id,labels)
  results <- matrix(0, N+1, N+1)
  results[2:(N+1),2:(N+1)] <- X

  for (i in 2:(N+1)){
    dist_temp <- 0.5*(X[i-1,index1]+X[i-1,index2]-X[index1,index2])
    results[1,i] <- dist_temp
    results[i,1] <- dist_temp
  }
  rownames(results) <- labels_new
  colnames(results) <- labels_new

  results <- tryCatch({
    results <- results[-c(index1+1,index2+1),]
    results <- results[,-c(index1+1,index2+1)]
    return(as.dist(results))
  }, error = function(e) {
    print(paste("Only one node left."))
    newick <- labels_new[1]
    return(newick)
  })
}
