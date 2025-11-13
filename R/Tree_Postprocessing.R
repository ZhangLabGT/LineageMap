
merge_2_trees <- function(subtree1,subtree2){
  # Step 1: create a dummy root with two tips
  dummy_root <- read.tree(text="(A:1,B:1);")

  # Step 2: bind subtree1 to tip A and subtree2 to tip B
  temp1 <- bind.tree(dummy_root, subtree1, where=which(dummy_root$tip.label=="A"))
  temp2 <- bind.tree(temp1, subtree2, where=which(temp1$tip.label=="B"))

  # Step 3: remove dummy tip labels (optional)
  temp2$tip.label[temp2$tip.label %in% c("A","B")] <- NULL

  temp2
}

reroot_tree <- function(tree){
  # --- 1️⃣ Get all internal nodes ---
  internal_nodes <- (Ntip(tree) + 1):max(tree$edge)
  half_count <- floor(Ntip(tree)/2)
  # --- 2️⃣ Compute subtree size balance for each internal node ---
  balance_info <- lapply(internal_nodes, function(node) {
    children <- Descendants(tree, node, "children")
    if (length(children) != 2) return(NULL)
    left_tips <- Descendants(tree, children[1], "tips")[[1]]
    right_tips <- Descendants(tree, children[2], "tips")[[1]]
    list(node = node,
         left_tips = left_tips,
         right_tips = right_tips,
         diff = min(abs(length(left_tips) - half_count),abs(length(right_tips) - half_count)))
  })
  balance_info <- Filter(Negate(is.null), balance_info)

  # --- 3️⃣ Pick node with smallest difference (most balanced split) ---
  best_info <- balance_info[[which.min(sapply(balance_info, `[[`, "diff"))]]


  # --- 4️⃣ Choose one side’s tip as outgroup ---
  # (doesn't matter which — it just determines orientation)
  outgroup_tip <-  Descendants(tree, best_info$node, "tips")[[1]]

  # --- 5️⃣ Reroot the tree using that tip ---
  tree_balanced_root <- root(tree, outgroup = outgroup_tip, resolve.root = FALSE)

  # --- 6️⃣ Ladderize for aesthetics ---
  tree_balanced_root$edge.length <- rep(1,length(tree_balanced_root$edge.length))

  root_children <- Descendants(tree_balanced_root, Ntip(tree_balanced_root) + 1, "children")

  if (length(root_children) > 2) {
    # Count tips under each child
    child_sizes <- sapply(root_children, function(ch) length(Descendants(tree_balanced_root, ch, "tips")[[1]]))
    # Pick the two smallest
    merge_children <- root_children[order(child_sizes)][1:2]

    # Extract subtrees
    subtree1 <- extract.clade(tree_balanced_root, merge_children[1])
    subtree2 <- extract.clade(tree_balanced_root, merge_children[2])

    # Merge them into a single subtree
    merged_subtree <- merge_2_trees(subtree1, subtree2)  # binary merge

    # Remove old children and reattach merged subtree
    #remaining_children <- setdiff(root_children, merge_children)
    remaining_subtree <- extract.clade(tree_balanced_root, root_children[order(child_sizes)][3])

    # Rebuild a new tree rooted at a trifurcation (now resolved to binary)
    tree_balanced_root <- merge_2_trees(merged_subtree, remaining_subtree)
  }

  # --- 7️⃣ Ladderize and return ---
  tree_balanced_root$edge.length <- rep(1,length(tree_balanced_root$edge.length))
  tree_balanced_root <- reorder(tree_balanced_root, index.only = FALSE)
  return(tree_balanced_root)
}

make_ultrametric <- function(tree) {
  depths <- node.depth.edgelength(tree)
  max_depth <- max(depths[1:Ntip(tree)])

  for (i in 1:Ntip(tree)) {
    edge_idx <- which(tree$edge[, 2] == i)
    parent <- tree$edge[edge_idx, 1]
    child_depth <- depths[i]
    # Extend this branch so all leaves reach max_depth
    tree$edge.length[edge_idx] <- tree$edge.length[edge_idx] + (max_depth - child_depth)
  }
  tree
}

get_root_to_leaf_paths <- function(tree) {
  # tree$edge is a matrix of parent-child relationships
  edges <- tree$edge

  # Identify the root node
  root <- setdiff(edges[, 1], edges[, 2])
  if (length(root) != 1) {
    stop("The tree must have a single root.")
  }
  root <- root[1]

  # Identify leaf nodes (tips)
  leaves <- setdiff(edges[, 2], edges[, 1])

  # Recursive helper function to find all paths
  find_paths <- function(current_node, current_path) {
    # Find children of the current node
    children <- edges[edges[, 1] == current_node, 2]

    # If no children → it's a leaf, return the completed path
    if (length(children) == 0) {
      return(list(current_path))
    }

    # Otherwise, recursively continue down each child
    paths <- lapply(children, function(child) {
      find_paths(child, c(current_path, child))
    })

    # Flatten nested list
    do.call(c, paths)
  }

  # Start recursion from the root
  paths <- find_paths(root, c(root))

  # Return as list of integer vectors
  return(paths)
}


#' Compute ancestral node coordinates based on child nodes
#'
#' @param tree phylo object
#' @param coords data.frame with columns: id (tip name or number), x, y
#' @param method how to combine child positions ("mean", "median", "weighted")
#' @return data.frame with x, y coordinates for all nodes (tips + internal)
compute_ancestral_coordinates <- function(tree, coords, method = c("mean", "median", "weighted")) {
  method <- match.arg(method)

  Ntip <- length(tree$tip.label)
  Nnode <- tree$Nnode
  total_nodes <- Ntip + Nnode

  # Initialize coordinate matrix
  xy <- matrix(NA, nrow = total_nodes, ncol = 2)
  rownames(xy) <- 1:total_nodes

  # Assign leaf coordinates
  tip_map <- match(coords$id, tree$tip.label)
  xy[tip_map, 1:2] <- as.matrix(coords[, c("x", "y")])

  # Compute subtree sizes for weighting if needed
  subtree_sizes <- rep(1, total_nodes)

  # Process edges in reverse order (postorder traversal)
  postorder <- rev(postorder(tree))

  for (node in postorder) {
    children <- tree$edge[tree$edge[, 1] == node, 2]

    # Skip tips
    if (length(children) == 0) next

    child_coords <- xy[children, , drop = FALSE]
    child_sizes <- subtree_sizes[children]

    # Skip if any child is NA
    if (any(is.na(child_coords))) next

    if (method == "mean") {
      xy[node, ] <- colMeans(child_coords)
    } else if (method == "median") {
      xy[node, ] <- apply(child_coords, 2, median)
    } else if (method == "weighted") {
      xy[node, ] <- colSums(child_coords * child_sizes) / sum(child_sizes)
    }

    # Update subtree size
    subtree_sizes[node] <- sum(child_sizes)
  }

  # Return results
  res <- data.frame(
    id = 1:total_nodes,
    x = xy[, 1],
    y = xy[, 2],
    type = c(rep("leaf", Ntip), rep("internal", Nnode))
  )

  return(res)
}
