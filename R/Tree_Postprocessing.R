
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

visualize_tree_3d <- function(tree_ultra,loc_lm_nj,cs_lm_nj){
  root_node <- Ntip(tree_ultra) + 1
  depths <- node.depth.edgelength(tree_ultra)
  depths <- depths / max(depths) * 10  # rescale to 0–10 range for visualization

  # --- 2️⃣ Combine location + z + inferred cell state ---
  df_nodes <- loc_lm_nj %>%
    mutate(
      node = as.numeric(loc_lm_nj$id),
      z = depths[node],
      state = as.factor(cs_lm_nj)   # directly assign inferred states
    )

  # --- 3️⃣ Build edge (parent-child) coordinate table ---
  edges <- as.data.frame(tree_ultra$edge)
  colnames(edges) <- c("parent", "child")

  edges_coords <- edges %>%
    left_join(df_nodes %>% rename(x_parent = x, y_parent = y, z_parent = z), by = c("parent" = "node")) %>%
    left_join(df_nodes %>% rename(x_child = x, y_child = y, z_child = z), by = c("child" = "node"))

  # --- 4️⃣ Define state colors ---
  num_states <- length(unique(df_nodes$state))
  palette_colors <- brewer.pal(max(3, min(num_states, 12)), "Set1")
  state_colors <- setNames(palette_colors[1:num_states], sort(unique(df_nodes$state)))

  # --- 5️⃣ Create interactive 3D plot ---
  fig <- plot_ly()

  # Add node markers
  fig <- fig %>%
    add_markers(
      data = df_nodes,
      x = ~x, y = ~y, z = ~z,
      color = ~state,
      colors = state_colors,
      marker = list(size = 6, opacity = 0.9, line = list(width = 1.5, color = "black")),
      text = ~paste("Node:", node, "<br>State:", state),
      hoverinfo = "text"
    )

  # Add connecting edges (ancestor → descendant)
  for (i in 1:nrow(edges_coords)) {
    fig <- fig %>% add_trace(
      type = "scatter3d",
      mode = "lines",
      x = c(edges_coords$x_parent[i], edges_coords$x_child[i]),
      y = c(edges_coords$y_parent[i], edges_coords$y_child[i]),
      z = c(edges_coords$z_parent[i], edges_coords$z_child[i]),
      line = list(color = 'black', width = 2.5),
      showlegend = FALSE
    )
  }

  # --- 6️⃣ Optionally: Add semi-transparent planes for generations ---
  for (z_level in unique(df_nodes$z)) {
    fig <- fig %>% add_trace(
      type = "mesh3d",
      x = c(min(df_nodes$x), max(df_nodes$x), max(df_nodes$x), min(df_nodes$x)),
      y = c(min(df_nodes$y), min(df_nodes$y), max(df_nodes$y), max(df_nodes$y)),
      z = rep(z_level, 4),
      color = I("gray"),
      opacity = 0.15,
      showscale = FALSE
    )
  }

  # --- 7️⃣ Final layout ---
  fig <- fig %>%
    layout(scene = list(
      xaxis = list(title = "Spatial X"),
      yaxis = list(title = "Spatial Y"),
      zaxis = list(title = "Lineage Depth"),
      camera = list(eye = list(x = 1.3, y = 1.3, z = 1))
    ))

  fig <- fig %>%
    layout(
      scene = list(
        xaxis = list(
          title = list(text = "Spatial X", font = list(size = 22)),
          tickfont = list(size = 18)
        ),
        yaxis = list(
          title = list(text = "Spatial Y", font = list(size = 22)),
          tickfont = list(size = 18)
        ),
        zaxis = list(
          title = list(text = "Lineage Depth", font = list(size = 22)),
          tickfont = list(size = 18)
        ),
        camera = list(eye = list(x = 1.3, y = 1.3, z = 1))
      ),
      margin = list(l = 0, r = 0, b = 0, t = 0),
      font = list(family = "Arial", size = 20, color = "black")
    )
  fig <- fig %>%
    layout(scene = list(
      xaxis = list(backgroundcolor = "rgb(245,245,245)"),
      yaxis = list(backgroundcolor = "rgb(245,245,245)"),
      zaxis = list(backgroundcolor = "rgb(245,245,245)")
    ))

  fig <- fig %>%
    layout(
      scene = list(
        camera = list(eye = list(x = 1.6, y = 1.6, z = 0.1))
      )
    )
  fig
}


visualize_tree_2d <- function(loc_data,cell_types,internal = FALSE){
  if (internal){
    loc_df <- loc_data %>%
      mutate(
        node = as.numeric(loc_data$id),
        state = state)

    # Extract edges as data frame (parent-child pairs)
    edges_df <- as.data.frame(tree_ultra$edge)
    colnames(edges_df) <- c("parent", "child")

    # Merge coordinates for edges
    edge_coords <- edges_df %>%
      left_join(loc_df, by = c("parent" = "node")) %>%
      rename(x_parent = x, y_parent = y, state_parent = state) %>%
      left_join(loc_df, by = c("child" = "node")) %>%
      rename(x_child = x, y_child = y, state_child = state)

    leaf_nodes <- setdiff(edge_coords$child, edge_coords$parent)
    loc_df <- loc_df %>%
      mutate(node_type = ifelse(node %in% leaf_nodes, "Leaf", "Internal"))

    # Plot
    p <- ggplot() +
      # Edges (light gray dashed)
      #geom_segment(
      #  data = edge_coords,
      #  aes(x = x_parent, y = y_parent, xend = x_child, yend = y_child),
      #  color = "gray70", linetype = "dashed",
      #  linewidth = 0.4, alpha = 0.8
      #) +

      # Nodes (leaf vs internal with different shapes)
      geom_point(
        data = loc_df,
        aes(
          x = x, y = y,
          color = factor(state),
          shape = node_type
        ),
        size = 4, stroke = 0.4, alpha = 0.9
      ) +

      # Manual color & shape scales
      scale_color_manual(values = colors, name = "Cell state") +
      scale_shape_manual(
        values = c("Internal" = 17, "Leaf" = 16),
        name = "Node type"
      ) +

      # Coordinate scaling
      coord_fixed() +

      # Labels
      labs(
        x = "x coordinate",
        y = "y coordinate",
        title = "2D Spatial Visualization of Lineage Tree"
      ) +

      # Publication-ready theme
      theme_minimal(base_family = "Helvetica", base_size = 16) +
      theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14, color = "black"),
        axis.ticks = element_line(size = 0.6, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 13),
        legend.key.size = unit(0.8, "cm"),
        panel.grid = element_blank(),
        plot.margin = margin(10, 10, 10, 10),
        legend.position = "right"
      )
    return(p)
  } else{
    loc_df <- loc_data %>%
      mutate(
        node = as.numeric(rownames(loc_data)),
        state = loc_data$state)


    # Define a color palette for states
    colors <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3",
                "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999")
    state_levels <- sort(unique(cs_lm_nj))
    colors <- setNames(palette_colors[1:num_states], sort(unique(df_nodes$state)))

    # Plot
    p <- ggplot() +
      # Edges (light gray dashed)
      # Nodes (colored by cell state)
      geom_point(data = loc_df,
                 aes(x = x, y = y, color = factor(state)),
                 size = 4, stroke = 0.4, alpha = 0.9) +
      scale_color_manual(values = colors, name = "Cell state") +
      theme_minimal(base_size = 14) +
      coord_fixed() +
      labs(x = "x coordinate", y = "y coordinate",
           title = "2D Spatial Visualization of Lineage Tree") +
      theme(
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "right"
      )
    p+

      # Publication-ready theme
      theme_minimal(base_family = "Helvetica", base_size = 16) +
      theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14, color = "black"),
        axis.ticks = element_line(size = 0.6, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 13),
        legend.key.size = unit(0.8, "cm"),
        panel.grid = element_blank(),
        plot.margin = margin(10, 10, 10, 10),
        legend.position = "right"
      )
    return(p)
  }

}
