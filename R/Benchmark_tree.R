get_descendant_tips <- function(tree, node) {
  Ntip <- Ntip(tree)
  # BFS stack
  stack <- node
  tips <- integer(0)
  while (length(stack) > 0) {
    cur <- stack[1]
    stack <- stack[-1]
    children <- tree$edge[tree$edge[,1] == cur, 2]
    if (length(children) == 0) next
    for (ch in children) {
      if (ch <= Ntip) {
        tips <- c(tips, ch)
      } else {
        stack <- c(stack, ch)
      }
    }
  }
  unique(tree$tip.label[sort(tips)])
}

# Greedy one-to-one matching by maximum overlap
match_internal_nodes_by_overlap <- function(tree_true, tree_infer) {
  Ntip_true <- Ntip(tree_true)
  Ntip_inf  <- Ntip(tree_infer)
  if (!all(tree_true$tip.label %in% tree_infer$tip.label)) {
    stop("Tip labels must match between true and inferred trees to compute ancestral MSE.")
  }
  # internal node ids
  int_true <- (Ntip_true + 1):(Ntip_true + tree_true$Nnode)
  int_inf  <- (Ntip_inf + 1):(Ntip_inf + tree_infer$Nnode)
  # clade lists (tip label vectors)
  clades_true <- lapply(int_true, function(n) get_descendant_tips(tree_true, n))
  clades_inf  <- lapply(int_inf,  function(n) get_descendant_tips(tree_infer, n))
  names(clades_true) <- as.character(int_true)
  names(clades_inf)  <- as.character(int_inf)
  # overlap matrix: rows true, cols infer
  M <- matrix(0, nrow = length(clades_true), ncol = length(clades_inf),
              dimnames = list(names(clades_true), names(clades_inf)))
  for (i in seq_along(clades_true)) {
    for (j in seq_along(clades_inf)) {
      M[i, j] <- length(intersect(clades_true[[i]], clades_inf[[j]]))
    }
  }
  # greedy one-to-one matching: pick largest overlap, then remove row & col, repeat
  matched <- list()
  Mcopy <- M
  while (max(Mcopy) > 0) {
    idx <- which(Mcopy == max(Mcopy), arr.ind = TRUE)[1, ]
    r <- rownames(Mcopy)[idx[1]]
    c <- colnames(Mcopy)[idx[2]]
    matched[[r]] <- c
    # zero out row and column
    Mcopy[idx[1], ] <- 0
    Mcopy[, idx[2]] <- 0
  }
  # return mapping as named character vector: names = true node ids, values = infer node ids
  unlist(matched)
}

# Main function: compute MSE of ancestral coordinates
ancestral_MSE <- function(tree_true, tree_infer, loc_true, loc_infer,
                          include_tips = FALSE, match_method = "overlap_greedy",
                          verbose = FALSE) {
  # tree_true, tree_infer: ape phylo objects with identical tip labels
  # loc_true, loc_infer: data.frame or matrix with rows for nodes (tips and internals) or with rownames
  #    Preferred: rownames = c(tree$tip.label, internal_node_names) OR you can pass rows in ape order:
  #    rows 1:Ntip = tips (in order of tree$tip.label), rows (Ntip+1):Nnode = internals in ape numeric order.
  # include_tips: if TRUE, also include tip coordinates in the MSE computation (tips usually already match)
  # Returns list: mse (numeric), per_node (data.frame), mapping (data.frame)

  # Ensure tip labels match
  if (!all(sort(tree_true$tip.label) == sort(tree_infer$tip.label))) {
    stop("Tip labels in true and inferred trees must match (same set).")
  }

  Ntip <- Ntip(tree_true)

  # prepare loc lookup: try to use rownames if present, else assume ape order (tips then internals)
  prepare_loc_lookup <- function(tree, loc) {
    if (is.null(rownames(loc))) {
      # expect rows in ape order: tips (1:Ntip) then internals (Ntip+1 : Ntip+Nnode)
      if (nrow(loc) != (Ntip + tree$Nnode)) {
        stop("loc has no rownames and number of rows doesn't match tree nodes.")
      }
      rownames(loc) <- c(tree$tip.label, as.character((Ntip + 1):(Ntip + tree$Nnode)))
    }
    return(as.matrix(loc))
  }

  loc_true <- prepare_loc_lookup(tree_true, loc_true)
  loc_infer <- prepare_loc_lookup(tree_infer, loc_infer)

  # Match internal nodes
  if (match_method == "overlap_greedy") {
    mapping <- match_internal_nodes_by_overlap(tree_true, tree_infer)
  } else {
    stop("Only 'overlap_greedy' match_method implemented.")
  }

  if (length(mapping) == 0) {
    warning("No internal node matches found (no overlaps). Returning NA.")
    return(list(mse = NA_real_, per_node = NULL, mapping = NULL))
  }

  # Build per-node errors
  rows <- list()
  for (true_node in names(mapping)) {
    inf_node <- mapping[[true_node]]
    # true node rowname in loc_true: if loc prepared used integer name for internals
    true_rowname <- true_node
    infer_rowname <- inf_node
    if (! (true_rowname %in% rownames(loc_true))) {
      # try with "nodeX" label if present in loc_true
      if (paste0("node", as.integer(true_node) - Ntip) %in% rownames(loc_true)) {
        true_rowname <- paste0("node", as.integer(true_node) - Ntip)
      }
    }
    if (! (infer_rowname %in% rownames(loc_infer))) {
      if (paste0("node", as.integer(inf_node) - Ntip) %in% rownames(loc_infer)) {
        infer_rowname <- paste0("node", as.integer(inf_node) - Ntip)
      }
    }
    if (!(true_rowname %in% rownames(loc_true)) || !(infer_rowname %in% rownames(loc_infer))) {
      # skip if coordinates not available
      next
    }
    z_true <- as.numeric(loc_true[true_rowname, ])
    z_inf  <- as.numeric(loc_infer[infer_rowname, ])
    se <- sum((z_true - z_inf)^2)
    rows[[length(rows) + 1]] <- data.frame(
      true_node = as.integer(true_node),
      infer_node = as.integer(infer_rowname),
      z_true_x = z_true[1], z_true_y = z_true[2],
      z_inf_x = z_inf[1], z_inf_y = z_inf[2],
      sq_error = se, stringsAsFactors = FALSE
    )
  }

  per_node_df <- do.call(rbind, rows)
  if (is.null(per_node_df) || nrow(per_node_df) == 0) {
    warning("No matched coordinates found for internals.")
    return(list(mse = NA_real_, per_node = NULL, mapping = mapping))
  }

  mse <- mean(per_node_df$sq_error)

  if (verbose) {
    message(sprintf("Matched %d internal nodes; MSE = %g", nrow(per_node_df), mse))
  }

  return(list(mse = mse, per_node = per_node_df, mapping = mapping))
}


define_clones <- function(tree, k) {
  n_tips <- Ntip(tree)

  # Simple approach: use hierarchical clustering on cophenetic distance
  coph_dist <- cophenetic(tree)
  hc <- hclust(as.dist(coph_dist))
  clone_assignments <- cutree(hc, k = min(k, n_tips))
  names(clone_assignments) <- tree$tip.label

  return(clone_assignments)
}

information_weighted_rf <- function(tree_true, tree_infer) {
  n_tips <- Ntip(tree_true)
  bp_true <- prop.part(tree_true)
  bp_infer <- prop.part(tree_infer)

  total_weighted_penalty <- 0
  max_possible_penalty <- 0

  for (split in bp_true) {
    weight <- min(length(split), n_tips - length(split))
    max_possible_penalty <- max_possible_penalty + weight
    if (!any(sapply(bp_infer, function(s) setequal(split, s)))) {
      total_weighted_penalty <- total_weighted_penalty + weight
    }
  }

  return(if (max_possible_penalty == 0) 0 else total_weighted_penalty / max_possible_penalty)
}

# -Spatio-Topological Consistency
spatio_topological_consistency <- function(tree_true, tree_infer, leaf_coords_true, leaf_coords_infer) {
  dist_pat_true <- cophenetic(tree_true)
  dist_euc_true <- as.matrix(dist(leaf_coords_true))
  dist_pat_infer <- cophenetic(tree_infer)
  dist_euc_infer <- as.matrix(dist(leaf_coords_infer))

  common_tips <- intersect(tree_true$tip.label, tree_infer$tip.label)
  common_tips <- as.numeric(common_tips)
  corr_true <- cor(as.vector(dist_pat_true[common_tips, common_tips]),
                   as.vector(dist_euc_true[common_tips, common_tips]),
                   method = "spearman")
  corr_inferred <- cor(as.vector(dist_pat_infer[common_tips, common_tips]),
                       as.vector(dist_euc_infer[common_tips, common_tips]),
                       method = "spearman")

  if (is.na(corr_true) || is.na(corr_inferred) || abs(corr_true) < 0.01) return(0)

  consistency <- abs(corr_inferred) / abs(corr_true)
  return(min(1.0, if (sign(corr_inferred) != sign(corr_true)) consistency * 0.5 else consistency))
}

#Chronological Divergence Score
chronological_divergence_score <- function(tree_true, tree_infer) {
  calculate_node_depths <- function(tree) {
    n_tips <- Ntip(tree)
    depths <- numeric(n_tips + tree$Nnode)
    tree_ordered <- reorder(tree, "pruningwise")
    for (i in seq_len(nrow(tree_ordered$edge))) {
      parent <- tree_ordered$edge[i, 1]
      child <- tree_ordered$edge[i, 2]
      depths[child] <- depths[parent] + tree_ordered$edge.length[i]
    }
    return(depths)
  }

  depths_true <- calculate_node_depths(tree_true)
  weight_fn <- function(d) 1 / (1 + d)

  bp_true <- prop.part(tree_true)
  bp_infer <- prop.part(tree_infer)
  n_tips <- Ntip(tree_true)

  total_weighted_penalty <- 0
  max_possible_penalty <- 0

  for (split in bp_true) {
    if (length(split) %in% c(0, n_tips)) next
    node_id <- getMRCA(tree_true, tree_true$tip.label[split])
    weight <- weight_fn(depths_true[node_id])
    max_possible_penalty <- max_possible_penalty + weight
    if (!any(sapply(bp_infer, function(s) setequal(split, s)))) {
      total_weighted_penalty <- total_weighted_penalty + weight
    }
  }

  return(if (max_possible_penalty == 0) 0 else total_weighted_penalty / max_possible_penalty)
}

# Generation-Aware Divergence (Temporal Trend Analysis)
generation_aware_divergence <- function(tree_true, tree_infer) {
  # Calculate generation (# edges from root) for each node
  calculate_node_generations <- function(tree) {
    n_tips <- Ntip(tree)
    n_nodes <- tree$Nnode
    generations <- rep(0, n_tips + n_nodes)
    root <- n_tips + 1

    tree_ordered <- reorder(tree, "pruningwise")
    for (i in rev(seq_len(nrow(tree_ordered$edge)))) {
      parent <- tree_ordered$edge[i, 1]
      child <- tree_ordered$edge[i, 2]
      generations[child] <- generations[parent] + 1
    }
    return(generations)
  }

  generations_true <- calculate_node_generations(tree_true)

  bp_true <- prop.part(tree_true)
  bp_infer <- prop.part(tree_infer)
  n_tips <- Ntip(tree_true)

  # Track errors by generation
  generation_errors <- list()
  generation_totals <- list()

  for (split in bp_true) {
    if (length(split) %in% c(0, n_tips)) next
    node_id <- getMRCA(tree_true, tree_true$tip.label[split])
    gen <- generations_true[node_id]
    gen_str <- as.character(gen)

    if (is.null(generation_totals[[gen_str]])) {
      generation_totals[[gen_str]] <- 0
      generation_errors[[gen_str]] <- 0
    }
    generation_totals[[gen_str]] <- generation_totals[[gen_str]] + 1

    # Check if split is missing in inferred tree
    if (!any(sapply(bp_infer, function(s) setequal(split, s)))) {
      generation_errors[[gen_str]] <- generation_errors[[gen_str]] + 1
    }
  }

  # Compute error rates per generation
  generations <- as.numeric(names(generation_totals))
  if (length(generations) < 2) return(list(slope = 0, pvalue = 1, mean_error = 0,
                                           early_error = 0, late_error = 0))

  error_rates <- sapply(names(generation_totals), function(g) {
    generation_errors[[g]] / generation_totals[[g]]
  })

  # Linear regression: error_rate ~ generation
  tryCatch({
    lm_fit <- lm(error_rates ~ generations)
    beta_1 <- coef(lm_fit)[2]
    pval <- summary(lm_fit)$coefficients[2, 4]

    # Early vs Late errors
    g_max <- max(generations)
    early_gens <- generations[generations <= floor(g_max / 3)]
    late_gens <- generations[generations >= ceiling(2 * g_max / 3)]

    early_error <- if (length(early_gens) > 0) {
      mean(error_rates[as.character(early_gens)])
    } else { 0 }

    late_error <- if (length(late_gens) > 0) {
      mean(error_rates[as.character(late_gens)])
    } else { 0 }

    return(list(
      slope = beta_1,
      pvalue = pval,
      mean_error = mean(error_rates),
      early_error = early_error,
      late_error = late_error
    ))
  }, error = function(e) {
    return(list(slope = 0, pvalue = 1, mean_error = mean(error_rates),
                early_error = 0, late_error = 0))
  })
}

# --- Cell-Type-Specific Metrics ---
celltype_specific_metrics <- function(tree_true, tree_infer, gt_tip_states, infer_tip_states) {
  # Get common tips
  common_tips <- intersect(tree_true$tip.label, tree_infer$tip.label)
  tree_true_pruned <- keep.tip(tree_true, common_tips)
  tree_infer_pruned <- keep.tip(tree_infer, common_tips)

  # Get cell types for common tips
  gt_states <- gt_tip_states[common_tips]
  unique_celltypes <- unique(gt_states[!is.na(gt_states)])

  results <- list()

  for (celltype in unique_celltypes) {
    # Get tips of this cell type
    celltype_tips <- names(gt_states)[gt_states == celltype & !is.na(gt_states)]

    if (length(celltype_tips) < 3) {
      # Not enough tips for meaningful metrics
      results[[paste0("CT_", celltype, "_RF")]] <- NA
      results[[paste0("CT_", celltype, "_PathCorr")]] <- NA
      results[[paste0("CT_", celltype, "_MonoScore")]] <- NA
      next
    }

    # 1. RF Distance for subtree containing only this cell type
    tryCatch({
      if (length(celltype_tips) > 2) {
        subtree_true <- keep.tip(tree_true_pruned, celltype_tips)
        subtree_infer <- keep.tip(tree_infer_pruned, celltype_tips)
        rf_dist <- RF.dist(subtree_true, subtree_infer, normalize = TRUE)
        results[[paste0("CT_", celltype, "_RF")]] <- rf_dist
      } else {
        results[[paste0("CT_", celltype, "_RF")]] <- NA
      }
    }, error = function(e) {
      results[[paste0("CT_", celltype, "_RF")]] <- NA
    })

    # 2. Path correlation for this cell type
    tryCatch({
      coph_true <- cophenetic(tree_true_pruned)
      coph_infer <- cophenetic(tree_infer_pruned)

      celltype_pairs_idx <- which(outer(rownames(coph_true), rownames(coph_true),
                                        function(x, y) x %in% celltype_tips & y %in% celltype_tips))

      if (length(celltype_pairs_idx) > 1) {
        corr <- cor(as.vector(coph_true)[celltype_pairs_idx],
                    as.vector(coph_infer)[celltype_pairs_idx],
                    method = "spearman", use = "complete.obs")
        results[[paste0("CT_", celltype, "_PathCorr")]] <- corr
      } else {
        results[[paste0("CT_", celltype, "_PathCorr")]] <- NA
      }
    }, error = function(e) {
      results[[paste0("CT_", celltype, "_PathCorr")]] <- NA
    })

    # 3. Monophyly Score: How well are cells of this type clustered together?
    tryCatch({
      # Check if celltype forms a monophyletic clade in true tree
      is_mono_true <- is.monophyletic(tree_true_pruned, celltype_tips)
      is_mono_infer <- is.monophyletic(tree_infer_pruned, celltype_tips)

      # Score: 1 if both monophyletic or both not, 0.5 if one is, 0 adjusted
      if (is_mono_true && is_mono_infer) {
        mono_score <- 1.0
      } else if (!is_mono_true && !is_mono_infer) {
        mono_score <- 1.0  # Both non-monophyletic is still correct
      } else {
        mono_score <- 0.0  # Mismatch
      }

      results[[paste0("CT_", celltype, "_MonoScore")]] <- mono_score
    }, error = function(e) {
      results[[paste0("CT_", celltype, "_MonoScore")]] <- NA
    })
  }

  # Overall cell-type accuracy: Average RF across all cell types
  ct_rf_values <- sapply(unique_celltypes, function(ct) {
    results[[paste0("CT_", ct, "_RF")]]
  })
  results$CellType_AvgRF <- mean(ct_rf_values, na.rm = TRUE)

  # Variance in cell-type performance
  results$CellType_RF_Variance <- var(ct_rf_values, na.rm = TRUE)

  return(results)
}
