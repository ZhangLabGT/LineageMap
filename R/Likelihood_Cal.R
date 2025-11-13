LikelihoodCal_ST_OU <- function(tree, muts, cell_state_labels, state_lineages, loc,
                             lambda1 = 10, lambda2 = 0.1, alpha = 1.0,
                             estimate_ou_params = TRUE,                 # fit OU params if TRUE
                             init_alpha = 1.0, init_sigma2 = 1.0) {    # initial guesses
  # Dependencies
  if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    stop("Please install the 'mvtnorm' package.")
  }

  N_char <- ncol(muts)
  nodes_internal <- sort(unique(tree$edge[,1]), decreasing = TRUE)

  # extend mutation table for internal nodes
  allele_unique <- unique(unlist(muts))
  muts <- rbind(muts, matrix(0, length(nodes_internal), N_char))
  internal_lookup <- data.frame(node = numeric(), index = numeric())

  # storage for state transitions
  state_transitions <- data.frame(from = character(),
                                  to = character(),
                                  dist = numeric(),
                                  stringsAsFactors = FALSE)

  # storage for losses
  ad_loss <- 0
  spatial_loss <- 0
  l_barcode <- 0
  l_expression <- 0

  # parameters for ad/asymmetric division heuristic (kept)
  p_a <- 0.8   # asymmetric division prob
  m_r <- 0.8   # migration rate (not used in OU)

  # -------------------------
  # 1) original upward pass: mutations & state aggregation (unchanged)
  # -------------------------
  for (node in nodes_internal) {
    children_index <- tree$edge[tree$edge[,1] == node, 2]
    children <- integer()
    for (cid in children_index) {
      if (cid > length(tree$tip.label)) {
        c_row <- internal_lookup$index[internal_lookup$node == cid]
      } else {
        c_row <- cid
      }
      children <- c(children, c_row)
    }

    ## --- MESSAGE PASSING: mutations ---
    barcodes_children <- muts[children, , drop = FALSE]
    for (i in seq_len(N_char)) {
      chars <- unique(barcodes_children[, i])
      if (length(chars) == 1) {
        muts[node, i] <- chars
      } else {
        muts[node, i] <- "0"  # unresolved
      }
    }

    ## --- MESSAGE PASSING: states ---
    cell_state_children <- cell_state_labels[children]
    state_parent <- "0"
    #browser()
    for (lineage in state_lineages) {
      #print(lineage)
      idx1 <- match(cell_state_children[1], lineage)
      idx2 <- match(cell_state_children[2], lineage)
      if (is.na(idx1) || is.na(idx2)) next
      state_parent <- lineage[min(idx1, idx2)]
      state_dist <- abs(idx1 - idx2)

      # record transition events (two child transitions)
      state_transitions <- rbind(state_transitions,
                                 data.frame(from = state_parent,
                                            to = cell_state_children[1],
                                            dist = abs(match(state_parent, lineage) - idx1),
                                            stringsAsFactors = FALSE))
      state_transitions <- rbind(state_transitions,
                                 data.frame(from = state_parent,
                                            to = cell_state_children[2],
                                            dist = abs(match(state_parent, lineage) - idx2),
                                            stringsAsFactors = FALSE))
      break
    }

    # update labels + lookup
    cell_state_labels <- c(cell_state_labels, state_parent)
    internal_lookup <- rbind(internal_lookup,
                             data.frame(node = node, index = length(cell_state_labels)))

    # asymmetric division penalty
    if (length(children) > 1) {
      if (cell_state_children[1] == cell_state_children[2]) {
        ad_loss <- ad_loss + log2(1 - p_a)
      } else {
        ad_loss <- ad_loss + log2(p_a)
      }
    }
  }

  # -------------------------
  # 2) Spatial: replace Brownian by state-dependent scalar OU + estimation
  # -------------------------
  # Helpers for scalar-OU (parent-state controlled) transitions:
  OU_edge_scalar <- function(t, alpha, sigma2, mu_vec) {
    # returns list(A (scalar), b (q-vector), Se (scalar))
    if (t <= 0) {
      A  <- 1.0
      b  <- rep(0, length(mu_vec))
      Se <- 0.0
    } else {
      A  <- exp(-alpha * t)
      b  <- (1 - A) * mu_vec
      # avoid divide-by-zero if alpha very small
      Se <- if (alpha > 1e-12) (sigma2 / (2 * alpha)) * (1 - exp(-2 * alpha * t)) else sigma2 * t
    }
    list(A = A, b = b, Se = Se)
  }

  # combine children contributions for scalar-A OU (information form)
  combine_children_scalarOU <- function(children_msgs, edge_pars, q, ridge = 1e-9) {
    Lambda <- matrix(0, q, q)
    eta    <- rep(0, q)
    cache  <- vector("list", length(children_msgs))
    for (j in seq_along(children_msgs)) {
      m_c <- children_msgs[[j]]$mean
      S_c <- children_msgs[[j]]$cov
      par <- edge_pars[[j]]
      A   <- par$A
      b   <- par$b
      Se  <- par$Se
      S_ctov <- S_c + diag(Se, q)
      # stable inverse
      cholS <- tryCatch(chol(S_ctov + diag(ridge, q)), error = function(e) NULL)
      if (is.null(cholS)) {
        S_inv <- solve(S_ctov + diag(1e-8, q))
      } else {
        S_inv <- chol2inv(cholS)
      }
      # Lambda contribution (A scalar)
      Lambda_c <- (A^2) * S_inv
      eta_c    <- A * (S_inv %*% (m_c - b))
      Lambda   <- Lambda + Lambda_c
      eta      <- eta + as.vector(eta_c)
      cache[[j]] <- list(m_c = m_c, S_ctov = S_ctov, A = A, b = b)
    }
    # posterior at parent
    cholL <- tryCatch(chol(Lambda + diag(ridge, q)), error = function(e) NULL)
    if (is.null(cholL)) {
      S_par <- solve(Lambda + diag(1e-8, q))
    } else {
      S_par <- chol2inv(cholL)
    }
    m_par <- as.numeric(S_par %*% eta)
    list(mean = m_par, cov = S_par, cache = cache)
  }

  # OU upward pass that returns loglik, loc, cov; it uses current cell_state_labels (including internals)
  OU_message_passing_state_scalar <- function(tree, Y, all_states,
                                              alpha_by_state, sigma2_by_state, mu_by_state,
                                              mu0 = NULL, Sigma0 = NULL, obs_cov = NULL) {
    n_tips  <- length(tree$tip.label)
    n_nodes <- n_tips + tree$Nnode
    q <- ncol(Y)
    if (is.null(mu0)) mu0 <- colMeans(Y)
    if (is.null(Sigma0)) Sigma0 <- diag(1e6, q)
    if (is.null(obs_cov)) obs_cov <- diag(1e-6, q)

    # storage
    m_mean <- vector("list", n_nodes)
    m_cov  <- vector("list", n_nodes)

    # initialize leaves
    for (i in seq_len(n_tips)) {
      m_mean[[i]] <- as.numeric(Y[i, ])
      m_cov [[i]] <- obs_cov
    }

    E <- tree$edge
    L <- if (is.null(tree$edge.length)) rep(1, nrow(E)) else tree$edge.length
    # children list
    ch_list <- vector("list", n_nodes)
    for (e in seq_len(nrow(E))) {
      parent <- E[e,1]; child <- E[e,2]
      ch_list[[parent]] <- c(ch_list[[parent]], child)
    }
    # internal nodes in postorder (decreasing parent indices works for ape trees)
    nodes_internal <- sort(unique(E[,1]), decreasing = TRUE)
    loglik <- 0.0

    for (u in nodes_internal) {
      ch <- ch_list[[u]]
      if (length(ch) == 0) next

      s_u <- as.character(all_states[u])
      alpha <- alpha_by_state[[s_u]]
      sigma2 <- sigma2_by_state[[s_u]]
      mu_s <- as.numeric(mu_by_state[[s_u]])

      # child messages and edge params
      child_msgs <- vector("list", length(ch))
      edge_pars  <- vector("list", length(ch))
      for (j in seq_along(ch)) {
        v <- ch[j]
        idx <- which(E[,1]==u & E[,2]==v)[1]
        t  <- if (length(idx)) L[idx] else 1.0
        par <- OU_edge_scalar(t, alpha, sigma2, mu_s)
        edge_pars[[j]]  <- par
        child_msgs[[j]] <- list(mean = m_mean[[v]], cov = m_cov[[v]])
      }

      comb <- combine_children_scalarOU(child_msgs, edge_pars, q)
      m_mean[[u]] <- comb$mean
      m_cov [[u]] <- comb$cov

      # accumulate per-child log-likelihood evaluated at parent mean
      for (entry in comb$cache) {
        mean_pred <- as.numeric(entry$A * m_mean[[u]] + entry$b)
        # use mvtnorm::dmvnorm for logdensity
        logdens <- mvtnorm::dmvnorm(entry$m_c, mean = mean_pred, sigma = entry$S_ctov, log = TRUE)
        loglik <- loglik + logdens
      }
    }

    # find root (node with no parent)
    parents <- unique(E[,1]); children_all <- unique(E[,2])
    root_candidates <- setdiff(parents, children_all)
    root <- if (length(root_candidates)) root_candidates[1] else max(nodes_internal)

    S_root_int <- m_cov[[root]] + Sigma0
    loglik <- loglik + mvtnorm::dmvnorm(m_mean[[root]], mean = mu0, sigma = S_root_int, log = TRUE)

    loc_mat <- matrix(NA_real_, nrow = n_nodes, ncol = q)
    for (i in seq_len(n_nodes)) loc_mat[i, ] <- as.numeric(m_mean[[i]])

    list(loglik = as.numeric(loglik), loc = loc_mat, cov = m_cov, states = all_states)
  }

  # Prepare Y (leaf coordinates) from input loc (assumes loc has x,y for tips in order)
  N_leaf <- length(tree$tip.label)
  if (!all(c("x","y") %in% colnames(loc))) {
    stop("loc must be a data.frame with columns 'x' and 'y' for the leaves.")
  }
  Y <- matrix(NA_real_, nrow = N_leaf, ncol = 2)
  for (i in seq_len(N_leaf)) {
    Y[i,1] <- loc$x[i]; Y[i,2] <- loc$y[i]
  }

  # Build list of states for all nodes: current cell_state_labels has leaves + internals appended earlier.
  # Ensure it's length equals number of rows in muts (nodes)
  n_nodes_total <- nrow(muts)
  if (length(cell_state_labels) != n_nodes_total) {
    stop("cell_state_labels length must equal number of rows in muts after extension.")
  }
  states_all <- factor(cell_state_labels)
  state_levels <- levels(states_all)

  # Initial per-state parameters or fit them by MLE
  # Parameter containers (named lists)
  alpha_by_state <- list(); sigma2_by_state <- list(); mu_by_state <- list()
  for (s in state_levels) {
    alpha_by_state[[s]] <- init_alpha
    sigma2_by_state[[s]] <- init_sigma2
    # initialize mu as mean of leaf coords in that state (fallback to global mean)
    idx <- which(as.character(states_all[seq_len(N_leaf)]) == s)
    if (length(idx) > 0) {
      mu_by_state[[s]] <- colMeans(Y[idx, , drop = FALSE])
    } else {
      mu_by_state[[s]] <- colMeans(Y)
    }
  }

  if (estimate_ou_params) {
    # pack/unpack helpers
    pack_params <- function(alpha_by_state, sigma2_by_state, mu_by_state, state_levels) {
      v <- numeric()
      for (s in state_levels) {
        v <- c(v, log(alpha_by_state[[s]]), log(sigma2_by_state[[s]]), as.numeric(mu_by_state[[s]]))
      }
      names(v) <- NULL
      v
    }
    unpack_params <- function(theta, state_levels, q) {
      res_alpha <- list(); res_sigma2 <- list(); res_mu <- list()
      off <- 1
      for (s in state_levels) {
        log_alpha <- theta[off]; log_sigma2 <- theta[off+1]
        mu_s <- theta[(off+2):(off+1+q)]
        res_alpha[[s]] <- exp(log_alpha)
        res_sigma2[[s]] <- exp(log_sigma2)
        res_mu[[s]] <- as.numeric(mu_s)
        off <- off + 2 + q
      }
      list(alpha = res_alpha, sigma2 = res_sigma2, mu = res_mu)
    }

    # objective: negative spatial log-likelihood (OU message passing)
    neg_spatial_nll <- function(theta) {
      pars <- unpack_params(theta, state_levels, q = 2)
      res <- OU_message_passing_state_scalar(tree, Y, states_all,
                                             pars$alpha, pars$sigma2, pars$mu,
                                             mu0 = NULL, Sigma0 = NULL, obs_cov = diag(1e-6,2))
      # return negative loglik
      return(-res$loglik)
    }

    theta0 <- pack_params(alpha_by_state, sigma2_by_state, mu_by_state, state_levels)
    # use L-BFGS-B with box constraints for numerical stability
    opt <- optim(theta0, neg_spatial_nll, method = "L-BFGS-B", control = list(maxit = 200))
    fitted <- unpack_params(opt$par, state_levels, q = 2)
    alpha_by_state <- fitted$alpha
    sigma2_by_state <- fitted$sigma2
    mu_by_state <- fitted$mu
    # optional: you can inspect opt$convergence / opt$value
  }

  # Run OU message passing with (fitted) params
  ou_res <- OU_message_passing_state_scalar(tree, Y, states_all,
                                            alpha_by_state, sigma2_by_state, mu_by_state,
                                            mu0 = colMeans(Y), Sigma0 = diag(1e3,2),
                                            obs_cov = diag(1e-6,2))
  spatial_loglik <- ou_res$loglik
  # base-2 loss like your earlier code
  spatial_loss <- -spatial_loglik / log(2)

  # reconstruct loc_data as a data.frame: rows are tips (1..N_leaf) then internals in increasing node index order
  loc_mat <- ou_res$loc  # rows 1..n_nodes (tips then internals by ape convention)
  loc_df <- data.frame(x = loc_mat[,1], y = loc_mat[,2])

  # -------------------------
  # 3) state transition (HMM-style) expression likelihood (unchanged)
  # -------------------------
  if (nrow(state_transitions) > 0) {
    for (lineage in state_lineages) {
      for (i in seq_len(nrow(state_transitions))) {
        if (!(state_transitions$from[i] %in% lineage &&
              state_transitions$to[i] %in% lineage)) next
        idx_from <- match(state_transitions$from[i], lineage)
        all_dists <- abs(seq_along(lineage) - idx_from)
        probs <- exp(-alpha * all_dists)
        probs <- probs / sum(probs)
        idx_to <- match(state_transitions$to[i], lineage)
        l_expression <- l_expression + log(probs[idx_to] + 1e-12)
      }
    }
  }

  # -------------------------
  # 4) barcode term (kept as previous placeholder l_barcode)
  #    If you have a barcode likelihood computation, plug it in here by filling l_barcode.
  # -------------------------
  # l_barcode remains as computed earlier (still zero unless you add code)

  # Final additive likelihood (loss)
  likelihood <- l_barcode + lambda1 * l_expression + lambda2 * spatial_loss

  return(list(likelihood = likelihood,
              loc_data = loc_df,
              cell_state_infer = cell_state_labels,
              muts_internal = muts,
              transitions = state_transitions,
              # also return OU parameters for inspection
              ou_params = list(alpha = alpha_by_state,
                               sigma2 = sigma2_by_state,
                               mu = mu_by_state),
              ou_loglik = spatial_loglik,
              ou_res = ou_res))
}

#' tree local search function, subtree swapping
#' @param tree lineage tree object
#' @param muts lineage barcode data
#' @param loc spatial location data
#' @param prob_internal probability to sample every internal node
#' @import TreeTools
#' @export
TreeLC2 <- function(tree,muts,loc,prob_internal = NULL){

  tree <- Preorder(tree)
  n_nodes <- tree$Nnode
  n_tips <- length(tree$tip.label)
  total_nodes <- n_nodes + n_tips
  tree$node.label <-paste0("node", seq(1:tree$Nnode))
  descendants <- Descendants(tree)

  if (is.null(prob_internal)){
    prob_internal <- c()
    for (i in 1:n_nodes){
      node <- n_tips + i
      #barcode_children <- muts[tree$edge[tree$edge[,1] == node,2],]
      #barcode_node <- muts[node,]
      #char_diff <- rowSums(barcode_children != barcode_node[col(barcode_children)])

      # spatial dispersion for descendants
      desc <- descendants[[node]]
      desc <- desc[desc <= n_tips]
      spatial_disp <- if (length(desc) > 1) mean(dist(loc[desc,])) else 0

      # combine (example: dispersion + mismatch)
      prob_internal <- c(prob_internal,spatial_disp)
    }

    prob_internal <- prob_internal / sum(prob_internal)

    prob_internal <- c(dunif(1:n_tips,1/n_tips,n_tips),prob_internal)
    #prob_internal <- dunif(1:(n_nodes + n_tips),1,n_nodes + n_tips)
  }

  subtree_swap <- sample(tree$Nnode + length(tree$tip.label),2,prob = prob_internal)
  subtree1 <- Subtree(tree,subtree_swap[1])
  subtree2 <- Subtree(tree,subtree_swap[2])

  repeat {
    subtree_swap <- sample(total_nodes, 2, prob = prob_internal)
    subtree1 <- Subtree(tree, subtree_swap[1])
    subtree2 <- Subtree(tree, subtree_swap[2])
    if (length(intersect(subtree1$tip.label, subtree2$tip.label)) == 0) break
  }

  subtree1$name <- subtree_swap[1]
  subtree1$edge.length <- rep(1, nrow(subtree1$edge))
  subtree2$name <- subtree_swap[2]
  subtree2$edge.length <- rep(1, nrow(subtree2$edge))

  root1 <- tree$edge[tree$edge[,2]==subtree1$name,1]
  root2 <- tree$edge[tree$edge[,2]==subtree2$name,1]

  if (root1 > root2) {
    reg <- subtree2
    subtree2 <- subtree1
    subtree1 <- reg
  }
  root1 <- tree$edge[tree$edge[,2]==subtree1$name,1]
  root2 <- tree$edge[tree$edge[,2]==subtree2$name,1]
  root1 <- paste0("node",root1-Ntip(tree))
  root2 <- paste0("node",root2-Ntip(tree))
  subtree1$root.edge <- 1
  subtree2$root.edge <- 1
  Ntip1 <- subtree1$Ntip
  Ntip2 <- subtree2$Ntip

  if (length(subtree1$tip.label)>1){
    tree_prune1 <- drop.tip(tree, subtree1$tip.label, trim.internal = FALSE, subtree = FALSE,root.edge = 1)
  } else{
    tree_prune1 <- tree
    tree_prune1$tip.label[tree_prune1$tip.label == subtree1$tip.label] <- "NA"
  }
  if (length(subtree2$tip.label)>1){
    tree_prune2 <- drop.tip(tree_prune1, subtree2$tip.label, trim.internal = FALSE, subtree = FALSE,root.edge = 1)
  } else{
    tree_prune2 <- tree_prune1
    tree_prune2$tip.label[tree_prune2$tip.label == subtree2$tip.label] <- "NA"
  }
  #tree_new1 <- bind.tree(tree_prune1, subtree1, where = subtree2$node.label[1]-Nnode1)
  #tree_prune2 <- drop.tip(tree_prune1, subtree2$tip.label, trim.internal = FALSE, subtree = FALSE,root.edge = 1)
  #print(regraft_loc)
  #N_na <- sum(tree_prune2$tip.label == 'NA')
  tree_new1 <- bind.tree(tree_prune2, subtree2, where = Ntip(tree_prune2)+match(root1,tree_prune2$node.label))
  tree_new2 <- bind.tree(tree_new1, subtree1, where = Ntip(tree_new1)+match(root2,tree_new1$node.label))
  tree_final <- drop.tip(tree_new2, tree_new2$tip.label[grepl("node", tree_new2$tip.label, fixed = TRUE)], trim.internal = TRUE)
  tree_final <- drop.tip(tree_final, "NA")
  tree_final$node.label <- NULL
  tree_final$edge.length <- rep(1, length(tree_final$edge.length))
  #plot(tree_final)
  return(tree_final)
}

#'
safe_solve <- function(A) {
  if (is.null(dim(A))) {
    A <- matrix(A, 1, 1)
  }
  return(solve(A))
}


barcode_loglik <- function(tree, muts, mu = 0.1) {
  N_sites <- ncol(muts)
  N_nodes <- nrow(muts)
  N_leaves <- length(tree$tip.label)

  edge_len <- if (!is.null(tree$edge.length)) tree$edge.length else rep(1, nrow(tree$edge))

  # Recursive likelihood function
  compute_likelihood <- function(node, site) {
    if (node <= N_leaves) {
      # Leaf: observed state
      obs <- muts[node, site]
      if (obs == 0) return(c(1, 0))
      else if (obs == 1) return(c(0, 1))
      else return(c(1, 1))  # unresolved / missing
    } else {
      # Internal: combine children
      children <- tree$edge[tree$edge[,1] == node, 2]
      L <- c(1, 1)  # initialize
      for (child in children) {
        e_idx <- which(tree$edge[,1] == node & tree$edge[,2] == child)[1]
        t <- edge_len[e_idx]

        # Transition probabilities
        P <- matrix(c(
          exp(-mu*t),        0,
          1 - exp(-mu*t),    1
        ), nrow = 2, byrow = TRUE)

        # Child likelihood
        L_child <- compute_likelihood(child, site)

        # Message passing
        L <- L * as.vector(P %*% L_child)
      }
      return(L)
    }
  }

  # Accumulate across sites
  loglik <- 0
  root <- max(tree$edge[,1])
  for (site in 1:N_sites) {
    L_root <- compute_likelihood(root, site)
    # Prior root = always unmutated at t=0
    prob_site <- L_root[1]
    loglik <- loglik + log(prob_site + 1e-12)
  }

  return(loglik / log(2))  # convert to log2 like others
}

LikelihoodCal_ST <- function(tree, muts, cell_state_labels, state_lineages, loc,
                             lambda1 = 10, lambda2 = 0.1, alpha = 1.0) {
  loc$id <- 1:nrow(loc)
  N_char <- ncol(muts)
  nodes_internal <- sort(unique(tree$edge[,1]), decreasing = TRUE)

  # extend mutation table for internal nodes
  allele_unique <- unique(unlist(muts))
  muts <- rbind(muts, matrix(0, length(nodes_internal), N_char))
  internal_lookup <- data.frame(node = numeric(), index = numeric())

  # storage for state transitions
  state_transitions <- data.frame(from = character(),
                                  to = character(),
                                  dist = numeric(),
                                  stringsAsFactors = FALSE)

  # storage for losses
  ad_loss <- 0
  spatial_loss <- 0
  l_barcode <- 0
  l_expression <- 0

  # parameters
  p_a <- 0.8   # asymmetric division prob
  m_r <- 0.8   # migration rate

  for (node in nodes_internal) {
    children_index <- tree$edge[tree$edge[,1] == node, 2]
    children <- integer()
    for (cid in children_index) {
      if (cid > length(tree$tip.label)) {
        c_row <- internal_lookup$index[internal_lookup$node == cid]
      } else {
        c_row <- cid
      }
      children <- c(children, c_row)
    }

    ## --- MESSAGE PASSING: mutations ---
    barcodes_children <- muts[children, ]
    for (i in seq_len(N_char)) {
      chars <- unique(barcodes_children[, i])
      if (length(chars) == 1) {
        muts[node, i] <- chars
      } else {
        muts[node, i] <- "0"  # unresolved
      }
    }
    l_barcode <- barcode_loglik(tree, muts, mu = 0.1)

    ## --- MESSAGE PASSING: states ---
    cell_state_children <- cell_state_labels[children]
    state_parent <- "0"

    for (lineage in state_lineages) {
      idx1 <- match(cell_state_children[1], lineage)
      idx2 <- match(cell_state_children[2], lineage)
      if (is.na(idx1) || is.na(idx2)) next
      state_parent <- lineage[min(idx1, idx2)]
      state_dist <- abs(idx1 - idx2)

      # record transition events
      state_transitions <- rbind(state_transitions,
                                 data.frame(from = state_parent,
                                            to = cell_state_children[1],
                                            dist = abs(match(state_parent, lineage) - idx1)))
      state_transitions <- rbind(state_transitions,
                                 data.frame(from = state_parent,
                                            to = cell_state_children[2],
                                            dist = abs(match(state_parent, lineage) - idx2)))
      break
    }

    # update labels + lookup
    cell_state_labels <- c(cell_state_labels, state_parent)
    internal_lookup <- rbind(internal_lookup,
                             data.frame(node = node, index = length(cell_state_labels)))

    # asymmetric division penalty
    if (length(children) > 1) {
      if (cell_state_children[1] == cell_state_children[2]) {
        ad_loss <- ad_loss + log2(1 - p_a)
      } else {
        ad_loss <- ad_loss + log2(p_a)
      }
    }
  }

  ## --- MESSAGE PASSING: locations ---
  # --- Parameters for Brownian spatial model ---
  sigma2 <- 1.0   # diffusion variance per unit branch length
  eps2   <- 1e-3  # observation noise variance at leaves
  mu0    <- c(mean(loc$x[seq_along(tree$tip.label)]),
              mean(loc$y[seq_along(tree$tip.label)])) # weak prior mean
  Sigma0 <- diag(1e3, 2)  # weak prior covariance

  # Pre-allocate per-node Gaussian messages (mean, cov)
  n_nodes <- nrow(muts)  # leaves + internals after your rbind
  m_mean  <- vector("list", n_nodes)  # each is length-2
  m_cov   <- vector("list", n_nodes)  # each is 2x2

  # Initialize leaf messages
  N_leaf <- length(tree$tip.label)
  for (i in seq_len(N_leaf)) {
    m_mean[[i]] <- c(loc$x[i], loc$y[i])
    m_cov[[i]]  <- diag(eps2, 2)
  }

  # Helper: log N(x; m, S)
  log_dmvnorm_iso <- function(x, m, S) {
    # S is 2x2; use Cholesky for stability
    L <- chol(S)
    z <- backsolve(L, forwardsolve(t(L), x - m, upper.tri = TRUE, transpose = TRUE))
    quad <- sum(z * z)
    logdet <- 2 * sum(log(diag(L)))
    return(-0.5 * (2 * log(2*pi) + logdet + quad))
  }

  # Upward Brownian message passing and log-likelihood accumulation
  spatial_loglik <- 0.0

  nodes_internal <- sort(unique(tree$edge[,1]), decreasing = TRUE)
  edge_len <- if (!is.null(tree$edge.length)) tree$edge.length else rep(1, nrow(tree$edge))

  for (node in nodes_internal) {
    child_rows <- tree$edge[tree$edge[,1] == node, , drop = FALSE]
    children <- child_rows[, 2]
    # Combine child Gaussians into parent via precision summation
    Lambda_v <- matrix(0, 2, 2)
    eta_v    <- c(0, 0)
    for (j in seq_along(children)) {
      cnode <- children[j]
      # t_vc from edge list
      e_idx <- which(tree$edge[,1] == node & tree$edge[,2] == cnode)[1]
      t_vc  <- edge_len[e_idx]
      S_c_to_v <- m_cov[[cnode]] + diag(sigma2 * t_vc, 2)
      # contribution to parent info
      #print(S_c_to_v)
      S_inv <- safe_solve(S_c_to_v)
      Lambda_v <- Lambda_v + S_inv
      eta_v    <- eta_v + S_inv %*% m_mean[[cnode]]
    }
    S_v <- solve(Lambda_v)
    m_v <- as.vector(S_v %*% eta_v)
    m_mean[[node]] <- m_v
    m_cov [[node]] <- S_v

    temp <- data.frame(x = m_v[1],y = m_v[2],id = node)
    loc <- rbind(loc,temp)

    # Add child likelihood contributions at the parent mean
    for (j in seq_along(children)) {
      cnode <- children[j]
      e_idx <- which(tree$edge[,1] == node & tree$edge[,2] == cnode)[1]
      t_vc  <- edge_len[e_idx]
      S_c_to_v <- m_cov[[cnode]] + diag(sigma2 * t_vc, 2)
      spatial_loglik <- spatial_loglik + log_dmvnorm_iso(m_mean[[cnode]], m_v, S_c_to_v)
    }
  }

  # Root prior integration
  root <- max(nodes_internal)
  S_root_int <- m_cov[[root]] + Sigma0
  spatial_loglik <- spatial_loglik + log_dmvnorm_iso(m_mean[[root]], mu0, S_root_int)

  # If you want base-2 logs to match the rest:
  spatial_loss <- spatial_loglik / log(2)

  ## --- Compute proper state transition likelihood ---
  if (nrow(state_transitions) > 0) {
    for (lineage in state_lineages) {
      for (i in seq_len(nrow(state_transitions))) {
        if (!(state_transitions$from[i] %in% lineage &&
              state_transitions$to[i] %in% lineage)) next
        # distances from "from" to all possible states
        idx_from <- match(state_transitions$from[i], lineage)
        all_dists <- abs(seq_along(lineage) - idx_from)
        probs <- exp(-alpha * all_dists)
        probs <- probs / sum(probs)

        idx_to <- match(state_transitions$to[i], lineage)
        l_expression <- l_expression + log(probs[idx_to] + 1e-12)
      }
    }
  }

  ## --- Final likelihood ---
  likelihood <- l_barcode + lambda1 * l_expression + lambda2 * spatial_loss

  return(list(likelihood = likelihood,
              cell_state_infer = cell_state_labels,
              loc_data = loc,
              muts_internal = muts,
              transitions = state_transitions))
}
