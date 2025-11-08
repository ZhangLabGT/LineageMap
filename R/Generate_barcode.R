GenerateBarcode <- function(ncells, p_d = 1, Sigma = 0.5,mu = 0.1, N_char = 16, N_ms = 100, unif_on = FALSE, T_cell = NULL){
  T_cell$edge.length <- rep(1, length(T_cell$edge[,1]))
  N_nodes <- length(T_cell$edge[,1])+1
  cell_edges <- cbind(T_cell$edge,T_cell$edge.length)
  cell_edges <- cbind(c(1:length(cell_edges[,1])),cell_edges)
  cell_connections <- table(c(cell_edges[,2],cell_edges[,3]))
  cell_root <- as.numeric(names(cell_connections)[cell_connections==2])
  Node_cell <- cell_root


  root_barcode <- rep(0,N_char)
  neutral <- SampleBarcode(Node_cell,0,cif_center, edges = cell_edges,p_a = p_a,p_d = p_d, mu = mu, barcode = root_barcode, N_ms = N_ms, unif_on = unif_on)
  muts <- neutral[,4:length(neutral[1,])]

  muts[muts == Inf] <- '-'

  return(muts)
}

SampleBarcode <- function(par,depth,edges, muts = NULL,p_d = 0.1, mu = 0.1, p_a = 0.8,barcode = NULL, N_ms = NULL, unif_on = FALSE){
  children <- edges[edges[,2]==par,3] # get the children of the current node
  result<-lapply(c(1:length(children)),function(j){
    edge<-edges[edges[,2]==par & edges[,3]==children[j],] # given the parent and child, find the edge
    if(sum(edges[,2]==children[j])==0){
      result <- SampleEdgeBarcode(edge,depth,edges,mu = mu,p_d = p_d, barcode = barcode,N_ms = N_ms, unif_on = unif_on)
      result <- result[c(1:(length(result[,1]-1))),]
    }else{
      result <- SampleEdgeBarcode(edge,depth,edges,mu = mu,p_d = p_d, barcode = barcode,N_ms = N_ms, unif_on = unif_on)
      barcode <- result[length(result[,1]),4:length(result[1,])]
      result <- result[c(1:(length(result[,1]-1))),]
      depth <- depth + edge[4]
      result1 <- SampleBarcode(children[j],depth,edges,mu = mu,p_d = p_d, barcode = barcode,N_ms = N_ms, unif_on = unif_on)
      result <- rbind(result,result1)
    }
    return(result)
  })
  #browser()
  result<-do.call(rbind,result)
  result <-result[!duplicated(result[,2],fromLast = TRUE),]
  result<- result[order(result[,2]),]
  return(result)
}


SampleEdgeBarcode <- function(edge,depth,edges,mu = 0.1,p_d = 0, barcode = NULL, N_ms = NULL, unif_on = FALSE){
  #browser()
  child_barcode <- barcode
  state_dist <- Mutated_state_dist(N_ms, cm)
  child_barcode <- generate_mutation(child_barcode,mu = mu,p_d = p_d,N_ms = N_ms, mutation_dist = state_dist, unif_on = unif_on)
  barcodes <- rbind(barcode,child_barcode)
  result<-cbind(c(depth,depth+1),barcodes)

  rownames(result) <- c()
  result <- cbind(rep(edge[2],length(result[,1])),rep(edge[3],length(result[,1])),result)
  return(result)
}

#' Generate mutated states from a real dataset
#' @param N_ms the number of mutated states required
#' @param cm character matrix of the real dataset
Mutated_state_dist <- function(N_ms,cm){
  cm_table <- table(unlist(cm))
  cm_table <- sort(cm_table,decreasing = TRUE)
  cm_table <- cm_table[dimnames(cm_table)[[1]]!='0']
  cm_table <- cm_table[dimnames(cm_table)[[1]]!='-']

  if (N_ms < length(dimnames(cm_table)[[1]])){
    cm_table <- cm_table[1:N_ms]
  }
  cm_prop <- prop.table(cm_table)
  states_dist <- cm_prop
  return(states_dist)
}

#' Generate CRISPR induced mutations
#' @param barcode input barcode (from the parent cell)
#' @param mu mutation rate per target, default is 1
#' @param N_ms the number of mutated states required
#' @param mutation_dist input distribution for the mutated states
#' @param p_d whether or not to simulate dropout effects, 0 or 1
#' @param unif_on distribution mode of mutated states. TRUE: uniform distribution; FALSE: real fitted distribution
generate_mutation <- function(barcode, mu=0.1 , N_ms = 100 ,mutation_dist = NULL, p_d = 0, unif_on = FALSE){
  if (unif_on){
    states <- as.list(seq (1,N_ms,1))
    prob_dist <- NULL
  }else{
    states <- dimnames(mutation_dist)[[1]]
    prob_dist <- as.vector(mutation_dist)
  }

  m <- length(barcode)
  child_barcode <- barcode
  mu_loc = runif(m) < mu
  mutation_cites = (child_barcode == 0) &  mu_loc
  n_mut = sum(mutation_cites)
  if (n_mut != 0) {
    child_barcode[mutation_cites] = as.integer(sample(states, n_mut, replace = T, prob = prob_dist))
    if ((n_mut >=2)&(p_d == 1)){
      child_barcode <- generate_dropout(child_barcode,mutation_cites)
    }
  }

  return (child_barcode)
}

#' Generate excision dropout
#' @param barcode input barcode with mutations
#' @param mutation_cites the mutated target sites
generate_dropout <- function(barcode,mutation_cites){
  dropout_between = sample(which(mutation_cites), 2 )
  barcode[dropout_between[1]:dropout_between[2]] <- Inf
  return(barcode)
}


#' Simulate capture dropout based on observed counts and barcode data
#' @param observed_counts Observed counts of gene expressions
#' @param muts Unprocessed mutation barcodes
#' @export
CaptureDrop <- function(observed_counts, muts) {
  #find the top 10 highly expressed genes in the observed counts
  #observed_counts <- t(observed_counts)
  n_char <- dim(muts)[2]
  observed_counts <-  log(observed_counts+1)
  gene_means <- rowMeans(observed_counts)
  gene_order <- order(gene_means,decreasing = TRUE)[1:10]
  reg <- sample(1:10,1)
  gene_reg <- observed_counts[,gene_order[reg]]
  drop_reg <- which(gene_reg == 0)
  N_drop <- length(drop_reg)
  for (i in drop_reg){
    muts[i,] <- rep('-',n_char)
  }
  print(sprintf("Barcodes of %d cells removed based on gene %d.", N_drop, gene_order[reg] ))
  return(muts)
}

#' Simulate technical dropout uniformly across rows, cols, and globally
#' @param mat Unprocessed mutation barcodes
#' @param drop_rate Dropout rate, 0-1
#' @param mode Dropout mode globally, across rows or columns
#' @param seed random seed
#' @export
drop_characters <- function(mat, drop_rate = 0.3, mode = c("global", "row", "col"), seed = NULL) {
  mode <- match.arg(mode)
  if (!is.null(seed)) set.seed(seed)

  mat_dropped <- mat
  n_row <- nrow(mat)
  n_col <- ncol(mat)

  if (mode == "global") {
    # Drop randomly across the entire matrix
    mask <- matrix(runif(length(mat)) < drop_rate, nrow = n_row)

  } else if (mode == "row") {
    # Ensure similar drop rate in each row
    mask <- matrix(FALSE, n_row, n_col)
    for (i in seq_len(n_row)) {
      mask[i, ] <- runif(n_col) < drop_rate
    }

  } else if (mode == "col") {
    # Ensure similar drop rate in each column
    mask <- matrix(FALSE, n_row, n_col)
    for (j in seq_len(n_col)) {
      mask[, j] <- runif(n_row) < drop_rate
    }
  }

  mat_dropped[mask] <- "-"
  return(mat_dropped)
}
