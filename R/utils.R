#' Normalize adjacency matrix from igraph graph object
#'
#' @param graph igraph graph object
#'
#' @return A column-wise normalized matrix
#' @export
#'
#' @importFrom igraph as_adjacency_matrix
#' @importFrom igraph sample_gnp
#' @importFrom Matrix Matrix
#'
#' @examples
#' node10
#' norm_W(node10)
norm_W <- function(graph) {
  adj_matrix <- as_adjacency_matrix(graph, sparse = FALSE, attr = "weight")

  ## calculate column sums
  D = diag(Cpp_colSums(adj_matrix))
  if (any(diag(D) == 0)) {
    warning("Zero-degree nodes present. Consider removing or adjusting them.")
  }

  ## Inverse of Degree Matrix (D^-1)
  D_inv = Arma_Inv(D)

  ##column-wise normalization: W = AD^-1
  W = adj_matrix %*% D_inv
  W_sparse = Matrix(W, sparse = TRUE)

  return(W_sparse)
}

#' Run random walk with restart (RWR)
#'
#' @param values a vector for initial scores, the orders should match node orders
#' @param adj_mat column-wise normalized adjacency (sparse) matrix
#' @param int_scores NULL or vector, intermediate scores if run the function on the fly
#' @param restart_prop 0.15, restart probability
#' @param num_iter maximum iterations to process
#' @param delta minimal differences for convergence
#'
#' @importFrom methods is
#'
#' @return a list include (1) whether RWR converges in defined iterations; (2) output vector after RWR
#' @export
#'
#' @examples
#' scores = runif(ncol(node10_norm))
#' runRWR_normAdj(values = scores, adj_mat = node10_norm)
runRWR_normAdj <- function(values, adj_mat, int_scores = NULL, restart_prop = 0.15, num_iter = 1000, delta = 1e-6){
  if(length(which(is(adj_mat)=="sparseMatrix"))!=0) {
    adj_mat = as.matrix(adj_mat)
  }
  tmp = rwr_cpp(adj_mat = adj_mat,
                init_scores = values,
                int_scores = int_scores,
                restart_prop = restart_prop,
                num_iter = num_iter,
                delta = delta)
  if (!tmp$isconverged) {
    warnings("RWR did not converge. Consider increasing the number of iterations.")
  }
  return(tmp)
}
