#' Sample graph with 10 nodes
#'
#' A random graph generated with 10 nodes.
#'
#' @format igraph object
"node10"

#' Normalized adjacency matrix using node10 object
#'
#' A column-wise normalized adjacency sparse matrix
#'
#' @format Matrix sparseMatrix
"node10_norm"

#' Example results from RWR
#'
#' Results generated from RWR with example graph and random initial scores
#'
#' @format a list include whether RWR converges within given iterations, and a vector for the post-scores for each node
"node10_rwr"

#' Normalized adjacency matrix for the largest connected community from human String database
#'
#' A sparse matrix for the largest connected community derived from human String database. The connections were pruned to only retain
#'  confidence scores greater than 900 (or 0.9). The connections were subsequently column-wisely normalized and converted into sparse matrix.
#'
#' @format sparse matrix
#' @source \url{https://stringdb-downloads.org/}
"stringdb_human_norm"

#' Example data for RWR input
#'
#' A dataframe to test the RWR function. In practice, only use the second column in the function. Genes are ordered to match the graph object.
#'
#' @format data frame: gene_name score
"testInput"
