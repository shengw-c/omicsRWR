## test graph: node10, node10_norm
set.seed(111)
node10 <- sample_gnp(n = 10, p = 0.5, directed = FALSE)
E(node10)$weight <- runif(ecount(node10), min = 1, max = 10)
node10_norm = norm_W(node10)

## test RWR
set.seed(222)
scores = runif(ncol(node10_norm))
node10_rwr = runRWR_normAdj(values = scores, adj_mat = node10_norm)

usethis::use_data(node10, node10_norm, node10_rwr, overwrite = TRUE)
