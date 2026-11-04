lfqpca_weighted_orthonormalize <- function(loadings, weights, eps = 1e-10) {
  loadings <- as.matrix(loadings)
  if (!ncol(loadings)) return(loadings)
  G <- lfqpca_weighted_gram(loadings, weights)
  loadings %*% lfqpca_sym_power(G, -0.5, eps = eps)
}
