lfqpca_sim_orthonormal_basis <- function(F, w, tol = 1e-9) {
  F <- as.matrix(F)
  if (!ncol(F)) return(matrix(0, nrow(F), 0))
  A <- sqrt(w) * F
  sv <- svd(A)
  r <- sum(sv$d > max(dim(A)) * max(sv$d) * tol)
  if (!r) return(matrix(0, nrow(F), 0))
  sweep(sv$u[, seq_len(r), drop = FALSE], 1, sqrt(w), "/")
}
