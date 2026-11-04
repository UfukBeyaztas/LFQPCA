lfqpca_sim_project <- function(F, basis, w, tol = 1e-10) {
  Q <- lfqpca_sim_orthonormal_basis(basis, w, tol)
  if (!ncol(Q)) stop("The projection basis has rank zero.", call. = FALSE)
  Q %*% crossprod(Q, sweep(as.matrix(F), 1, w, "*"))
}