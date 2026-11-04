lfqpca_safe_independent_columns <- function(X, tol = 1e-9) {
  X <- as.matrix(X)
  if (!ncol(X)) return(integer(0))
  qx <- qr(X, tol = tol, LAPACK = FALSE)
  if (qx$rank <= 0) return(integer(0))
  sort(qx$pivot[seq_len(qx$rank)])
}