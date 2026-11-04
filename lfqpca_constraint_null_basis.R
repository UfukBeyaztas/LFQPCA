lfqpca_constraint_null_basis <- function(C_i, tol = 1e-10) {
  C_i <- as.matrix(C_i)
  J <- nrow(C_i)
  if (J < 1) return(matrix(0, 0, 0))
  if (!ncol(C_i)) return(diag(J))
  
  qC <- qr(C_i, tol = tol, LAPACK = FALSE)
  rank_C <- qC$rank
  
  if (rank_C < min(2, ncol(C_i))) {
    stop(
      "Subject time design is rank deficient; distinct visit times are needed.",
      call. = FALSE
    )
  }
  if (J <= rank_C) return(matrix(0, J, 0))
  
  Q <- qr.Q(qC, complete = TRUE)
  N <- Q[, (rank_C + 1):J, drop = FALSE]
  
  err <- max(abs(crossprod(C_i, N)))
  if (!is.finite(err) || err > 1e-7) {
    stop("Failed to construct the within-subject null-space basis.",
         call. = FALSE)
  }
  N
}