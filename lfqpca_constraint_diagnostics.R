lfqpca_constraint_diagnostics <- function(zeta, subject_id, H) {
  zeta <- as.matrix(zeta)
  H <- as.matrix(H)
  I <- max(subject_id)
  out <- matrix(0, I, 3)
  colnames(out) <- c(
    "max_abs_intercept_sum",
    "max_abs_linear_time_sum",
    "max_abs_all_time_basis"
  )
  
  for (i in seq_len(I)) {
    rows <- which(subject_id == i)
    if (!ncol(zeta)) next
    violation <- crossprod(H[rows, , drop = FALSE], zeta[rows, , drop = FALSE])
    out[i, 1] <- max(abs(violation[1, , drop = FALSE]))
    out[i, 2] <- if (nrow(violation) >= 2) {
      max(abs(violation[2, , drop = FALSE]))
    } else 0
    out[i, 3] <- max(abs(violation))
  }
  out
}