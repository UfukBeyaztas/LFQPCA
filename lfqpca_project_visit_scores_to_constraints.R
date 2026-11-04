lfqpca_project_visit_scores_to_constraints <- function(zeta, subject_id, H) {
  zeta <- as.matrix(zeta)
  if (!ncol(zeta)) return(zeta)
  
  for (i in seq_len(max(subject_id))) {
    rows_i <- which(subject_id == i)
    N_i <- lfqpca_constraint_null_basis(H[rows_i, , drop = FALSE])
    if (!ncol(N_i)) {
      zeta[rows_i, ] <- 0
    } else {
      zeta[rows_i, ] <- N_i %*% crossprod(N_i, zeta[rows_i, , drop = FALSE])
    }
  }
  zeta
}