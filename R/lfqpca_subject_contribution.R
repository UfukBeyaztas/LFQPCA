lfqpca_subject_contribution <- function(subject_id, time_std, xi, phi0, phi1) {
  n <- length(subject_id)
  D <- nrow(phi0)
  Kx <- ncol(phi0)
  if (!Kx) return(matrix(0, n, D))
  
  Xi_full <- xi[subject_id, , drop = FALSE]
  Xi_full %*% t(phi0) + (Xi_full * time_std) %*% t(phi1)
}
