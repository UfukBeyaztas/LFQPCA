lfqpca_visit_contribution <- function(zeta, psi) {
  if (!ncol(psi)) return(matrix(0, nrow(zeta), nrow(psi)))
  as.matrix(zeta) %*% t(as.matrix(psi))
}
