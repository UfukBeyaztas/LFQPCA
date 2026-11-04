lfqpca_mean_surface <- function(time_basis, eta_functions) {
  as.matrix(time_basis) %*% t(as.matrix(eta_functions))
}