lfqpca_sym_power <- function(A, power, eps = 1e-10) {
  A <- (as.matrix(A) + t(as.matrix(A))) / 2
  eg <- eigen(A, symmetric = TRUE)
  values <- pmax(eg$values, eps)^power
  eg$vectors %*% diag(values, length(values), length(values)) %*% t(eg$vectors)
}
