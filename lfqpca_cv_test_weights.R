lfqpca_cv_test_weights <- function(
    Y,
    held_out,
    subject_id,
    quadrature_weights) {
  
  Y <- as.matrix(Y)
  held_out <- as.matrix(held_out)
  
  n <- nrow(Y)
  D <- ncol(Y)
  
  if (!all(dim(held_out) == c(n, D))) {
    stop("held_out must have the same dimensions as Y.", call. = FALSE)
  }
  
  I <- max(subject_id)
  visit_counts <- tabulate(subject_id, nbins = I)
  
  objective_weights <- matrix(0, nrow = n, ncol = D)
  
  for (r in seq_len(n)) {
    
    test_points <- which(held_out[r,] & is.finite(Y[r,]))
    
    if (!length(test_points)) {
      next
    }
    
    denominator <- sum(quadrature_weights[test_points])
    
    if (!is.finite(denominator) || denominator <= 0) {
      stop(
        paste0(
          "Held-out quadrature weights have nonpositive total for curve ",
          r,
          "."
        ),
        call. = FALSE
      )
    }
    
    omega_test <- quadrature_weights[test_points] / denominator
    
    i <- subject_id[r]
    J_i <- visit_counts[i]
    
    objective_weights[r, test_points] <- (1 / I) * (1 / J_i) * omega_test
  }
  
  objective_weights
}