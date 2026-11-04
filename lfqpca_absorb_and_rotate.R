lfqpca_absorb_and_rotate <- function(
    eta_functions,
    xi,
    phi0,
    phi1,
    zeta,
    psi,
    weights) {
  
  D <- nrow(phi0)
  Kx <- ncol(phi0)
  Ku <- ncol(psi)
  
  if (Kx) {
    xi_mean <- colMeans(xi)
    shift0 <- as.numeric(phi0 %*% xi_mean)
    shift1 <- as.numeric(phi1 %*% xi_mean)
    
    rot_x <- lfqpca_rotate_block(
      scores = xi,
      loadings = rbind(phi0, phi1),
      weights = rep(weights, 2),
      center_scores = TRUE
    )
    
    pair <- rot_x$loadings
    phi0 <- pair[seq_len(D), , drop = FALSE]
    phi1 <- pair[D + seq_len(D), , drop = FALSE]
    xi <- rot_x$scores
    
    eta_functions[, 1] <- eta_functions[, 1] + shift0
    eta_functions[, 2] <- eta_functions[, 2] + shift1
    lambda <- rot_x$variances
  } else {
    lambda <- numeric(0)
  }
  
  if (Ku) {
    rot_u <- lfqpca_rotate_block(
      scores = zeta,
      loadings = psi,
      weights = weights,
      center_scores = FALSE
    )
    zeta <- rot_u$scores
    psi <- rot_u$loadings
    nu <- rot_u$variances
  } else {
    nu <- numeric(0)
  }
  
  list(
    eta_functions = eta_functions,
    xi = xi,
    phi0 = phi0,
    phi1 = phi1,
    zeta = zeta,
    psi = psi,
    lambda = lambda,
    nu = nu
  )
}