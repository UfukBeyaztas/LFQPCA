lfqpcr_build_predictor_spaces <- function(lfqpca_fit, tol = 1e-10) {
  lfqpcr_validate_lfqpca_fit(lfqpca_fit)
  
  w_tilde <- lfqpca_trapezoid_weights(lfqpca_fit$grid, normalize = FALSE)
  
  F_B <- cbind(as.matrix(lfqpca_fit$phi0), as.matrix(lfqpca_fit$phi1))
  
  F_U <- as.matrix(lfqpca_fit$psi)
  
  B_space <- lfqpcr_orthonormal_space_basis(
    generating_functions = F_B,
    integration_weights = w_tilde,
    tol = tol
  )
  
  U_space <- lfqpcr_orthonormal_space_basis(
    generating_functions = F_U,
    integration_weights = w_tilde,
    tol = tol
  )
  
  B_error <- B_space$orthonormality_error
  U_error <- U_space$orthonormality_error
  
  list(
    chi_B = B_space$basis,
    chi_U = U_space$basis,
    r_B = B_space$rank,
    r_U = U_space$rank,
    integration_weights = w_tilde,
    B_orthonormality_error = if (is.null(B_error)) 0 else B_error,
    U_orthonormality_error = if (is.null(U_error)) 0 else U_error
  )
}
