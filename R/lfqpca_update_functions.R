lfqpca_update_functions <- function(
    Y,
    B,
    H,
    subject_id,
    time_std,
    xi,
    zeta,
    tau,
    objective_weights,
    rq_method = "fnb") {
  
  mask <- is.finite(Y)
  
  design <- lfqpca_build_function_design(
    B = B,
    H = H,
    subject_id = subject_id,
    time_std = time_std,
    xi = xi,
    zeta = zeta,
    mask = mask,
    objective_weights = objective_weights
  )
  
  y <- Y[design$index]
  
  beta <- lfqpca_quantile_fit(
    X = design$X,
    y = y,
    tau = tau,
    case_weights = design$case_weights,
    method = rq_method
  )
  
  q <- ncol(B)
  M <- ncol(H)
  Kx <- ncol(xi)
  Ku <- ncol(zeta)
  
  Theta <- matrix(beta, nrow = q, ncol = M + 2 * Kx + Ku)
  
  eta_index <- seq_len(M)
  phi0_index <- M + seq_len(Kx)
  phi1_index <- M + Kx + seq_len(Kx)
  psi_index <- M + 2 * Kx + seq_len(Ku)
  
  theta_eta <- Theta[, eta_index, drop = FALSE]
  
  theta_phi0 <-
    if (Kx) {
      Theta[, phi0_index, drop = FALSE]
    } else {
      matrix(0, q, 0)
    }
  
  theta_phi1 <-
    if (Kx) {
      Theta[, phi1_index, drop = FALSE]
    } else {
      matrix(0, q, 0)
    }
  
  theta_psi <-
    if (Ku) {
      Theta[, psi_index, drop = FALSE]
    } else {
      matrix(0, q, 0)
    }
  
  list(
    eta_functions = B %*% theta_eta,
    phi0 = B %*% theta_phi0,
    phi1 = B %*% theta_phi1,
    psi = B %*% theta_psi,
    theta_eta = theta_eta,
    theta_phi0 = theta_phi0,
    theta_phi1 = theta_phi1,
    theta_psi = theta_psi
  )
}
