lfqpca_state_objective <- function(
    Y,
    H,
    eta_functions,
    subject_id,
    time_std,
    xi,
    phi0,
    phi1,
    zeta,
    psi,
    tau,
    objective_weights) {
  
  fitted <- lfqpca_reconstruct(
    time_basis = H,
    eta_functions = eta_functions,
    subject_id = subject_id,
    time_std = time_std,
    xi = xi,
    phi0 = phi0,
    phi1 = phi1,
    zeta = zeta,
    psi = psi
  )
  
  residual <- Y - fitted
  
  objective <- lfqpca_weighted_check_loss(
    residual = residual,
    tau = tau,
    objective_weights = objective_weights
  )
  
  list(fitted = fitted, objective = objective)
}
