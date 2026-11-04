lfqpcr_reconstruct_components <- function(lfqpca_fit) {
  lfqpcr_validate_lfqpca_fit(lfqpca_fit)
  
  B_hat <- lfqpca_subject_contribution(
    subject_id = lfqpca_fit$subject_id,
    time_std = lfqpca_fit$time_std,
    xi = lfqpca_fit$xi,
    phi0 = lfqpca_fit$phi0,
    phi1 = lfqpca_fit$phi1
  )
  
  U_hat <- lfqpca_visit_contribution(
    zeta = lfqpca_fit$zeta,
    psi = lfqpca_fit$psi
  )
  
  list(B = B_hat, U = U_hat)
}