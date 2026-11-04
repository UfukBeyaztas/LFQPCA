lfqpcr_coordinates_from_scores <- function(
    subject_id,
    time_std,
    xi,
    zeta,
    phi0,
    phi1,
    psi,
    chi_B,
    chi_U,
    integration_weights) {
  
  B_curves <- lfqpca_subject_contribution(
    subject_id = subject_id,
    time_std = time_std,
    xi = xi,
    phi0 = phi0,
    phi1 = phi1
  )
  
  U_curves <- lfqpca_visit_contribution(zeta = zeta, psi = psi)
  
  coords <- lfqpcr_coordinates_from_curves(
    B_curves = B_curves,
    U_curves = U_curves,
    chi_B = chi_B,
    chi_U = chi_U,
    integration_weights = integration_weights
  )
  
  coords$B <- B_curves
  coords$U <- U_curves
  coords
}