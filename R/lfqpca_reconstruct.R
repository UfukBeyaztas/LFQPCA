lfqpca_reconstruct <- function(time_basis, eta_functions,
                               subject_id, time_std,
                               xi, phi0, phi1,
                               zeta, psi) {
  lfqpca_mean_surface(time_basis, eta_functions) +
    lfqpca_subject_contribution(subject_id, time_std, xi, phi0, phi1) +
    lfqpca_visit_contribution(zeta, psi)
}
