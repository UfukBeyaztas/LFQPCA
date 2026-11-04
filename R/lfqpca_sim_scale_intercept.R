lfqpca_sim_scale_intercept <- function(phi0, phi1, psi, subject_coef,
                                       visit_coef, J_design, settings) {
  J_use <- settings$bound_visits
  if (J_design > J_use) {
    stop("The design has more than ", J_use, " visits per subject.",
         call. = FALSE)
  }
  leverage <- sqrt(1 - 1 / J_use)
  subject_profile <- abs(phi0) + settings$max_abs_time * abs(phi1)
  visit_bound <- visit_coef * sqrt(J_use) * leverage
  settings$min_scale +
    max(subject_coef * subject_profile + visit_bound * abs(psi))
}
