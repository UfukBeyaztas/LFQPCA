lfqpca_sim_truth <- function(data, tau, settings = lfqpca_sim_settings()) {
  b <- data$truth
  q_shift <- qnorm(tau) - qnorm(settings$reference_tau)
  phi0 <- b$phi0
  phi1 <- b$phi1
  psi <- b$psi
  xi <- b$xi
  zeta <- b$zeta
  if (b$heterogeneity == "score") {
    xi <- b$xi + q_shift * b$xi_scale
    zeta <- b$zeta + q_shift * b$zeta_scale
  } else if (b$heterogeneity == "loading") {
    phi0 <- b$phi0 + q_shift * b$scale_phi0
    phi1 <- b$phi1 + q_shift * b$scale_phi1
    psi <- b$psi + q_shift * b$scale_psi
  }
  population <- b$location_population + q_shift * b$scale_population
  B_intercept_slope <- b$location_subject_intercept_slope +
    q_shift * b$scale_subject
  B <- B_intercept_slope + b$location_subject_quadratic
  U <- b$location_visit + q_shift * b$scale_visit
  list(
    tau = tau,
    q_shift = q_shift,
    eta = b$eta + q_shift * b$scale_eta,
    phi0 = phi0,
    phi1 = phi1,
    psi = psi,
    xi = xi,
    zeta = zeta,
    population = population,
    B_intercept_slope = B_intercept_slope,
    B = B,
    U = U,
    quantile_surface = population + B + U
  )
}
