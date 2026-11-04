lfqpca_sim_score_multiplier <- function(settings) {
  q_cal <- qnorm(settings$calibration_tau) - qnorm(settings$reference_tau)
  sqrt(1 / settings$score_correlation^2 - 1) / abs(q_cal)
}