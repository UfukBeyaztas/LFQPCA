summary.lfqpca_fit <- function(object, ...) {
  output <- list(
    tau = object$tau,
    converged = object$converged,
    termination = object$termination,
    nondegenerate = object$nondegenerate,
    rank_diagnostics = object$rank_diagnostics,
    iterations = object$iterations,
    objective = object$objective,
    empirical_coverage = lfqpca_weighted_coverage(
      Y = object$Y,
      fitted = object$fitted,
      objective_weights = object$objective_weights
    ),
    subject_score_variances = object$lambda,
    visit_score_variances = object$nu,
    subject_pve = object$pve_subject,
    visit_pve = object$pve_visit,
    maximum_intercept_constraint_error =
      max(object$constraint_diagnostics[, 1]),
    maximum_slope_constraint_error = max(object$constraint_diagnostics[, 2]),
    maximum_all_basis_constraint_error =
      max(object$constraint_diagnostics[, 3]),
    all_start_objectives = object$all_start_objectives,
    best_start = object$best_start
  )
  class(output) <- "summary.lfqpca_fit"
  output
}