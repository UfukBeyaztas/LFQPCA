print.summary.lfqpca_fit <- function(x, ...) {
  cat("LFQPCA summary\n")
  cat("  Quantile level:       ", x$tau, "\n", sep = "")
  cat("  Converged:            ", x$converged, "\n", sep = "")
  cat("  Termination:          ", x$termination, "\n", sep = "")
  cat("  Nondegenerate blocks: ", x$nondegenerate, "\n", sep = "")
  cat("  Iterations:           ", x$iterations, "\n", sep = "")
  cat("  LFQPCA objective:     ",
      format(x$objective, digits = 8), "\n", sep = "")
  cat("  Empirical P(Y<=Qhat): ", format(x$empirical_coverage, digits = 6),
      "\n", sep = "")
  cat("  Maximum C'Z errors:   ",
      format(x$maximum_intercept_constraint_error, scientific = TRUE), ", ",
      format(x$maximum_slope_constraint_error, scientific = TRUE), ", ",
      format(x$maximum_all_basis_constraint_error, scientific = TRUE), "\n",
      sep = "")
  cat("  Subject score variances:\n")
  print(x$subject_score_variances)
  cat("  Visit score variances:\n")
  print(x$visit_score_variances)
  invisible(x)
}