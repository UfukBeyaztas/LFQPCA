print.lfqpca_fit <- function(x, ...) {
  cat("Longitudinal Functional Quantile PCA fit\n")
  cat("  Quantile level:             ", x$tau, "\n", sep = "")
  cat("  Subject components:         ", x$K_subject, "\n", sep = "")
  cat("  Visit components:           ", x$K_visit, "\n", sep = "")
  cat("  Iterations:                 ", x$iterations, "\n", sep = "")
  cat("  Converged:                  ", x$converged, "\n", sep = "")
  cat("  Termination:                ", x$termination, "\n", sep = "")
  cat("  Nondegenerate blocks:       ", x$nondegenerate, "\n", sep = "")
  cat("  Rejected block updates:     ",
      x$n_rejected_function_updates + x$n_rejected_score_updates, "\n",
      sep = "")
  cat("  LFQPCA objective:           ",
      format(x$objective, digits = 8), "\n", sep = "")
  cat("  Best start:                 ", x$best_start, "\n", sep = "")
  cat("  Max intercept constraint:   ",
      format(max(x$constraint_diagnostics[, 1]), scientific = TRUE), "\n",
      sep = "")
  cat("  Max slope constraint:       ",
      format(max(x$constraint_diagnostics[, 2]), scientific = TRUE), "\n",
      sep = "")
  cat("  Max all-basis constraint:   ",
      format(max(x$constraint_diagnostics[, 3]), scientific = TRUE), "\n",
      sep = "")
  invisible(x)
}