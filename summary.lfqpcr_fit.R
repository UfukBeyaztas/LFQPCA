summary.lfqpcr_fit <- function(object, ...) {
  eligible_residuals <- object$residuals[object$eligible]
  eligible_weights <- object$response_weights[object$eligible]
  
  weighted_coverage <- sum(
    eligible_weights *
      as.numeric(object$V[object$eligible] <= object$fitted[object$eligible])
  ) / sum(eligible_weights)
  
  out <- list(
    tau_F = object$tau_F,
    tau_R = object$tau_R,
    I_R = object$I_R,
    r_B = object$r_B,
    r_U = object$r_U,
    design_rank = object$design_rank,
    n_stage2 = sum(object$eligible),
    objective = object$objective,
    weighted_coverage = weighted_coverage,
    B_orthonormality_error = object$B_orthonormality_error,
    U_orthonormality_error = object$U_orthonormality_error,
    dropped_columns = if (length(object$dropped_columns)) {
      colnames(object$design)[object$dropped_columns]
    } else {
      character(0)
    },
    coefficients = object$coefficients
  )
  
  class(out) <- "summary.lfqpcr_fit"
  out
}