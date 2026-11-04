lfqpca_weighted_check_loss <- function(residual, tau, objective_weights) {
  residual <- as.matrix(residual)
  objective_weights <- as.matrix(objective_weights)
  
  if (!all(dim(residual) == dim(objective_weights))) {
    stop(
      "residual and objective_weights must have identical dimensions.",
      call. = FALSE
    )
  }
  
  valid <- is.finite(residual) & is.finite(objective_weights) &
    objective_weights > 0
  
  if (!any(valid)) return(NA_real_)
  
  sum(objective_weights[valid] * lfqpca_check_loss_values(residual[valid], tau))
}
