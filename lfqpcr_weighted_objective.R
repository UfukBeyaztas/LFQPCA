lfqpcr_weighted_objective <- function(residual, tau_R, weights) {
  residual <- as.numeric(residual)
  weights <- as.numeric(weights)
  
  if (length(residual) != length(weights)) {
    stop("residual and weights must have the same length.", call. = FALSE)
  }
  
  valid <- is.finite(residual) & is.finite(weights) & weights > 0
  
  if (!any(valid)) {
    return(NA_real_)
  }
  
  sum(weights[valid] * lfqpca_check_loss_values(residual[valid], tau_R))
}