lfqpca_weighted_coverage <- function(
    Y,
    fitted,
    objective_weights) {
  
  Y <- as.matrix(Y)
  fitted <- as.matrix(fitted)
  objective_weights <- as.matrix(objective_weights)
  
  if (!all(dim(Y) == dim(fitted)) || !all(dim(Y) == dim(objective_weights))) {
    stop(
      "Y, fitted, and objective_weights must have identical dimensions.",
      call. = FALSE
    )
  }
  
  valid <-
    is.finite(Y) & is.finite(fitted) & is.finite(objective_weights) &
    objective_weights > 0
  
  if (!any(valid)) {
    return(NA_real_)
  }
  
  weights_valid <- objective_weights[valid]
  
  indicators <- as.numeric(Y[valid] <= fitted[valid])
  
  sum(weights_valid * indicators) / sum(weights_valid)
}
