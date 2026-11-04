residuals.lfqpca_fit <- function(object, ...) {
  object$Y - object$fitted
}
