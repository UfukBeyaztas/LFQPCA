lfqpcr_coefficient_functions <- function(object) {
  if (!inherits(object, "lfqpcr_fit")) {
    stop("object must inherit from 'lfqpcr_fit'.", call. = FALSE)
  }
  
  data.frame(
    grid = object$lfqpca_fit$grid,
    gamma_B = object$gamma_B,
    gamma_U = object$gamma_U
  )
}
