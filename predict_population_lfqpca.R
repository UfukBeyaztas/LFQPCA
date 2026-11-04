predict_population_lfqpca <- function(object, new_visit_time) {
  if (!inherits(object, "lfqpca_fit")) {
    stop("object must inherit from 'lfqpca_fit'.", call. = FALSE)
  }
  u <- (as.numeric(new_visit_time) - object$time_center) / object$time_scale
  H <- lfqpca_time_basis(u, degree = object$mean_time_degree)
  lfqpca_mean_surface(H, object$eta_functions)
}