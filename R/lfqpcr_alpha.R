lfqpcr_alpha <- function(object, visit_time) {
  if (!inherits(object, "lfqpcr_fit")) {
    stop("object must inherit from 'lfqpcr_fit'.", call. = FALSE)
  }
  
  u <-
    (as.numeric(visit_time) - object$lfqpca_fit$time_center) /
    object$lfqpca_fit$time_scale
  
  G <- lfqpcr_response_time_basis(
    time_std = u,
    degree = object$response_time_degree,
    minimum_degree = object$lfqpca_fit$mean_time_degree
  )
  
  as.numeric(G %*% object$a_alpha)
}
