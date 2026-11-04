predict.lfqpcr_fit <- function(
    object,
    newY,
    new_subject,
    new_visit_time,
    new_covariates = NULL,
    trajectory_only = FALSE,
    ...) {
  
  if (!inherits(object, "lfqpcr_fit")) {
    stop("object must inherit from 'lfqpcr_fit'.", call. = FALSE)
  }
  
  pred_F <- predict(
    object$lfqpca_fit,
    newY = newY,
    new_subject = new_subject,
    new_visit_time = new_visit_time
  )
  
  coords <- lfqpcr_coordinates_from_scores(
    subject_id = pred_F$subject_id,
    time_std = pred_F$time_std,
    xi = pred_F$xi,
    zeta = pred_F$zeta,
    phi0 = object$lfqpca_fit$phi0,
    phi1 = object$lfqpca_fit$phi1,
    psi = object$lfqpca_fit$psi,
    chi_B = object$chi_B,
    chi_U = object$chi_U,
    integration_weights = object$integration_weights
  )
  
  if (isTRUE(trajectory_only) && ncol(coords$x_U)) {
    coords$x_U[,] <- 0
  }
  
  design <- lfqpcr_build_stage2_design(
    time_std = pred_F$time_std,
    x_B = coords$x_B,
    x_U = coords$x_U,
    covariates = new_covariates,
    response_time_degree = object$response_time_degree,
    minimum_time_degree = object$lfqpca_fit$mean_time_degree
  )
  
  if (ncol(design$X) != length(object$coefficients)) {
    stop(
      "The new LFQPCR design is incompatible with the fitted model. ",
      "Check the scalar-covariate columns and their order.",
      call. = FALSE
    )
  }
  
  fitted <- as.numeric(design$X %*% object$coefficients)
  
  list(
    fitted = fitted,
    x_B = coords$x_B,
    x_U = coords$x_U,
    B_hat = coords$B,
    U_hat = if (trajectory_only) {
      matrix(0, nrow(coords$U), ncol(coords$U))
    } else {
      coords$U
    },
    xi = pred_F$xi,
    zeta = pred_F$zeta,
    time_std = pred_F$time_std,
    subject_id = pred_F$subject_id,
    trajectory_only = isTRUE(trajectory_only)
  )
}
